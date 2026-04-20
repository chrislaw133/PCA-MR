#!/bin/bash
#SBATCH --array=1-100
#SBATCH --mem=4G
#SBATCH --time=00-1:00:00

set -euo pipefail

module load apps/r/4.3.3

Rscript --vanilla - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(MASS)
})

############################################################
## PATHS
############################################################
gene_file <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/cumexpr_stuff/snpcorrelations/genesandinstruments.txt"
beta_dir  <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/cumexpr_stuff/betavectors"
ld_dir    <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/cumexpr_stuff/snpcorrelations"

out_dir   <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/results/cum_pca_methods"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

############################################################
## READ + STANDARDIZE GENES TABLE
############################################################
genes <- fread(gene_file)
setDT(genes)

.pick_col <- function(dt, candidates) {
  nms <- names(dt)
  low <- tolower(nms)
  cand_low <- tolower(candidates)
  hit <- match(cand_low, low)
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) return(NULL)
  nms[hit[1]]
}

gene_col <- .pick_col(genes, c("gene", "genes", "symbol", "gene_symbol", "hgnc", "id"))
chr_col  <- .pick_col(genes, c("chr", "chrom", "chromosome", "chromosome_number"))
nx_col   <- .pick_col(genes, c("n_x", "nx", "sample_size_x", "nexp", "n_exp", "nexposure", "n_exposure", "n"))

if (is.null(gene_col) || is.null(chr_col)) {
  stop(
    paste0(
      "Could not find required columns in genes file.\n",
      "Found columns: ", paste(names(genes), collapse = ", "), "\n",
      "Need something like gene + chr."
    )
  )
}

if (gene_col != "gene") setnames(genes, gene_col, "gene")
if (chr_col  != "chr")  setnames(genes, chr_col,  "chr")
if (!is.null(nx_col) && nx_col != "N_X") setnames(genes, nx_col, "N_X")

genes[, gene := as.character(gene)]
genes[, chr := as.integer(gsub("[^0-9]+", "", as.character(chr)))]
if ("N_X" %in% names(genes)) genes[, N_X := as.integer(N_X)]

genes <- genes[is.finite(chr) & !is.na(chr) & nzchar(gene)]
genes <- unique(genes, by = c("gene", "chr"))
setorder(genes, chr, gene)

############################################################
## SLURM ARRAY INFO
############################################################
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("SLURM_ARRAY_TASK_ID not set")

n_tasks <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
if (is.na(n_tasks) || n_tasks < 1) n_tasks <- 100

gene_idx <- seq_len(nrow(genes))
assigned_task <- ((gene_idx - 1) %% n_tasks) + 1
genes_task <- genes[assigned_task == task_id]

cat("Array task", task_id, "of", n_tasks, "processing", nrow(genes_task), "genes\n")

job_id <- Sys.getenv("SLURM_JOB_ID")
if (is.na(job_id) || job_id == "") job_id <- "nojobid"
out_file <- file.path(out_dir, paste0("pca_methods_", job_id, "_part_", task_id, ".tsv"))

############################################################
## TUNABLES
## kept aligned with your pasted PCA simulation code
############################################################
KEEP_VAR_IVW   <- 0.999
KEEP_VAR_EGGER <- 0.999
KEEP_VAR_ML    <- 0.999

F_MIN <- 10

MIN_MODES_IVW   <- 2
MIN_MODES_EGGER <- 3
MIN_MODES_ML    <- 2

KAPPA_MAX <- 100
ML_MODEL  <- "random"
IVW_MODEL <- "random"

N_X_DEFAULT <- 117

############################################################
## HELPERS
############################################################
load_ld_matrix <- function(gene, chr) {
  f <- file.path(ld_dir, paste0(gene, "_chr", chr, ".unphased.vcor1.zst"))
  if (!file.exists(f)) return(NULL)

  tmp <- tempfile(fileext = ".vars")
  on.exit(unlink(tmp), add = TRUE)

  system(paste("unzstd -c", shQuote(f), ">", shQuote(tmp)))
  ld <- tryCatch(as.matrix(fread(tmp)), error = function(e) NULL)
  if (is.null(ld)) return(NULL)

  ld <- (ld + t(ld)) / 2
  diag(ld) <- 1
  ld
}

compute_zeta <- function(LD) {
  LD <- as.matrix(LD)
  eig <- pmax(svd(LD)$d, 1e-12)

  total <- sum(eig)
  if (total <= 0) return(NA_real_)

  p_i <- eig / total
  H <- -sum(p_i * log(p_i))
  r_eff <- exp(H)
  r_eff / ncol(LD)
}

condnum <- function(M) {
  ev <- tryCatch(
    svd(M, nu = 0, nv = 0)$d,
    error = function(e) NA_real_
  )
  if (any(!is.finite(ev))) return(Inf)
  ev <- pmax(ev, 1e-12)
  max(ev) / min(ev)
}

############################################################
## EXACT PCA/BURGESS METHOD LOGIC FROM YOUR SIM SCRIPT
############################################################
gls_egger <- function(bx_tilde, by_tilde, Omega, lambda) {
  n <- length(bx_tilde)
  p <- 2L

  if (length(by_tilde) != n) stop("gls_egger: length mismatch")
  if (!is.matrix(Omega) || any(dim(Omega) != c(n, n))) {
    stop("gls_egger: Omega has wrong dimension")
  }

  if (n < 3) {
    return(list(
      intercept = NA_real_, slope = NA_real_,
      se_intercept = NA_real_, se_slope = NA_real_,
      p_intercept = NA_real_, p_slope = NA_real_,
      Q = NA_real_, df = NA_real_, rse = NA_real_, Q_pval = NA_real_
    ))
  }

  sgn <- ifelse(bx_tilde >= 0, 1, -1)
  bx_tilde <- abs(bx_tilde)
  by_tilde <- sgn * by_tilde
  S <- diag(sgn, nrow = n, ncol = n)
  Omega <- S %*% Omega %*% S

  X <- cbind(1, bx_tilde)
  y <- matrix(by_tilde, ncol = 1)

  Omega_inv <- tryCatch(solve(Omega), error = function(e) NULL)
  if (is.null(Omega_inv)) {
    return(list(
      intercept = NA_real_, slope = NA_real_,
      se_intercept = NA_real_, se_slope = NA_real_,
      p_intercept = NA_real_, p_slope = NA_real_,
      Q = NA_real_, df = NA_real_, rse = NA_real_, Q_pval = NA_real_
    ))
  }
  
  Lambda <- diag(lambda, nrow = length(lambda), ncol = length(lambda))
  Omega_inv <- Lambda %*% Omega_inv

  XtOiX <- t(X) %*% Omega_inv %*% X
  bread <- tryCatch(solve(XtOiX), error = function(e) NULL)
  if (is.null(bread)) {
    return(list(
      intercept = NA_real_, slope = NA_real_,
      se_intercept = NA_real_, se_slope = NA_real_,
      p_intercept = NA_real_, p_slope = NA_real_,
      Q = NA_real_, df = NA_real_, rse = NA_real_, Q_pval = NA_real_
    ))
  }

  beta_hat <- bread %*% t(X) %*% Omega_inv %*% y
  resid <- y - X %*% beta_hat

  Q <- as.numeric(t(resid) %*% Omega_inv %*% resid)
  df <- n - 2

  if (!is.finite(df) || df <= 0) {
    return(list(
      intercept = NA_real_, slope = NA_real_,
      se_intercept = NA_real_, se_slope = NA_real_,
      p_intercept = NA_real_, p_slope = NA_real_,
      Q = Q, df = df, rse = NA_real_, Q_pval = NA_real_
    ))
  }

  rse <- sqrt(Q / df)
  infl <- max(1, rse)

  se_intercept <- sqrt(bread[1, 1]) * infl
  se_slope     <- sqrt(bread[2, 2]) * infl

  intercept_hat <- as.numeric(beta_hat[1, 1])
  slope_hat     <- as.numeric(beta_hat[2, 1])

  t_int   <- intercept_hat / se_intercept
  t_slope <- slope_hat / se_slope

  p_intercept <- 2 * pt(abs(t_int), df = df, lower.tail = FALSE)
  p_slope     <- 2 * pt(abs(t_slope), df = df, lower.tail = FALSE)
  Q_pval <- pchisq(Q, df = df, lower.tail = FALSE)

  list(
    intercept = intercept_hat,
    slope = slope_hat,
    se_intercept = se_intercept,
    se_slope = se_slope,
    p_intercept = p_intercept,
    p_slope = p_slope,
    Q = Q,
    df = df,
    rse = rse,
    Q_pval = Q_pval
  )
}

burgess_pca_ivw <- function(bx, by, byse, ld, var_expl = 0.999) {
  p <- length(bx)
  if (p < 2) return(NULL)

  z <- bx / byse
  z[!is.finite(z)] <- 0
  Phi <- (z %o% z) * ld

  pca <- prcomp(Phi, scale. = FALSE)
  eigvals <- pca$sdev^2

  if (!all(is.finite(eigvals)) || sum(eigvals) <= 0) return(NULL)

  K <- which(cumsum(eigvals) / sum(eigvals) >= var_expl)[1]
  if (is.na(K) || K < 1) K <- 1

  R <- pca$rotation[, 1:K, drop = FALSE]

  bx0 <- as.numeric(t(R) %*% bx)
  by0 <- as.numeric(t(R) %*% by)

  Omega <- (byse %o% byse) * ld
  pcOmega <- t(R) %*% Omega %*% R

  Omega_inv <- tryCatch(solve(pcOmega), error = function(e) NULL)
  if (is.null(Omega_inv)) return(NULL)

  denom <- as.numeric(t(bx0) %*% Omega_inv %*% bx0)
  if (!is.finite(denom) || denom <= 0) return(NULL)

  slope <- as.numeric((t(bx0) %*% Omega_inv %*% by0) / denom)
  se <- sqrt(1 / denom)
  pval <- 2 * pnorm(-abs(slope / se))

  list(
    slope = slope,
    se = se,
    p = pval,
    K = K
  )
}

pca_prepare <- function(bx, by, sex, sey, ld, N_X,
                        keep_var = KEEP_VAR_IVW,
                        F_min = F_MIN,
                        min_modes = MIN_MODES_IVW,
                        kappa_max = KAPPA_MAX) {

  eig <- svd(ld)
  Q <- eig$u
  lambda <- eig$d

  bx_pc <- as.numeric(crossprod(Q, bx))
  by_pc <- as.numeric(crossprod(Q, by))

  SigX <- (sex %o% sex) * ld
  SigY <- (sey %o% sey) * ld
  SigX_pc <- t(Q) %*% SigX %*% Q
  SigY_pc <- t(Q) %*% SigY %*% Q

  cumvar <- cumsum(lambda) / sum(lambda)
  K0 <- which(cumvar >= keep_var)[1]
  if (is.na(K0)) K0 <- length(lambda)
  keep <- seq_len(K0)

  if (length(keep) < min_modes) {
    keep <- seq_len(min(min_modes, length(lambda)))
  }
  if (length(keep) < min_modes) return(NULL)

  SigmaX <- SigX_pc[keep, keep, drop = FALSE]
  SigmaY <- SigY_pc[keep, keep, drop = FALSE]
  se_bx <- sqrt(diag(SigmaX))
  se_by <- sqrt(diag(SigmaY))

  R2 <- bx_pc[keep]^2 / (bx_pc[keep]^2 + N_X * se_bx^2)
  Fpc <- (N_X - 2) * R2 / (1 - R2)

  keep <- keep[Fpc > F_min]
  if (length(keep) < min_modes) return(NULL)
  if (length(keep) < min_modes) return(NULL)

  Omega <- SigY_pc[keep, keep, drop = FALSE]
  kappa <- condnum(Omega)

  while (length(keep) > min_modes && is.finite(kappa) && kappa > kappa_max) {
    keep <- keep[-length(keep)]
    Omega <- SigY_pc[keep, keep, drop = FALSE]
    kappa <- condnum(Omega)
  }

  if (length(keep) < min_modes) return(NULL)
  if (!is.finite(kappa) || kappa > kappa_max) return(NULL)

  SigmaX <- SigX_pc[keep, keep, drop = FALSE]
  SigmaY <- SigY_pc[keep, keep, drop = FALSE]

  se_bx <- sqrt(diag(SigmaX))
  se_by <- sqrt(diag(SigmaY))

  R2_first_stage <- bx^2 / (bx^2 + N_X * sex^2)
  F_first_stage <- (N_X - 2) * R2_first_stage / pmax(1 - R2_first_stage, 1e-12)
  F_first_stage[!is.finite(F_first_stage)] <- NA_real_
  mean_F_first_stage <- mean(F_first_stage, na.rm = TRUE)

  list(
    bx = bx_pc[keep],
    by = by_pc[keep],
    bxse = se_bx,
    byse = se_by,
    Omega = Omega,
    SigmaX = SigmaX,
    SigmaY = SigmaY,
    keep = keep,
    n_modes = length(keep),
    mean_F = mean(Fpc[match(keep, seq_len(K0))], na.rm = TRUE),
    mean_F_marginal = mean_F_first_stage,
    kappa = kappa,
    lambda = lambda[keep]
  )
}

pca_gls <- function(pc, model = IVW_MODEL, phi_floor = 1) {
  if (is.null(pc)) return(NULL)
  if (!model %in% c("fixed", "random")) stop("model must be 'fixed' or 'random'")

  Omega_inv <- tryCatch(solve(pc$Omega), error = function(e) NULL)
  if (is.null(Omega_inv)) return(NULL)

  Lambda <- diag(pc$lambda, nrow = length(pc$lambda), ncol = length(pc$lambda))
  W <- Lambda %*% Omega_inv

  denom <- as.numeric(t(pc$bx) %*% W %*% pc$bx)
  if (!is.finite(denom) || denom <= 0) return(NULL)

  slope <- as.numeric((t(pc$bx) %*% W %*% pc$by) / denom)

  resid <- as.numeric(pc$by - slope * pc$bx)
  Q <- as.numeric(t(resid) %*% W %*% resid)
  df <- length(pc$bx) - 1

  if (!is.finite(df) || df <= 0) return(NULL)

  phi <- if (model == "fixed") 1 else max(phi_floor, Q / df)
  se <- sqrt(phi / denom)
  p <- 2 * pnorm(-abs(slope / se))
  Q_pvalue <- pchisq(Q, df = df, lower.tail = FALSE)
  rse <- sqrt(Q / df)

  list(
    slope = slope,
    se = se,
    p = p,
    Q = Q,
    df = df,
    phi = phi,
    rse = rse,
    model = model,
    Q_pvalue = Q_pvalue
  )
}

pca_maxlik <- function(pc, model = ML_MODEL, alpha = 0.05) {
  if (is.null(pc)) return(NULL)

  loglikelihoodcorrel <- function(param, x, Taux, y, Tauy) {
    k <- length(x)
    b <- param[1:k]
    theta <- param[k + 1]

    as.numeric(
      0.5 * t(x - b) %*% Taux %*% (x - b) +
      0.5 * t(y - theta * b) %*% Tauy %*% (y - theta * b)
    )
  }

  bx_pc <- as.numeric(pc$bx)
  by_pc <- as.numeric(pc$by)
  SigmaX <- pc$SigmaX
  SigmaY <- pc$SigmaY

  k <- length(bx_pc)
  if (k < 1) return(NULL)
  if (model == "random" && k < 2) return(NULL)

  Lambda <- diag(pc$lambda, nrow = length(pc$lambda), ncol = length(pc$lambda))

  Taux <- tryCatch(solve(SigmaX), error = function(e) NULL)
  Taux <- Lambda %*% Taux
  Tauy <- tryCatch(solve(SigmaY), error = function(e) NULL)
  Tauy <- Lambda %*% Tauy
  if (is.null(Taux) || is.null(Tauy)) return(NULL)

  init_gls <- tryCatch(pca_gls(pc, model = model), error = function(e) NULL)
  theta0 <- if (!is.null(init_gls) && is.finite(init_gls$slope)) init_gls$slope else 0

  par0 <- c(bx_pc, theta0)

  opt <- tryCatch(
    optim(
      par = par0,
      fn = loglikelihoodcorrel,
      x = bx_pc,
      Taux = Taux,
      y = by_pc,
      Tauy = Tauy,
      hessian = TRUE,
      method = "BFGS",
      control = list(maxit = 25000, reltol = 1e-10)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || is.null(opt$par) || is.null(opt$hessian)) return(NULL)
  if (!is.finite(opt$value)) return(NULL)

  H <- opt$hessian
  if (is.null(H)) return(NULL)

  V <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V)) return(NULL)

  theta <- as.numeric(opt$par[k + 1])
  theta_var_fixed <- as.numeric(V[k + 1, k + 1])
  if (!is.finite(theta_var_fixed) || theta_var_fixed < 0) return(NULL)

  if (model == "fixed") {
    theta_se <- sqrt(theta_var_fixed)
  } else {
    phi <- max(2 * opt$value / (k - 1), 1)
    theta_se <- sqrt(theta_var_fixed * phi)
  }

  if (!is.finite(theta_se) || theta_se <= 0) return(NULL)

  zval <- theta / theta_se
  pval <- 2 * pnorm(-abs(zval))
  ci_lo <- theta - qnorm(1 - alpha / 2) * theta_se
  ci_hi <- theta + qnorm(1 - alpha / 2) * theta_se

  heter_stat <- 2 * opt$value
  heter_df <- k - 1
  heter_p <- if (heter_df > 0) pchisq(heter_stat, df = heter_df, lower.tail = FALSE) else NA_real_
  rse <- if (heter_df > 0) sqrt(heter_stat / heter_df) else NA_real_

  list(
    slope = theta,
    se = theta_se,
    p = pval,
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    K = k,
    model = model,
    rse = rse,
    heter_stat = heter_stat,
    heter_df = heter_df,
    heter_p = heter_p,
    conv = opt$convergence,
    value = opt$value
  )
}

pca_mr_gegger_gene <- function(bx, by, sex, sey, ld, N_X,
                               keep_var = KEEP_VAR_EGGER,
                               F_min = F_MIN,
                               kappa_max = KAPPA_MAX,
                               min_modes = MIN_MODES_EGGER) {

  eig <- svd(ld)
  Q <- eig$u
  Lambda <- eig$d

  var_frac <- Lambda / sum(Lambda)
  K0 <- which(cumsum(var_frac) >= keep_var)[1]
  if (is.na(K0) || K0 < 1) K0 <- 1
  keep0 <- seq_len(K0)

  if (length(keep0) < min_modes) {
    keep0 <- seq_len(min(min_modes, length(Lambda)))
  }
  if (length(keep0) < min_modes) return(NULL)

  bx_all <- as.numeric(crossprod(Q, bx))
  by_all <- as.numeric(crossprod(Q, by))

  Sigma_y <- (sey %o% sey) * ld
  Sigma_x <- (sex %o% sex) * ld
  Sigma_y_tilde <- t(Q) %*% Sigma_y %*% Q
  Sigma_x_tilde <- t(Q) %*% Sigma_x %*% Q

  bx0 <- bx_all[keep0]
  SigX0 <- Sigma_x_tilde[keep0, keep0, drop = FALSE]
  se_bx0 <- sqrt(diag(SigX0))

  R2_mode0 <- (bx0^2) / (bx0^2 + N_X * (se_bx0^2))
  F_mode0 <- (N_X - 2) * (R2_mode0 / (1 - R2_mode0))

  keepF <- keep0[which(F_mode0 > F_min)]
  n_after_F <- length(keepF)

  if (length(keepF) < min_modes) {
    ord <- order(F_mode0, decreasing = TRUE)
    take <- ord[seq_len(min(min_modes, length(ord)))]
    keepF <- keep0[take]
  }

  keepF <- sort(unique(keepF))
  if (length(keepF) < min_modes) return(NULL)

  bx1 <- bx_all[keepF]
  by1 <- by_all[keepF]
  Omega1 <- Sigma_y_tilde[keepF, keepF, drop = FALSE]
  SigX1 <- Sigma_x_tilde[keepF, keepF, drop = FALSE]
  se_bx1 <- sqrt(diag(SigX1))

  R2_mode1 <- (bx1^2) / (bx1^2 + N_X * (se_bx1^2))
  F_mode1 <- (N_X - 2) * (R2_mode1 / (1 - R2_mode1))

  kappa <- condnum(Omega1)

  while (length(keepF) > min_modes && is.finite(kappa) && kappa > kappa_max) {
    keepF <- keepF[-length(keepF)]
    bx1 <- bx_all[keepF]
    by1 <- by_all[keepF]
    Omega1 <- Sigma_y_tilde[keepF, keepF, drop = FALSE]
    SigX1 <- Sigma_x_tilde[keepF, keepF, drop = FALSE]
    se_bx1 <- sqrt(diag(SigX1))
    R2_mode1 <- (bx1^2) / (bx1^2 + N_X * (se_bx1^2))
    F_mode1 <- (N_X - 2) * (R2_mode1 / (1 - R2_mode1))
    kappa <- condnum(Omega1)
  }

  if (length(keepF) < min_modes) return(NULL)
  if (!is.finite(kappa) || kappa > kappa_max) return(NULL)

  R2_total <- sum(R2_mode1, na.rm = TRUE)
  k <- length(keepF)

  F_joint <- if (k >= 1 && (N_X - k - 1) > 0 && R2_total > 0) {
    ((N_X - k - 1) / k) * (R2_total / (1 - R2_total))
  } else NA_real_

  D <- sqrt(diag(Omega1))
  Corr <- Omega1 / (D %o% D)
  offdiag <- Corr[upper.tri(Corr)]
  mean_abs_corr <- if (length(offdiag)) mean(abs(offdiag)) else NA_real_

  out <- list(
    gegger_intercept = NA_real_,
    gegger_slope = NA_real_,
    gegger_se_intercept = NA_real_,
    gegger_se_slope = NA_real_,
    gegger_p_intercept = NA_real_,
    gegger_p_slope = NA_real_,
    gegger_phi = NA_real_,
    n_modes_initial = length(keep0),
    n_modes_after_F = n_after_F,
    n_modes_used = k,
    cum_var_expl_initial = 100 * sum(var_frac[keep0]),
    cum_var_expl_used = 100 * sum(var_frac[keepF]),
    kappa_used = kappa,
    mean_abs_corr = mean_abs_corr,
    mean_F_mode = mean(F_mode1, na.rm = TRUE),
    median_F_mode = median(F_mode1, na.rm = TRUE),
    min_F_mode = min(F_mode1, na.rm = TRUE),
    mean_R2_mode = mean(R2_mode1, na.rm = TRUE),
    median_R2_mode = median(R2_mode1, na.rm = TRUE),
    min_R2_mode = min(R2_mode1, na.rm = TRUE),
    R2_total_RSS = R2_total,
    F_RSS_joint = F_joint,
    N_X_used = N_X
  )

  if (k < 3) return(out)

  fit <- gls_egger(bx1, by1, Omega = Omega1, lambda = Lambda[keepF])

  out$gegger_intercept <- fit$intercept
  out$gegger_slope <- fit$slope
  out$gegger_se_intercept <- fit$se_intercept
  out$gegger_se_slope <- fit$se_slope
  out$gegger_p_intercept <- fit$p_intercept
  out$gegger_p_slope <- fit$p_slope

  out
}

############################################################
## MAIN LOOP
############################################################
results <- list()

for (i in seq_len(nrow(genes_task))) {
  gene <- genes_task$gene[i]
  chr  <- genes_task$chr[i]

  cat("\n>>>", gene, "(chr", chr, ")\n")

  paths <- file.path(
    beta_dir,
    paste0(gene, c(".eqtl.beta", ".eqtl.se", ".gwas.beta", ".gwas.se"))
  )

  if (!all(file.exists(paths))) {
    cat("  ✗ missing beta/se files\n")
    next
  }

  bx  <- as.numeric(readLines(paths[1]))
  sex <- as.numeric(readLines(paths[2]))
  by  <- as.numeric(readLines(paths[3]))
  sey <- as.numeric(readLines(paths[4]))

  ld <- load_ld_matrix(gene, chr)
  if (is.null(ld) || length(bx) != nrow(ld) || nrow(ld) != ncol(ld)) {
    cat("  ✗ LD missing or dimension mismatch\n")
    next
  }

  if (length(by) != length(bx) || length(sex) != length(bx) || length(sey) != length(bx)) {
    cat("  ✗ beta/se length mismatch\n")
    next
  }

  sx_med <- median(sex[is.finite(sex) & sex > 0], na.rm = TRUE)
  sy_med <- median(sey[is.finite(sey) & sey > 0], na.rm = TRUE)
  if (!is.finite(sx_med)) sx_med <- 1
  if (!is.finite(sy_med)) sy_med <- 1
  sex[!is.finite(sex) | sex <= 0] <- sx_med
  sey[!is.finite(sey) | sey <= 0] <- sy_med

  N_X <- if ("N_X" %in% names(genes_task) &&
             is.finite(genes_task$N_X[i]) &&
             genes_task$N_X[i] > 10) {
    as.integer(genes_task$N_X[i])
  } else {
    N_X_DEFAULT
  }

  zeta <- compute_zeta(ld)

  fit_burgess <- tryCatch(
    burgess_pca_ivw(
      bx = bx, by = by, byse = sey, ld = ld, var_expl = 0.99
    ),
    error = function(e) NULL
  )

  pc_ivw <- tryCatch(
    pca_prepare(
      bx = bx, by = by, sex = sex, sey = sey, ld = ld, N_X = N_X,
      keep_var = KEEP_VAR_IVW, F_min = F_MIN,
      min_modes = MIN_MODES_IVW, kappa_max = KAPPA_MAX
    ),
    error = function(e) NULL
  )

  fit_ivw <- if (!is.null(pc_ivw)) {
    tryCatch(pca_gls(pc_ivw), error = function(e) NULL)
  } else NULL

  pc_ml <- tryCatch(
    pca_prepare(
      bx = bx, by = by, sex = sex, sey = sey, ld = ld, N_X = N_X,
      keep_var = KEEP_VAR_ML, F_min = F_MIN,
      min_modes = MIN_MODES_ML, kappa_max = KAPPA_MAX
    ),
    error = function(e) NULL
  )

  fit_ml <- if (!is.null(pc_ml)) {
    tmp <- tryCatch(pca_maxlik(pc_ml, model = ML_MODEL, alpha = 0.05), error = function(e) NULL)
    if (!is.null(tmp) && tmp$conv != 0) NULL else tmp
  } else NULL

  fit_gegger <- tryCatch(
    pca_mr_gegger_gene(
      bx = bx, by = by, sex = sex, sey = sey, ld = ld, N_X = N_X,
      keep_var = KEEP_VAR_EGGER, F_min = F_MIN,
      kappa_max = KAPPA_MAX, min_modes = MIN_MODES_EGGER
    ),
    error = function(e) NULL
  )

  cat(sprintf(
    "  ✓ Burgess=%s | PCA-IVW=%s | PCA-ML=%s | PCA-Egger=%s\n",
    ifelse(!is.null(fit_burgess) && is.finite(fit_burgess$slope), sprintf("%.4f", fit_burgess$slope), "NA"),
    ifelse(!is.null(fit_ivw)     && is.finite(fit_ivw$slope),     sprintf("%.4f", fit_ivw$slope), "NA"),
    ifelse(!is.null(fit_ml)      && is.finite(fit_ml$slope),      sprintf("%.4f", fit_ml$slope), "NA"),
    ifelse(!is.null(fit_gegger)  && is.finite(fit_gegger$gegger_slope), sprintf("%.4f", fit_gegger$gegger_slope), "NA")
  ))

  results[[length(results) + 1]] <- data.frame(
    gene = gene,
    chr = chr,
    N_X = N_X,
    p = length(bx),
    zeta = zeta,

    burgess_beta = if (!is.null(fit_burgess)) fit_burgess$slope else NA_real_,
    burgess_se   = if (!is.null(fit_burgess)) fit_burgess$se else NA_real_,
    burgess_p    = if (!is.null(fit_burgess)) fit_burgess$p else NA_real_,
    burgess_K    = if (!is.null(fit_burgess)) fit_burgess$K else NA_integer_,

    pca_ivw_beta = if (!is.null(fit_ivw)) fit_ivw$slope else NA_real_,
    pca_ivw_se   = if (!is.null(fit_ivw)) fit_ivw$se else NA_real_,
    pca_ivw_p    = if (!is.null(fit_ivw)) fit_ivw$p else NA_real_,
    pca_ivw_Q    = if (!is.null(fit_ivw)) fit_ivw$Q else NA_real_,
    pca_ivw_df   = if (!is.null(fit_ivw)) fit_ivw$df else NA_real_,
    pca_ivw_Q_p  = if (!is.null(fit_ivw)) fit_ivw$Q_pvalue else NA_real_,
    pca_ivw_phi  = if (!is.null(fit_ivw)) fit_ivw$phi else NA_real_,
    pca_ivw_rse  = if (!is.null(fit_ivw)) fit_ivw$rse else NA_real_,
    pca_ivw_npcs = if (!is.null(pc_ivw)) pc_ivw$n_modes else NA_integer_,
    pca_ivw_mean_F = if (!is.null(pc_ivw)) pc_ivw$mean_F else NA_real_,
    pca_ivw_mean_F_marginal = if (!is.null(pc_ivw)) pc_ivw$mean_F_marginal else NA_real_,
    pca_ivw_kappa = if (!is.null(pc_ivw)) pc_ivw$kappa else NA_real_,

    pca_ml_beta = if (!is.null(fit_ml)) fit_ml$slope else NA_real_,
    pca_ml_se   = if (!is.null(fit_ml)) fit_ml$se else NA_real_,
    pca_ml_p    = if (!is.null(fit_ml)) fit_ml$p else NA_real_,
    pca_ml_ci_lo = if (!is.null(fit_ml)) fit_ml$ci_lo else NA_real_,
    pca_ml_ci_hi = if (!is.null(fit_ml)) fit_ml$ci_hi else NA_real_,
    pca_ml_npcs = if (!is.null(fit_ml)) fit_ml$K else NA_integer_,
    pca_ml_conv = if (!is.null(fit_ml)) fit_ml$conv else NA_integer_,
    pca_ml_rse  = if (!is.null(fit_ml)) fit_ml$rse else NA_real_,
    pca_ml_heter_stat = if (!is.null(fit_ml)) fit_ml$heter_stat else NA_real_,
    pca_ml_heter_df   = if (!is.null(fit_ml)) fit_ml$heter_df else NA_real_,
    pca_ml_heter_p    = if (!is.null(fit_ml)) fit_ml$heter_p else NA_real_,
    pca_ml_mean_F = if (!is.null(pc_ml)) pc_ml$mean_F else NA_real_,
    pca_ml_mean_F_marginal = if (!is.null(pc_ml)) pc_ml$mean_F_marginal else NA_real_,
    pca_ml_kappa = if (!is.null(pc_ml)) pc_ml$kappa else NA_real_,

    pca_egger_beta = if (!is.null(fit_gegger)) fit_gegger$gegger_slope else NA_real_,
    pca_egger_se   = if (!is.null(fit_gegger)) fit_gegger$gegger_se_slope else NA_real_,
    pca_egger_p    = if (!is.null(fit_gegger)) fit_gegger$gegger_p_slope else NA_real_,
    pca_egger_intercept = if (!is.null(fit_gegger)) fit_gegger$gegger_intercept else NA_real_,
    pca_egger_intercept_se = if (!is.null(fit_gegger)) fit_gegger$gegger_se_intercept else NA_real_,
    pca_egger_intercept_p  = if (!is.null(fit_gegger)) fit_gegger$gegger_p_intercept else NA_real_,
    pca_egger_npcs = if (!is.null(fit_gegger)) fit_gegger$n_modes_used else NA_integer_,
    pca_egger_mean_F = if (!is.null(fit_gegger)) fit_gegger$mean_F_mode else NA_real_,
    pca_egger_kappa = if (!is.null(fit_gegger)) fit_gegger$kappa_used else NA_real_,

    stringsAsFactors = FALSE
  )
}

final <- rbindlist(results, fill = TRUE)
fwrite(final, out_file, sep = "\t")

cat("\n=== DONE (array", task_id, ") ===\n")
cat("Genes analyzed:", nrow(final), "\n")
cat("Output:", out_file, "\n")
EOF
