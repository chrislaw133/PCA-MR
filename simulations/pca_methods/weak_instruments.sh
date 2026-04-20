#!/bin/bash
#SBATCH --array=1-1000
#SBATCH --mem=2G
#SBATCH --time=00-1:00:00

set -euo pipefail

module load apps/r/4.3.3

Rscript --vanilla - <<'EOF'
suppressPackageStartupMessages({
  library(simmrd)
  library(Matrix)
  library(MASS)
  library(data.table)
})

REPS_PER_SCENARIO <- 1000

outdir <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/null_simulations/weak_pca"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ar1_vals <- c(0.2, 0.4, 0.6, 0.8, 0.9)
causal_vals <- c(0, 0.01, 0.03, 0.05, 0.075, 0.1)
pleio_vals  <- c(0, 0.01, 0.03, 0.05, 0.075, 0.1)
F_stats <- c(5, 10, 25)

# Main power / null grid
grid_main <- CJ(
  scenario_class = "main",
  ar1 = ar1_vals,
  causal_effect = causal_vals,
  fixed_Fstat = F_stats
)[
  , `:=`(
    sample_size = 250L,
    pleio_type = "none",
    pleio_effect = 0,
    n_uhp = 0L,
    n_chp = 0L,
    y_var_uhp = 0,
    u_var_chp = 0
  )
]

grid_uhp <- CJ(
  pleio_effect = pleio_vals,
  fixed_Fstat = F_stats
)[
  , `:=`(
    scenario_class = "pleio",
    sample_size = 250L,
    ar1 = 0.9,
    causal_effect = 0,
    pleio_type = "UHP",
    n_uhp = 4L,
    n_chp = 0L,
    y_var_uhp = pleio_effect,
    u_var_chp = 0
  )
]

grid_chp <- CJ(
  pleio_effect = pleio_vals,
  fixed_Fstat = F_stats
)[
  , `:=`(
    scenario_class = "pleio",
    sample_size = 250L,
    ar1 = 0.9,
    causal_effect = 0,
    pleio_type = "CHP",
    n_uhp = 0L,
    n_chp = 4L,
    y_var_uhp = 0,
    u_var_chp = pleio_effect
  )
]

scenario_grid <- rbindlist(list(grid_main, grid_uhp, grid_chp), fill = TRUE)
scenario_grid[, scenario_id := .I]
stopifnot(nrow(scenario_grid) == 126)

array_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (array_id < 1L || array_id > REPS_PER_SCENARIO) {
  stop("array_id out of bounds")
}

rep_id <- array_id
scenario_ids <- seq_len(nrow(scenario_grid))


make_params <- function(sc, rep_seed) {
  list(
    sample_size_Xs = sc$sample_size,
    sample_size_Y  = 50000,
    prop_gwas_overlap_Xs_and_Y = 0,
    number_of_exposures = 1,
    phenotypic_correlation_Xs = 0,
    genetic_correlation_Xs = 0,
    Xs_variance_explained_by_U = sc$u_var_chp,

    Y_variance_explained_by_Xs = sc$causal_effect^2,
    signs_of_causal_effects = 1,

    Y_variance_explained_by_U = sc$u_var_chp,
    number_of_causal_SNPs = 10,
    mafs_of_causal_SNPs = stats::runif(10, 0.1, 0.5),

    Xs_variance_explained_by_g = 0.3,

    number_of_UHP_causal_SNPs = sc$n_uhp,
    number_of_CHP_causal_SNPs = sc$n_chp,
    Y_variance_explained_by_UHP = sc$y_var_uhp,
    U_variance_explained_by_CHP = sc$u_var_chp,

    LD_causal_SNPs = sprintf("ar1(%s)", sc$ar1),
    number_of_LD_blocks = 2,
    MR_standardization = "none",
    simtype = "weak",
    MVMR_IV_selection_type = "joint",
    IV_Pvalue_threshold = 1e-4,
    LD_pruning_r2 = 1,
    N_of_LD_ref = 500,
    fix_Fstatistic_at = sc$fixed_Fstat
  )
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

effective_rank <- function(vals) {
  vals <- as.numeric(vals)
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) return(NA_real_)
  p_i <- vals / sum(vals)
  H <- -sum(p_i * log(p_i))
  exp(H)
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

gls_egger <- function(bx_tilde, by_tilde, Omega) {
  n <- length(bx_tilde)
  p <- 2L

  if (length(by_tilde) != n) stop("gls_egger: length mismatch")
  if (!is.matrix(Omega) || any(dim(Omega) != c(n, n))) {
    stop("gls_egger: Omega has wrong dimension")
  }

  # actual Egger requires at least 3 instruments / modes
  if (n < 3) {
    return(list(
      intercept = NA_real_, slope = NA_real_,
      se_intercept = NA_real_, se_slope = NA_real_,
      p_intercept = NA_real_, p_slope = NA_real_,
      Q = NA_real_, df = NA_real_, rse = NA_real_, Q_pval = NA_real_
    ))
  }

  # source-style sign harmonization
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
    K = K,
    erank = effective_rank(eigvals[seq_len(K)])
  )
}

KEEP_VAR_IVW <- 0.999
KEEP_VAR_EGGER <- 0.999
KEEP_VAR_ML <- 0.999
F_MIN <- 0
MIN_MODES_IVW <- 1
MIN_MODES_EGGER <- 3
MIN_MODES_ML <- 1
KAPPA_MAX <- 100
ML_MODEL <- "random" 
IVW_MODEL <- "random"

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
  
  SigmaX = SigX_pc[keep, keep, drop = FALSE]
  SigmaY = SigY_pc[keep, keep, drop = FALSE]
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
  
  SigmaX = SigX_pc[keep, keep, drop = FALSE]
  SigmaY = SigY_pc[keep, keep, drop = FALSE]
  
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
    lambda = lambda[keep],
    erank = effective_rank(lambda[keep])
  )
}

pca_gls <- function(pc, model = IVW_MODEL, phi_floor = 1) {
  if (is.null(pc)) return(NULL)
  if (!model %in% c("fixed", "random")) stop("model must be 'fixed' or 'random'")

  Omega_inv <- tryCatch(solve(pc$Omega), error = function(e) NULL)
  if (is.null(Omega_inv)) return(NULL)

  Lambda <- diag(pc$lambda, nrow = length(pc$lambda), ncol = length(pc$lambda))
  W <- Lambda %*% Omega_inv
  #W <- Omega_inv
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
    N_X_used = N_X,
    erank = effective_rank(Lambda[keepF])
  )

  if (k < 3) return(out)

  fit <- gls_egger(bx1, by1, Omega = Omega1)

  out$gegger_intercept <- fit$intercept
  out$gegger_slope <- fit$slope
  out$gegger_se_intercept <- fit$se_intercept
  out$gegger_se_slope <- fit$se_slope
  out$gegger_p_intercept <- fit$p_intercept
  out$gegger_p_slope <- fit$p_slope

  out
}

res_list <- vector("list", length(scenario_ids))

for (ii in seq_along(scenario_ids)) {
  scenario_id <- scenario_ids[ii]
  sc <- scenario_grid[scenario_id]

  seed_i <- 1e6 * scenario_id + rep_id
  set.seed(seed_i)

  individual_params <- make_params(sc, seed_i)

  gen_err <- NULL
gwas_data <- tryCatch(
  generate_individual(params = individual_params, seed = seed_i),
  error = function(e) {
    gen_err <<- conditionMessage(e)
    NULL
  }
)

if (is.null(gwas_data)) {
  res_list[[ii]] <- data.frame(
    array_id = array_id,
    rep = rep_id,
    scenario_id = scenario_id,

    scenario_class = sc$scenario_class,
    sample_size = sc$sample_size,
    ar1 = sc$ar1,
    fixed_Fstat = sc$fixed_Fstat,
    causal_effect = sc$causal_effect,
    causal_var_expl = sc$causal_effect^2,
    pleio_type = sc$pleio_type,
    pleio_effect = sc$pleio_effect,
    uhp_effect = sc$y_var_uhp,
    chp_effect = sc$u_var_chp,

    seed = seed_i,
    zeta = NA_real_,
    mean_F_marginal = NA_real_,

    burgess_beta = NA_real_,
    burgess_se   = NA_real_,
    burgess_p    = NA_real_,
    burgess_K    = NA_integer_,
    burgess_erank = NA_real_,

    pca_ivw_beta = NA_real_,
    pca_ivw_se   = NA_real_,
    pca_ivw_p    = NA_real_,
    pca_ivw_Q_p  = NA_real_,
    pca_ivw_npcs = NA_integer_,
    pca_ivw_mean_F = NA_real_,
    pca_ivw_kappa = NA_real_,
    pca_ivw_erank = NA_real_,

    pca_ml_beta = NA_real_,
    pca_ml_se   = NA_real_,
    pca_ml_p    = NA_real_,
    pca_ml_ci_lo = NA_real_,
    pca_ml_ci_hi = NA_real_,
    pca_ml_het_p = NA_real_,
    pca_ml_npcs = NA_integer_,
    pca_ml_conv = NA_integer_,
    pca_ml_mean_F = NA_real_,
    pca_ml_kappa  = NA_real_,
    pca_ml_erank = NA_real_,

    pca_egger_beta = NA_real_,
    pca_egger_se   = NA_real_,
    pca_egger_p    = NA_real_,
    pca_egger_intercept = NA_real_,
    pca_egger_intercept_se = NA_real_,
    pca_egger_intercept_p = NA_real_,
    pca_egger_npcs = NA_integer_,
    pca_egger_mean_F = NA_real_,
    pca_egger_kappa = NA_real_,
    pca_egger_erank = NA_real_,

    true_K = NA_integer_,
    true_K_variants = NA_character_,
    sim_error = gen_err,
    stringsAsFactors = FALSE
  )
  next
}

  bx  <- as.numeric(gwas_data$bx)
  by  <- as.numeric(gwas_data$by)
  sex <- as.numeric(gwas_data$bxse)
  sey <- as.numeric(gwas_data$byse)
  ld  <- as.matrix(gwas_data$LDhatMatrix)

  zeta <- compute_zeta(ld)
  N_X <- individual_params$sample_size_Xs

  fit_burgess <- tryCatch(
    burgess_pca_ivw(bx = bx, by = by, byse = sey, ld = ld, var_expl = 0.999),
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

  fit_ivw <- if (!is.null(pc_ivw)) tryCatch(pca_gls(pc_ivw), error = function(e) NULL) else NULL

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

  p <- length(bx)
  snp_id <- as.character(seq_len(p))

  fmt_vars <- function(idx) {
    idx <- sort(unique(as.integer(idx)))
    idx <- idx[!is.na(idx) & idx >= 1 & idx <= p]
    if (length(idx) == 0) return("none")
    paste(snp_id[idx], collapse = ",")
  }

  TRUE_invalid_idx <- which(gwas_data$IVtype != "valid")
  true_K <- length(TRUE_invalid_idx)
  true_K_variants <- fmt_vars(TRUE_invalid_idx)

  res_list[[ii]] <- data.frame(
    array_id = array_id,
    rep = rep_id,
    scenario_id = scenario_id,

    scenario_class = sc$scenario_class,
    sample_size = sc$sample_size,
    ar1 = sc$ar1,
    fixed_Fstat = sc$fixed_Fstat,
    causal_effect = sc$causal_effect,
    causal_var_expl = sc$causal_effect^2,
    pleio_type = sc$pleio_type,
    pleio_effect = sc$pleio_effect,
    uhp_effect = sc$y_var_uhp,
    chp_effect = sc$u_var_chp,

    seed = seed_i,
    zeta = zeta,
    
    mean_F_marginal = if (!is.null(pc_ivw)) pc_ivw$mean_F_marginal else if (!is.null(pc_ml)) pc_ml$mean_F_marginal else NA_real_,
    burgess_beta = if (!is.null(fit_burgess)) fit_burgess$slope else NA_real_,
    burgess_se   = if (!is.null(fit_burgess)) fit_burgess$se else NA_real_,
    burgess_p    = if (!is.null(fit_burgess)) fit_burgess$p else NA_real_,
    burgess_K    = if (!is.null(fit_burgess)) fit_burgess$K else NA_integer_,
    burgess_erank = if (!is.null(fit_burgess)) fit_burgess$erank else NA_real_,

    pca_ivw_beta = if (!is.null(fit_ivw)) fit_ivw$slope else NA_real_,
    pca_ivw_se   = if (!is.null(fit_ivw)) fit_ivw$se else NA_real_,
    pca_ivw_p    = if (!is.null(fit_ivw)) fit_ivw$p else NA_real_,
    pca_ivw_Q_p  = if (!is.null(fit_ivw)) fit_ivw$Q_pvalue else NA_real_,
    pca_ivw_npcs = if (!is.null(pc_ivw)) pc_ivw$n_modes else NA_integer_,
    pca_ivw_mean_F = if (!is.null(pc_ivw)) pc_ivw$mean_F else NA_real_,
    pca_ivw_kappa = if (!is.null(pc_ivw)) pc_ivw$kappa else NA_real_,
    pca_ivw_erank = if (!is.null(pc_ivw)) pc_ivw$erank else NA_real_,

    pca_ml_beta = if (!is.null(fit_ml)) fit_ml$slope else NA_real_,
    pca_ml_se   = if (!is.null(fit_ml)) fit_ml$se else NA_real_,
    pca_ml_p    = if (!is.null(fit_ml)) fit_ml$p else NA_real_,
    pca_ml_ci_lo = if (!is.null(fit_ml)) fit_ml$ci_lo else NA_real_,
    pca_ml_ci_hi = if (!is.null(fit_ml)) fit_ml$ci_hi else NA_real_,
    pca_ml_het_p = if (!is.null(fit_ml)) fit_ml$heter_p else NA_real_,
    pca_ml_npcs = if (!is.null(fit_ml)) fit_ml$K else NA_integer_,
    pca_ml_conv = if (!is.null(fit_ml)) fit_ml$conv else NA_integer_,
    pca_ml_mean_F = if (!is.null(pc_ml)) pc_ml$mean_F else NA_real_,
    pca_ml_kappa  = if (!is.null(pc_ml)) pc_ml$kappa else NA_real_,
    pca_ml_erank = if (!is.null(pc_ml)) pc_ml$erank else NA_real_,

    pca_egger_beta = if (!is.null(fit_gegger)) fit_gegger$gegger_slope else NA_real_,
    pca_egger_se   = if (!is.null(fit_gegger)) fit_gegger$gegger_se_slope else NA_real_,
    pca_egger_p    = if (!is.null(fit_gegger)) fit_gegger$gegger_p_slope else NA_real_,
    pca_egger_intercept = if (!is.null(fit_gegger)) fit_gegger$gegger_intercept else NA_real_,
    pca_egger_intercept_se = if (!is.null(fit_gegger)) fit_gegger$gegger_se_intercept else NA_real_,
    pca_egger_intercept_p = if (!is.null(fit_gegger)) fit_gegger$gegger_p_intercept else NA_real_,
    pca_egger_npcs = if (!is.null(fit_gegger)) fit_gegger$n_modes_used else NA_integer_,
    pca_egger_mean_F = if (!is.null(fit_gegger)) fit_gegger$mean_F_mode else NA_real_,
    pca_egger_kappa = if (!is.null(fit_gegger)) fit_gegger$kappa_used else NA_real_,
    pca_egger_erank = if (!is.null(fit_gegger)) fit_gegger$erank else NA_real_,
    
    sim_error = NA_character_,
    true_K = true_K,
    true_K_variants = true_K_variants,
    stringsAsFactors = FALSE
  )
}

outfile <- file.path(
  outdir,
  sprintf("rep_%d.tsv", array_id)
)

res <- data.table::rbindlist(res_list, fill = TRUE)
data.table::fwrite(res, outfile, sep = "\t")

EOF