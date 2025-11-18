#!/bin/bash

Rscript --vanilla - <<'EOF'
suppressPackageStartupMessages({
  library(simmrd)
  library(MendelianRandomization)
  library(Matrix)
  library(data.table)
})

#Replicate ID
rep_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(rep_id)) rep_id <- 1
cat(sprintf("🧬 Running replicate %d\n", rep_id))

#Output directory
outdir <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/null_simulations"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Invidual-level simulation
individual_params <- list(
  sample_size_Xs = 2e4,
  sample_size_Y = 2e4,
  prop_gwas_overlap_Xs_and_Y = 0,
  number_of_exposures = 1,
  phenotypic_correlation_Xs = 0,
  genetic_correlation_Xs = c(0, 0),
  Xs_variance_explained_by_U = 0.001, #Sqrt for causal effect
  Y_variance_explained_by_Xs = 0,
  signs_of_causal_effects = 1,
  Y_variance_explained_by_U = 0.0025, #Sqrt for causal effect
  number_of_causal_SNPs = 10,
  mafs_of_causal_SNPs = stats::runif(100, 0.1, 0.5),
  Xs_variance_explained_by_g = 0.1,
  number_of_UHP_causal_SNPs = 2,
  number_of_CHP_causal_SNPs = 2,
  Y_variance_explained_by_UHP = 0.0025, #Sqrt for causal effect (Effect of confounding is split between UHP and CHP for outcome)
  U_variance_explained_by_CHP = 1,
  LD_causal_SNPs = 'ar1(0.4)',
  number_of_LD_blocks = 3,
  MR_standardization = 'Z',
  simtype = 'winners',
  MVMR_IV_selection_type = 'joint',
  IV_Pvalue_threshold = 1e-4,
  LD_pruning_r2 = 1,
  N_of_LD_ref = Inf,
  fix_Fstatistic_at = 10
)

#PCA-MR Function
pca_mr <- function(bx, by, se_y, ld) {
  eig <- eigen(ld, symmetric = TRUE)
  Q <- eig$vectors
  Lambda <- pmax(eig$values, 1e-8)
  bx_tilde <- as.numeric(crossprod(Q, bx))
  by_tilde <- as.numeric(crossprod(Q, by))
  cov_gwas <- (se_y %o% se_y) * ld
  cov_tilde <- t(Q) %*% cov_gwas %*% Q
  se_tilde <- sqrt(diag(cov_tilde))
  var_frac <- Lambda / sum(Lambda)
  keep <- which(cumsum(var_frac) <= 0.99)
  if (length(keep) == 0) keep <- which.max(var_frac)
  bx_tilde <- bx_tilde[keep]; by_tilde <- by_tilde[keep]; se_tilde <- se_tilde[keep]
  Lambda <- Lambda[keep]
  w <- Lambda / (se_tilde^2)
  w_se <- 1 / se_tilde^2
  slope <- sum(w * bx_tilde * by_tilde) / sum(w * bx_tilde^2)
  resid <- by_tilde - slope * bx_tilde
  Qstat <- sum(w * resid^2)
  Q_pval <- pchisq(Qstat, df = (length(bx_tilde) -1), lower.tail = FALSE)
  phi <- max(1, Qstat / (length(bx_tilde) - 1))
  se <- sqrt(1 / sum(w_se * bx_tilde^2)) * sqrt(phi)
  p <- 2 * pnorm(-abs(slope / se))
  list(slope = slope, se = se, p = p, Q_pval = Q_pval)
}

#One replicate
run_one_sim <- function(individual_params, r) {
    gwas_data <- generate_individual(params = individual_params, seed = sample.int(1e8, 1))
      
      bx <- as.numeric(gwas_data$bx)
      by <- as.numeric(gwas_data$by)
      bxse <- as.numeric(gwas_data$bxse)
      byse <- as.numeric(gwas_data$byse)
      LD <- as.matrix(gwas_data$LDMatrix)
      
      input <- mr_input(bx, bxse, by, byse)
      
      # Run each MR
      pcgmm <- tryCatch(mr_pcgmm(mr_input(bx, bxse, by, byse, correlation = LD),
                                 nx = 2e4, ny = 2e4, thres = 0.99, robust = TRUE),
                        error = function(e) NULL)
      pca   <- tryCatch(pca_mr(bx, by, byse, LD), error = function(e) NULL)
      
      data.frame(
        rep = r,
        slope_pca   = if (!is.null(pca)) pca$slope else NA,
        p_pca       = if (!is.null(pca)) pca$p else NA,
        Q_pval      = if (!is.null(pca)) pca$Q_pval else NA,
        slope_pcgmm = if (!is.null(pcgmm)) pcgmm@Estimate else NA,
        p_pcgmm     = if (!is.null(pcgmm)) pcgmm@Pvalue else NA
      )
    }

#Run replicate
set.seed(rep_id)
res <- run_one_sim(individual_params, rep_id)

# === Write output ===
outfile <- file.path(outdir, sprintf("null_rep_%04d.csv", rep_id))
fwrite(res, outfile, sep = "\t")

EOF
