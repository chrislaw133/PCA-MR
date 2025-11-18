#!/bin/bash

Rscript --vanilla - <<'EOF'
suppressPackageStartupMessages({
  library(simmrd)
  library(MendelianRandomization)
  library(Matrix)
  library(data.table)
})

#Get replicate ID
rep_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(rep_id)) rep_id <- 1
cat(sprintf("🧬 Running replicate %d\n", rep_id))

#Output directory
outdir <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/null_simulations"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

Individual-level simulation
individual_params <- list(
  sample_size_Xs = 2e4,
  sample_size_Y = 2e4,
  prop_gwas_overlap_Xs_and_Y = 0,
  number_of_exposures = 1,
  phenotypic_correlation_Xs = 0,
  genetic_correlation_Xs = c(0, 0),
  Xs_variance_explained_by_U = 0.01,
  Y_variance_explained_by_Xs = 0,
  signs_of_causal_effects = 1,
  Y_variance_explained_by_U = 0.0025,
  number_of_causal_SNPs = 10,
  mafs_of_causal_SNPs = stats::runif(100, 0.1, 0.5),
  Xs_variance_explained_by_g = 0.1,
  number_of_UHP_causal_SNPs = 2,
  number_of_CHP_causal_SNPs = 2,
  Y_variance_explained_by_UHP = 0.0025,
  U_variance_explained_by_CHP = 1,
  LD_causal_SNPs = 'ar1(0.4)',
  number_of_LD_blocks = 3,
  MR_standardization = 'Z',
  simtype = 'winners',
  MVMR_IV_selection_type = 'joint',
  IV_Pvalue_threshold = 1e-4,
  LD_pruning_r2 = 0.1,
  N_of_LD_ref = Inf,
  fix_Fstatistic_at = 10
)

#One replicate function
run_one_sim <- function(individual_params, r) {
  gwas_data <- generate_individual(individual_params, seed = sample.int(1e8, 1))
  bx <- as.numeric(gwas_data$bx)
  by <- as.numeric(gwas_data$by)
  bxse <- as.numeric(gwas_data$bxse)
  byse <- as.numeric(gwas_data$byse)
  LD <- as.matrix(gwas_data$LDMatrix)

  input <- mr_input(bx, bxse, by, byse)
  ivw   <- tryCatch(mr_ivw(input), error = function(e) NULL)
  egger <- tryCatch(mr_egger(input), error = function(e) NULL)
  ml    <- tryCatch(mr_maxlik(input), error = function(e) NULL)
  wm    <- tryCatch(mr_median(input), error = function(e) NULL)
  lasso <- tryCatch(mr_lasso(input), error = function(e) NULL)
  cml   <- tryCatch(
    mr_cML(
      input,
      MA = TRUE,
      DP = TRUE,
      K_vec = 0:(length(bx) - 2),
      random_start = 0,
      num_pert = 200,
      random_start_pert = 0,
      maxit = 100,
      random_seed = 314,
      n = individual_params$sample_size_Y,
      Alpha = 0.05
    ),
    error = function(e) NULL
  )

  data.frame(
    rep = r,
    slope_ivw   = if (!is.null(ivw)) ivw@Estimate else NA,
    p_ivw       = if (!is.null(ivw)) ivw@Pvalue else NA,
    slope_egger = if (!is.null(egger)) egger@Estimate else NA,
    p_egger     = if (!is.null(egger)) egger@Pvalue.Est else NA,
    slope_ml    = if (!is.null(ml)) ml@Estimate else NA,
    p_ml        = if (!is.null(ml)) ml@Pvalue else NA,
    slope_wm    = if (!is.null(wm)) wm@Estimate else NA,
    p_wm        = if (!is.null(wm)) wm@Pvalue else NA,
    slope_lasso = if (!is.null(lasso)) lasso@Estimate else NA,
    p_lasso     = if (!is.null(lasso)) lasso@Pvalue else NA,
    slope_cml = if (!is.null(cml)) cml@Estimate else NA,
    p_cml     = if (!is.null(cml)) cml@Pvalue else NA
  )
}

#Run one replicate
set.seed(rep_id)
res <- run_one_sim(individual_params, rep_id)

#Write output
outfile <- file.path(outdir, sprintf("null_rep_%04d.csv", rep_id))
fwrite(res, outfile, sep = "\t")

EOF
