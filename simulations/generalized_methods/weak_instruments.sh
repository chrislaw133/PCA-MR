#!/bin/bash
#SBATCH --array=1-1000
#SBATCH --mem=2G
#SBATCH --time=00-1:00:00


set -euo pipefail

module load apps/r/4.3.3

Rscript --vanilla - <<'EOF'
suppressPackageStartupMessages({
  library(simmrd)
  library(MendelianRandomization)
  library(Matrix)
  library(data.table)
  library(MASS)
})

REPS_PER_SCENARIO <- 1000

outdir <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/null_simulations/weak_generalized"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################################################
## Change this in one place for all three methods
############################################################
MODEL_TYPE <- "random"  

ar1_vals      <- c(0.2, 0.4, 0.6, 0.8, 0.9)
causal_vals   <- c(0, 0.01, 0.03, 0.05, 0.075, 0.1)
pleio_vals    <- c(0, 0.01, 0.03, 0.05, 0.075, 0.1)
F_stats <- c(5, 10, 25)
pruning_vals  <- c(0.1, 0.4, 0.6, 0.8)

############################################################
## Scenario grid
############################################################

grid_main <- CJ(
  scenario_class = "main",
  ar1 = ar1_vals,
  causal_effect = causal_vals,
  fixed_Fstat = F_stats,
  pruning_r2 = pruning_vals
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
  fixed_Fstat = F_stats,
  pleio_effect = pleio_vals,
  pruning_r2 = pruning_vals
)[
  , `:=`(
    scenario_class = "pleio",
    sample_size = 250L,
    ar1 = 0.6,
    causal_effect = 0,
    pleio_type = "UHP",
    n_uhp = 4L,
    n_chp = 0L,
    y_var_uhp = pleio_effect,
    u_var_chp = 0
  )
]

grid_chp <- CJ(
  fixed_Fstat = F_stats,
  pleio_effect = pleio_vals,
  pruning_r2 = pruning_vals
)[
  , `:=`(
    scenario_class = "pleio",
    sample_size = 250L,
    ar1 = 0.6,
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

# 360 + 72 + 72 = 504
stopifnot(nrow(scenario_grid) == 504)

array_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(array_id) || array_id < 1L || array_id > REPS_PER_SCENARIO) {
  stop("array_id out of bounds")
}

rep_id <- array_id
scenario_ids <- seq_len(nrow(scenario_grid))

############################################################
## Helpers
############################################################

make_params <- function(sc, rep_seed) {
  list(
    sample_size_Xs = sc$sample_size,
    sample_size_Y  = 50000,
    prop_gwas_overlap_Xs_and_Y = 0,
    number_of_exposures = 1,
    phenotypic_correlation_Xs = 0,
    genetic_correlation_Xs = 0,
    Xs_variance_explained_by_U = 0,

    Y_variance_explained_by_Xs = sc$causal_effect^2,
    signs_of_causal_effects = 1,

    Y_variance_explained_by_U = 0,
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
    LD_pruning_r2 = sc$pruning_r2,
    N_of_LD_ref = 500,
    fix_Fstatistic_at = sc$fixed_Fstat
  )
}

compute_zeta <- function(LD) {
  LD <- as.matrix(LD)
  LD <- (LD + t(LD)) / 2

  eig <- eigen(LD, symmetric = TRUE, only.values = TRUE)$values
  eig[eig < 0] <- 0

  total <- sum(eig)
  if (total <= 0) {
    return(list(r_eff = NA_real_, zeta = NA_real_))
  }

  p_i <- eig / total
  p_i <- p_i[p_i > 0]

  H <- -sum(p_i * log(p_i))
  r_eff <- exp(H)
  zeta <- r_eff / ncol(LD)

  list(r_eff = r_eff, zeta = zeta)
}

pick_first <- function(x, n) {
  if (!is.null(x) && length(x) == n) x else NULL
}

############################################################
## Main simulation loop
############################################################

res_list <- vector("list", length(scenario_ids))

for (ii in seq_along(scenario_ids)) {
  scenario_id <- scenario_ids[ii]
  sc <- scenario_grid[scenario_id]

  seed_i <- 1000000L * scenario_id + rep_id
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
      pruning_r2 = sc$pruning_r2,
      fixed_Fstat = sc$fixed_Fstat,
      causal_effect = sc$causal_effect,
      causal_var_expl = sc$causal_effect^2,
      pleio_type = sc$pleio_type,
      pleio_effect = sc$pleio_effect,
      uhp_effect = sc$y_var_uhp,
      chp_effect = sc$u_var_chp,

      seed = seed_i,
      r_eff = NA_real_,
      zeta = NA_real_,

      beta_givw = NA_real_,
      se_givw   = NA_real_,
      p_givw    = NA_real_,

      beta_gegger = NA_real_,
      se_gegger   = NA_real_,
      p_gegger    = NA_real_,
      intercept_gegger = NA_real_,
      se_intercept_gegger = NA_real_,
      p_intercept_gegger = NA_real_,

      beta_gml = NA_real_,
      se_gml   = NA_real_,
      p_gml    = NA_real_,

      true_K = NA_integer_,
      true_K_variants = NA_character_,
      sim_error = gen_err,
      stringsAsFactors = FALSE
    )
    next
  }

  bx   <- as.numeric(gwas_data$bx)
  by   <- as.numeric(gwas_data$by)
  bxse <- as.numeric(gwas_data$bxse)
  byse <- as.numeric(gwas_data$byse)

  LD <- NULL
  if (!is.null(gwas_data$LDMatrix)) {
    LD <- as.matrix(gwas_data$LDMatrix)
  } else if (!is.null(gwas_data$LDhatMatrix)) {
    LD <- as.matrix(gwas_data$LDhatMatrix)
  }

  if (is.null(LD)) {
    res_list[[ii]] <- data.frame(
      array_id = array_id,
      rep = rep_id,
      scenario_id = scenario_id,

      scenario_class = sc$scenario_class,
      sample_size = sc$sample_size,
      ar1 = sc$ar1,
      pruning_r2 = sc$pruning_r2,
      fixed_Fstat = sc$fixed_Fstat,
      causal_effect = sc$causal_effect,
      causal_var_expl = sc$causal_effect^2,
      pleio_type = sc$pleio_type,
      pleio_effect = sc$pleio_effect,
      uhp_effect = sc$y_var_uhp,
      chp_effect = sc$u_var_chp,

      seed = seed_i,
      r_eff = NA_real_,
      zeta = NA_real_,

      beta_givw = NA_real_,
      se_givw   = NA_real_,
      p_givw    = NA_real_,

      beta_gegger = NA_real_,
      se_gegger   = NA_real_,
      p_gegger    = NA_real_,
      intercept_gegger = NA_real_,
      se_intercept_gegger = NA_real_,
      p_intercept_gegger = NA_real_,

      beta_gml = NA_real_,
      se_gml   = NA_real_,
      p_gml    = NA_real_,

      true_K = NA_integer_,
      true_K_variants = NA_character_,
      sim_error = "LD matrix missing from generate_individual output",
      stringsAsFactors = FALSE
    )
    next
  }

  p <- length(bx)
  eff_rank <- compute_zeta(LD)
  r_eff <- eff_rank$r_eff
  zeta  <- eff_rank$zeta

  snp_id <- pick_first(gwas_data$snp_id, p)
  if (is.null(snp_id)) snp_id <- pick_first(gwas_data$rsid, p)
  if (is.null(snp_id)) snp_id <- pick_first(gwas_data$SNP, p)
  if (is.null(snp_id)) snp_id <- as.character(seq_len(p))

  fmt_vars <- function(idx) {
    idx <- sort(unique(as.integer(idx)))
    idx <- idx[!is.na(idx) & idx >= 1 & idx <= p]
    if (length(idx) == 0) return("none")
    paste(snp_id[idx], collapse = ",")
  }

  TRUE_invalid_idx <- which(gwas_data$IVtype != "valid")
  true_K <- length(TRUE_invalid_idx)
  true_K_variants <- fmt_vars(TRUE_invalid_idx)

  ############################################################
  ## Generalized estimators
  ############################################################
  input <- tryCatch(
    mr_input(bx = bx, bxse = bxse, by = by, byse = byse, correlation = LD),
    error = function(e) NULL
  )

  givw <- if (!is.null(input)) {
    tryCatch(
      mr_ivw(input, model = MODEL_TYPE, correl = TRUE),
      error = function(e) NULL
    )
  } else NULL

  gegg <- if (!is.null(input)) {
    tryCatch(
      mr_egger(input, model = MODEL_TYPE, correl = TRUE),
      error = function(e) NULL
    )
  } else NULL

  gml <- if (!is.null(input)) {
    tryCatch(
      mr_maxlik(input, model = MODEL_TYPE, correl = TRUE),
      error = function(e) NULL
    )
  } else NULL

  ############################################################
  ## Save row
  ############################################################
  res_list[[ii]] <- data.frame(
    array_id = array_id,
    rep = rep_id,
    scenario_id = scenario_id,

    scenario_class = sc$scenario_class,
    sample_size = sc$sample_size,
    ar1 = sc$ar1,
    pruning_r2 = sc$pruning_r2,
    fixed_Fstat = sc$fixed_Fstat,
    causal_effect = sc$causal_effect,
    causal_var_expl = sc$causal_effect^2,
    pleio_type = sc$pleio_type,
    pleio_effect = sc$pleio_effect,
    uhp_effect = sc$y_var_uhp,
    chp_effect = sc$u_var_chp,

    seed = seed_i,
    model_type = MODEL_TYPE,
    r_eff = r_eff,
    zeta = zeta,

    beta_givw = if (!is.null(givw)) givw@Estimate else NA_real_,
    se_givw   = if (!is.null(givw)) givw@StdError else NA_real_,
    p_givw    = if (!is.null(givw)) givw@Pvalue else NA_real_,

    beta_gegger = if (!is.null(gegg)) gegg@Estimate else NA_real_,
    se_gegger   = if (!is.null(gegg)) gegg@StdError.Est else NA_real_,
    p_gegger    = if (!is.null(gegg)) gegg@Pvalue.Est else NA_real_,
    intercept_gegger = if (!is.null(gegg)) gegg@Intercept else NA_real_,
    se_intercept_gegger = if (!is.null(gegg)) gegg@StdError.Int else NA_real_,
    p_intercept_gegger = if (!is.null(gegg)) gegg@Pvalue.Int else NA_real_,

    beta_gml = if (!is.null(gml)) gml@Estimate else NA_real_,
    se_gml   = if (!is.null(gml)) gml@StdError else NA_real_,
    p_gml    = if (!is.null(gml)) gml@Pvalue else NA_real_,

    true_K = true_K,
    true_K_variants = true_K_variants,

    sim_error = NA_character_,
    stringsAsFactors = FALSE
  )
}

outfile <- file.path(outdir, sprintf("rep_%d.tsv", array_id))
res <- data.table::rbindlist(res_list, fill = TRUE)
data.table::fwrite(res, outfile, sep = "\t")

EOF