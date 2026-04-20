#!/bin/bash
#SBATCH --array=1-100
#SBATCH --mem=4G
#SBATCH --time=00-24:00:00

set -euo pipefail

module load apps/r/4.3.3

Rscript --vanilla - <<'EOF'
suppressPackageStartupMessages({
  library(MendelianRandomization)
  library(data.table)
  library(Matrix)
})

############################################################
## PATHS
############################################################
gene_file <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/cumexpr_stuff/snpcorrelations/genesandinstruments.txt"
beta_dir  <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/GLS_0.9/betavectors"
ld_dir    <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/GLS_0.9/snpcorr"

out_dir   <- "/deac/bio/lackGrp/lawrcm22/serotonin_eqtl/results/cum_generalized"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

MODEL_TYPE <- "random"

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
out_file <- file.path(out_dir, paste0("generalized_methods_", job_id, "_part_", task_id, ".tsv"))

############################################################
## HELPERS
############################################################
load_ld_matrix <- function(gene, chr) {
  f <- file.path(ld_dir, paste0(gene, "_chr", chr, "_pruned.unphased.vcor1.zst"))
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
  LD <- (LD + t(LD)) / 2

  eig <- eigen(LD, symmetric = TRUE, only.values = TRUE)$values
  eig[eig < 0] <- 0

  total <- sum(eig)
  if (total <= 0) return(NA_real_)

  p_i <- eig / total
  p_i <- p_i[p_i > 0]

  H <- -sum(p_i * log(p_i))
  r_eff <- exp(H)
  r_eff / ncol(LD)
}

read_numeric_file <- function(path) {
  as.numeric(readLines(path))
}

############################################################
## MAIN LOOP
############################################################
results <- vector("list", nrow(genes_task))

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
    results[[i]] <- data.frame(
      gene = gene,
      chr = chr,
      N_X = if ("N_X" %in% names(genes_task)) genes_task$N_X[i] else NA_integer_,
      p = NA_integer_,
      model_type = MODEL_TYPE,
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

      error = "missing beta/se files",
      stringsAsFactors = FALSE
    )
    next
  }

  bx   <- read_numeric_file(paths[1])
  bxse <- read_numeric_file(paths[2])
  by   <- read_numeric_file(paths[3])
  byse <- read_numeric_file(paths[4])

  ld <- load_ld_matrix(gene, chr)
  if (is.null(ld)) {
    cat("  ✗ LD missing\n")
    results[[i]] <- data.frame(
      gene = gene,
      chr = chr,
      N_X = if ("N_X" %in% names(genes_task)) genes_task$N_X[i] else NA_integer_,
      p = length(bx),
      model_type = MODEL_TYPE,
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

      error = "LD missing",
      stringsAsFactors = FALSE
    )
    next
  }

  if (length(by) != length(bx) || length(bxse) != length(bx) || length(byse) != length(bx)) {
    cat("  ✗ beta/se length mismatch\n")
    results[[i]] <- data.frame(
      gene = gene,
      chr = chr,
      N_X = if ("N_X" %in% names(genes_task)) genes_task$N_X[i] else NA_integer_,
      p = length(bx),
      model_type = MODEL_TYPE,
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

      error = "beta/se length mismatch",
      stringsAsFactors = FALSE
    )
    next
  }

  if (nrow(ld) != length(bx) || ncol(ld) != length(bx)) {
    cat("  ✗ LD dimension mismatch\n")
    results[[i]] <- data.frame(
      gene = gene,
      chr = chr,
      N_X = if ("N_X" %in% names(genes_task)) genes_task$N_X[i] else NA_integer_,
      p = length(bx),
      model_type = MODEL_TYPE,
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

      error = "LD dimension mismatch",
      stringsAsFactors = FALSE
    )
    next
  }

  bxse_med <- median(bxse[is.finite(bxse) & bxse > 0], na.rm = TRUE)
  byse_med <- median(byse[is.finite(byse) & byse > 0], na.rm = TRUE)
  if (!is.finite(bxse_med)) bxse_med <- 1
  if (!is.finite(byse_med)) byse_med <- 1
  bxse[!is.finite(bxse) | bxse <= 0] <- bxse_med
  byse[!is.finite(byse) | byse <= 0] <- byse_med

  zeta <- compute_zeta(ld)

  input <- tryCatch(
    mr_input(bx = bx, bxse = bxse, by = by, byse = byse, correlation = ld),
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

  cat(sprintf(
    "  ✓ GIVW=%s | GEgger=%s | GML=%s\n",
    ifelse(!is.null(givw) && is.finite(givw@Estimate), sprintf("%.4f", givw@Estimate), "NA"),
    ifelse(!is.null(gegg) && is.finite(gegg@Estimate), sprintf("%.4f", gegg@Estimate), "NA"),
    ifelse(!is.null(gml)  && is.finite(gml@Estimate),  sprintf("%.4f", gml@Estimate), "NA")
  ))

  results[[i]] <- data.frame(
    gene = gene,
    chr = chr,
    N_X = if ("N_X" %in% names(genes_task)) genes_task$N_X[i] else NA_integer_,
    p = length(bx),
    model_type = MODEL_TYPE,
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

    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

final <- rbindlist(results, fill = TRUE)
fwrite(final, out_file, sep = "\t")

cat("\n=== DONE (array", task_id, ") ===\n")
cat("Genes analyzed:", nrow(final), "\n")
cat("Output:", out_file, "\n")
EOF
