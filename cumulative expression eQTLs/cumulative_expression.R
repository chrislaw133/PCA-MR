#Fetch protein coding genes
library(biomaRt)
library(data.table)
library(tidyverse)
library(fdapace)
library(pracma)
library(sceasy)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)

# Connect to Ensembl GRCh37 (hg19) archive
ensembl37 <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  host = "grch37.ensembl.org"
)

# Get all protein-coding genes
genes_hg19 <- getBM(
  attributes = c("ensembl_gene_id",    # Ensembl gene ID
                 "hgnc_symbol",        # HGNC gene symbol
                 "chromosome_name",    # Chromosome
                 "start_position",     # Genomic start
                 "end_position",       # Genomic end
                 "strand"),            # Strand (+1 / -1)
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl37
)

genes_hg19 <- genes_hg19[genes_hg19$hgnc_symbol != "", ]

#Read in RNA expression
use_condaenv("atacseq-env",
             conda="/path/to/anaconda3/bin/conda",
             required=TRUE)

rds <- sceasy::convertFormat("/path/to/FPP_SERT.h5ad", from="anndata", to="seurat")

#Quality control
mat <- rds@assays$RNA@data
rds$percent.mt <- Matrix::colSums(mat[mito.genes, , drop = FALSE]) / Matrix::colSums(mat) * 100
rds$nFeature_RNA <- Matrix::colSums(rds@assays$RNA@data > 0)
rds <- subset(rds, subset = nFeature_RNA > 200 & percent.mt < 10)

#Intersect genes in RNA expression dataset and genes_hg19
genelist <- rownames(rds)
genelist <- genelist[genelist %in% genes_hg19$hgnc_symbol]

rds[["RNA"]]@meta.features <- data.frame(row.names = rownames(rds))

VariableFeatures(rds) <- genelist

#Scaling
rds <- ScaleData(rds)

#PCA
rds <- RunPCA(rds)

#Batch correction
library(harmony)
rds <- RunHarmony(rds,"batch")

#Dimensionality reduction
rds <- RunUMAP(rds, reduction="harmony", dims=1:10)

#Clustering
rds <- FindNeighbors(rds, dims = 1:10)
rds <- FindClusters(rds, resolution = 0.5)
SetIdent(rds,value="seurat_clusters")

#Read in slingshot trajectories derived from HVGs
dfti <- fread("/path/to/dfti.txt")
#setnames(dfti, "rn", "cellid")

rds@meta.data$cellid <- rownames(rds@meta.data)

#first extract the processed RNA abundance levels
df_rna <- rds@assays$RNA$scale.data %>% t() %>% as_tibble()
df_rna$cellid <- colnames(rds)

# combine with cell information
df_rna <- df_rna %>%
  left_join(rds@meta.data, by = "cellid") %>%
  dplyr::select(donor_id, cellid, any_of(genelist)) %>%
  left_join(dfti %>% dplyr::select(cellid, slingPseudotime_1), by = "cellid")

# transfer pseudotime points to time period (by rounding the pseudotime)
df_rna <- df_rna %>% mutate(pseudotime=round(slingPseudotime_1))
setnames(df_rna, "donor_id", "id")

cum_expression <- function(exposureDat, genelist, id_var = "id", pseudotime_var = "pseudotime") {
  library(tidyverse)
  library(fdapace)
  library(pracma)
  
  all_donors <- unique(exposureDat[[id_var]])
  
  # collapse replicates per donor/pseudotime
  exposureDat <- exposureDat %>%
    group_by(.data[[id_var]], .data[[pseudotime_var]]) %>%
    summarise(across(all_of(genelist), mean), .groups = "drop") %>%
    drop_na()
  
  fullPACE <- function(geneName) {
  s <- exposureDat %>%
    dplyr::select(all_of(c(id_var, pseudotime_var, geneName))) %>%
    dplyr::rename(expression = !!sym(geneName)) %>%
    group_by(.data[[id_var]]) %>%
    arrange(.data[[pseudotime_var]]) %>%
    group_split()

  # keep only donors with >1 timepoint
  s <- s[lengths(s) > 1]
  if (length(s) == 0) {
    return(setNames(rep(NA_real_, length(all_donors)), all_donors))
  }

  ly <- map(s, ~ .x$expression)
  lt <- map(s, ~ .x[[pseudotime_var]])

  res <- tryCatch(
    FPCA(ly, lt, optns = list(dataType = "Sparse")),
    error = function(e) NULL
  )
  if (is.null(res) || is.null(res$xiEst) || is.null(res$phi)) {
    return(setNames(rep(NA_real_, length(all_donors)), all_donors))
  }

  # reconstruct donor curves and integrate
  Ti_est <- res$workGrid
  x_it <- t(matrix(replicate(length(s), res$mu), ncol = length(s))) +
    res$xiEst %*% t(res$phi)

  xcum <- map_dbl(1:length(s), function(i) {
    trapz(x = Ti_est, y = x_it[i, ])
  })

  donor_ids <- map_chr(s, ~ as.character(unique(.x[[id_var]])[1]))
  out <- setNames(rep(NA_real_, length(all_donors)), as.character(all_donors))
  out[donor_ids] <- xcum

  return(out[as.character(all_donors)])
}

  # build donor × gene matrix
  cum_list <- lapply(genelist, fullPACE)
  cum_mat <- do.call(cbind, cum_list)
  rownames(cum_mat) <- all_donors
  colnames(cum_mat) <- genelist
  
  return(cum_mat)
}

cum_mat <- cum_expression(df_rna, genelist, id_var="id", pseudotime_var="pseudotime")
