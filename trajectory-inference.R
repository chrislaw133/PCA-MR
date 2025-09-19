library(sceasy)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)

#Load in rna seq
use_condaenv("atacseq-env",
             conda="path/to/conda",
             required=TRUE)

rds <- sceasy::convertFormat("path/to/combined.h5ad", from="anndata", to="seurat")

```

rds[["RNA"]]@meta.features <- data.frame(row.names = rownames(rds))

#High variable features
rds <- FindVariableFeatures(rds, selection.method = "vst", nfeatures = 2000)

#Scaling
rds <- ScaleData(rds)

#PCA
rds <- RunPCA(rds)

#Batch correction
library(harmony)
rds <- RunHarmony(rds,"batch")

#Dimension reduction
rds <- RunUMAP(rds, reduction="harmony", dims=1:10)

#Clustering
rds <- FindNeighbors(rds, dims = 1:10)
rds <- FindClusters(rds, resolution = 0.5)

#Subset
SetIdent(rds,value="seurat_clusters")
UMAPPlot(rds)

#Visualization
DimPlot(rds, group.by = 'celltype')
DimPlot(rds, group.by = 'time_point')

#Slingshot - Trajectory inference
library(SingleCellExperiment)
library(slingshot)
library(data.table)
sce <- as.SingleCellExperiment(rds)
sce <- slingshot(sce,reducedDim="UMAP")
dfti <- sce@colData
dfti$slingshot <- NULL
dfti <- as.data.table(dfti,keep.rownames = T) 
setnames(dfti, "rn", "cellid")


rds <- SetIdent(rds, value = "time_point")
levels(Idents(rds)) 

#Pairwise markers
markers_D11_D30 <- FindMarkers(rds, ident.1 = "D11", ident.2 = "D30") %>%
  rownames_to_column("gene") %>%
  mutate

markers_D11_D52 <- FindMarkers(rds, ident.1 = "D11", ident.2 = "D52") %>%
  rownames_to_column("gene") %>%
  mutate(contrast = "D11_vs_D52")

markers_D30_D52 <- FindMarkers(rds, ident.1 = "D30", ident.2 = "D52") %>%
  rownames_to_column("gene") %>%
  mutate(contrast = "D30_vs_D52")

markers_all <- bind_rows(markers_D11_D30,
                         markers_D11_D52,
                         markers_D30_D52)

genelist <- markers_all %>%
  filter(p_val_adj < 0.05) %>%
  distinct(gene) %>%
  pull(gene)

#Smoothed, per gene trajectories
df_rna <- rds@assays$RNA$scale.data %>% 
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("cellid")

rds@meta.data$cellid <- rownames(rds@meta.data)

df_rna <- df_rna %>%
  left_join(rds@meta.data, by = "cellid", suffix = c("", ".meta")) %>%
  dplyr::select(donor_id, cellid, any_of(genelist)) %>%
  left_join(dfti %>% dplyr::select(cellid, slingPseudotime_1), by = "cellid")


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
valid_genes <- intersect(genelist, colnames(df_rna))


#Donor x cumulative expression matrix
cum_expression <- function(exposureDat, genelist, id_var="id", pseudotime_var="pseudotime") {
  library(tidyverse)
  library(doParallel)
  library(foreach)
  library(fdapace)
  library(pracma)
  registerDoParallel()
  loaded_packages <- sub("package:", "", search()[grepl("package:", search())])
  
  all_donors <- unique(exposureDat[[id_var]])
  
  # dynamically use id_var and pseudotime_var
  exposureDat <- exposureDat %>%
    group_by(id, pseudotime) %>%
    summarise(across(all_of(genelist), mean), .groups = "drop") %>%  # <--- note `.groups="drop"`
    drop_na() %>%
    dplyr::rename(
      id = !!sym(id_var),
      pseudotime = !!sym(pseudotime_var)
    )
  
  fullPACE <- function(geneName) {
  s <- exposureDat %>%
    dplyr::select(all_of(c(id_var, pseudotime_var, geneName))) %>%
    dplyr::rename(expression = !!sym(geneName)) %>%
    group_by(.data[[id_var]]) %>%
    arrange(.data[[pseudotime_var]]) %>%
    group_split()
  
  # filter out donors with <=1 point
  lens <- map_dbl(s, nrow)
  keep <- which(lens > 1)
  if (length(keep) == 0) {
    # no valid donors → return NA for all donors
    all_donors <- unique(exposureDat[[id_var]])
    return(setNames(rep(NA_real_, length(all_donors)), all_donors))
  }
  
  s <- s[keep]
  ly <- map(s, ~ .x$expression)
  lt <- map(s, ~ .x[[pseudotime_var]])
  
  res <- FPCA(ly, lt, optns = list(dataType = "Sparse"))
  Ti_est <- res$workGrid
  x_it <- t(matrix(replicate(length(s), res$mu), ncol=length(s))) + res$xiEst %*% t(res$phi)
  
  vals <- map_dbl(seq_along(s), function(i) trapz(x=Ti_est, y=x_it[i,]))
  donor_ids <- map_chr(s, ~ unique(.x[[id_var]])[1])
  out <- setNames(rep(NA_real_, length(all_donors)), all_donors)
  out[donor_ids] <- vals
  return(out)
}

  
  cum_mat <- foreach(geneName = genelist, .combine = 'cbind', .packages = loaded_packages) %dopar% {
    fullPACE(geneName)
  }
  
  colnames(cum_mat) <- genelist
  cum_mat <- as_tibble(cum_mat)
  
  # add donor IDs back
  cum_mat[[id_var]] <- all_donors
  
  cum_mat <- as.data.frame(cum_mat)

  # attach donors properly
  rownames(cum_mat) <- all_donors     # donors as rownames
  expr <- as.data.frame(t(cum_mat))   # transpose: donors → rows, genes → columns
  expr$donor_id <- rownames(expr)     # bring donor IDs back as column
  expr <- expr[, c("donor_id", setdiff(names(expr), "donor_id"))]
  
  return(cum_mat)
}

cum_mat <- cum_expression(df_rna, valid_genes, id_var="id",  pseudotime_var="pseudotime")

expr <- as.data.frame(t(cum_mat))
rownames(expr) <- colnames(cum_mat_clean)
expr <- expr[-959, ]
GE_pace <- cbind(id = rownames(expr), expr)


Get PCs to use as fixed covariates in LMM or in matrix eqtl
library(data.table)
library(FactoMineR)

#Remove "id"
cum_mat_clean <- cum_mat[, colnames(cum_mat) != "id"]

# Run PCA
pcs <- prcomp(cum_mat_clean, center = TRUE, scale. = TRUE)

# First 15 PCs
covariates <- as.data.frame(pcs$x[, 1:15])

# Add donor IDs as first column (assuming rownames(cum_mat) are donor IDs)
covariates$donor_id <- rownames(cum_mat)

# Reorder so donor_id is the first column
covariates <- covariates[, c("donor_id", paste0("PC", 1:15))]



