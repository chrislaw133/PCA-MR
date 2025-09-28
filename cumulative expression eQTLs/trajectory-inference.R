library(sceasy)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)

#Load in rna seq
use_condaenv("atacseq-env",
             conda="path/to/conda",
             required=TRUE)

rds <- sceasy::convertFormat("path/to/combined.h5ad", from="anndata", to="seurat")

#Quality control
mat <- rds@assays$RNA@data
rds$percent.mt <- Matrix::colSums(mat[mito.genes, , drop = FALSE]) / Matrix::colSums(mat) * 100
rds$percent.rb <- Matrix::colSums(mat[ribo.genes, , drop = FALSE]) / Matrix::colSums(mat) * 100
rds$nFeature_RNA <- Matrix::colSums(rds@assays$RNA@data > 0)
rds <- subset(rds, subset = nFeature_RNA > 200 & percent.mt < 20 & percent.rb < 20)

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
fwrite(dfti, "")

