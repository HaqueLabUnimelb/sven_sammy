# Standard Seurat processing of the Cell Ranger output into GEX_HTO_Processed.rds, the starting object of
# liver_scRNAseq_processed_NI.R: hashtag normalisation and demultiplexing, then the standard RNA workflow
# and a first clustering.
#
# Written for this repository (20 November 2026). The calls and settings are the ones recorded in
# GEX_HTO_Processed.rds (object@commands, in the recorded order); HTODemux uses the standard
# positive.quantile = 0.99. No cell or gene filter: the object keeps all genes and the hashtag doublets
# and negatives, which the later scripts remove by selecting hashtags.
#
# Usage: Rscript 0_make_GEX_HTO_Processed.R <Cell Ranger outs/filtered_feature_bc_matrix> <output .rds>

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
matrix_dir <- if (length(args) >= 1) args[1] else "new_cellranger/outs/filtered_feature_bc_matrix"
output_rds <- if (length(args) >= 2) args[2] else "GEX_HTO_Processed.rds"

counts <- Read10X(data.dir = matrix_dir)   # list with "Gene Expression" and "Antibody Capture"
seurat_obj <- CreateSeuratObject(counts = counts[["Gene Expression"]])
seurat_obj[["HTO"]] <- CreateAssayObject(counts = counts[["Antibody Capture"]])

# Hashtags: CLR normalisation and demultiplexing
seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")
seurat_obj <- HTODemux(seurat_obj, assay = "HTO", positive.quantile = 0.99)

# RNA: standard workflow and first clustering
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

print(table(seurat_obj$HTO_classification.global))
saveRDS(seurat_obj, file = output_rds)
