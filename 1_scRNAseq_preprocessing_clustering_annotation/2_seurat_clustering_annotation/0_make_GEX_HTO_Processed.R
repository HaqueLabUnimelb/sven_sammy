# Standard Seurat processing of the Cell Ranger output into GEX_HTO_Processed.rds, the starting object of
# liver_scRNAseq_processed_NI.R: hashtag normalisation and demultiplexing, cell filter, then the standard RNA
# workflow and a first clustering.
#
# Written for this repository (20 November 2026; settings aligned with the Methods of the manuscript on
# 21 November 2026). The steps and settings are the ones given in the Methods: HTODemux with
# positive.quantile = 0.95 after CLR normalisation of the hashtag counts; cells with more than 20 %
# mitochondrial counts, more than 60,000 UMIs or more than 7,500 detected genes are removed; 2,000 variable
# genes, 30 principal components, 20 dimensions, resolution 0.3. The RNA calls are also the ones recorded in
# GEX_HTO_Processed.rds (object@commands, in the recorded order). The object keeps all genes and the hashtag
# doublets and negatives, which the later scripts remove by selecting hashtags.
#
# Usage: Rscript 0_make_GEX_HTO_Processed.R <Cell Ranger outs/filtered_feature_bc_matrix> <output .rds>

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
matrix_dir <- if (length(args) >= 1) args[1] else "new_cellranger/outs/filtered_feature_bc_matrix"
output_rds <- if (length(args) >= 2) args[2] else "GEX_HTO_Processed.rds"

counts <- Read10X(data.dir = matrix_dir)   # list with "Gene Expression" and "Antibody Capture"
gex <- counts[["Gene Expression"]]
hto <- counts[["Antibody Capture"]]

# Cell barcodes present in both matrices
joint_barcodes <- intersect(colnames(gex), colnames(hto))
seurat_obj <- CreateSeuratObject(counts = gex[, joint_barcodes])
seurat_obj[["HTO"]] <- CreateAssayObject(counts = hto[, joint_barcodes])

# Hashtags: CLR normalisation and demultiplexing
seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")
seurat_obj <- HTODemux(seurat_obj, assay = "HTO", positive.quantile = 0.95)

# Cell filter: mitochondrial content, UMIs and detected genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj <- subset(seurat_obj, subset = percent.mt <= 20 & nCount_RNA <= 60000 & nFeature_RNA <= 7500)

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
