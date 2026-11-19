setwd('/data/gpfs/projects/punim1236/liver_curio')

library(cowplot)
library(tidyverse)
library(spacexr)
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)

## prepare scRNA-seq reference

## get reference dataset
load(file = '/data/gpfs/projects/punim1236/liver_curio/tabula_muris/droplet_Liver_seurat_tiss.Robj')
tiss <- UpdateSeuratObject(tiss)

table(tiss@meta.data$cell_ontology_class)

## remove "leukocytes", "duct epithelial cell" and "endothelial cell of hepatic sinusoid", jsut left with "hepatocytes"
tiss_v2 <- subset(tiss, cells = tiss@meta.data %>% filter(cell_ontology_class %in% c('hepatocyte')) %>% rownames())
tiss_v2 <- SetIdent(tiss_v2, value = tiss_v2$cell_ontology_class)

table(tiss_v2@meta.data$cell_ontology_class)

## load immune data
immune <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS")

# Check immune object structure and extract counts properly
cat("Immune object version:", as.character(immune@version), "\n")
cat("Immune assay class:", class(immune@assays$RNA), "\n")

# For Seurat v5 Assay5 objects, use GetAssayData with layer parameter
immune_counts <- GetAssayData(immune, assay = "RNA", layer = "counts")
immune_meta <- immune@meta.data

# Recreate immune object with current Seurat version
immune <- CreateSeuratObject(counts = immune_counts, meta.data = immune_meta)
immune <- SetIdent(immune, value = immune$annotation)

table(immune@meta.data$annotation)

tiss <- merge(tiss_v2, immune)

table(tiss@active.ident)

rm(immune, immune_counts, immune_meta)
rm(tiss_v2)
gc()

## save in cell_type slot
tiss$cell_type <- tiss@active.ident

table(tiss$cell_type)

### do RCTD ###

## import data
puck1_loc <- "/data/projects/punim1236/liver_curio/2026714RCTDnew/bin20/Y01453J8_tissue_bin20_raw_counts.rds"

rctd_save_loc <- "/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/rctd_MNP1_2_Tcell_splitup_no_sinusoids.rds"

## import puck
puck1 <- readRDS(file=puck1_loc)
puck1 <- UpdateSeuratObject(puck1)

# For Seurat v5 Assay5 objects, need to join layers first
if (class(puck1@assays$RNA)[1] == "Assay5") {
  puck1 <- JoinLayers(puck1, assay = "RNA")
}

## get materials for converting puck to SpatialRNA
counts <- GetAssayData(puck1, assay = "RNA", layer = "counts")

# Extract coordinates - check multiple sources
coords <- NULL

# Method 1: Try to get from metadata (most reliable)
if ("x" %in% colnames(puck1@meta.data) && "y" %in% colnames(puck1@meta.data)) {
  coords <- puck1@meta.data[, c("x", "y")]
  coords <- as.data.frame(coords)
}

# Method 2: Try reduction named "spatial"
if (is.null(coords) && "spatial" %in% names(puck1@reductions)) {
  coords <- as.data.frame(puck1@reductions$spatial@cell.embeddings)
  colnames(coords) <- c("x", "y")
}

# Method 3: Try reduction named "map"
if (is.null(coords) && "map" %in% names(puck1@reductions)) {
  coords <- as.data.frame(puck1@reductions$map@cell.embeddings)
  colnames(coords) <- c("x", "y")
}

if (is.null(coords)) {
  stop("Cannot find spatial coordinates in puck object")
}

colnames(coords) <- c("x","y")

nUMI <- as.numeric(puck1@meta.data$nCount_RNA)
names(nUMI) <- colnames(counts)

## convert reference to SummarizedExperiment for spacexr
cell_types <- as.factor(tiss@meta.data$cell_type)
names(cell_types) <- rownames(tiss@meta.data)

# For Seurat v5 Assay5 objects, need to join layers first
if (class(tiss@assays$RNA)[1] == "Assay5") {
  tiss <- JoinLayers(tiss, assay = "RNA")
}

# Get reference counts
ref_counts <- GetAssayData(tiss, assay = "RNA", layer = "counts")

# Create SummarizedExperiment for reference
reference <- SummarizedExperiment(
  assays = list(counts = ref_counts),
  colData = DataFrame(cell_type = cell_types)
)

## Create SpatialExperiment for puck
spatial_experiment <- SpatialExperiment(
  assay = counts,
  spatialCoords = as.matrix(coords),
  colData = DataFrame(nUMI = nUMI)
)

## Create RCTD object using spacexr API
rctd_obj <- createRctd(
  spatial_experiment = spatial_experiment,
  reference_experiment = reference
)

## remove
rm(puck1, reference, spatial_experiment, coords, counts, nUMI)
gc()

## run RCTD - this is the step that runs RCTD, takes around 1 day for a typical puck of 30-40k spots.
results_spe <- runRctd(rctd_obj, rctd_mode = "full", max_cores = 4)
saveRDS(results_spe, file=rctd_save_loc)
