#!/usr/bin/env Rscript

# =============================================================================
# Run RCTD on a single Y01453J8 Stereo-seq bin20 dataset
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(spacexr)
  library(Seurat)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(S4Vectors)
})


# =============================================================================
# Input and output
# =============================================================================

work_dir <- paste0(
  "/data/gpfs/projects/punim1236/",
  "liver_curio/2026714RCTDnew/bin20"
)

puck_file <- file.path(
  work_dir,
  "Y01453J8_tissue_bin20_raw_counts.rds"
)

rctd_output <- file.path(
  work_dir,
  "Y01453J8_tissue_bin20_rctd_results.rds"
)

tabula_muris_file <- paste0(
  "/data/gpfs/projects/punim1236/",
  "liver_curio/tabula_muris/",
  "droplet_Liver_seurat_tiss.Robj"
)

immune_file <- paste0(
  "/data/projects/punim1236/",
  "liver_scRNAseq/",
  "sven_liver_new_cytotoxic_cluster_forRCTD.RDS"
)

setwd(work_dir)


# =============================================================================
# Check files
# =============================================================================

required_files <- c(
  puck_file,
  tabula_muris_file,
  immune_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "The following input files do not exist:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

cat("============================================================\n")
cat("RCTD analysis: Y01453J8 bin20\n")
cat("============================================================\n\n")

cat(sprintf(
  "spacexr version: %s\n",
  as.character(packageVersion("spacexr"))
))

cat(sprintf(
  "Spatial input:\n%s\n\n",
  puck_file
))


# =============================================================================
# 1. Prepare Tabula Muris liver reference
# =============================================================================

cat("Loading Tabula Muris liver reference...\n")

load(
  file = tabula_muris_file
)

if (!exists("tiss")) {
  stop(
    "The Tabula Muris Robj file does not contain an object named 'tiss'."
  )
}

tiss <- UpdateSeuratObject(tiss)

if (
  !"cell_ontology_class" %in%
  colnames(tiss@meta.data)
) {
  stop(
    "Tabula Muris metadata does not contain cell_ontology_class."
  )
}

cat("\nOriginal Tabula Muris cell types:\n")
print(
  table(
    tiss$cell_ontology_class,
    useNA = "ifany"
  )
)

selected_liver_types <- c(
  "duct epithelial cell",
  "endothelial cell of hepatic sinusoid",
  "hepatocyte"
)

selected_cells <- rownames(tiss@meta.data)[
  tiss$cell_ontology_class %in%
  selected_liver_types
]

if (length(selected_cells) == 0) {
  stop(
    "No selected liver cell types were found."
  )
}

tiss_v2 <- subset(
  tiss,
  cells = selected_cells
)

tiss_v2$cell_type <- as.character(
  tiss_v2$cell_ontology_class
)

Idents(tiss_v2) <- "cell_type"

cat("\nSelected liver reference cell types:\n")
print(table(tiss_v2$cell_type))

rm(tiss)
gc()


# =============================================================================
# 2. Prepare immune reference
# =============================================================================

cat("\nLoading immune reference...\n")

immune <- readRDS(
  immune_file
)

immune <- UpdateSeuratObject(
  immune
)

if (
  !"annotation" %in%
  colnames(immune@meta.data)
) {
  stop(
    "Immune reference metadata does not contain annotation."
  )
}

cat(
  "Immune object version:",
  as.character(immune@version),
  "\n"
)

cat(
  "Immune RNA assay class:",
  class(immune[["RNA"]])[1],
  "\n"
)

if (inherits(immune[["RNA"]], "Assay5")) {
  immune <- JoinLayers(
    immune,
    assay = "RNA"
  )
}

immune_counts <- GetAssayData(
  immune,
  assay = "RNA",
  layer = "counts"
)

immune_meta <- immune@meta.data[
  colnames(immune_counts),
  ,
  drop = FALSE
]

valid_immune_cells <- (
  !is.na(immune_meta$annotation) &
  immune_meta$annotation != ""
)

immune_counts <- immune_counts[
  ,
  valid_immune_cells,
  drop = FALSE
]

immune_meta <- immune_meta[
  valid_immune_cells,
  ,
  drop = FALSE
]

immune <- CreateSeuratObject(
  counts = immune_counts,
  meta.data = immune_meta
)

immune$cell_type <- as.character(
  immune$annotation
)

Idents(immune) <- "cell_type"

cat("\nImmune cell types:\n")
print(table(immune$cell_type))

rm(
  immune_counts,
  immune_meta
)

gc()


# =============================================================================
# 3. Merge references
# =============================================================================

cat("\nMerging liver and immune references...\n")

reference_seurat <- merge(
  x = tiss_v2,
  y = immune,
  merge.data = FALSE
)

rm(
  tiss_v2,
  immune
)

gc()

if (inherits(reference_seurat[["RNA"]], "Assay5")) {
  reference_seurat <- JoinLayers(
    reference_seurat,
    assay = "RNA"
  )
}

ref_counts <- GetAssayData(
  reference_seurat,
  assay = "RNA",
  layer = "counts"
)

ref_counts <- methods::as(
  ref_counts,
  "dgCMatrix"
)

cell_types <- factor(
  reference_seurat$cell_type
)

names(cell_types) <- colnames(
  reference_seurat
)

if (anyNA(cell_types)) {
  stop(
    "Some reference cells do not have a cell_type annotation."
  )
}

cat("\nMerged reference cell types:\n")
print(table(cell_types))

reference_coldata <- S4Vectors::DataFrame(
  cell_type = cell_types,
  row.names = colnames(ref_counts)
)

reference <- SummarizedExperiment(
  assays = list(
    counts = ref_counts
  ),
  colData = reference_coldata
)

cat(sprintf(
  "\nReference dimensions: %s genes x %s cells\n",
  format(nrow(reference), big.mark = ","),
  format(ncol(reference), big.mark = ",")
))

rm(
  reference_seurat,
  ref_counts,
  reference_coldata,
  cell_types
)

gc()


# =============================================================================
# 4. Load Y01453J8 bin20 spatial Seurat object
# =============================================================================

cat("\nLoading Y01453J8 bin20 spatial data...\n")

puck <- readRDS(
  puck_file
)

puck <- UpdateSeuratObject(
  puck
)

if (inherits(puck[["RNA"]], "Assay5")) {
  puck <- JoinLayers(
    puck,
    assay = "RNA"
  )
}

counts <- GetAssayData(
  puck,
  assay = "RNA",
  layer = "counts"
)

counts <- methods::as(
  counts,
  "dgCMatrix"
)


# =============================================================================
# 5. Extract spatial coordinates
# =============================================================================

coords <- NULL

# Preferred source: Seurat metadata
if (
  all(
    c("x", "y") %in%
    colnames(puck@meta.data)
  )
) {
  coords <- as.data.frame(
    puck@meta.data[
      ,
      c("x", "y"),
      drop = FALSE
    ]
  )

  cat(
    "Coordinates extracted from Seurat metadata.\n"
  )
}

# Alternative: spatial reduction
if (
  is.null(coords) &&
  "spatial" %in%
  names(puck@reductions)
) {
  coords <- as.data.frame(
    Embeddings(
      puck,
      reduction = "spatial"
    )
  )

  coords <- coords[, 1:2, drop = FALSE]
  colnames(coords) <- c("x", "y")

  cat(
    "Coordinates extracted from spatial reduction.\n"
  )
}

# Alternative: map reduction
if (
  is.null(coords) &&
  "map" %in%
  names(puck@reductions)
) {
  coords <- as.data.frame(
    Embeddings(
      puck,
      reduction = "map"
    )
  )

  coords <- coords[, 1:2, drop = FALSE]
  colnames(coords) <- c("x", "y")

  cat(
    "Coordinates extracted from map reduction.\n"
  )
}

if (is.null(coords)) {
  stop(
    "Cannot find spatial coordinates in the Seurat object."
  )
}

missing_coordinates <- setdiff(
  colnames(counts),
  rownames(coords)
)

if (length(missing_coordinates) > 0) {
  stop(
    sprintf(
      "%d bins do not have spatial coordinates.",
      length(missing_coordinates)
    )
  )
}

# Explicitly reorder coordinates to match counts
coords <- coords[
  colnames(counts),
  c("x", "y"),
  drop = FALSE
]

if (!identical(
  rownames(coords),
  colnames(counts)
)) {
  stop(
    "Coordinate order does not match count matrix order."
  )
}

if (any(!is.finite(as.matrix(coords)))) {
  stop(
    "Coordinates contain NA or infinite values."
  )
}


# =============================================================================
# 6. Calculate nUMI directly from raw counts
# =============================================================================

nUMI <- Matrix::colSums(
  counts
)

names(nUMI) <- colnames(
  counts
)

cat(sprintf(
  "\nSpatial data dimensions: %s genes x %s bins\n",
  format(nrow(counts), big.mark = ","),
  format(ncol(counts), big.mark = ",")
))

cat(sprintf(
  "Mean UMI per bin: %.2f\n",
  mean(nUMI)
))

cat(sprintf(
  "Median UMI per bin: %.2f\n",
  median(nUMI)
))

cat(sprintf(
  "Bins with at least 100 UMI: %s / %s\n",
  format(sum(nUMI >= 100), big.mark = ","),
  format(length(nUMI), big.mark = ",")
))


# =============================================================================
# 7. Check shared genes
# =============================================================================

common_genes <- intersect(
  rownames(counts),
  rownames(reference)
)

cat(sprintf(
  "\nShared genes between spatial data and reference: %s\n",
  format(length(common_genes), big.mark = ",")
))

if (length(common_genes) < 100) {
  stop(
    paste0(
      "Too few shared genes between spatial data and reference. ",
      "Check whether both datasets use the same gene-name format."
    )
  )
}


# =============================================================================
# 8. Create SpatialExperiment
# =============================================================================

spatial_coldata <- S4Vectors::DataFrame(
  nUMI = as.numeric(nUMI),
  row.names = colnames(counts)
)

spatial_experiment <- SpatialExperiment(
  assays = list(
    counts = counts
  ),
  colData = spatial_coldata,
  spatialCoords = as.matrix(coords)
)

cat("\nSpatialExperiment created successfully.\n")


# =============================================================================
# 9. Create and run RCTD
# =============================================================================

cat("\nCreating RCTD object...\n")

rctd_obj <- createRctd(
  spatial_experiment = spatial_experiment,
  reference_experiment = reference
)

rm(
  puck,
  counts,
  coords,
  nUMI,
  spatial_coldata,
  spatial_experiment,
  reference
)

gc()

max_cores <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "SLURM_CPUS_PER_TASK",
      unset = "4"
    )
  )
)

if (
  is.na(max_cores) ||
  max_cores < 1
) {
  max_cores <- 4L
}

cat(sprintf(
  "\nRunning RCTD in doublet mode with %d cores...\n",
  max_cores
))

results_spe <- runRctd(
  rctd_obj,
  rctd_mode = "doublet",
  max_cores = max_cores
)


# =============================================================================
# 10. Save result
# =============================================================================

cat(sprintf(
  "\nSaving RCTD result:\n%s\n",
  rctd_output
))

saveRDS(
  results_spe,
  file = rctd_output,
  compress = FALSE
)

cat("\nRCTD analysis completed successfully.\n")