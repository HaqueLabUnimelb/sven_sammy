#!/usr/bin/env Rscript

# =============================================================================
# Convert Y01453J8 bin20 raw-count H5AD into a Seurat RDS object
# for downstream RCTD analysis.
#
# This version reads the AnnData dataframe index name from the "_index"
# attribute. It therefore supports named indices such as:
#   obs/bin_id
#   var/gene
# as well as the default:
#   obs/_index
#   var/_index
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(hdf5r)
})

# -----------------------------------------------------------------------------
# Input and output
# -----------------------------------------------------------------------------

h5ad_file <- paste0(
  "/data/gpfs/projects/punim1236/",
  "liver_curio/2026714RCTDnew/bin20/",
  "Y01453J8_tissue_bin20_raw_counts.h5ad"
)

output_rds <- paste0(
  "/data/gpfs/projects/punim1236/",
  "liver_curio/2026714RCTDnew/bin20/",
  "Y01453J8_tissue_bin20_raw_counts.rds"
)

if (!file.exists(h5ad_file)) {
  stop(sprintf("Input H5AD file does not exist:\n%s", h5ad_file))
}

cat("============================================================\n")
cat("Convert Y01453J8 bin20 H5AD to Seurat RDS\n")
cat("============================================================\n\n")
cat(sprintf("Input:\n%s\n\n", h5ad_file))

# -----------------------------------------------------------------------------
# Helper: read obs/var index according to AnnData dataframe specification
# -----------------------------------------------------------------------------

read_anndata_index <- function(h5, group_name) {

  if (!h5$exists(group_name)) {
    stop(sprintf("H5AD group '%s' does not exist.", group_name))
  }

  group <- h5[[group_name]]
  attrs <- hdf5r::h5attributes(group)

  index_name <- NULL

  if ("_index" %in% names(attrs)) {
    index_name <- as.character(attrs[["_index"]])[1]
  }

  # Fallback for older/default H5AD files
  candidates <- unique(c(
    index_name,
    "_index",
    "index",
    if (group_name == "obs") "bin_id" else "gene"
  ))

  candidates <- candidates[
    !is.na(candidates) &
      nzchar(candidates)
  ]

  for (candidate in candidates) {
    dataset_path <- paste0(group_name, "/", candidate)

    if (h5$exists(dataset_path)) {
      values <- h5[[dataset_path]]$read()

      cat(sprintf(
        "Using %s index dataset: %s\n",
        group_name,
        dataset_path
      ))

      return(as.character(values))
    }
  }

  cat(sprintf(
    "\nContents of H5AD group '%s':\n",
    group_name
  ))
  print(group$ls(recursive = FALSE))

  stop(sprintf(
    paste0(
      "Could not find the dataframe index for group '%s'. ",
      "The group _index attribute was: %s"
    ),
    group_name,
    ifelse(is.null(index_name), "NULL", index_name)
  ))
}

# -----------------------------------------------------------------------------
# Open H5AD
# -----------------------------------------------------------------------------

h5 <- H5File$new(h5ad_file, mode = "r")

on.exit(
  try(h5$close_all(), silent = TRUE),
  add = TRUE
)

# -----------------------------------------------------------------------------
# Read bin IDs and gene names
# -----------------------------------------------------------------------------

cell_ids <- read_anndata_index(h5, "obs")
gene_names <- read_anndata_index(h5, "var")

n_bins <- length(cell_ids)
n_genes <- length(gene_names)

cat(sprintf(
  "\nH5AD dimensions: %s bins x %s genes\n",
  format(n_bins, big.mark = ","),
  format(n_genes, big.mark = ",")
))

if (anyNA(cell_ids) || any(cell_ids == "")) {
  stop("Missing or empty bin IDs were found.")
}

if (anyDuplicated(cell_ids)) {
  stop("Duplicated bin IDs were found.")
}

if (anyNA(gene_names) || any(gene_names == "")) {
  stop("Missing or empty gene names were found.")
}

if (anyDuplicated(gene_names)) {
  cat("Warning: duplicated gene names found; applying make.unique().\n")
  gene_names <- make.unique(gene_names)
}

# -----------------------------------------------------------------------------
# Read count matrix
# -----------------------------------------------------------------------------

if (
  h5$exists("X/data") &&
  h5$exists("X/indices") &&
  h5$exists("X/indptr")
) {

  x_attrs <- hdf5r::h5attributes(h5[["X"]])
  encoding_type <- ""

  if ("encoding-type" %in% names(x_attrs)) {
    encoding_type <- as.character(x_attrs[["encoding-type"]])[1]
  }

  cat(sprintf(
    "Sparse X encoding: %s\n",
    ifelse(nzchar(encoding_type), encoding_type, "not recorded")
  ))

  data_values <- as.numeric(h5[["X/data"]]$read())
  indices <- as.numeric(h5[["X/indices"]]$read())
  indptr <- as.numeric(h5[["X/indptr"]]$read())

  cat(sprintf(
    "Non-zero entries: %s\n",
    format(length(data_values), big.mark = ",")
  ))

  if (length(data_values) != length(indices)) {
    stop("X/data and X/indices have different lengths.")
  }

  if (
    grepl("csr", encoding_type, ignore.case = TRUE) ||
    length(indptr) == n_bins + 1
  ) {

    cat("Reading X as CSR: bins x genes.\n")

    row_indices <- rep(
      seq_len(n_bins),
      times = as.integer(diff(indptr))
    )

    if (length(row_indices) != length(data_values)) {
      stop(
        "CSR reconstruction failed: row-index length does not match X/data."
      )
    }

    bin_by_gene <- sparseMatrix(
      i = row_indices,
      j = indices + 1,
      x = data_values,
      dims = c(n_bins, n_genes),
      giveCsparse = TRUE
    )

  } else if (
    grepl("csc", encoding_type, ignore.case = TRUE) ||
    length(indptr) == n_genes + 1
  ) {

    cat("Reading X as CSC: bins x genes.\n")

    column_indices <- rep(
      seq_len(n_genes),
      times = as.integer(diff(indptr))
    )

    if (length(column_indices) != length(data_values)) {
      stop(
        "CSC reconstruction failed: column-index length does not match X/data."
      )
    }

    bin_by_gene <- sparseMatrix(
      i = indices + 1,
      j = column_indices,
      x = data_values,
      dims = c(n_bins, n_genes),
      giveCsparse = TRUE
    )

  } else {
    stop(sprintf(
      paste0(
        "Cannot determine sparse X format. ",
        "length(indptr)=%d; expected %d for CSR or %d for CSC."
      ),
      length(indptr),
      n_bins + 1,
      n_genes + 1
    ))
  }

  # Seurat expects genes x bins
  counts <- t(bin_by_gene)
  counts <- methods::as(counts, "dgCMatrix")

} else if (h5$exists("X")) {

  cat("Reading dense X matrix...\n")

  dense_x <- as.matrix(h5[["X"]]$read())

  if (
    nrow(dense_x) == n_bins &&
    ncol(dense_x) == n_genes
  ) {
    counts <- Matrix::Matrix(t(dense_x), sparse = TRUE)

  } else if (
    nrow(dense_x) == n_genes &&
    ncol(dense_x) == n_bins
  ) {
    counts <- Matrix::Matrix(dense_x, sparse = TRUE)

  } else {
    stop(sprintf(
      paste0(
        "Unexpected dense X dimensions: %d x %d. ",
        "Expected %d x %d."
      ),
      nrow(dense_x),
      ncol(dense_x),
      n_bins,
      n_genes
    ))
  }

  counts <- methods::as(counts, "dgCMatrix")

} else {
  stop("Cannot find the X count matrix in the H5AD file.")
}

rownames(counts) <- gene_names
colnames(counts) <- cell_ids

# -----------------------------------------------------------------------------
# Validate raw counts
# -----------------------------------------------------------------------------

if (length(counts@x) > 0) {

  if (any(!is.finite(counts@x))) {
    stop("Count matrix contains NA, NaN or infinite values.")
  }

  if (any(counts@x < 0)) {
    stop("Count matrix contains negative values.")
  }

  non_integer <- abs(counts@x - round(counts@x)) > 1e-6

  if (any(non_integer)) {
    examples <- head(counts@x[non_integer], 10)

    stop(paste0(
      "Count matrix contains non-integer values and may not be raw counts. ",
      "Examples: ",
      paste(examples, collapse = ", ")
    ))
  }

  counts@x <- round(counts@x)
}

counts <- Matrix::drop0(counts)

cat(sprintf(
  "Count matrix: %s genes x %s bins\n",
  format(nrow(counts), big.mark = ","),
  format(ncol(counts), big.mark = ",")
))

# -----------------------------------------------------------------------------
# Read coordinates
# -----------------------------------------------------------------------------

required_coordinate_paths <- c("obs/x", "obs/y")

missing_coordinate_paths <- required_coordinate_paths[
  !vapply(
    required_coordinate_paths,
    function(path) h5$exists(path),
    logical(1)
  )
]

if (length(missing_coordinate_paths) > 0) {
  stop(sprintf(
    "Missing coordinate fields: %s",
    paste(missing_coordinate_paths, collapse = ", ")
  ))
}

coords <- data.frame(
  x = as.numeric(h5[["obs/x"]]$read()),
  y = as.numeric(h5[["obs/y"]]$read()),
  row.names = cell_ids,
  check.names = FALSE
)

if (nrow(coords) != ncol(counts)) {
  stop("Coordinate count does not match the number of bins.")
}

if (any(!is.finite(as.matrix(coords)))) {
  stop("Coordinates contain NA, NaN or infinite values.")
}

# Close the HDF5 file before building the Seurat object
h5$close_all()

# -----------------------------------------------------------------------------
# Recalculate QC directly from raw counts
# -----------------------------------------------------------------------------

nCount_RNA <- Matrix::colSums(counts)
nFeature_RNA <- Matrix::colSums(counts > 0)

names(nCount_RNA) <- colnames(counts)
names(nFeature_RNA) <- colnames(counts)

cat(sprintf(
  "Mean UMI per bin: %.2f\n",
  mean(nCount_RNA)
))

cat(sprintf(
  "Median UMI per bin: %.2f\n",
  median(nCount_RNA)
))

cat(sprintf(
  "Mean genes per bin: %.2f\n",
  mean(nFeature_RNA)
))

cat(sprintf(
  "Median genes per bin: %.2f\n",
  median(nFeature_RNA)
))

# -----------------------------------------------------------------------------
# Create Seurat object
# -----------------------------------------------------------------------------

cat("\nCreating Seurat object...\n")

metadata <- data.frame(
  x = coords$x,
  y = coords$y,
  row.names = cell_ids,
  check.names = FALSE
)

seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "Y01453J8_bin20",
  meta.data = metadata,
  min.cells = 0,
  min.features = 0
)

coords <- coords[
  colnames(seurat_obj),
  ,
  drop = FALSE
]

spatial_coords <- as.matrix(coords)
colnames(spatial_coords) <- c("SPATIAL_1", "SPATIAL_2")

seurat_obj[["spatial"]] <- CreateDimReducObject(
  embeddings = spatial_coords,
  key = "SPATIAL_",
  assay = "RNA"
)

cat("\nSeurat object summary:\n")
print(seurat_obj)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

cat(sprintf(
  "\nSaving RDS:\n%s\n",
  output_rds
))

saveRDS(
  seurat_obj,
  file = output_rds,
  compress = FALSE
)

if (!file.exists(output_rds) || file.info(output_rds)$size == 0) {
  stop("RDS saving failed or produced an empty file.")
}

cat(sprintf(
  "Saved RDS size: %.2f GB\n",
  file.info(output_rds)$size / 1024^3
))

cat("\nDone.\n")
