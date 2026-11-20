#!/usr/bin/env Rscript
# =============================================================================
# test_rctd_seurat.R: helper of build_and_test_rctd_seurat.sbatch.
#
#   Rscript test_rctd_seurat.R bin20-copy    <archive_dir> <h5ad> <out_dir> <copy.R>
#   Rscript test_rctd_seurat.R bin20-compare <new.rds> <archived.rds> <report.txt>
#   Rscript test_rctd_seurat.R rctd-copy     <archive_dir> <tabula_muris.Robj> <immune.RDS> <bin20.rds> <out_dir> <copy.R>
#   Rscript test_rctd_seurat.R rctd-compare  <new.rds> <archived.rds> <report.txt>
#   Rscript test_rctd_seurat.R seurat        <archive_dir> <GEX_HTO_Processed.rds> <sven_liver_for_paper.RDS> <out_dir>
#
# bin20-copy, rctd-copy: copies of 02_Y01453J8_bin20_h5ad_to_rds.R and rctd_script_NI.R in which only
#   paths are changed (input files, output file, working directory); the copies are run with Rscript.
# seurat (part 1, test D): loads the packages of the Seurat scripts, then runs lines 1-3, 6, 20, 65,
#   98 and 101 of a copy of liver_scRNAseq_processed_NI.R (paths changed) and compares the
#   resolution-0.7 clusters with those stored in sven_liver_for_paper.RDS.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
source(file.path(here, "test_helpers.R"))
mode <- args[1]

report <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); report <<- c(report, s) }
max_abs <- function(a, b) if (length(a) == 0) 0 else max(abs(a - b))

if (mode == "bin20-copy") {
  archive <- args[2]; h5ad <- args[3]; out <- args[4]; copy <- args[5]
  subs <- data.frame(
    line = c(27L, 28L, 29L, 33L, 34L),
    old = c('"/data/gpfs/projects/punim1236/",', '"liver_curio/2026714RCTDnew/bin20/",', '"Y01453J8_tissue_bin20_raw_counts.h5ad"',
            '"/data/gpfs/projects/punim1236/",', '"liver_curio/2026714RCTDnew/bin20/",'),
    new = c(sprintf('"%s/",', dirname(h5ad)), '"",', sprintf('"%s"', basename(h5ad)),
            sprintf('"%s/",', out), '"",'),
    stringsAsFactors = FALSE)
  copy_with_paths(file.path(archive, "3_STOmics_liver_v1.3/2_bin20_matrix/02_Y01453J8_bin20_h5ad_to_rds.R"), copy, subs)

} else if (mode == "bin20-compare") {
  suppressPackageStartupMessages(library(Seurat))
  a <- readRDS(args[2]); b <- readRDS(args[3])
  ca <- GetAssayData(a, assay = "RNA", layer = "counts"); cb <- GetAssayData(b, assay = "RNA", layer = "counts")
  say("dimensions: new %s, archived %s", paste(dim(a), collapse = " x "), paste(dim(b), collapse = " x "))
  say("counts matrix identical(): %s; same gene and bin names: %s", identical(ca, cb),
      identical(dimnames(ca), dimnames(cb)))
  say("meta.data identical(): %s", identical(a@meta.data, b@meta.data))
  say("spatial embedding identical(): %s", identical(Embeddings(a, "spatial"), Embeddings(b, "spatial")))
  say("whole object identical(): %s (version stamps and command logs may differ)", identical(a, b))
  say("SeuratObject version stamp: new %s, archived %s", as.character(a@version), as.character(b@version))
  writeLines(report, args[4])

} else if (mode == "rctd-copy") {
  archive <- args[2]; tm <- args[3]; immune <- args[4]; puck <- args[5]; out <- args[6]; copy <- args[7]
  subs <- data.frame(
    line = c(1L, 13L, 25L, 57L, 59L),
    old = c("setwd('/data/gpfs/projects/punim1236/liver_curio')",
            "'/data/gpfs/projects/punim1236/liver_curio/tabula_muris/droplet_Liver_seurat_tiss.Robj'",
            '"/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS"',
            '"/data/projects/punim1236/liver_curio/2026714RCTDnew/bin20/Y01453J8_tissue_bin20_raw_counts.rds"',
            '"/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/rctd_MNP1_2_Tcell_splitup_no_sinusoids.rds"'),
    new = c(sprintf("setwd('%s')", out), sprintf("'%s'", tm), sprintf('"%s"', immune), sprintf('"%s"', puck),
            sprintf('"%s"', file.path(out, "rctd_MNP1_2_Tcell_combined_no_sinusoids_TEST.rds"))),
    stringsAsFactors = FALSE)
  copy_with_paths(file.path(archive, "3_STOmics_liver_v1.3/3_RCTD/rctd_script_NI.R"), copy, subs)

} else if (mode == "rctd-compare") {
  suppressPackageStartupMessages(library(SpatialExperiment))
  a <- readRDS(args[2]); b <- readRDS(args[3])
  say("class: new %s, archived %s", class(a)[1], class(b)[1])
  say("dimensions (cell types x bins): new %s, archived %s", paste(dim(a), collapse = " x "), paste(dim(b), collapse = " x "))
  say("assays: new %s; archived %s", paste(assayNames(a), collapse = ", "), paste(assayNames(b), collapse = ", "))
  say("same cell types (rownames, in order): %s; same bins (colnames, in order): %s",
      identical(rownames(a), rownames(b)), identical(colnames(a), colnames(b)))
  for (an in intersect(assayNames(a), assayNames(b))) {
    x <- as.matrix(assay(a, an)); y <- as.matrix(assay(b, an))
    if (identical(dim(x), dim(y))) {
      d <- abs(x - y)
      say("assay %s: identical() %s; max |new - archived| = %.3g; mean |diff| = %.3g; bins with any |diff| > 1e-6: %d of %d",
          an, identical(x, y), max(d), mean(d), sum(colSums(d > 1e-6) > 0), ncol(d))
      if (an == "weights") {
        say("weights per cell type, max |diff|: %s",
            paste(sprintf("%s=%.3g", rownames(d), apply(d, 1, max)), collapse = "; "))
        say("weights column sums: new range %.6f-%.6f, archived range %.6f-%.6f",
            min(colSums(x)), max(colSums(x)), min(colSums(y)), max(colSums(y)))
      }
    } else say("assay %s: different dimensions", an)
  }
  say("spatialCoords identical(): %s", identical(spatialCoords(a), spatialCoords(b)))
  ma <- metadata(a); mb <- metadata(b)
  say("metadata entries: new %s; archived %s", paste(names(ma), collapse = ", "), paste(names(mb), collapse = ", "))
  writeLines(report, args[4])

} else if (mode == "seurat") {
  archive <- args[2]; gex <- args[3]; ref <- args[4]; out <- args[5]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  # 1. install-and-load test of the packages that the three part 1 scripts load
  t0 <- Sys.time()
  for (p in c("Seurat", "ggplot2", "tidyverse", "dplyr", "ggrepel")) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
    say("library(%s): loaded, version %s, from %s", p, as.character(packageVersion(p)), dirname(find.package(p)))
  }
  # 2. clustering check: copy of the revision script, paths changed, selected lines
  src <- file.path(archive, "1_scRNAseq_preprocessing_clustering_annotation/2_seurat_clustering_annotation/liver_scRNAseq_processed_NI.R")
  copy <- file.path(out, "liver_scRNAseq_processed_NI_copy.R")
  subs <- data.frame(
    line = c(NA, NA, 6L),
    old = c("/data/projects/punim1236/liver_scRNAseq/", "/data/gpfs/projects/punim1236/liver_scRNAseq/",
            sprintf('"%s/GEX_HTO_Processed.rds"', out)),
    new = c(paste0(out, "/"), paste0(out, "/"), sprintf('"%s"', gex)),
    stringsAsFactors = FALSE)
  copy_with_paths(src, copy, subs)
  setwd(out)
  run_lines(copy, "1-3,6,20,65,98,101")
  new_cl <- seurat_obj$seurat_clusters
  ref_obj <- readRDS(ref)
  ref_cl <- ref_obj$seurat_clusters
  say("cells: new %d, sven_liver_for_paper.RDS %d, same cell names in the same order: %s",
      length(new_cl), length(ref_cl), identical(names(new_cl), names(ref_cl)))
  common <- intersect(names(new_cl), names(ref_cl))
  a <- as.character(new_cl[common]); b <- as.character(ref_cl[common])
  say("clusters: new %d (%s), stored %d", nlevels(droplevels(new_cl)), paste(levels(droplevels(new_cl)), collapse = ","),
      nlevels(droplevels(ref_cl)))
  say("cells with the same cluster number: %d of %d", sum(a == b), length(common))
  tab_new <- table(factor(a, levels = levels(new_cl))); tab_ref <- table(factor(b, levels = levels(ref_cl)))
  say("cluster sizes, new:    %s", paste(names(tab_new), tab_new, sep = ":", collapse = " "))
  say("cluster sizes, stored: %s", paste(names(tab_ref), tab_ref, sep = ":", collapse = " "))
  # adjusted Rand index (base R)
  ct <- table(a, b); nij <- sum(choose(ct, 2)); ai <- sum(choose(rowSums(ct), 2)); bj <- sum(choose(colSums(ct), 2))
  n2 <- choose(length(a), 2); ari <- (nij - ai * bj / n2) / (0.5 * (ai + bj) - ai * bj / n2)
  say("adjusted Rand index new vs stored: %.6f", ari)
  if ("RNA_snn_res.0.7" %in% colnames(ref_obj@meta.data)) {
    say("RNA_snn_res.0.7 identical to stored: %s",
        identical(as.character(seurat_obj$RNA_snn_res.0.7), as.character(ref_obj$RNA_snn_res.0.7[colnames(seurat_obj)])))
  }
  say("SeuratObject version stamp of GEX_HTO_Processed.rds: %s; of sven_liver_for_paper.RDS: %s",
      as.character(sven@version), as.character(ref_obj@version))
  say("graph used by FindClusters: %s", paste(names(seurat_obj@graphs), collapse = ", "))
  say("run time %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  writeLines(report, file.path(out, "compare_clusters.txt"))
  write_session_info(file.path(out, "sessionInfo_seurat_scripts.txt"))
} else stop("unknown mode ", mode)

if (mode %in% c("bin20-compare", "rctd-compare")) write_session_info(sub("\\.txt$", "_sessionInfo.txt", args[4]))
