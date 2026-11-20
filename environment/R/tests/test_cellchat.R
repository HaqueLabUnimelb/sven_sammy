#!/usr/bin/env Rscript
# =============================================================================
# test_cellchat.R: helper of build_and_test_cellchat.sbatch.
#
#   Rscript test_cellchat.R copy    <run> <archive_dir> <input_rds> <out_dir> <copy.R>
#   Rscript test_cellchat.R compare <run> <archive_dir> <out_dir> <report.tsv>
#
# <run> is "infected" (run 1) or "cytotoxic" (run 2).
# copy:    writes a copy of the archived CellChat script in which only input_rds (line 19) and
#          outdir (line 20) are changed.
# compare: compares the tables 01 to 04 written by the copy with the archived tables of that run.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
source(file.path(here, "test_helpers.R"))

runs <- list(
  infected = list(
    script = "2_CellChat/run_cellchat_MNP2_infected_vs_naive_CD4_latest.R",
    archived = "2_CellChat/MNP2_CD4_20261023",
    input = "/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS",
    outdir = "/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/cellchat_latest_MNP2_CD4_20261023_matched_original_params",
    tables = c("01_input_annotation_counts.csv", "02_input_annotation_by_HTO.csv",
               "03_differential_LR_CD4_to_MNP2_infected_minus_naive.csv",
               "04_differential_LR_MNP2_to_CD4_infected_minus_naive.csv")),
  cytotoxic = list(
    script = "2_CellChat/run_cellchat_MNP2_cytotoxic_vs_naive_CD4_latest.R",
    archived = "2_CellChat/MNP2_cytotoxicCD4_vs_naiveCD4_20261023",
    input = "/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS",
    outdir = "/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/cellchat_latest_MNP2_cytotoxicCD4_vs_naiveCD4_20261023_matched_original_params",
    tables = c("01_input_annotation_counts.csv", "02_input_annotation_by_HTO.csv",
               "03_differential_LR_CD4_to_MNP2_cytotoxic_minus_naive.csv",
               "04_differential_LR_MNP2_to_CD4_cytotoxic_minus_naive.csv"))
)

mode <- args[1]; run <- runs[[args[2]]]; archive <- args[3]
if (is.null(run)) stop("unknown run: ", args[2])

if (mode == "copy") {
  input_rds <- args[4]; out_dir <- args[5]; copy <- args[6]
  subs <- data.frame(
    line = c(19L, 20L),
    old = c(sprintf('input_rds <- "%s"', run$input), sprintf('outdir <- "%s"', run$outdir)),
    new = c(sprintf('input_rds <- "%s"', input_rds), sprintf('outdir <- "%s"', out_dir)),
    stringsAsFactors = FALSE)
  copy_with_paths(file.path(archive, run$script), copy, subs)
} else if (mode == "compare") {
  out_dir <- args[4]; report <- args[5]
  res <- do.call(rbind, lapply(run$tables, function(f)
    compare_csv(file.path(out_dir, f), file.path(archive, run$archived, f), label = f)))
  print(res, row.names = FALSE)
  utils::write.table(res, report, sep = "\t", quote = FALSE, row.names = FALSE)
  # 00_run_config.txt: show which lines differ (paths and versions are expected to differ)
  a <- readLines(file.path(out_dir, "00_run_config.txt")); b <- readLines(file.path(archive, run$archived, "00_run_config.txt"))
  cat("\n00_run_config.txt lines that differ from the archived file:\n")
  for (i in seq_len(max(length(a), length(b)))) if (!identical(a[i], b[i])) cat(sprintf("  new: %s\n  old: %s\n", a[i], b[i]))
  # number of pairs per direction and sign, as given in the README
  for (f in run$tables[3:4]) {
    d <- utils::read.csv(file.path(out_dir, f), check.names = FALSE)
    cat(sprintf("%s: %d pairs, %d positive, %d negative\n", f, nrow(d), sum(d$signed_log_prob > 0), sum(d$signed_log_prob < 0)))
  }
} else stop("mode must be copy or compare")
