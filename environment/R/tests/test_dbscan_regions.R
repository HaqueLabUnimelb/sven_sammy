#!/usr/bin/env Rscript
# =============================================================================
# test_dbscan_regions.R: test driver for part 3, step 5 (DBSCAN regions, enrichment, colocalisation).
#
#   Rscript test_dbscan_regions.R <archive_dir> <rctd_rds> <hotspot_csv> <archived_output_dir> <out_dir>
#
# 1. Writes copies of liverstomics_hotspot_dbscan_NI.R and colocalisation_script_NI.R to <out_dir>,
#    in which only paths are changed: the two input files, the output folder of the figures and
#    tables (-> <out_dir>/), the folder of setwd() (-> <out_dir>/) and, because of that, the RCTD
#    file name that was relative to it (line 77, -> <rctd_rds>).
# 2. Sources the colocalisation functions, as the README says, and runs these line ranges of the
#    copy, in this order:
#      1-40     read the RCTD result and the Hotspot scores, plotting function
#      64-253   long table of RCTD weights, selection of the bins (module_5 > 4), DBSCAN, region
#               outlines and region membership (hit, sub_df), check figures
#      307-615  save hs_out and sub_df, colocalisation per region, pooled, outside the regions and
#               in 12 random windows, enrichment (all regions together and per region), figures
#    Not run: lines 42-61 (module maps and a save/reload of the RCTD object, not used later) and
#    258-304 (cell-type maps; the cell type is set by hand, see the README).
# 3. Compares with the archived results: enrichment_of_celltypes.csv, sub_df_containsdbscan_info.rds,
#    hs_out.rds and the colocalisation tables rho_region_40um<k>.csv and rho_rand_<k>.csv.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
archive <- normalizePath(args[1]); rctd_rds <- args[2]; hotspot_csv <- args[3]
orig_out <- sub("/+$", "", args[4]); out <- sub("/+$", "", args[5])
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
source(file.path(here, "test_helpers.R"))
t_start <- Sys.time()

for (d in c("zoom", "background_regions", "enrichment")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
step5 <- file.path(archive, "3_STOmics_liver_v1.3/5_DBSCAN_regions_enrichment")
copy_main <- file.path(out, "liverstomics_hotspot_dbscan_NI_copy.R")
copy_coloc <- file.path(out, "colocalisation_script_NI_copy.R")

orig_dir <- "/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/"
subs <- data.frame(
  line = c(5L, 7L, NA, NA, 75L, 77L),
  old = c("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/rctd_MNP1_2_Tcell_combined_no_sinusoids.rds",
          "/data/projects/punim1236/liver_stomicsv1.3/hotspot/hotspot_results/hotspot_spot_scores.csv",
          orig_dir,
          "/data/projects/punim1236/liver_stomicsv1.3/hotspot/hotspot_results/rctd_res_hotspot_modules.RDS",
          "setwd('/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/')",
          'rctd_loc <- "rctd_MNP1_2_Tcell_combined_no_sinusoids.rds"'),
  new = c(rctd_rds, hotspot_csv, paste0(out, "/"), file.path(out, "rctd_res_hotspot_modules.RDS"),
          sprintf("setwd('%s/')", out), sprintf('rctd_loc <- "%s"', rctd_rds)),
  stringsAsFactors = FALSE)
copy_with_paths(file.path(step5, "liverstomics_hotspot_dbscan_NI.R"), copy_main, subs)
file.copy(file.path(step5, "colocalisation_script_NI.R"), copy_coloc, overwrite = TRUE)
left <- grep(orig_dir, readLines(copy_main), fixed = TRUE, value = TRUE)
if (length(left)) stop("the original output folder is left in the copy:\n", paste(left, collapse = "\n"))

# ----------------------------------------------------------------------------- run
source(copy_coloc, encoding = "UTF-8")
run_lines(copy_main, "1-40,64-253,307-615")

# per-region enrichment is not written by the script; keep it for the record (driver addition)
utils::write.csv(enrich_reg, file.path(out, "enrichment", "enrichment_per_region_TEST_DRIVER.csv"), row.names = FALSE)
t_run <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

# ----------------------------------------------------------------------------- compare
cat("\n################ comparison with the archived results\n")
rep_lines <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); rep_lines <<- c(rep_lines, s) }

say("bins with module_5 > 4.0: %d; DBSCAN clusters kept (>= 20 bins): %d", nrow(sel), length(keep))

tab <- compare_csv(file.path(out, "enrichment/enrichment_of_celltypes.csv"),
                   file.path(step5, "result_tables/enrichment_of_celltypes.csv"), label = "enrichment_of_celltypes.csv")

old_sub <- readRDS(file.path(orig_out, "sub_df_containsdbscan_info.rds"))
new_sub <- readRDS(file.path(out, "sub_df_containsdbscan_info.rds"))
say("sub_df_containsdbscan_info.rds: identical() = %s; all.equal() = %s; rows new %d, archived %d",
    identical(new_sub, old_sub), paste(isTRUE(all.equal(new_sub, old_sub)), collapse = " "), nrow(new_sub), nrow(old_sub))
region_bins <- function(s) { u <- unique(s[, c("x", "y", "cluster")]); table(factor(u$cluster, levels = sort(unique(as.integer(as.character(u$cluster))))) ) }
rb_new <- region_bins(new_sub); rb_old <- region_bins(old_sub)
say("regions: new %d, archived %d; bins per region identical: %s", length(rb_new), length(rb_old),
    identical(as.vector(rb_new), as.vector(rb_old)) && identical(names(rb_new), names(rb_old)))
say("bins per region (new):      %s", paste(names(rb_new), as.vector(rb_new), sep = ":", collapse = " "))
say("bins per region (archived): %s", paste(names(rb_old), as.vector(rb_old), sep = ":", collapse = " "))
say("bins inside any region: new %d, archived %d", nrow(unique(new_sub[, c("x", "y")])), nrow(unique(old_sub[, c("x", "y")])))
w_diff <- if (nrow(new_sub) == nrow(old_sub)) max(abs(new_sub$weight - old_sub$weight)) else NA
say("sub_df weight column: max |new - archived| = %s", format(w_diff))

old_hs <- readRDS(file.path(orig_out, "hs_out.rds")); new_hs <- readRDS(file.path(out, "hs_out.rds"))
say("hs_out.rds: identical() = %s; all.equal() = %s", identical(new_hs, old_hs), paste(isTRUE(all.equal(new_hs, old_hs)), collapse = " "))

rho <- list()
for (nm in names(outs_regions)) {
  f <- paste0("rho_region_40um", nm, ".csv")
  rho[[f]] <- cbind(compare_csv(file.path(out, f), file.path(orig_out, f), label = f),
                    bins = nrow(region_list[[nm]]), sampled_bins = nrow(outs_regions[[nm]]$pooled_means))
}
for (nm in names(outs_rand)) {
  f <- paste0("rho_", nm, ".csv")
  rho[[f]] <- cbind(compare_csv(file.path(out, "background_regions", f), file.path(orig_out, "background_regions", f), label = f),
                    bins = nrow(rand_list[[nm]]), sampled_bins = nrow(outs_rand[[nm]]$pooled_means))
}
rho <- do.call(rbind, rho)
tab$bins <- NA; tab$sampled_bins <- NA
all_tab <- rbind(tab, rho)
print(all_tab, row.names = FALSE)
utils::write.table(all_tab, file.path(out, "compare_tables.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(rep_lines, file.path(out, "compare_objects.txt"))
say("driver run time (script sections): %.1f min", t_run)
write_session_info(file.path(out, "sessionInfo_dbscan_regions.txt"))
cat("DONE\n")
