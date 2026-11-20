# =============================================================================
# sessioninfo_at_exit.R: R profile that writes sessionInfo() when an Rscript session ends normally.
#
#   R_PROFILE_USER=sessioninfo_at_exit.R SESSIONINFO_FILE=<out.txt> Rscript <script.R>
#
# Used for the copies of archived scripts that do not write sessionInfo() themselves
# (02_Y01453J8_bin20_h5ad_to_rds.R, rctd_script_NI.R). The script itself is run unchanged with
# Rscript; this profile only defines .Last(), which R calls at the normal end of the session.
# (Sourcing the script from a wrapper instead would change its behaviour: a top-level on.exit(),
# as in 02_Y01453J8_bin20_h5ad_to_rds.R, then runs at once instead of at the end.)
# =============================================================================
.Last <- function() {
  out <- Sys.getenv("SESSIONINFO_FILE")
  if (nzchar(out)) {
    writeLines(c(capture.output(utils::sessionInfo()), "", ".libPaths():", paste0("  ", .libPaths())), out)
  }
}
