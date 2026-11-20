#!/usr/bin/env Rscript
# =============================================================================
# install_dbscan_regions_env.R
#
# Builds an R package library for the DBSCAN regions, enrichment and colocalisation step of part 3
# (3_STOmics_liver_v1.3/5_DBSCAN_regions_enrichment/liverstomics_hotspot_dbscan_NI.R and
# colocalisation_script_NI.R).
#
# PROPOSAL: the original R session ran R 4.3.1 (Producer field of its PDF files, version stamp of its .rds
# files); its package versions are not recorded.
# This lock table takes every package that the two scripts load, and their dependencies, at the
# versions of the RCTD library (rctd_seurat_packages.tsv, R 4.5.0, installed 14 July 2026), which
# wrote the RCTD result that the scripts read. dbscan is not in that library; the lock table uses
# dbscan 1.2.7, the version on CRAN when the archived enrichment table was written (27 October 2026;
# dbscan 1.2.7 was published on 5 October 2026 and replaced by 1.3.0 later).
#
# usage:   Rscript install_dbscan_regions_env.R <library_dir> [n_cpus]
#
# R:             4.5.0 (proposed; another 4.5.x is accepted with a warning)
# Bioconductor:  3.22
# Lock table:    dbscan_regions_packages.tsv, in the folder of this script (columns package, version,
#                source, origin)
#
# What the script does:
#   1. reads the lock table;
#   2. downloads the source tarball of every package at exactly the locked version:
#        CRAN          <CRAN_MIRROR>/src/contrib/<pkg>_<ver>.tar.gz or .../src/contrib/Archive/<pkg>/
#        Bioconductor  <BIOC_MIRROR>/packages/3.22/bioc/src/contrib/<pkg>_<ver>.tar.gz or its Archive/
#        GitHub        https://codeload.github.com/<user>/<repo>/tar.gz/<commit>
#   3. puts the CRAN and Bioconductor tarballs into a local package repository and installs them all
#      with install.packages() from that repository only, so that no other version can be chosen and
#      R works out the installation order; packages are compiled from source;
#   4. installs the GitHub packages from the downloaded commit and records the commit in their
#      DESCRIPTION (RemoteSha), as remotes::install_github() does;
#   5. compares every installed version with the lock table, writes <env>_lock_check.tsv to the
#      working directory and exits with status 1 if a package is missing or has another version.
#
# The library is built from <library_dir> and the base R library only: R_LIBS, R_LIBS_USER and
# R_LIBS_SITE are set to <library_dir> for the duration of the script, so packages in a personal or
# site library are neither used nor changed.
#
# Optional environment variables:
#   CRAN_MIRROR      default https://cloud.r-project.org
#   BIOC_MIRROR      default https://bioconductor.org
#   LOCK_WORK_DIR    folder for downloads and build logs, default <tempdir>/lock_build
#
# Packages that need system libraries (compiled packages): a C, C++ and Fortran compiler, make,
# pkg-config, and the development files of zlib, libcurl, openssl and ImageMagick (magick, needed by
# SpatialExperiment). On Spartan the module R/4.5.0 loads all of them
# (module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0).
# =============================================================================

ENV_NAME       <- "dbscan_regions"
LOCK_FILE      <- "dbscan_regions_packages.tsv"
R_RECORDED     <- "4.5.0"
BIOC_VERSION   <- "3.22"

# ----------------------------------------------------------------------------- arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: Rscript install_", ENV_NAME, "_env.R <library_dir> [n_cpus]", call. = FALSE)
}
lib <- args[1]
ncpus <- if (length(args) >= 2) as.integer(args[2]) else {
  n <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  if (is.na(n)) parallel::detectCores() else n
}
if (is.na(ncpus) || ncpus < 1) ncpus <- 1L

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
lock_path <- file.path(script_dir, LOCK_FILE)
if (!file.exists(lock_path)) stop("lock table not found: ", lock_path, call. = FALSE)

cran_mirror <- sub("/+$", "", Sys.getenv("CRAN_MIRROR", "https://cloud.r-project.org"))
bioc_mirror <- sub("/+$", "", Sys.getenv("BIOC_MIRROR", "https://bioconductor.org"))
work <- Sys.getenv("LOCK_WORK_DIR", file.path(tempdir(), "lock_build"))

log_msg <- function(...) { cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n"); flush(stdout()) }

# ----------------------------------------------------------------------------- R version
r_now <- as.character(getRversion())
r_minor_now <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
r_minor_rec <- sub("^([0-9]+\\.[0-9]+).*$", "\\1", R_RECORDED)
if (r_minor_now != r_minor_rec) {
  stop(sprintf("this lock table needs R %s.x (recorded: R %s); this is R %s", r_minor_rec, R_RECORDED, r_now),
       call. = FALSE)
}
if (r_now != R_RECORDED) {
  warning(sprintf("the archived results were produced with R %s; this is R %s", R_RECORDED, r_now),
          call. = FALSE, immediate. = TRUE)
}

# ----------------------------------------------------------------------------- library paths
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
lib <- normalizePath(lib)
Sys.setenv(R_LIBS = lib, R_LIBS_USER = lib, R_LIBS_SITE = lib)
.libPaths(c(lib, .Library))
options(timeout = max(600, getOption("timeout")), warn = 1)

log_msg("environment:", ENV_NAME, "| R", r_now, "| library", lib, "| CPUs", ncpus)
log_msg("lock table:", lock_path)

lock <- read.delim(lock_path, colClasses = "character", check.names = FALSE)
stopifnot(all(c("package", "version", "source") %in% colnames(lock)), !anyDuplicated(lock$package))
lock$kind <- ifelse(lock$source == "CRAN", "CRAN",
             ifelse(grepl("^Bioconductor", lock$source), "Bioconductor",
             ifelse(grepl("^GitHub ", lock$source), "GitHub", NA)))
if (anyNA(lock$kind)) stop("unknown source in lock table: ", paste(lock$package[is.na(lock$kind)], collapse = ", "))
bioc_rows <- lock$kind == "Bioconductor"
if (any(bioc_rows & lock$source != paste("Bioconductor", BIOC_VERSION))) {
  stop("lock table names a Bioconductor release other than ", BIOC_VERSION)
}
log_msg(sprintf("packages: %d (CRAN %d, Bioconductor %s %d, GitHub %d)", nrow(lock),
                sum(lock$kind == "CRAN"), BIOC_VERSION, sum(bioc_rows), sum(lock$kind == "GitHub")))

# ----------------------------------------------------------------------------- download
repo_dir <- file.path(work, "repo")
contrib <- file.path(repo_dir, "src", "contrib")
log_dir <- file.path(work, "install_logs")
dir.create(contrib, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

fetch <- function(urls, dest) {
  for (u in urls) {
    ok <- tryCatch(suppressWarnings(utils::download.file(u, dest, mode = "wb", quiet = TRUE)) == 0,
                   error = function(e) FALSE)
    if (isTRUE(ok) && file.exists(dest) && file.size(dest) > 0) return(u)
    if (file.exists(dest)) unlink(dest)
  }
  NA_character_
}

lock$url <- NA_character_
t0 <- Sys.time()
for (i in which(lock$kind != "GitHub")) {
  p <- lock$package[i]; v <- lock$version[i]
  tarball <- paste0(p, "_", v, ".tar.gz")
  dest <- file.path(contrib, tarball)
  if (file.exists(dest) && file.size(dest) > 0) { lock$url[i] <- "(already downloaded)"; next }
  urls <- if (lock$kind[i] == "CRAN") {
    c(file.path(cran_mirror, "src/contrib", tarball),
      file.path(cran_mirror, "src/contrib/Archive", p, tarball))
  } else {
    unlist(lapply(c("bioc", "data/annotation", "data/experiment", "workflows"), function(r) {
      base <- file.path(bioc_mirror, "packages", BIOC_VERSION, r, "src/contrib")
      c(file.path(base, tarball), file.path(base, "Archive", p, tarball))
    }))
  }
  lock$url[i] <- fetch(urls, dest)
  if (is.na(lock$url[i])) log_msg("NOT FOUND:", p, v, "(", lock$source[i], ")")
}
missing_src <- lock$package[lock$kind != "GitHub" & is.na(lock$url)]
log_msg(sprintf("downloads done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (length(missing_src)) stop("source tarballs not found: ", paste(missing_src, collapse = ", "), call. = FALSE)

# ----------------------------------------------------------------------------- install CRAN + Bioconductor
tools::write_PACKAGES(contrib, type = "source")
local_repo <- paste0("file://", normalizePath(repo_dir))
avail <- utils::available.packages(repos = local_repo, type = "source", filters = list())
pkgs <- lock$package[lock$kind != "GitHub"]
stopifnot(all(pkgs %in% rownames(avail)))

t1 <- Sys.time()
log_msg("installing", length(pkgs), "packages from the local repository with", ncpus, "parallel jobs")
withCallingHandlers(
  utils::install.packages(pkgs, lib = lib, repos = local_repo, type = "source", Ncpus = ncpus,
                          dependencies = c("Depends", "Imports", "LinkingTo"),
                          keep_outputs = log_dir, quiet = FALSE),
  warning = function(w) log_msg("WARNING:", conditionMessage(w))
)
log_msg(sprintf("install.packages() finished in %.1f min", as.numeric(difftime(Sys.time(), t1, units = "mins"))))

# ----------------------------------------------------------------------------- install GitHub packages
for (i in which(lock$kind == "GitHub")) {
  spec <- sub("^GitHub ", "", lock$source[i])              # user/repo@sha
  user <- sub("/.*$", "", spec)
  repo <- sub("^[^/]+/([^@]+)@.*$", "\\1", spec)
  sha  <- sub("^.*@", "", spec)
  tb <- file.path(work, paste0(repo, "_", sha, ".tar.gz"))
  lock$url[i] <- fetch(sprintf("https://codeload.github.com/%s/%s/tar.gz/%s", user, repo, sha), tb)
  if (is.na(lock$url[i])) { log_msg("NOT FOUND on GitHub:", spec); next }
  exdir <- file.path(work, "github", repo)
  unlink(exdir, recursive = TRUE)
  dir.create(exdir, recursive = TRUE)
  utils::untar(tb, exdir = exdir)
  pkgdir <- list.dirs(exdir, recursive = FALSE)[1]
  desc_file <- file.path(pkgdir, "DESCRIPTION")
  desc <- readLines(desc_file)
  desc <- desc[!grepl("^(Remote|Github)[A-Za-z0-9]*:", desc)]
  writeLines(c(desc,
               "RemoteType: github", "RemoteHost: api.github.com", paste0("RemoteRepo: ", repo),
               paste0("RemoteUsername: ", user), paste0("RemoteRef: ", sha), paste0("RemoteSha: ", sha),
               paste0("GithubRepo: ", repo), paste0("GithubUsername: ", user), paste0("GithubRef: ", sha),
               paste0("GithubSHA1: ", sha)), desc_file)
  log_msg("installing", lock$package[i], "from GitHub", spec)
  out <- file.path(log_dir, paste0(lock$package[i], ".github.out"))
  status <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-multiarch", paste0("--library=", shQuote(lib)), shQuote(pkgdir)),
                    stdout = out, stderr = out)
  if (status != 0) log_msg("FAILED:", lock$package[i], "- see", out)
}

# ----------------------------------------------------------------------------- check
ip <- utils::installed.packages(lib.loc = lib, fields = "RemoteSha", noCache = TRUE)
lock$installed_version <- ifelse(lock$package %in% rownames(ip), ip[match(lock$package, rownames(ip)), "Version"], NA)
lock$installed_sha <- ifelse(lock$package %in% rownames(ip), ip[match(lock$package, rownames(ip)), "RemoteSha"], NA)
lock$ok <- !is.na(lock$installed_version) & lock$installed_version == lock$version
gh <- lock$kind == "GitHub"
lock$ok[gh] <- lock$ok[gh] & lock$installed_sha[gh] == sub("^.*@", "", lock$source[gh])
lock$ok[is.na(lock$ok)] <- FALSE
check_file <- file.path(getwd(), paste0(ENV_NAME, "_lock_check.tsv"))
utils::write.table(lock[, c("package", "version", "source", "installed_version", "ok", "url")], check_file,
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "")
extra <- setdiff(rownames(ip), lock$package)
log_msg(sprintf("installed and matching the lock table: %d of %d; check table: %s", sum(lock$ok), nrow(lock), check_file))
if (length(extra)) log_msg("installed but not in the lock table:", paste(extra, collapse = ", "))
log_msg(sprintf("total time %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (!all(lock$ok)) {
  bad <- lock[!lock$ok, c("package", "version", "installed_version")]
  print(bad, row.names = FALSE)
  log_msg("build logs of the failed packages:", log_dir)
  quit(save = "no", status = 1)
}
log_msg("DONE: all", nrow(lock), "packages installed at the locked versions")
