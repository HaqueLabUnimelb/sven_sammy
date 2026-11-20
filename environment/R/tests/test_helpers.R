# =============================================================================
# test_helpers.R: functions shared by the test drivers in this folder.
#
#   copy_with_paths(src, dst, subs)   copy an archived script and change only path strings
#   run_lines(file, ranges)           run selected line ranges of a script (interactive scripts)
#   compare_csv(new, old)             byte equality, else cell-by-cell comparison of two CSV tables
#   write_session_info(file)          sessionInfo() and the library paths, as text
#
# Plain base R; no package is needed.
# =============================================================================

# subs: data.frame with columns line (integer, NA = every line), old, new (fixed strings).
# Every substitution must match at least once, otherwise the copy is not written.
copy_with_paths <- function(src, dst, subs) {
  x <- readLines(src, warn = FALSE, encoding = "UTF-8")
  orig <- x
  for (k in seq_len(nrow(subs))) {
    ln <- subs$line[k]
    idx <- if (is.na(ln)) seq_along(x) else ln
    hit <- idx[grepl(subs$old[k], x[idx], fixed = TRUE)]
    if (!length(hit)) stop(sprintf("substitution %d not found in %s (line %s): %s", k, src, ln, subs$old[k]))
    x[hit] <- gsub(subs$old[k], subs$new[k], x[hit], fixed = TRUE)
  }
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, dst, useBytes = TRUE)
  changed <- which(x != orig)
  cat(sprintf("copy_with_paths: %s -> %s; %d line(s) changed\n", src, dst, length(changed)))
  for (i in changed) cat(sprintf("  line %d\n    - %s\n    + %s\n", i, orig[i], x[i]))
  invisible(changed)
}

# ranges: character such as "1-18,64-104,108"; each range must be a complete set of expressions.
run_lines <- function(file, ranges, envir = globalenv()) {
  x <- readLines(file, warn = FALSE, encoding = "UTF-8")
  parts <- strsplit(strsplit(ranges, ",", fixed = TRUE)[[1]], "-", fixed = TRUE)
  for (p in parts) {
    a <- as.integer(p[1]); b <- as.integer(p[length(p)])
    cat(sprintf("\n######## %s lines %d-%d (%s)\n", basename(file), a, b, format(Sys.time(), "%H:%M:%S")))
    exprs <- parse(text = x[a:b], keep.source = TRUE, encoding = "UTF-8")
    source(exprs = exprs, local = envir, echo = TRUE, print.eval = TRUE, max.deparse.length = 500,
           spaced = FALSE)
  }
  invisible(TRUE)
}

# Compare two CSV files. Returns a one-row data.frame.
compare_csv <- function(new, old, label = basename(new), row.names = FALSE) {
  res <- data.frame(table = label, bytes_identical = NA, nrow_new = NA, nrow_old = NA, ncol_new = NA,
                    ncol_old = NA, same_columns = NA, text_cells_differing = NA, numeric_columns = NA,
                    max_abs_diff = NA, max_rel_diff = NA, status = NA, stringsAsFactors = FALSE)
  if (!file.exists(new)) { res$status <- "new file missing"; return(res) }
  if (!file.exists(old)) { res$status <- "archived file missing"; return(res) }
  res$bytes_identical <- unname(tools::md5sum(new) == tools::md5sum(old))
  a <- utils::read.csv(new, check.names = FALSE, stringsAsFactors = FALSE)
  b <- utils::read.csv(old, check.names = FALSE, stringsAsFactors = FALSE)
  res$nrow_new <- nrow(a); res$nrow_old <- nrow(b); res$ncol_new <- ncol(a); res$ncol_old <- ncol(b)
  res$same_columns <- identical(colnames(a), colnames(b))
  if (!res$same_columns || nrow(a) != nrow(b)) {
    res$status <- if (isTRUE(res$bytes_identical)) "identical" else "different shape"
    return(res)
  }
  ntext <- 0L; nnum <- 0L; mad <- 0; mrd <- 0
  for (j in seq_along(a)) {
    u <- a[[j]]; v <- b[[j]]
    if (is.numeric(u) && is.numeric(v)) {
      nnum <- nnum + 1L
      na_diff <- sum(is.na(u) != is.na(v))
      ntext <- ntext + na_diff
      ok <- !is.na(u) & !is.na(v)
      d <- abs(u[ok] - v[ok])
      if (length(d)) {
        mad <- max(mad, d)
        den <- pmax(abs(u[ok]), abs(v[ok]))
        mrd <- max(mrd, ifelse(den > 0, d / den, 0))
      }
    } else {
      ntext <- ntext + sum(as.character(u) != as.character(v) | is.na(u) != is.na(v), na.rm = TRUE)
    }
  }
  res$text_cells_differing <- ntext; res$numeric_columns <- nnum
  res$max_abs_diff <- mad; res$max_rel_diff <- mrd
  res$status <- if (isTRUE(res$bytes_identical)) "identical (bytes)" else
    if (ntext == 0 && mad == 0) "same values, different bytes" else
      if (ntext == 0) "numeric differences only" else "different"
  res
}

write_session_info <- function(file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(capture.output(sessionInfo()), "", ".libPaths():", paste0("  ", .libPaths())), file)
  invisible(file)
}
