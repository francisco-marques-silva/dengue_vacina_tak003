# =============================================================================
# 42_manifest.R — record of what produced the results
#
#   outputs/data_manifest.csv  md5, size and date of each data file
#                                (DATASUS cache, raw, interim, processed):
#                                this is what allows saying, months later, that
#                                the analysis ran on exactly these files
#   outputs/R_session.txt         version of R, of the system and of each package
#   renv.lock                    package versions (renv::snapshot from
#                                the library in use; whoever reproduces runs
#                                renv::restore() and gets the same versions)
# MANIFEST_NO_HASH=1 skips the md5 (the cache is ~2 GB; takes a few minutes).
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })

# ---- 1. data manifest --------------------------------------------------------
fpaths <- c(list.files(DIR$cache, recursive = TRUE, full.names = TRUE, pattern = "\\.parquet$"),
          list.files(DIR$raw, recursive = TRUE, full.names = TRUE),
          list.files(DIR$interim, full.names = TRUE), list.files(DIR$processed, full.names = TRUE))
fpaths <- fpaths[!grepl("/sipni/tmp/", fpaths)]
m <- data.table(path = substring(gsub("\\\\", "/", fpaths), nchar(ROOT) + 2L),
                bytes = file.size(fpaths), modified_em = format(file.mtime(fpaths), "%Y-%m-%d %H:%M:%S"))
m[, folder := sub("^data/([^/]+)/.*$", "\\1", path)]
if (Sys.getenv("MANIFEST_NO_HASH") != "1") {
  log_msg("Computing md5 of ", nrow(m), " files (", round(sum(m$bytes) / 2^30, 2), " GB)...")
  m[, md5 := unname(tools::md5sum(fpaths))]
} else m[, md5 := NA_character_]
setorder(m, path)
fwrite(m, file.path(ROOT, "outputs", "data_manifest.csv"))
print(m[, .(fnames = .N, GB = round(sum(bytes) / 2^30, 2)), by = folder])

# SIH months present and SINAN/SIM years (what is in the cache)
cache <- m[folder == "cache_datasus"]
cache[, system := sub("^data/cache_datasus/([^/]+)/.*$", "\\1", path)]
cache[, year := as.integer(sub(".*?_(\\d{4})(_\\d{2})?(_\\w+)?\\.parquet$", "\\1", path))]
print(cache[, .(fnames = .N, years = paste(range(year), collapse = "-")), by = system])

# ---- 2. R session ---------------------------------------------------------------
pk <- c("data.table", "fixest", "ggplot2", "patchwork", "arrow", "curl", "jsonlite", "stringi",
        "flextable", "officer", "sf", "geobr", "microdatasus", "brpop")
vers <- vapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "missing"), "")
writeLines(c(R.version.string, paste("Platform:", R.version$platform), paste("SO:", Sys.info()[["sysname"]], Sys.info()[["release"]]),
             paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "", "Packages:",
             sprintf("  %-14s %s", pk, vers), "", capture.output(sessionInfo())),
           file.path(ROOT, "outputs", "R_session.txt"))

# ---- 3. renv.lock ----------------------------------------------------------------
# DESCRIPTION lists the dependencies; the "explicit" snapshot records only them (and
# what they require), from the library in use — without needing renv::init()
if (requireNamespace("renv", quietly = TRUE)) {
  ok <- tryCatch({
    renv::snapshot(project = ROOT, library = .libPaths(), type = "explicit", prompt = FALSE)
    TRUE }, error = function(e) { log_msg("renv::snapshot failed: ", conditionMessage(e)); FALSE })
  if (ok) log_msg("renv.lock updated")
} else log_msg("renv missing: renv.lock not generated")
log_msg("42 done")
