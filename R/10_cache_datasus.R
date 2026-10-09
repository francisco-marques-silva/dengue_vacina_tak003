# =============================================================================
# 10_cache_datasus.R — the raw microdata the project uses
#
# The project reads DATASUS microdata from data/cache_datasus/ (Parquet). This
# script only COPIES there the files still missing, from an already existing
# cache (CACHE_SOURCE in .Renviron, or the sibling folder of the old repository,
# detected automatically; CACHE_SOURCE=none disables copying). Whatever is not
# found will be DOWNLOADED from DATASUS by scripts 11-13 (hours; see README).
#
# Required files (11-13 request exactly these):
#   SINAN-DENGUE_<year>.parquet and SINAN-CHIKUNGUNYA_<year>.parquet, 2017-2026
#   SIM-DO_<UF>_<year>.parquet, 27 states × 2017-2025
#   SIH-RD_<UF>_<year>_<mm>_full.parquet, 27 states × Jan/2017-Mar/2026
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })

# ---- 1. list of what the project needs ---------------------------------------
years <- PARAM$year_start:PARAM$year_end
ult  <- PARAM$sih_last_billmonth
billmonths <- CJ(year = PARAM$year_start:ult[["year"]], month = 1:12)[year < ult[["year"]] | month <= ult[["month"]]]
required <- c(
  unlist(lapply(c(years, PARAM$year_end + 1L), function(a)
    c(cache_file("SINAN-DENGUE", a), cache_file("SINAN-CHIKUNGUNYA", a)))),
  unlist(lapply(years, function(a) vapply(UFS, function(u) cache_file("SIM-DO", a, u), ""))),
  unlist(lapply(UFS, function(u) billmonths[, cache_file("SIH-RD", year, u, month), by = .(year, month)]$V1)))
missing <- required[!file.exists(required)]
log_msg(length(required), " files needed; ", length(missing), " not yet in ", DIR$cache)

# ---- 2. where to copy from ---------------------------------------------------
origin <- Sys.getenv("CACHE_SOURCE")
if (tolower(origin) == "none") {
  origin <- ""                                   # rebuild from scratch: no copying, everything downloaded
  log_msg("CACHE_SOURCE=none: nothing copied; scripts 11-13 will download everything from DATASUS")
} else if (!nzchar(origin)) {
  # old repository as a sibling folder of this project, or this project inside it
  cand <- c(file.path(dirname(ROOT), "20260921 - Microdatasus", "dados_comuns", "cache_datasus"),
            file.path(dirname(ROOT), "dados_comuns", "cache_datasus"))
  cand <- cand[dir.exists(cand)]
  if (length(cand)) origin <- cand[1]
}

# ---- 3. copy whatever exists at the source -----------------------------------
if (length(missing) && nzchar(origin) && dir.exists(origin)) {
  rel <- substring(missing, nchar(DIR$cache) + 2L)            # "<system>/<file>"
  de  <- file.path(origin, rel)
  has <- file.exists(de)
  log_msg("Copying ", sum(has), " files from ", origin, " (", round(sum(file.size(de[has])) / 2^30, 2), " GB)")
  for (d in unique(dirname(missing[has]))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(de[has], missing[has], copy.date = TRUE)
  log_msg(sum(ok), " copied; ", sum(!has), " not in the source and will be downloaded by scripts 11-13")
} else if (length(missing) && tolower(Sys.getenv("CACHE_SOURCE")) != "none") {
  log_msg("No copy source (set CACHE_SOURCE); scripts 11-13 will download from DATASUS")
}

missing <- required[!file.exists(required)]
fwrite(data.table(fname = substring(required, nchar(DIR$cache) + 2L),
                  presente = file.exists(required)),
       file.path(DIR$logs, "cache_datasus_inventario.csv"))
log_msg("Cache: ", length(required) - length(missing), "/", length(required), " files present",
        if (length(missing)) paste0(" (", length(missing), " to download)") else "")
