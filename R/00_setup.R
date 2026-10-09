# =============================================================================
# 00_setup.R — packages, parameters and project folders
#
# Every script starts by loading this file. It:
#   1. locates the project root (the folder containing dengue_vaccine_tak003.Rproj),
#      so any script runs from any working directory;
#   2. installs whatever is missing (skip with SKIP_INSTALL=1);
#   3. defines PARAM (study parameters) and DIR (folders);
#   4. loads R/01_functions.R.
#
# Recognized environment variables (all optional; set them in .Renviron):
#   SKIP_INSTALL=1       does not check/install packages
#   DATASUS_UPDATE=1  allows re-downloading recent DATASUS files (by
#                        default the cache is FROZEN: a file that exists
#                        is never replaced — this guarantees the same result)
#   CACHE_SOURCE=<dir>   folder of an old cache from which 10_cache_datasus.R
#                        copies the files instead of downloading them
#   N_PERM, N_PERM_DOSE  number of permutation draws (default 400 and 300)
# =============================================================================

ROOT <- local({
  mark <- "dengue_vacina_tak003.Rproj"
  cand <- character()
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) cand <- c(cand, dirname(normalizePath(f[1], mustWork = FALSE)))
  for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) {
    cand <- c(cand, dirname(normalizePath(fr$ofile, mustWork = FALSE))); break
  }
  cand <- c(cand, getwd())
  for (d in cand) {
    while (nzchar(d) && d != dirname(d)) {
      if (file.exists(file.path(d, mark))) return(normalizePath(d, winslash = "/"))
      d <- dirname(d)
    }
  }
  stop("Project root not found (looked for '", mark, "')")
})

# ---- packages ---------------------------------------------------------------
if (Sys.getenv("SKIP_INSTALL") != "1") {
  # in non-interactive Rscript CRAN has no mirror ("@CRAN@") and install.packages fails
  if (is.null(getOption("repos")) || identical(unname(getOption("repos")["CRAN"]), "@CRAN@"))
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  cran <- c("data.table", "fixest", "ggplot2", "patchwork", "arrow", "curl", "jsonlite",
            "stringi", "flextable", "officer", "sf", "geobr", "remotes", "renv")
  missing <- setdiff(cran, rownames(installed.packages()))
  if (length(missing)) install.packages(missing)
  # microdatasus >= 3.0 (reads DBC without the read.dbc package) and brpop (populations)
  if (!requireNamespace("microdatasus", quietly = TRUE) ||
      utils::packageVersion("microdatasus") < "3.0.0")
    remotes::install_github("rfsaldanha/microdatasus")
  if (!requireNamespace("brpop", quietly = TRUE)) remotes::install_github("rfsaldanha/brpop")
}
suppressPackageStartupMessages(library(data.table))

# ---- study parameters -------------------------------------------------------
PARAM <- list(
  year_start    = 2017L,                     # 7 years without vaccine before 2024
  year_end       = 2025L,                     # last year of hospitalization/symptoms
  vaccine_start = as.Date("2024-02-01"),     # first phase of the strategy (521 municipalities)
  # SIH: billing months up to 3 months after the end (late-recorded hospitalizations)
  sih_last_billmonth = c(year = 2026L, month = 3L),
  # threshold of cumulative D1 coverage in the target age group defining exposure
  # ALTERNATIVE used only in the sensitivity analysis (treat_cov); high on purpose: at 1%,
  # 2,581 municipalities outside the official list would become "adopters" (residents
  # vaccinated outside the municipality) and the control group would vanish
  threshold_coverage = 0.20,
  timeout = 600
)

# ---- folders ----------------------------------------------------------------
DIR <- list(
  cache     = file.path(ROOT, "data", "cache_datasus"),   # raw microdata (Parquet)
  raw       = file.path(ROOT, "data", "raw"),             # manual/downloaded inputs (lists, SI-PNI)
  interim   = file.path(ROOT, "data", "interim"),         # aggregates by state-year (the frozen snapshot)
  processed = file.path(ROOT, "data", "processed"),       # panels read by the analysis
  mod       = file.path(ROOT, "outputs", "models"),      # .rds models and intermediate tables
  tab       = file.path(ROOT, "outputs", "tables"),      # final tables (docx + csv)
  fig       = file.path(ROOT, "outputs", "figures"),      # final figures
  logs      = file.path(ROOT, "logs")
)
# fixed row order of the group and region tables
GROUP_ORDER  <- c("first phase", "expanded", "never included")
REGION_ORDER <- c("Central-West", "Northeast", "North", "Southeast", "South")

invisible(lapply(c(DIR, file.path(DIR$raw, c("sipni", "lists")),
                   DIR$fig),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

options(timeout = max(PARAM$timeout, getOption("timeout")), encoding = "UTF-8")
# scripts and outputs are UTF-8 (accents in tables and figures)
if (!isTRUE(l10n_info()[["UTF-8"]])) suppressWarnings(try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE))
if (requireNamespace("fixest", quietly = TRUE))
  fixest::setFixest_nthreads(max(1, parallel::detectCores() - 1))

source(file.path(ROOT, "R", "01_functions.R"))
