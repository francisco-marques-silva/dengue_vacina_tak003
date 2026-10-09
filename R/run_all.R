# =============================================================================
# run_all.R — runs the whole project, in order, logging to logs/
#
#   source("R/run_all.R")                 # all: data -> analyses -> outputs
#   Sys.setenv(ANALYSIS_ONLY = "1"); source(...) # skips the data phase (10-16):
#                                             # uses data/interim and data/raw as they are
#   Sys.setenv(FROM = "31_validity"); source(...)  # resumes from a given script
#   Sys.setenv(CACHE_SOURCE = "none"); source(...)  # rebuild from scratch: 10 copies no
#       cache from anywhere and 11-15 download everything (DATASUS and SI-PNI); requires empty data/
#
# Phases and approximate time (ordinary laptop, 8 cores):
#   data      10-16   hours the first time without cache (DATASUS and SI-PNI download);
#                     minutes if data/cache_datasus/ already exists; seconds with
#                     data/interim/ ready (the scripts skip what is already done).
#                     The data scripts are resumable; if one fails (unstable DATASUS
#                     FTP), it is retried up to 3 times.
#   panel     20      ~2 min
#   analyses  30-35   ~2 h (the spatial permutation in 32 takes ~1 h with N_PERM=400)
#   outputs   40-44   ~5 min (42 computes the cache md5; 44 downloads the IBGE boundaries the first time)
#   check     99      seconds
# Each script runs in its own process (Rscript), so an error in one does not
# affect the following ones, and each one's log goes to logs/<script>.log.
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })

dat    <- c("10_cache_datasus", "11_sinan", "12_sih", "13_sim", "14_population", "15_sipni",
              "16_municipality_list")
analyses <- c("20_panel", "30_main_model", "31_validity", "32_spatial_permutation",
              "33_dose_response", "34_birth_cohort", "35_sensitivity")
outputs   <- c("40_tables", "41_figures", "42_manifest", "43_manuscript_numbers", "44_maps", "99_check_reproduction")
scripts <- c(if (Sys.getenv("ANALYSIS_ONLY") != "1") dat, analyses, outputs)
# FROM=<script name> resumes from it (e.g. Sys.setenv(FROM = "16_municipality_list"))
since <- match(Sys.getenv("FROM"), scripts)
if (!is.na(since)) scripts <- scripts[since:length(scripts)]

# ---- initial state: what already exists decides how much gets downloaded ----
n_cache <- length(list.files(DIR$cache, "\\.parquet$", recursive = TRUE))
n_int   <- length(list.files(DIR$interim, "\\.rds$"))
n_red   <- length(list.files(file.path(DIR$raw, "sipni", "reduced"), "\\.rds$"))
log_msg("Initial state: DATASUS cache ", n_cache, " files; interim ", n_int, " files; reduced SI-PNI ",
        n_red, " months -> ",
        if (n_int > 0 && n_red > 0) "reproduction from the existing aggregates (downloads only what is missing)"
        else if (n_cache > 0) "aggregation from the existing cache"
        else "full rebuild: download from DATASUS and SI-PNI (hours)")

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
skip_before <- Sys.getenv("SKIP_INSTALL", unset = NA)
Sys.setenv(SKIP_INSTALL = "1")            # packages were already checked here; child processes inherit
summ <- data.table(script = scripts, start = as.POSIXct(NA), end = as.POSIXct(NA), status = NA_character_,
                     attempts = NA_integer_)
for (i in seq_along(scripts)) {
  s <- scripts[i]
  fpath <- file.path(ROOT, "R", paste0(s, ".R")); fpath_log <- file.path(DIR$logs, paste0(s, ".log"))
  summ[i, start := Sys.time()]
  log_msg("== ", s)
  max_tent <- if (s %in% dat) 3L else 1L          # the data scripts are resumable
  for (tent in seq_len(max_tent)) {
    st <- system2(rscript, shQuote(fpath), stdout = fpath_log, stderr = fpath_log)
    if (st == 0) break
    if (tent < max_tent) { log_msg("   failed (attempt ", tent, "); retrying in 60 s"); Sys.sleep(60) }
  }
  summ[i, `:=`(end = Sys.time(), status = if (st == 0) "ok" else paste("ERROR", st), attempts = tent)]
  if (st != 0) {
    log_msg("ERROR in ", s, " — see ", fpath_log, ". Stopping here.")
    break
  }
}
if (is.na(skip_before)) Sys.unsetenv("SKIP_INSTALL") else Sys.setenv(SKIP_INSTALL = skip_before)
summ[, minutes := round(as.numeric(difftime(end, start, units = "mins")), 1)]
print(summ[!is.na(status)])
fwrite(summ, file.path(DIR$logs, "run_all_summ.csv"))
