# =============================================================================
# 11_sinan.R — SINAN: dengue (secondary outcomes) and chikungunya (negative
# control), aggregated by municipality of residence × age group × month of
# symptom onset
#
# Input:  data/cache_datasus/SINAN-*/SINAN-*_<year>.parquet (downloaded if missing)
# Output: data/interim/sinan_<disease>_<year>.rds (one per year; resumable)
#          data/processed/sinan_aggregate.rds
# Probable case = notified and not discarded (CLASSI_FIN != 5); for dengue,
# also excludes reclassification as chikungunya (13). One extra year (2026)
# is included to capture notifications of December 2025 symptoms.
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })

VARS <- c("ID_MN_RESI", "DT_SIN_PRI", "DT_NOTIFIC", "NU_IDADE_N", "CLASSI_FIN", "CRITERIO",
          "HOSPITALIZ", "EVOLUCAO")
SYSTEMS <- c(dengue = "SINAN-DENGUE", chik = "SINAN-CHIKUNGUNYA")

for (ag in names(SYSTEMS)) for (year in PARAM$year_start:(PARAM$year_end + 1L)) {
  fpath <- file.path(DIR$interim, sprintf("sinan_%s_%d.rds", ag, year))
  if (file.exists(fpath)) next
  log_msg(SYSTEMS[[ag]], " ", year)
  d <- get_datasus(SYSTEMS[[ag]], year, vars = VARS, timeout = PARAM$timeout)
  if (is.null(d) || !nrow(d)) { log_msg("no data for ", year); next }

  for (v in c("CLASSI_FIN", "CRITERIO", "HOSPITALIZ", "EVOLUCAO"))
    set(d, j = v, value = trimws(as.character(d[[v]])))
  d[, onset_date := as_date(DT_SIN_PRI)]
  d[is.na(onset_date), onset_date := as_date(DT_NOTIFIC)]
  d[, `:=`(code_muni = mun6(ID_MN_RESI), month = month_of(onset_date),
           agegroup = to_agegroup(age_sinan(NU_IDADE_N)))]
  d <- d[!is.na(agegroup) & !is.na(month) & nchar(code_muni) == 6 & !grepl("0000$", code_muni)]
  d[, probable := is.na(CLASSI_FIN) | CLASSI_FIN != "5"]
  if (ag == "dengue") d[CLASSI_FIN == "13", probable := FALSE]

  agg <- d[probable == TRUE, .(
    cases         = .N,
    cases_lab     = sum(CRITERIO == "1" & CLASSI_FIN %in% c("10", "11", "12"), na.rm = TRUE),
    cases_severe  = sum(CLASSI_FIN %in% c("11", "12", "2", "3", "4"), na.rm = TRUE),
    cases_hosp    = sum(HOSPITALIZ == "1", na.rm = TRUE),
    deaths_disease = sum(EVOLUCAO == "2", na.rm = TRUE)), by = .(code_muni, agegroup, month)]
  agg[, disease := ag]
  saveRDS(agg, fpath)
  log_msg("ok: ", nrow(d), " notifications -> ", nrow(agg), " cells")
  rm(d, agg); gc()
}

fpaths  <- list.files(DIR$interim, "^sinan_(dengue|chik)_\\d{4}\\.rds$", full.names = TRUE)
sinan <- rbindlist(lapply(fpaths, readRDS))
sinan <- sinan[, lapply(.SD, sum), by = .(disease, code_muni, agegroup, month),
               .SDcols = c("cases", "cases_lab", "cases_severe", "cases_hosp", "deaths_disease")]
sinan <- sinan[month >= as.Date(sprintf("%d-01-01", PARAM$year_start)) &
               month <= as.Date(sprintf("%d-12-01", PARAM$year_end))]
saveRDS(sinan, file.path(DIR$processed, "sinan_aggregate.rds"))
log_msg("SINAN consolidated: ", nrow(sinan), " rows; ", sinan[disease == "dengue", sum(cases)],
        " probable dengue cases")
