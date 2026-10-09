# =============================================================================
# 13_sim.R — SIM: dengue deaths (underlying cause A90/A91/A97), descriptive
#
# Deaths at ages 10-14 are too rare for a municipal model (89 in the whole
# series); they appear only in the description (note to Table 1).
# Input: SIM-DO by state-year, 2017-2025. Output: data/interim/sim_<year>.rds and
# data/processed/sim_aggregate.rds. Resumable.
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })

VARS <- c("CODMUNRES", "DTOBITO", "IDADE", "CAUSABAS")

for (year in PARAM$year_start:PARAM$year_end) {
  fpath <- file.path(DIR$interim, sprintf("sim_%d.rds", year))
  if (file.exists(fpath)) next
  log_msg("SIM-DO ", year)
  d <- get_sim(year, vars = VARS, timeout = PARAM$timeout)
  if (is.null(d) || !nrow(d)) next
  d <- d[grepl(ICD_DENGUE, CAUSABAS)]
  d[, `:=`(code_muni = mun6(CODMUNRES), month = month_of(as_date(DTOBITO, "dmy")),
           agegroup = to_agegroup(age_sim(IDADE)))]
  saveRDS(d[!is.na(agegroup), .(deaths_sim = .N), by = .(code_muni, agegroup, month)], fpath)
  rm(d); gc()
}

sim <- rbindlist(lapply(list.files(DIR$interim, "^sim_\\d{4}\\.rds$", full.names = TRUE), readRDS))
saveRDS(sim, file.path(DIR$processed, "sim_aggregate.rds"))
log_msg("SIM: ", sum(sim$deaths_sim), " dengue deaths (ages 5-59, 2017-2025)")
