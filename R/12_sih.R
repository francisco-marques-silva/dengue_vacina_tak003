# =============================================================================
# 12_sih.R — SIH/SUS: dengue hospitalizations (PRIMARY outcome)
#
# Principal diagnosis A90, A91 or A97.x; excludes continuation AIHs (IDENT = 5)
# so the same hospitalization is not counted twice. Two aggregations in the same
# pass over the microdata:
#   by age group      -> data/interim/sih_<UF>_<year>.rds       -> sih_aggregate.rds
#   by single age     -> data/interim/sih_age_<UF>_<year>.rds -> sih_age_aggregate.rds
# (the second feeds the birth-cohort comparison, script 34)
# Input: SIH-RD by state-month, Jan/2017 to Mar/2026 (billing months up to 3 months
# after the end of the series capture late-recorded hospitalizations). Resumable.
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })

VARS <- c("MUNIC_RES", "DIAG_PRINC", "IDADE", "COD_IDADE", "DT_INTER", "MORTE", "IDENT")
year_final <- PARAM$sih_last_billmonth[["year"]]

for (uf in UFS) for (year in PARAM$year_start:year_final) {
  fpath_agegroup <- file.path(DIR$interim, sprintf("sih_%s_%d.rds", uf, year))
  fpath_age <- file.path(DIR$interim, sprintf("sih_age_%s_%d.rds", uf, year))
  if (file.exists(fpath_agegroup) && file.exists(fpath_age)) next
  month_end <- if (year == year_final) PARAM$sih_last_billmonth[["month"]] else 12L
  log_msg("SIH-RD ", uf, " ", year)
  d <- get_sih(uf, year, months = 1:month_end, vars = VARS, timeout = PARAM$timeout)
  if (is.null(d) || !nrow(d)) next

  d <- d[grepl(ICD_DENGUE, DIAG_PRINC) & !(as.character(IDENT) %in% "5")]
  d[, `:=`(code_muni = mun6(MUNIC_RES), month = month_of(as_date(DT_INTER, "ymd")),
           age = age_sih(IDADE, COD_IDADE))]
  d[, agegroup := to_agegroup(age)]
  saveRDS(d[!is.na(agegroup) & !is.na(month),
            .(admissions = .N, deaths_hosp = sum(as.character(MORTE) == "1", na.rm = TRUE)),
            by = .(code_muni, agegroup, month)], fpath_agegroup)
  saveRDS(d[!is.na(age) & age >= 5L & age <= 59L & !is.na(month),
            .(admissions = .N), by = .(code_muni, age, month)], fpath_age)
  rm(d); gc()
}

# ---- consolidate ------------------------------------------------------------
within <- function(x) x[month >= as.Date(sprintf("%d-01-01", PARAM$year_start)) &
                        month <= as.Date(sprintf("%d-12-01", PARAM$year_end))]
sih <- rbindlist(lapply(list.files(DIR$interim, "^sih_[A-Z]{2}_\\d{4}\\.rds$", full.names = TRUE),
                        readRDS), fill = TRUE)
sih <- within(sih[, .(admissions = sum(admissions), deaths_hosp = sum(deaths_hosp)),
                  by = .(code_muni, agegroup, month)])
saveRDS(sih, file.path(DIR$processed, "sih_aggregate.rds"))

sih_age <- rbindlist(lapply(list.files(DIR$interim, "^sih_age_[A-Z]{2}_\\d{4}\\.rds$", full.names = TRUE),
                              readRDS), fill = TRUE)
sih_age <- within(sih_age[, .(admissions = sum(admissions)), by = .(code_muni, age, month)])
saveRDS(sih_age, file.path(DIR$processed, "sih_age_aggregate.rds"))

log_msg("SIH consolidated: ", sum(sih$admissions), " admissions by age group; ",
        sum(sih_age$admissions), " by single year of age (should be equal)")

# =============================================================================
# Read from the parquet cache, without changing the aggregates above:
#   sih_birth_cohort_aggregate.rds  municipality × birth year × month × dengue hospitalizations
#   sih_controls_aggregate.rds    municipality × age group × month × hospitalizations for 4 negative controls
# Same filter as the block above: continuation AIHs out (IDENT = 5), month by
# admission date, Jan/2017-Dec/2025.
# Cohort: from the NASC field when the file has it. This project's cache does NOT have NASC
# (CACHE_COLUMNS never included it), and in that case the cohort is EXPECTED: a hospitalization at
# age a (completed years) on date d of year Y belongs to Y-a if the birthday has passed,
# which, with birthdays uniform over the year, has probability p = day of year of d /
# days in the year; the hospitalization enters with weight p in Y-a and 1-p in Y-a-1 (fractional counts).
# Resumable: if both files exist, the block does not re-read the cache.
# =============================================================================
FILE_COHORT <- file.path(DIR$processed, "sih_birth_cohort_aggregate.rds")
FILE_CTRL   <- file.path(DIR$processed, "sih_controls_aggregate.rds")
fs_sih <- list.files(file.path(DIR$cache, "SIH-RD"), "_full\\.parquet$", full.names = TRUE)
if (!(file.exists(FILE_COHORT) && file.exists(FILE_CTRL)) && length(fs_sih)) {
  log_msg("reading ", length(fs_sih), " SIH cache files for birth cohort and negative controls")
  CTRL <- c(ctrl_all = "all causes except A90/A91 and chapters O, P, Z", ctrl_a00a09 = "A00-A09",
            ctrl_j00j22 = "J00-J22", ctrl_a920 = "A92.0")
  ini <- as.Date(sprintf("%d-01-01", PARAM$year_start)); end <- as.Date(sprintf("%d-12-01", PARAM$year_end))
  parts <- lapply(fs_sih, function(f) {
    cols <- intersect(c("MUNIC_RES", "DIAG_PRINC", "IDENT", "DT_INTER", "IDADE", "COD_IDADE", "NASC"),
                      names(arrow::open_dataset(f)$schema))
    d <- as.data.table(arrow::read_parquet(f, col_select = tidyselect::all_of(cols)))
    d <- d[!(as.character(IDENT) %in% "5")]
    d[, `:=`(data = as_date(DT_INTER, "ymd"), dg = as.character(DIAG_PRINC))]
    d[, `:=`(month = month_of(data), code_muni = mun6(MUNIC_RES), age = age_sih(IDADE, COD_IDADE))]
    d <- d[!is.na(month) & month >= ini & month <= end]
    # negative controls, by panel age group
    d[, agegroup := to_agegroup(age)]
    ctrl <- d[!is.na(agegroup), .(ctrl_all = sum(!grepl(ICD_DENGUE, dg) & !grepl("^[OPZ]", dg)),
                               ctrl_a00a09 = sum(grepl("^A0[0-9]", dg)),
                               ctrl_j00j22 = sum(grepl("^J(0[0-9]|1[0-9]|2[0-2])", dg)),
                               ctrl_a920 = sum(grepl("^A920", dg))), by = .(code_muni, agegroup, month)]
    # dengue by birth cohort
    den <- d[grepl(ICD_DENGUE, dg)]
    has_birth <- "NASC" %in% names(den)
    if (has_birth) den[, birth := as_date(NASC, "ymd")] else den[, birth := as.Date(NA)]
    tally <- den[, .(dengue = .N, with_birth = sum(!is.na(birth)), without_birth_with_age = sum(is.na(birth) & !is.na(age)),
                    without_birth_without_age = sum(is.na(birth) & is.na(age)))]
    exact <- den[!is.na(birth), .(admissions = .N), by = .(code_muni, cohort = as.integer(format(birth, "%Y")), month)]
    expd <- den[is.na(birth) & !is.na(age)]
    expd[, `:=`(year = as.integer(format(data, "%Y")),
               p = as.integer(format(data, "%j")) / fifelse(as.integer(format(data, "%Y")) %% 4L == 0L, 366, 365))]
    expd <- rbind(expd[, .(code_muni, month, cohort = year - age, w = p)], expd[, .(code_muni, month, cohort = year - age - 1L, w = 1 - p)])
    expd <- expd[, .(admissions = sum(w)), by = .(code_muni, cohort, month)]
    list(ctrl = ctrl, cohort = rbind(exact[, method := "NASC"], expd[, method := "expected (age and admission date)"]),
         tally = tally)
  })
  ctrl <- rbindlist(lapply(parts, `[[`, "ctrl"))[, lapply(.SD, sum), by = .(code_muni, agegroup, month)]
  saveRDS(ctrl, FILE_CTRL)
  cohort <- rbindlist(lapply(parts, `[[`, "cohort"))[, .(admissions = sum(admissions)), by = .(code_muni, cohort, month, method)]
  saveRDS(cohort, FILE_COHORT)
  tally <- rbindlist(lapply(parts, `[[`, "tally"))[, lapply(.SD, sum)]
  fwrite(tally, file.path(DIR$processed, "sih_birth_cohort_count.csv"))
  log_msg("dengue ", tally$dengue, " admissions: ", tally$with_birth, " with NASC; ", tally$without_birth_with_age,
          " without NASC (cohort expected from age); ", tally$without_birth_without_age, " without NASC or age (left out of the aggregate)")
  log_msg("negative controls: ", paste(sprintf("%s = %s", names(CTRL), ctrl[, sapply(.SD, sum), .SDcols = names(CTRL)]),
                                                    collapse = "; "))
} else if (!(file.exists(FILE_COHORT) && file.exists(FILE_CTRL))) {
  log_msg("WARNING: no SIH cache; birth cohort and negative controls not generated")
}
