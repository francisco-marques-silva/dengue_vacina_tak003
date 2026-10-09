# =============================================================================
# 20_panel.R — analysis panels
#
# panel.rds        municipality × age group × month, Jan/2017-Dec/2025 (5,570
#                   municipalities × 7 age groups × 108 months), with outcomes, population,
#                   vaccination coverage and the design variables
# panel_cohort.rds municipality × birth cohort (2000-2018) × month, with
#                   hospitalizations only (no denominator), for script 34
#
# Main-model exposure = OFFICIAL LIST (intention to treat): the 521
# first-phase municipalities are treated from Feb/2024; the 2,230 expansion ones
# get expn = 1 and are left OUT of the main model (neither first-phase treated nor
# clean controls); the remaining 2,819 are the controls. SI-PNI coverage
# does not define treatment (it records the vaccinee's municipality of RESIDENCE, and
# municipalities outside the list show doses); it enters the dose-response and,
# with a high threshold, the sensitivity analysis (treat_cov).
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })

pop   <- readRDS(file.path(DIR$processed, "population.rds"))
sih   <- readRDS(file.path(DIR$processed, "sih_aggregate.rds"))
sinan <- readRDS(file.path(DIR$processed, "sinan_aggregate.rds"))
sim   <- readRDS(file.path(DIR$processed, "sim_aggregate.rds"))
vac   <- readRDS(file.path(DIR$processed, "dengue_vaccination.rds"))
lst <- fread(file.path(DIR$raw, "lists", "dengue_vaccine_municipalities_full.csv"), colClasses = "character")

months <- seq(as.Date(sprintf("%d-01-01", PARAM$year_start)), as.Date(sprintf("%d-12-01", PARAM$year_end)), by = "month")
munis <- sort(unique(pop$code_muni)); munis <- munis[!grepl("0000$", munis)]
NEVER <- 10000L                                   # "g" of never-treated municipalities
t_de  <- function(d) (as.integer(format(d, "%Y")) - PARAM$year_start) * 12L + as.integer(format(d, "%m"))

# ---- 1. skeleton + population + outcomes ------------------------------------
panel <- CJ(code_muni = munis, agegroup = AGEGROUP_LABELS, month = months)
panel[, year := as.integer(format(month, "%Y"))]
panel <- merge(panel, pop, by = c("code_muni", "year", "agegroup"), all.x = TRUE)

den <- sinan[disease == "dengue", .(code_muni, agegroup, month, den_cases = cases, den_lab = cases_lab,
                                   den_severe = cases_severe, den_hosp = cases_hosp, den_deaths = deaths_disease)]
chk <- sinan[disease == "chik", .(code_muni, agegroup, month, chik_cases = cases)]
for (d in list(sih, den, chk, sim)) {
  panel <- merge(panel, d, by = c("code_muni", "agegroup", "month"), all.x = TRUE)
  for (v in setdiff(names(d), c("code_muni", "agegroup", "month"))) set(panel, which(is.na(panel[[v]])), v, 0L)
}

# ---- 2. cumulative D1/D2 coverage in the target age group, lagged 1 month ---
# fixed denominator = 2024 population aged 10-14; truncated at 100%
cov <- vac[agegroup == TARGET_AGEGROUP & month %in% months, .(d1 = sum(d1), d2 = sum(d2)), by = .(code_muni, month)]
cov <- merge(CJ(code_muni = munis, month = months), cov, by = c("code_muni", "month"), all.x = TRUE)
cov[is.na(d1), d1 := 0L][is.na(d2), d2 := 0L]
cov <- merge(cov, pop[year == 2024L & agegroup == TARGET_AGEGROUP, .(code_muni, pop_target = pop)], by = "code_muni", all.x = TRUE)
setorder(cov, code_muni, month)
cov[, `:=`(cov_d1 = pmin(cumsum(d1) / pop_target, 1), cov_d2 = pmin(cumsum(d2) / pop_target, 1)), by = code_muni]
cov[, `:=`(cov_d1_l1 = shift(cov_d1, 1, fill = 0), cov_d2_l1 = shift(cov_d2, 1, fill = 0)), by = code_muni]
panel <- merge(panel, cov[, .(code_muni, month, cov_d1_l1, cov_d2_l1)], by = c("code_muni", "month"), all.x = TRUE)
panel[is.na(cov_d1_l1), `:=`(cov_d1_l1 = 0, cov_d2_l1 = 0)]

# ---- 3. exposure ----------------------------------------------------------------
adoption <- data.table(code_muni = munis)
adoption <- merge(adoption, lst[phase == "1", .(g = t_de(month_of(min(as.Date(start_date))))), by = .(code_muni = mun6(code_muni))],
                by = "code_muni", all.x = TRUE)
adoption <- merge(adoption, lst[phase != "1", .(g_expn = t_de(month_of(min(as.Date(start_date))))), by = .(code_muni = mun6(code_muni))],
                by = "code_muni", all.x = TRUE)
adoption <- merge(adoption, cov[cov_d1 >= PARAM$threshold_coverage, .(g_cov = t_de(min(month))), by = code_muni],
                by = "code_muni", all.x = TRUE)
for (v in c("g", "g_expn", "g_cov")) set(adoption, which(is.na(adoption[[v]])), v, NEVER)
panel <- merge(panel, adoption, by = "code_muni", all.x = TRUE)

panel[, `:=`(
  t      = t_de(month),
  target   = as.integer(agegroup == TARGET_AGEGROUP),
  expn    = as.integer(g_expn < NEVER),
  uf     = substr(code_muni, 1, 2),
  region = c("1" = "North", "2" = "Northeast", "3" = "Southeast", "4" = "South", "5" = "Central-West")[substr(code_muni, 1, 1)])]
panel[, cohort := fifelse(target == 1L, g, NEVER)]                    # adoption cohort (event study)
panel[, treat     := as.integer(target == 1L & t >= g)]                # main model
panel[, treat_expn := as.integer(target == 1L & (t >= g | t >= g_expn))] # expansion as treated (sensitivity)
panel[, treat_cov := as.integer(target == 1L & t >= g_cov)]            # exposure by coverage (sensitivity)
panel[, window_contamination := month >= as.Date("2025-02-01") & month <= as.Date("2025-04-01")]  # doses outside the age group
panel <- panel[!is.na(pop) & pop > 0]
saveRDS(panel, file.path(DIR$processed, "panel.rds"))
log_msg("panel.rds: ", nrow(panel), " rows; ", uniqueN(panel$code_muni), " municipalities; ",
        panel[g < NEVER, uniqueN(code_muni)], " in the first phase; ", panel[expn == 1L, uniqueN(code_muni)], " in the expansion")

# ---- 4. birth-cohort panel (main-model municipalities) ---------------------
sih_age <- readRDS(file.path(DIR$processed, "sih_age_aggregate.rds"))
sih_age[, cohort := as.integer(format(month, "%Y")) - age]         # approximate birth year
COHORTS <- 2000:2018
mun <- unique(panel[, .(code_muni, g, expn, uf, region)])
pc <- CJ(code_muni = mun[expn == 0L, code_muni], cohort = COHORTS, month = months)
pc <- merge(pc, sih_age[cohort %in% COHORTS, .(admissions = sum(admissions)), by = .(code_muni, cohort, month)],
            by = c("code_muni", "cohort", "month"), all.x = TRUE)
pc[is.na(admissions), admissions := 0L]
pc <- merge(pc, mun[, .(code_muni, g, uf, region)], by = "code_muni")
pc[, `:=`(year = as.integer(format(month, "%Y")), t = t_de(month), treated = as.integer(g < NEVER))]
# exposure: cohorts 2010-2014 from Feb/2024; 2015 from 2025 (enters the target age group)
pc[, exposed := as.integer((cohort %in% 2010:2014 & t >= t_de(PARAM$vaccine_start)) | (cohort == 2015L & year >= 2025L))]
pc[, `:=`(treat = treated * exposed, cohort_exposed = as.integer(cohort %in% 2010:2015))]
saveRDS(pc, file.path(DIR$processed, "panel_cohort.rds"))
log_msg("panel_cohort.rds: ", nrow(pc), " rows; ", uniqueN(pc$code_muni), " municipalities")
