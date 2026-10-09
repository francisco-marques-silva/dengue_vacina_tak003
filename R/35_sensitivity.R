# =============================================================================
# 35_sensitivity.R — negative controls, robustness and specifications
#
# 1. Sensitivity of the pooled effect (Table S1): negative control
#    (chikungunya), alternative comparison groups, exclusion of Feb-Apr/2025
#    (near-expiry doses given outside the age group), exposure by
#    SI-PNI coverage, temporal placebo (adoption 24 months earlier), expansion
#    as treated / as control, negative binomial, municipalities >= 100k,
#    aggregation by health region, by macro-region.
# 2. Additional regional fixed effects (age×UF×year; age×region×month; both).
# 3. Serotypes in SINAN: typing completeness and national composition by year.
# 4. Placebos by season (Feb-Jan), each season estimated against the three
#    previous ones — the earlier version of the validity analysis (Table S2), kept
#    to show why the 2019/20 placebo looked significant.
# 5. Contamination of the fixed age groups: doses at ages 15-19, ageing of the
#    vaccinated cohorts and dilution of the target age group in 2025 (cited in the Discussion).
# Outputs (outputs/models/): table3_sensitivity.csv,
#   table_specifications_regional.csv, table_serotype_completeness_year.csv,
#   table_serotype_composition_national.csv, table_placebos_season.csv,
#   table_placebos_permutation.csv, table_diagnostic_cohort.csv
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(fixest)

NEVER <- 10000L
panel_full <- readRDS(file.path(DIR$processed, "panel.rds"))
panel <- panel_full[expn == 0L]                                       # main model sample
fml <- function(rhs, y = "admissions", fe = FE_MAIN) as.formula(sprintf("%s ~ %s | %s", y, rhs, fe))

# ---- 1. sensitivity of the pooled effect -----------------------------------------
sens <- list()
SPEC <- list(
  list(lbl = "Negative control: chikungunya (SINAN)", d = quote(panel), y = "chik_cases", rhs = "treat"),
  list(lbl = "Comparison with adults 20-59 only", d = quote(panel[agegroup %in% c("10-14", "20-29", "30-39", "40-49", "50-59")]), y = "admissions", rhs = "treat"),
  list(lbl = "Comparison with 5-9 and 15-19 only", d = quote(panel[agegroup %in% c("05-09", "10-14", "15-19")]), y = "admissions", rhs = "treat"),
  list(lbl = "Excluding Feb-Apr/2025", d = quote(panel[window_contamination == FALSE]), y = "admissions", rhs = "treat"),
  list(lbl = sprintf("Adoption by D1 coverage >= %.0f%% (SI-PNI)", 100 * PARAM$threshold_coverage), d = quote(panel), y = "admissions", rhs = "treat_cov"),
  list(lbl = "Placebo: adoption 24 months earlier", d = quote(pre), y = "admissions", rhs = "treat_placebo"),
  list(lbl = "Including expanded municipalities (treated at entry)", d = quote(panel_full), y = "admissions", rhs = "treat_expn"),
  list(lbl = "Expanded municipalities kept as controls (conservative bias)", d = quote(panel_full), y = "admissions", rhs = "treat"),
  list(lbl = "Negative binomial", d = quote(panel), y = "admissions", rhs = "treat", fun = fenegbin),
  list(lbl = "Municipalities ≥ 100 thousand inhab.", d = quote(panel[code_muni %in% large]), y = "admissions", rhs = "treat"))
pre <- panel[month < PARAM$vaccine_start]
pre[, treat_placebo := as.integer(target == 1L & g < NEVER & t >= (g - 24L))]
large <- panel[year == 2022L, .(pt = sum(pop)), by = code_muni][pt >= 1e5, code_muni]
for (e in SPEC) {
  fun <- if (is.null(e$fun)) fepois else e$fun
  m <- tryCatch(fun(fml(e$rhs, e$y), data = eval(e$d), offset = ~log(pop), cluster = ~code_muni),
                error = function(err) { warning(e$lbl, ": ", conditionMessage(err)); NULL })
  sens[[e$lbl]] <- if (is.null(m)) data.table(analysis = e$lbl) else rr_row(m, e$rhs, e$lbl)
}

# aggregation by health region (region "treated" if >= 50% of the target population is in treated municipalities)
pr <- merge(panel, health_region(), by = "code_muni")
ag <- pr[, .(admissions = sum(admissions), pop = sum(pop), frac_treat = sum(treat * pop) / sum(pop)),
         by = .(health_reg, agegroup, t)]
ag[, treat := as.integer(frac_treat >= 0.5)]
m_rs <- fepois(admissions ~ treat | health_reg^t + health_reg^agegroup + agegroup^t, data = ag, offset = ~log(pop), cluster = ~health_reg)
sens$healthreg <- rr_row(m_rs, "treat", "Aggregated by health region (449)")
# main model standard error clustered by health region (computed in 31)
cl <- fread(file.path(DIR$mod, "table_cluster_alternative.csv"))[grepl("^health region", clustering)]
sens$cluster_rs <- data.table(analysis = "Main model, SE clustered by health region", RR = cl$RRR,
                              CI95_lwr = cl$CI95_lwr, CI95_upr = cl$CI95_upr, p = cl$p, se_log = cl$se_log, n_obs = cl$n_obs)
# by macro-region
for (r in intersect(REGION_ORDER, unique(panel$region))) {
  m <- fepois(fml("treat"), data = panel[region == r], offset = ~log(pop), cluster = ~code_muni)
  sens[[paste("Region", r)]] <- rr_row(m, "treat", paste("Region", r))
}
sens <- rbindlist(sens, fill = TRUE)
setnames(sens, "RR", "RRR")
fwrite(sens[, .(analysis, RRR, CI95_lwr, CI95_upr, p, n_obs)], file.path(DIR$mod, "table3_sensitivity.csv"))
cat("\n== sensitivity\n"); print(sens[, .(analysis, RRR = round(RRR, 3), CI95_lwr = round(CI95_lwr, 3), CI95_upr = round(CI95_upr, 3))])

# ---- 2. additional regional fixed effects -----------------------------------------
FE_REG <- list("base (main model)" = FE_MAIN,
               "+ agegroup^uf^year" = paste(FE_MAIN, "+ agegroup^uf^year"),
               "+ agegroup^region^month" = paste(FE_MAIN, "+ agegroup^region^t"),
               "+ agegroup^uf^year and agegroup^region^month" = paste(FE_MAIN, "+ agegroup^uf^year + agegroup^region^t"))
expd <- rbindlist(lapply(names(FE_REG), function(nm) {
  m <- fepois(fml("treat", fe = FE_REG[[nm]]), data = panel, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
  rr_row(m, "treat", nm)
}))
setnames(expd, c("analysis", "RR"), c("specification", "RRR"))
fwrite(expd[, .(specification, RRR, CI95_lwr, CI95_upr, p, n_obs)], file.path(DIR$mod, "table_specifications_regional.csv"))
cat("\n== regional specifications\n"); print(expd[, .(specification, RRR = round(RRR, 3), CI95_lwr = round(CI95_lwr, 3), CI95_upr = round(CI95_upr, 3))])

# ---- 3. serotypes in SINAN --------------------------------------------------------
FILE <- file.path(DIR$interim, "serotype_uf_month.rds")
if (!file.exists(FILE)) {
  parts <- list()
  for (a in PARAM$year_start:PARAM$year_end) {
    d <- get_datasus("SINAN-DENGUE", a, vars = c("ID_MN_RESI", "DT_SIN_PRI", "DT_NOTIFIC", "SOROTIPO"), timeout = PARAM$timeout)
    if (is.null(d) || !nrow(d)) next
    d[, onset_date := as_date(DT_SIN_PRI)][is.na(onset_date), onset_date := as_date(DT_NOTIFIC)]
    d[, `:=`(uf = substr(mun6(ID_MN_RESI), 1, 2), month = month_of(onset_date), serotype = trimws(as.character(SOROTIPO)))]
    d <- d[grepl("^[1-5][0-9]$", uf) & !is.na(month)]
    parts[[as.character(a)]] <- d[, .(cases = .N, typed = sum(serotype %in% c("1", "2", "3", "4")),
                                       denv1 = sum(serotype == "1", na.rm = TRUE), denv2 = sum(serotype == "2", na.rm = TRUE),
                                       denv3 = sum(serotype == "3", na.rm = TRUE), denv4 = sum(serotype == "4", na.rm = TRUE)),
                                   by = .(uf, month)]
    log_msg("SEROTYPE ", a, ": ", nrow(d), " notifications"); rm(d); gc()
  }
  sor <- rbindlist(parts)[, lapply(.SD, sum), by = .(uf, month), .SDcols = c("cases", "typed", paste0("denv", 1:4))]
  saveRDS(sor, FILE)
}
sor <- readRDS(FILE)
sor <- sor[month >= as.Date(sprintf("%d-01-01", PARAM$year_start)) & month <= as.Date(sprintf("%d-12-01", PARAM$year_end))]
sor[, year := as.integer(format(month, "%Y"))]
cmpr_uf <- sor[, .(cases = sum(cases), typed = sum(typed)), by = .(year, uf)]
cmpr_year <- cmpr_uf[, .(ufs = .N, ufs_with_30_typed = sum(typed >= 30), cases = sum(cases), typed = sum(typed),
                        pct_typed = round(100 * sum(typed) / sum(cases), 2), median_typed_by_uf = as.integer(median(typed))), by = year]
cmpr_year <- merge(cmpr_year, sor[, .(ufmes = .N, ufmes_with_min = sum(typed >= 20L)), by = year], by = "year")
fwrite(cmpr_year, file.path(DIR$mod, "table_serotype_completeness_year.csv"))
cmpr_nat <- sor[, .(denv1 = sum(denv1), denv2 = sum(denv2), denv3 = sum(denv3), denv4 = sum(denv4), typed = sum(typed)), by = year]
for (v in paste0("denv", 1:4)) cmpr_nat[, (paste0("f_", v)) := round(fifelse(typed > 0, get(v) / typed, NA_real_), 3)]
fwrite(cmpr_nat[order(year)], file.path(DIR$mod, "table_serotype_composition_national.csv"))
cat("\n== serotypes\n"); print(cmpr_nat[order(year), .(year, typed, f_denv1, f_denv2, f_denv3, f_denv4)])

# ---- 4. placebos by season (earlier version; short, moving pre-period) ---------
SEASON_REAL <- as.integer(format(PARAM$vaccine_start, "%Y"))
MONTH_T <- as.integer(format(PARAM$vaccine_start, "%m"))
plc <- list()
for (window in c("3 seasons of pre-period (main)", "all available pre-period (sensitivity)")) {
  n_pre <- if (startsWith(window, "3")) 3L else Inf
  for (tp in (SEASON_REAL - 5L):SEASON_REAL) {
    t0  <- as.Date(sprintf("%d-%02d-01", tp, MONTH_T))
    ini <- if (is.finite(n_pre)) as.Date(sprintf("%d-%02d-01", tp - n_pre, MONTH_T)) else as.Date("1900-01-01")
    end <- seq(t0, by = "month", length.out = 12L)[12L]
    d <- panel[month >= ini & month <= end]
    d[, treat_pl := as.integer(target == 1L & g < NEVER & month >= t0)]
    m <- fepois(fml("treat_pl"), data = d, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
    r <- rr_row(m, "treat_pl", window)
    plc[[paste(window, tp)]] <- data.table(window = window, season = sprintf("%d/%02d", tp, (tp + 1L) %% 100L),
                                           vaccine = tp == SEASON_REAL, months = uniqueN(d$month), months_pre = uniqueN(d[month < t0, month]),
                                           admissions_target_post = d[treat_pl == 1L, sum(admissions)],
                                           RRT = r$RR, CI95_lwr = r$CI95_lwr, CI95_upr = r$CI95_upr, p = r$p,
                                           log_RRT = log(r$RR), n_obs = r$n_obs)
  }
}
plc <- rbindlist(plc)
perm <- plc[, { real <- .SD[vaccine == TRUE]; pl <- .SD[vaccine == FALSE]
  data.table(n_placebos = nrow(pl), RRT_real = real$RRT, placebos_smaller = sum(pl$RRT <= real$RRT),
             placebos_more_extreme = sum(abs(pl$log_RRT) >= abs(real$log_RRT)),
             p_permutation_unilateral = (1 + sum(pl$RRT <= real$RRT)) / (1 + nrow(pl)),
             p_permutation_twosided = (1 + sum(abs(pl$log_RRT) >= abs(real$log_RRT))) / (1 + nrow(pl)),
             placebos_significant = sum(pl$p < 0.05), p_minimum_possible = 1 / (1 + nrow(pl))) }, by = window]
fwrite(plc[, .(window, season, vaccine, months, months_pre, admissions_target_post, RRT, CI95_lwr, CI95_upr, p, n_obs)],
       file.path(DIR$mod, "table_placebos_season.csv"))
fwrite(perm, file.path(DIR$mod, "table_placebos_permutation.csv"))
cat("\n== placebos by season\n"); print(plc[, .(window, season, RRT = round(RRT, 3), CI95 = sprintf("%.2f-%.2f", CI95_lwr, CI95_upr), p = signif(p, 2))])

# ---- 5. contamination of the fixed age groups (first phase) --------------------------------------
trs <- panel[g < NEVER, unique(code_muni)]
vac <- readRDS(file.path(DIR$processed, "dengue_vaccination.rds"))[code_muni %in% trs][, year := as.integer(format(month, "%Y"))]
pop <- readRDS(file.path(DIR$processed, "population.rds"))[code_muni %in% trs & year %in% 2024:2025, .(pop = sum(pop)), by = .(year, agegroup)]
p1519_25 <- pop[year == 2025L & agegroup == "15-19", pop]
d1_1519  <- vac[agegroup == "15-19", sum(d1)]                          # direct doses at 15-19 (Feb-Apr/2025)
env      <- vac[agegroup == TARGET_AGEGROUP & year == 2024L, sum(d1)] / 5   # vaccinated at 14 in 2024 -> 15-19 in 2025 (uniform)
dil      <- pop[year == 2024L & agegroup == "05-09", pop] / 5           # aged 9 in 2024 -> 10-14 in 2025 without eligibility
cf <- panel[code_muni %in% trs & agegroup == TARGET_AGEGROUP & month == max(month), .(code_muni, cov_d1_l1)]
diag_c <- data.table(
  item = c("D1 given at ages 15-19 in the 521 municipalities", "2024 D1 recipients aging into 15-19 in 2025 (approx.)",
           "Total contamination of the 15-19 comparison group in 2025 (approx.)",
           "Entering 10-14 in 2025 without eligibility in 2024 (approx.)",
           "Municipalities with D1 coverage = 100% (Dec/2025)", "Municipalities with D1 coverage >= 90% (Dec/2025)"),
  val = c(d1_1519, round(env), NA, round(dil), cf[cov_d1_l1 >= 0.999, .N], cf[cov_d1_l1 >= 0.90, .N]),
  pct = c(100 * d1_1519 / p1519_25, 100 * env / p1519_25, 100 * (d1_1519 + env) / p1519_25,
          100 * dil / pop[year == 2025L & agegroup == TARGET_AGEGROUP, pop], NA, NA),
  observation = c("% of the 2025 population aged 15-19", "assuming a uniform distribution within 10-14",
                 "bias toward the null: part of the control group is protected", "% of the 2025 target population",
                 sprintf("of %d; median coverage %.1f%%", nrow(cf), 100 * median(cf$cov_d1_l1)), ""))
fwrite(diag_c, file.path(DIR$mod, "table_diagnostic_cohort.csv"))
cat("\n== contamination of the fixed age groups\n"); print(diag_c)

# =====================================================================================
# Blocks 6-9: coverage, monthly doses, hospitalizations by ICD code and serotype by UF.
# =====================================================================================
grp <- unique(panel_full[, .(code_muni, group = fifelse(g < NEVER, "first phase", fifelse(expn == 1L, "expanded", "never included")))])
vac_all <- readRDS(file.path(DIR$processed, "dengue_vaccination.rds"))
pop_all <- readRDS(file.path(DIR$processed, "population.rds"))
REF <- as.Date(c("2024-12-01", "2025-12-01"))                  # doses through the end of these months

# ---- 6. D1 and D2 coverage by municipality (10-14) and D1 by single year of age and cohort ----
# 6a. same method as 20: cumulative doses / 2024 population aged 10-14, capped at 100%
pa24 <- pop_all[year == 2024L & agegroup == TARGET_AGEGROUP, .(code_muni, pop_target = pop)]
cm <- rbindlist(lapply(REF, function(r) {
  x <- vac_all[agegroup == TARGET_AGEGROUP & month <= r, .(d1 = sum(d1), d2 = sum(d2)), by = code_muni]
  x <- merge(grp, x, by = "code_muni", all.x = TRUE)
  x[is.na(d1), d1 := 0L][is.na(d2), d2 := 0L]
  x <- merge(x, pa24, by = "code_muni")[pop_target > 0]
  x[, `:=`(month_ref = format(r, "%Y-%m"), cov_d1 = pmin(d1 / pop_target, 1), cov_d2 = pmin(d2 / pop_target, 1))]
}))
cov_mun <- cm[, .(municipalities = .N,
                  d1_weighted_pct = 100 * sum(cov_d1 * pop_target) / sum(pop_target), d1_mean_pct = 100 * mean(cov_d1),
                  d1_median_pct = 100 * median(cov_d1),
                  d2_weighted_pct = 100 * sum(cov_d2 * pop_target) / sum(pop_target), d2_mean_pct = 100 * mean(cov_d2),
                  d2_median_pct = 100 * median(cov_d2)), by = .(group, month_ref)][order(match(group, GROUP_ORDER), month_ref)]
fwrite(cov_mun, file.path(DIR$mod, "table_coverage_d1_d2_municipal.csv"))
cat("\n== D1 and D2 coverage at ages 10-14 by group\n"); print(cov_mun)

# 6b. D1 by single year of age (age aggregate saved by 15). Cohort = year of dose - age
# at dose (same approximation as 20). Population by single age: brpop only has five-year
# age groups, so each age gets the group's pop / group width (uniform split).
# Ages 5-19 and 20 (so the 2005 cohort has a value in Dec/2025).
FILE_AGE <- file.path(DIR$processed, "dengue_vaccination_age.rds")
if (file.exists(FILE_AGE)) {
  vi <- merge(readRDS(FILE_AGE)[!is.na(age)], grp, by = "code_muni")
  vi[, cohort := as.integer(format(month, "%Y")) - age]
  pg <- merge(pop_all[year %in% 2024:2025], grp, by = "code_muni")[, .(pop_agegroup = sum(pop)), by = .(group, year, agegroup)]
  width <- function(fx) as.integer(sub(".*-", "", fx)) - as.integer(sub("-.*", "", fx)) + 1L
  ci <- rbindlist(lapply(REF, function(r) {
    a <- as.integer(format(r, "%Y"))
    x <- CJ(group = unique(grp$group), age = 5:20)[, `:=`(year_ref = a, cohort = a - age)]
    x <- merge(x, vi[month <= r, .(d1 = sum(d1)), by = .(group, cohort)], by = c("group", "cohort"), all.x = TRUE)
    x[is.na(d1), d1 := 0L][, agegroup := to_agegroup(age)]
    x <- merge(x, pg[year == a, .(group, agegroup, pop_agegroup)], by = c("group", "agegroup"))
    x[, pop := pop_agegroup / width(agegroup)][, coverage_d1 := d1 / pop]
  }))
  ci <- ci[order(match(group, GROUP_ORDER), year_ref, age), .(group, year_ref, age, cohort, d1, pop, coverage_d1)]
  fwrite(ci, file.path(DIR$mod, "table_coverage_single_age.csv"))
  SET_COV <- list("2005-2007" = 2005:2007, "2007-2009" = 2007:2009, "2016-2018" = 2016:2018, "2017-2019" = 2017:2019,
                   "2009-2015 (exposed, clean set)" = 2009:2015, "2010-2015 (exposed, script 20 set)" = 2010:2015)
  cc <- rbindlist(lapply(names(SET_COV), function(k)
    ci[cohort %in% SET_COV[[k]], .(cohorts = k, ages = paste(range(age), collapse = "-"), d1 = sum(d1), pop = sum(pop)),
       by = .(group, year_ref)]))
  cc[, coverage_d1 := d1 / pop]
  fwrite(cc, file.path(DIR$mod, "table_coverage_cohorts.csv"))
  cat("\n== D1 coverage by set of cohorts\n"); print(cc)
} else log_msg("WARNING: ", FILE_AGE, " missing; coverage by single year of age not computed")

# ---- 7. monthly doses: first phase and expansion, and target-age hospitalizations in the first phase -----
f1 <- grp[group == "first phase", code_muni]; fa <- grp[group == "expanded", code_muni]
va <- vac_all[agegroup == TARGET_AGEGROUP]
dm <- merge(va[code_muni %in% f1, .(d1_first_phase = sum(d1), d2_first_phase = sum(d2)), by = month],
            va[code_muni %in% fa, .(d1_expanded = sum(d1), d2_expanded = sum(d2)), by = month], by = "month", all = TRUE)
dm <- merge(panel_full[g < NEVER & agegroup == TARGET_AGEGROUP & month >= as.Date("2023-01-01"),
                        .(admissions_target_first_phase = sum(admissions)), by = month], dm, by = "month", all.x = TRUE)
for (v in setdiff(names(dm), "month")) set(dm, which(is.na(dm[[v]])), v, 0L)
fwrite(dm[order(month)], file.path(DIR$mod, "table_doses_monthly.csv"))
d1f <- va[code_muni %in% f1]
doses_res <- data.table(
  item = c("D1 at ages 10-14 in the 521 first-phase municipalities (Feb/2024-Dec/2025)", "D1 given through May/2024", "D1 given from Jun/2024 onward",
           "% of D1 given after May 2024 (Jun/2024 onward)", "% of D1 given after May 2024 and still in 2024 (Jun-Dec/2024)",
           "D2 at ages 10-14 in the 521 first-phase municipalities", "% of D2 given after May 2024"),
  val = c(d1f[, sum(d1)], d1f[month < as.Date("2024-06-01"), sum(d1)], d1f[month >= as.Date("2024-06-01"), sum(d1)],
            100 * d1f[month >= as.Date("2024-06-01"), sum(d1)] / d1f[, sum(d1)],
            100 * d1f[month >= as.Date("2024-06-01") & month <= as.Date("2024-12-01"), sum(d1)] / d1f[, sum(d1)],
            d1f[, sum(d2)], 100 * d1f[month >= as.Date("2024-06-01"), sum(d2)] / d1f[, sum(d2)]))
fwrite(doses_res, file.path(DIR$mod, "table_doses_summ.csv"))
cat("\n== doses in the first phase\n"); print(doses_res)

# ---- 8. SIH hospitalizations by code (A90, A91, A97) and year -----------------------------
# the aggregates from 12 do not keep the ICD code: reads only the needed columns from the SIH
# cache, with the same filter as 12 (principal diagnosis ^A9(0|1|7), no continuation AIH, year by
# admission date). Result saved in data/interim (like the serotype in block 3) only if the
# cache is complete.
FILE_ICD <- file.path(DIR$interim, "sih_codes_year.rds")
if (file.exists(FILE_ICD)) icd <- readRDS(FILE_ICD) else {
  fs_sih <- list.files(file.path(DIR$cache, "SIH-RD"), "_full\\.parquet$", full.names = TRUE)
  ult <- PARAM$sih_last_billmonth
  n_expected <- length(UFS) * ((ult[["year"]] - PARAM$year_start) * 12L + ult[["month"]])
  cols_icd <- c("DIAG_PRINC", "IDENT", "DT_INTER", "IDADE", "COD_IDADE")
  icd <- if (length(fs_sih)) rbindlist(lapply(fs_sih, function(f) {
    d <- as.data.table(arrow::read_parquet(f, col_select = tidyselect::any_of(cols_icd)))
    for (v in setdiff(cols_icd, names(d))) set(d, j = v, value = NA_character_)
    d <- d[grepl(ICD_DENGUE, DIAG_PRINC) & !(as.character(IDENT) %in% "5")]
    if (!nrow(d)) return(NULL)
    d[, .(admissions = .N), by = .(year = as.integer(format(as_date(DT_INTER, "ymd"), "%Y")),
                                   code = substr(DIAG_PRINC, 1, 3), age = age_sih(IDADE, COD_IDADE))]
  }))[, .(admissions = sum(admissions)), by = .(year, code, age)] else NULL
  if (length(fs_sih) >= n_expected) saveRDS(icd, FILE_ICD)
  else log_msg("WARNING: SIH cache with ", length(fs_sih), " of ", n_expected, " files; table by code not saved in interim")
}
if (!is.null(icd) && nrow(icd)) {
  icd <- icd[!is.na(year) & year >= PARAM$year_start & year <= PARAM$year_end]
  tc <- dcast(icd[, .(n = sum(admissions)), by = .(year, code)], year ~ code, value.var = "n", fill = 0)
  for (v in setdiff(c("A90", "A91", "A97"), names(tc))) set(tc, j = v, value = 0)
  tc[, total := A90 + A91 + A97]
  tc <- merge(tc, icd[!is.na(age) & age >= 5L & age <= 59L, .(total_5a59 = sum(admissions)), by = year], by = "year", all.x = TRUE)
  tc <- merge(tc, icd[!is.na(age) & age >= 10L & age <= 14L, .(total_10a14 = sum(admissions)), by = year], by = "year", all.x = TRUE)
  tc[, `:=`(pct_A97 = 100 * A97 / total)]
  fwrite(tc[order(year), .(year, A90, A91, A97, total, pct_A97, total_5a59, total_10a14)], file.path(DIR$mod, "table_sih_codes_year.csv"))
  cat("\n== admissions by ICD code and year\n"); print(tc)
} else log_msg("WARNING: no SIH cache; table by ICD code not generated")

# ---- 9. predominant serotype and composition by UF, 2023 and 2024 -----------------------------
ABBR <- c("11" = "RO", "12" = "AC", "13" = "AM", "14" = "RR", "15" = "PA", "16" = "AP", "17" = "TO", "21" = "MA",
           "22" = "PI", "23" = "CE", "24" = "RN", "25" = "PB", "26" = "PE", "27" = "AL", "28" = "SE", "29" = "BA",
           "31" = "MG", "32" = "ES", "33" = "RJ", "35" = "SP", "41" = "PR", "42" = "SC", "43" = "RS", "50" = "MS",
           "51" = "MT", "52" = "GO", "53" = "DF")
su <- sor[year %in% 2023:2024, lapply(.SD, sum), by = .(year, uf), .SDcols = c("cases", "typed", paste0("denv", 1:4))]
su[, abbr := ABBR[uf]]
mx <- as.matrix(su[, paste0("denv", 1:4), with = FALSE])
su[, `:=`(predominant = fifelse(typed > 0, paste0("DENV-", max.col(mx, ties.method = "first")), NA_character_),
          f_predominant = fifelse(typed > 0, apply(mx, 1, max) / typed, NA_real_))]
for (v in paste0("denv", 1:4)) su[, (paste0("f_", v)) := fifelse(typed > 0, get(v) / typed, NA_real_)]
su[, highlight := abbr %in% c("DF", "RJ", "BA")]
su <- su[order(year, abbr), .(year, uf, abbr, highlight, cases, typed, pct_typed = 100 * typed / cases,
                              predominant, f_predominant, f_denv1, f_denv2, f_denv3, f_denv4)]
fwrite(su, file.path(DIR$mod, "table_serotype_uf.csv"))
cat("\n== serotype by state, highlights\n"); print(su[highlight == TRUE])

# =====================================================================================
# Blocks 10-14. Block 10 appends a row to the end of table3_sensitivity.csv; script 43
# publishes it under its own key (sensitivity_march_start).
# =====================================================================================
pa2 <- copy(panel)[, tot_mt := sum(admissions), by = .(code_muni, t)][tot_mt > 0]
pa2[, ta := as.integer(g < NEVER & target == 1L)]
years <- PARAM$year_start:PARAM$year_end
fml_year <- as.formula(sprintf("admissions ~ i(year, ta, ref = %d) | %s", years[1], FE_MAIN))
fit_year <- function(d) {                                   # 2024 and 2025 re-anchored on 2017-2023
  m <- fepois(fml_year, data = d, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
  reanchor_by_year(m, "ta")[year >= 2024L]
}

# ---- 10. treatment starting March 2024 (February stays in the pre-period) -------------
panel[, treat_mar := as.integer(target == 1L & g < NEVER & month >= as.Date("2024-03-01"))]
m_mar <- fepois(fml("treat_mar"), data = panel, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
l_mar <- rr_row(m_mar, "treat_mar", "Treatment from March 2024")
t3_current <- fread(file.path(DIR$mod, "table3_sensitivity.csv"))
t3_current <- rbind(t3_current[!grepl("^Treatment from March", analysis)], l_mar[, .(analysis, RRR = RR, CI95_lwr, CI95_upr, p, n_obs)])
fwrite(t3_current, file.path(DIR$mod, "table3_sensitivity.csv"))
cat("\n== start in March 2024\n"); print(l_mar)

# ---- 11. exclusion of DF, RJ and BA (together and each one) ----------------------------------
EXCL <- list("DF, RJ and BA" = c("53", "33", "29"), "DF" = "53", "RJ" = "33", "BA" = "29")
exc <- rbindlist(lapply(names(EXCL), function(nm) {
  d <- panel[!uf %in% EXCL[[nm]]]
  mg <- fepois(fml("treat"), data = d, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
  lg <- rr_row(mg, "treat", nm)
  an <- fit_year(pa2[!uf %in% EXCL[[nm]]])
  rbind(data.table(exclusion = nm, specification = "pooled", year = NA_integer_, RR = lg$RR, CI95_lwr = lg$CI95_lwr,
                   CI95_upr = lg$CI95_upr, p = lg$p, p_calibrated_unilateral = NA_real_,
                   municipalities_treated = d[g < NEVER, uniqueN(code_muni)]),
        an[, .(exclusion = nm, specification = "by year (ref. 2017-2023)", year, RR, CI95_lwr, CI95_upr, p = NA_real_,
               p_calibrated_unilateral, municipalities_treated = d[g < NEVER, uniqueN(code_muni)])])
}))
fwrite(exc, file.path(DIR$mod, "table_exclusion_uf.csv"))
cat("\n== state exclusion\n"); print(exc[, .(exclusion, specification, year, RR = round(RR, 3), CI95 = sprintf("%.3f-%.3f", CI95_lwr, CI95_upr))])

# ---- 12. leaving out one treated health region at a time --------------------
rsm2 <- merge(unique(panel[, .(code_muni, g)]), health_region(), by = "code_muni")
reg_treat <- rsm2[g < NEVER, unique(health_reg)]
nms_rs <- fread(file.path(DIR$raw, "lists", "dengue_vaccine_municipalities_full.csv"), colClasses = "character")
nms_rs <- merge(rsm2[, .(code_muni, health_reg)], nms_rs[phase == "1", .(code_muni = mun6(code_muni), uf, name_rs = health_region)],
                  by = "code_muni")[, .(name_rs = paste0(name_rs[1], " (", uf[1], ")")), by = health_reg]
log_msg("leave-one-out: ", length(reg_treat), " health regions with a first-phase municipality")
full <- fit_year(pa2)
loo <- rbindlist(lapply(reg_treat, function(r) {
  out <- rsm2[health_reg == r, code_muni]
  x <- fit_year(pa2[!code_muni %in% out])
  x[, .(health_reg = r, municipalities_excluded = length(out), treated_excluded = rsm2[health_reg == r & g < NEVER, .N],
        year, log_b, RR, CI95_lwr, CI95_upr)]
}))
loo <- merge(loo, nms_rs, by = "health_reg", all.x = TRUE)
loo <- merge(loo, full[, .(year, log_b_full = log_b)], by = "year")[, change_log := log_b - log_b_full]
fwrite(loo[order(year, health_reg)], file.path(DIR$mod, "table_loo_health_region.csv"))
loo_res <- loo[, .(regions = .N, RR_full = exp(log_b_full[1]), RR_min = min(RR), RR_median = median(RR), RR_max = max(RR),
                   region_largest_change = name_rs[which.max(abs(change_log))], RR_without_it = RR[which.max(abs(change_log))]), by = year]
fwrite(loo_res, file.path(DIR$mod, "table_loo_health_region_summ.csv"))
cat("\n== leave-one-out by health region\n"); print(loo_res)

# ---- 13. SIH negative controls: main binary model ------------------------
FILE_CTRL <- file.path(DIR$processed, "sih_controls_aggregate.rds")
CTRL <- c(ctrl_all = "All causes except A90/A91 and chapters O, P, Z", ctrl_a00a09 = "A00-A09 (intestinal infectious diseases)",
          ctrl_j00j22 = "J00-J22 (acute respiratory infections)", ctrl_a920 = "A92.0 (chikungunya)")
if (file.exists(FILE_CTRL)) {
  pcb <- merge(panel, readRDS(FILE_CTRL), by = c("code_muni", "agegroup", "month"), all.x = TRUE)
  for (v in names(CTRL)) set(pcb, which(is.na(pcb[[v]])), v, 0L)
  pcb[, `:=`(treat_2024 = treat * (year == 2024L), treat_2025 = treat * (year == 2025L))]
  bin <- rbindlist(lapply(names(CTRL), function(y) {
    d <- pcb[, tot := sum(get(y)), by = .(code_muni, t)][tot > 0]
    m1 <- tryCatch(fepois(fml("treat", y), data = d, offset = ~log(pop), cluster = ~code_muni, notes = FALSE), error = function(e) NULL)
    m2 <- tryCatch(fepois(fml("treat_2024 + treat_2025", y), data = d, offset = ~log(pop), cluster = ~code_muni, notes = FALSE), error = function(e) NULL)
    lin <- function(m, term, per) if (is.null(m)) data.table(period = per) else {
      r <- rr_row(m, term, per); data.table(period = per, RRR = r$RR, CI95_lwr = r$CI95_lwr, CI95_upr = r$CI95_upr, p = r$p, n_obs = r$n_obs) }
    rbind(lin(m1, "treat", "pooled Feb/2024-Dec/2025"), lin(m2, "treat_2024", "2024"), lin(m2, "treat_2025", "2025"),
          fill = TRUE)[, `:=`(outcome = CTRL[[y]], var = y, events_target_treated_post = pcb[treat == 1L, sum(get(y))])]
  }), fill = TRUE)
  fwrite(bin, file.path(DIR$mod, "table_controls_sih_binary.csv"))
  cat("\n== SIH negative controls: binary model\n")
  print(bin[, .(outcome, period, RRR = round(RRR, 3), CI95 = sprintf("%.3f-%.3f", CI95_lwr, CI95_upr))])
} else log_msg("WARNING: ", FILE_CTRL, " missing; SIH negative controls not estimated")

# ---- 14. single table: estimate × plausibility ceiling --------------------------
# ceiling = smallest RR compatible with coverage: 1 - VE × coverage (binary, months, first-phase
# cohort); gradient per 10 p.p.: 1 - VE × 0.10 × effective coverage / final coverage;
# coverage category (RR relative to < 10%): (1 - VE × c_cat) / (1 - VE × c_ref), with the
# literal formula 1 - VE × c_cat in a separate column. "violates" = the (point) estimate or the
# whole CI falls below the ceiling.
VE <- c("100" = 1, "746" = 0.746, "841" = 0.841, "879" = 0.879)
lm_ <- function(nm) { f <- file.path(DIR$mod, nm); if (file.exists(f)) fread(f) else NULL }
pl1 <- lm_("table_plausibility_coverage.csv"); pl2 <- lm_("table_plausibility_v2.csv"); seas <- lm_("table_event_monthly_seasonal.csv")
cvm <- lm_("table_inforce_coverage_agegroup.csv"); dta <- lm_("table_crosssec_dose_annual.csv")
rob <- lm_("table_robustness_crosssec_dose.csv"); cat33 <- lm_("table_crosssec_dose_categories.csv")
cnx <- lm_("table_birth_cohort_exact.csv")
row_t <- function(analysis, period, est, lo, hi, cov_used, cov, type = "direto", ratio_final = NA_real_, c_ref = NA_real_) {
  x <- data.table(analysis, period, estimate = est, cint_lwr = lo, cint_upr = hi, coverage_used = cov_used, coverage_effective = cov)
  for (v in names(VE)) {
    ceiling <- switch(type, direto = 1 - VE[[v]] * cov, gradient = 1 - VE[[v]] * 0.10 * ratio_final,
                   category = (1 - VE[[v]] * cov) / (1 - VE[[v]] * c_ref))
    set(x, j = paste0("ceiling_ve", v), value = ceiling)
  }
  x[, `:=`(ceiling_literal_ve100 = if (type == "category") 1 - cov else NA_real_,
           ceiling_literal_ve841 = if (type == "category") 1 - 0.841 * cov else NA_real_)]
  x[, `:=`(violates_point_ve100 = estimate < ceiling_ve100, violates_cint_ve100 = cint_upr < ceiling_ve100,
           violates_point_ve841 = estimate < ceiling_ve841, violates_cint_ve841 = cint_upr < ceiling_ve841)]
}
tt <- list()
PERS <- c("2024 (Feb-Dec)", "2025", "Whole period (Feb/2024-Dec/2025)")
if (!is.null(pl1)) for (p in PERS) {
  a <- pl1[level == "aggregate" & label_txt == p]
  tt[[length(tt) + 1]] <- row_t("Binary first phase x never (Table 1B)", p, a$obs_RR, a$obs_lwr, a$obs_upr,
                                  "current, weighted by admissions: observed", a$coverage_d1)
  if (!is.null(pl2)) { v2 <- pl2[group == "first phase" & period == p]
    for (cc in c("current_expected", "inforce_observed", "inforce_expected"))
      tt[[length(tt) + 1]] <- row_t("Binary first phase x never (Table 1B)", p, a$obs_RR, a$obs_lwr, a$obs_upr,
                                      sub("^inforce", "in-force", sub("_", ", weighted by admissions: ", cc)),
                                      v2[[cc]]) }
}
if (!is.null(seas) && !is.null(cvm)) for (k_ in 1:5) {
  s <- seas[type == "monthly" & k == k_]; cf <- cvm[group == "first phase" & as.Date(month) == as.Date(s$month)]
  for (cc in c("cov_current", "cov_inforce"))
    tt[[length(tt) + 1]] <- row_t("Month, seasonally matched re-anchoring", format(as.Date(s$month), "%Y-%m"), s$RR, s$CI95_lwr, s$CI95_upr,
                                    sub("cov_", "", cc), cf[[cc]])
}
cov_final_expn <- if (!is.null(rob)) rob[grepl("^final coverage Dec/2025, expanded", item), val] else NA_real_
if (!is.null(dta) && !is.null(pl2)) {
  g25 <- dta[samp == "expanded (2,230)" & year == 2025L]; v2 <- pl2[group == "expanded" & period == "2025"]
  for (cc in c("current_observed", "current_expected", "inforce_observed", "inforce_expected"))
    tt[[length(tt) + 1]] <- row_t("2025 gradient in the expanded municipalities, per 10 p.p.", "2025", g25$RR_10pp, g25$CI95_lwr, g25$CI95_upr,
                                    paste(cc, "/ final coverage"), v2[[cc]], "gradient", v2[[cc]] / cov_final_expn)
}
# final coverage categories (expansion, 2025): effective current coverage in 2025 for each category
cvmun <- file.path(DIR$processed, "inforce_coverage_municipality.rds")
if (!is.null(cat33) && file.exists(cvmun)) {
  cvx <- readRDS(cvmun)[group == "expanded" & month >= as.Date("2025-01-01")]
  cf_mun <- panel_full[expn == 1L, .(cov_final = max(cov_d1_l1)), by = code_muni]
  cf_mun[, category := as.character(cut(cov_final, c(-Inf, 0.10, 0.25, 0.40, 0.55, Inf),
                                         labels = c("<10%", "10-25%", "25-40%", "40-55%", ">=55%"), right = FALSE))]
  cvx <- merge(cvx, cf_mun[, .(code_muni, category)], by = "code_muni")
  cc_cat <- cvx[, .(c = sum(cov_inforce * admissions) / sum(admissions)), by = category]
  c_ref <- cc_cat[category == "<10%", c]
  for (ct in cat33[year == 2025L, category])
    tt[[length(tt) + 1]] <- with(cat33[year == 2025L & category == ct],
      row_t(paste("Final coverage category", ct, "vs <10% (expanded)"), "2025", RR, CI95_lwr, CI95_upr,
              "in-force, weighted by admissions: observed (relative to <10%)", cc_cat[category == ct, c], "category", c_ref = c_ref))
}
# eligible cohorts (2010-2015): current coverage in 2025 = cumulative D1 through the previous month /
# population of these cohorts (ages 10-14 + 1/5 of 15-19 in 2025), mean over the months of 2025
if (!is.null(cnx) && file.exists(FILE_AGE)) {
  vc <- merge(readRDS(FILE_AGE)[!is.na(age)], grp, by = "code_muni")
  vc[, cohort := as.integer(format(month, "%Y")) - age]
  vc <- vc[cohort %in% 2010:2015, .(d1 = sum(d1)), by = .(group, month)][order(match(group, GROUP_ORDER), month)]
  vc[, cum_l1 := shift(cumsum(d1), 1, fill = 0), by = group]
  pc25 <- merge(pop_all[year == 2025L], grp, by = "code_muni")[, .(pop = sum(pop[agegroup == "10-14"]) + sum(pop[agegroup == "15-19"]) / 5), by = group]
  cov_cohorts <- merge(vc[month >= as.Date("2025-01-01")], pc25, by = "group")[, .(c = mean(pmin(cum_l1 / pop, 1))), by = group]
  for (nm in unique(cnx$cohortset)) {
    x <- cnx[cohortset == nm & specification == "by year, reference 2017-2023" & year == 2025L]
    f1 <- x[grepl("^first phase", contrast)]; am <- x[grepl("^Expanded", contrast)]
    tt[[length(tt) + 1]] <- row_t(paste("Cohort first phase x never,", nm), "2025", f1$RR, f1$CI95_lwr, f1$CI95_upr,
                                    "in-force, cohorts 2010-2015 (2025 mean)", cov_cohorts[group == "first phase", c])
    cam <- cov_cohorts[group == "expanded", c]
    tt[[length(tt) + 1]] <- row_t(paste("Cohort gradient in the expanded municipalities, per 10 p.p.,", nm), "2025", am$RR, am$CI95_lwr, am$CI95_upr,
                                    "in-force, cohorts 2010-2015 / final coverage", cam, "gradient", cam / cov_final_expn)
  }
}
tt <- rbindlist(tt, fill = TRUE)
fwrite(tt, file.path(DIR$mod, "table_estimate_vs_ceiling.csv"))
cat("\n== estimate vs bound\n")
print(tt[, .(analysis = substr(analysis, 1, 40), period, est = round(estimate, 3), cov = round(coverage_effective, 3),
             ceiling100 = round(ceiling_ve100, 3), ceiling841 = round(ceiling_ve841, 3), violates_point_ve100, violates_cint_ve100)])
log_msg("35 done")
