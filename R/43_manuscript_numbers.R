# =============================================================================
# 43_manuscript_numbers.R — every number cited in the text, in a single JSON
#
# Reads the tables in outputs/models/ (and the panel, for the descriptive counts)
# and writes outputs/manuscript_numbers.json with a STABLE NAME for each number.
# It is the file that 99_check_reproduction.R compares with the frozen version in
# docs/reference/ — and what the manuscript cites.
# Conventions: ratios and fractions are numbers; estimates come as
# {val, cint_lwr, cint_upr, p}; `_meta` holds the date, R version and missing tables.
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(jsonlite)

NEVER <- 10000L
missing <- character()
tab <- function(name) {                                   # NULL (and logged) if the table does not exist
  f <- file.path(DIR$mod, name)
  if (!file.exists(f)) { missing <<- c(missing, name); return(NULL) }
  fread(f)
}
est <- function(d, val, lwr, upr, pp = NULL) {          # {val, cint_lwr, cint_upr, p}
  if (is.null(d) || !nrow(d)) return(NULL)
  o <- list(val = d[[val]][1], cint_lwr = d[[lwr]][1], cint_upr = d[[upr]][1])
  if (!is.null(pp)) o$p <- d[[pp]][1]
  o
}
J <- list()

# ---- 1. design and sample ----------------------------------------------------------
panel <- readRDS(file.path(DIR$processed, "panel.rds"))
pa  <- panel[expn == 0L]
trs <- pa[g < NEVER, unique(code_muni)]
J$design <- list(year_start = PARAM$year_start, year_end = PARAM$year_end, vaccine_start = format(PARAM$vaccine_start, "%Y-%m-%d"),
                  months_series = uniqueN(panel$month), agegroup_target = TARGET_AGEGROUP, agegroups_comparison = setdiff(AGEGROUP_LABELS, TARGET_AGEGROUP),
                  threshold_coverage_sensitivity = PARAM$threshold_coverage)
J$samp <- list(
  municipalities_sample_main = uniqueN(pa$code_muni), municipalities_first_phase = length(trs),
  municipalities_never_included = uniqueN(pa$code_muni) - length(trs), municipalities_expanded = panel[expn == 1L, uniqueN(code_muni)],
  municipalities_included_total = length(trs) + panel[expn == 1L, uniqueN(code_muni)],
  pop_target_first_phase_2024 = pa[code_muni %in% trs & agegroup == TARGET_AGEGROUP & year == 2024L, sum(pop) / uniqueN(month)],
  pop_target_comparison_2024 = pa[!code_muni %in% trs & agegroup == TARGET_AGEGROUP & year == 2024L, sum(pop) / uniqueN(month)],
  admissions_5a59 = pa[, sum(admissions)], admissions_agegroup_target = pa[agegroup == TARGET_AGEGROUP, sum(admissions)],
  cases_probable = pa[, sum(den_cases)], admissions_target_treated_post = pa[treat == 1L, sum(admissions)],
  admissions_target_comparison_post = pa[agegroup == TARGET_AGEGROUP & !code_muni %in% trs & month >= PARAM$vaccine_start, sum(admissions)],
  deaths_target_series = pa[agegroup == TARGET_AGEGROUP, sum(deaths_sim)],
  deaths_target_treated_post = pa[agegroup == TARGET_AGEGROUP & code_muni %in% trs & month >= PARAM$vaccine_start, sum(deaths_sim)],
  deaths_target_comparison_post = pa[agegroup == TARGET_AGEGROUP & !code_muni %in% trs & month >= PARAM$vaccine_start, sum(deaths_sim)])

# ---- 2. vaccination coverage (SI-PNI) -------------------------------------------------
vac <- readRDS(file.path(DIR$processed, "dengue_vaccination.rds"))[month >= PARAM$vaccine_start & month <= max(panel$month)]
peak   <- vac[, .(d1 = sum(d1)), by = month][order(-d1)][1]
peak_a <- vac[agegroup == TARGET_AGEGROUP, .(d1 = sum(d1)), by = month][order(-d1)][1]
lst_official <- panel[g < NEVER | expn == 1L, unique(code_muni)]
ult <- pa[agegroup == TARGET_AGEGROUP & month == max(month)]
# population-weighted, simple mean and median across municipalities: cited in
# different contexts and NOT interchangeable
summ_cov <- function(d, v) list(weighted_by_population = 100 * sum(d[[v]] * d$pop) / sum(d$pop),
                                  mean_between_municipalities = 100 * mean(d[[v]]), median_between_municipalities = 100 * median(d[[v]]))
J$coverage <- list(
  d1_total = sum(vac$d1), d2_total = sum(vac$d2),
  pct_d1_agegroup_target = 100 * vac[agegroup == TARGET_AGEGROUP, sum(d1)] / sum(vac$d1), pct_d2_agegroup_target = 100 * vac[agegroup == TARGET_AGEGROUP, sum(d2)] / sum(vac$d2),
  peak_d1_month = format(peak$month, "%Y-%m"), peak_d1_doses_all_ages = peak$d1, peak_d1_doses_agegroup_target = peak_a$d1,
  pct_d1_em_municipalities_da_lst_official = 100 * vac[code_muni %in% lst_official, sum(d1)] / sum(vac$d1),
  pct_d1_em_municipalities_first_phase = 100 * vac[code_muni %in% trs, sum(d1)] / sum(vac$d1),
  cov_d1_dec2025_first_phase = summ_cov(ult[code_muni %in% trs], "cov_d1_l1"), cov_d1_dec2025_never = summ_cov(ult[!code_muni %in% trs], "cov_d1_l1"),
  cov_d2_dec2025_first_phase = summ_cov(ult[code_muni %in% trs], "cov_d2_l1"), cov_d2_dec2025_never = summ_cov(ult[!code_muni %in% trs], "cov_d2_l1"),
  municipalities_with_coverage_100pct = ult[code_muni %in% trs & cov_d1_l1 >= 0.999, .N],
  cov_final_expanded_percentiles_10_50_90 = as.list(setNames(100 * quantile(panel[expn == 1L, .(cf = max(cov_d1_l1)), by = code_muni]$cf, c(.1, .5, .9)),
                                                            c("p10", "p50", "p90"))))

# ---- 3. main effect, by vaccination year, event study ---------------------------------
t2 <- tab("table2_main_effect.csv")
if (!is.null(t2)) {
  key <- c("Admissions (SIH)" = "admissions_sih", "Probable cases (SINAN)" = "cases_probable_sinan",
             "Dengue with warning signs/severe (SINAN)" = "cases_severe_sinan", "Hospitalized (SINAN)" = "hospitalized_sinan")
  J$main_effect <- setNames(lapply(names(key), function(d) {
    x <- t2[outcome == d]; c(est(x, "RRR", "CI95_lwr", "CI95_upr", "p"), list(reduction_pct = x$reduction_pct[1], n_obs = x$n_obs[1])) }), unname(key))
}
ev <- tab("averted_admissions.csv")
if (!is.null(ev)) J$averted_admissions <- list(observed = ev$observed, averted = ev$averted)
het <- tab("heterogeneity_season.csv")
if (!is.null(het)) {
  setnames(het, c("Estimate", "Std. Error"), c("b", "se"))
  J$effect_by_season <- setNames(lapply(het$term, function(tm) { x <- het[term == tm]
    list(val = exp(x$b), cint_lwr = exp(x$b - 1.96 * x$se), cint_upr = exp(x$b + 1.96 * x$se), p = x[[5]]) }),
    sub("^treat_", "season_", het$term))
}
pre <- tab("table_pretrend.csv")
if (!is.null(pre)) {
  v <- function(m) pre[measure == m, val]
  J$event_study <- list(
    att_reanchored = list(val = v("Re-anchored ATT (post mean - pre mean)"), cint_lwr = v("CI95_lwr"), cint_upr = v("CI95_upr"), p = v("p ATT")),
    wald_pre_chi2 = v("Wald pre (X2)"), wald_pre_gl = v("gl"), wald_pre_p = v("p Wald pre"),
    flatness_chi2 = v("Wald pre flatness (k<=-3, no anchor)"), flatness_gl = v("df flatness"), flatness_p = v("p flatness"),
    trend_linear_pre = v("Linear pre-trend (log-RR/month)"), trend_linear_p = v("p linear trend"))
}
ee <- tab("table_event_study_reanchored.csv")
if (!is.null(ee)) J$event_study$coefficients_reanchored <- ee[, .(k, RR = RR_reanch, cint_lwr = lwr_reanch, cint_upr = upr_reanch,
                                                                        RR_ref_m1, cint_lwr_ref_m1 = lwr_ref_m1, cint_upr_ref_m1 = upr_ref_m1)]

# ---- 4. validity: plausibility, periods, permutation, dose-response, cohorts ----------
pl <- tab("table_plausibility_coverage.csv")
if (!is.null(pl)) {
  ag <- pl[level == "aggregate"]
  J$plausibility <- list(
    efficacies = list(ef675 = "Ranzani et al. 2025 (ages 10-14, SP, against hospitalization)", ef746 = "Krug Mareto et al. 2026, 1 dose, 1 year",
                     ef879 = "Krug Mareto et al. 2026, 2 doses, 1 year", ef841 = "Tricou et al. 2024, phase 3 trial, 54 months"),
    aggregates = setNames(lapply(seq_len(nrow(ag)), function(i) { x <- ag[i]
      list(coverage_effective_pct = 100 * x$coverage_d1, coverage_wtd_time_pct = 100 * x$coverage_d1_wtd_time, admissions = x$admissions,
           bound_efficacy_100 = x$bound_ef100, bound_efficacy_675 = x$bound_ef675, bound_efficacy_746 = x$bound_ef746,
           bound_efficacy_841 = x$bound_ef841, bound_efficacy_879 = x$bound_ef879,
           observed = x$obs_RR, observed_cint_lwr = x$obs_lwr, observed_cint_upr = x$obs_upr, violates_bound_efficacy_100 = isTRUE(x$violates_ef100)) }),
      c("season_2024", "season_2025", "period_all")),
    monthly = pl[level == "monthly", .(month = label_txt, coverage_pct = 100 * coverage_d1, admissions, bound_efficacy_100 = bound_ef100,
                                     bound_efficacy_675 = bound_ef675, observed_reanchored = obs_RR)])
}
ea <- tab("table_event_annual.csv");     if (!is.null(ea)) J$event_annual <- ea[, .(period, RRR, cint_lwr = CI95_lwr, cint_upr = CI95_upr)]
et <- tab("table_event_season.csv"); if (!is.null(et)) J$event_season <- et[, .(period, RRR, cint_lwr = CI95_lwr, cint_upr = CI95_upr)]
cp <- tab("table_calibration_placebo.csv"); if (!is.null(cp)) J$calibration_placebo <- cp[, .(subset, period, sd_log_placebos, t, gl, p_unilateral, p_twosided)]
ca <- tab("table_cluster_alternative.csv"); if (!is.null(ca)) J$cluster_alternative <- ca[, .(clustering, RRR, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p, se_log)]
ar <- tab("table_att_reanchorings.csv"); if (!is.null(ar)) J$att_reanchorings <- ar[, .(reanchoring, months_pre, RRR, cint_lwr = CI95_lwr, cint_upr = CI95_upr)]
pe <- tab("table_spatial_permutation.csv"); if (!is.null(pe)) J$spatial_permutation <- pe
dt <- tab("table_crosssec_dose_annual.csv")
if (!is.null(dt)) J$crosssec_dose_annual <- dt[, .(samp, year, RR_10pp, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p_calibrated_unilateral, municipalities, admissions_target_2025)]
dd <- tab("table_crosssec_dose_2024_2025.csv")
if (!is.null(dd)) J$crosssec_dose_outcomes <- dd[, .(samp, outcome, year, RR_10pp, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p)]
dcat <- tab("table_crosssec_dose_categories.csv"); if (!is.null(dcat)) J$crosssec_dose_categories <- dcat
dpm <- tab("table_crosssec_dose_permutation.csv"); if (!is.null(dpm)) J$crosssec_dose_permutation <- dpm
rb <- tab("table_robustness_crosssec_dose.csv"); if (!is.null(rb)) J$robustness_crosssec_dose <- rb
dri <- tab("table_dose_response_withinseason.csv")
if (!is.null(dri)) J$dose_response_withinseason <- setNames(lapply(seq_len(nrow(dri)), function(i) { x <- dri[i]
  list(analysis = x$analysis, RR_by_10pp = x$RR_10pp, cint_lwr = x$CI95_lwr, cint_upr = x$CI95_upr, p = x$p, municipalities = x$municipalities,
       admissions_agegroup_target = x$admissions_target, coverage_p25 = x$cov_p25, coverage_p50 = x$cov_p50, coverage_p75 = x$cov_p75) }),
  c("season_2025", "season_2024", "seasons_joint", "falsification_2019", "falsification_2023")[seq_len(nrow(dri))])
dr <- tab("table7_dose_response.csv"); if (!is.null(dr)) J$dose_response_longitudinal <- dr[, .(model, term, RR, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p)]
cn <- tab("table_birth_cohort.csv"); if (!is.null(cn)) J$birth_cohort <- cn
cna <- tab("table_birth_cohort_annual.csv"); if (!is.null(cna)) J$birth_cohort_annual <- cna[, .(year, RRR, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p_calibrated_unilateral)]
cg <- tab("table_cohort_group.csv")
if (!is.null(cg)) J$cohort_group <- cg[, .(analysis, year, RR, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p_calibrated_unilateral, municipalities, admissions_exposed_2025)]

# ---- 5. sensitivity, specifications, serotype, placebos, contamination -------------------
t3 <- tab("table3_sensitivity.csv")
# the March-start row (block 10 of script 35) stays out of this section and goes to the key
# sensitivity_march_start
t3_mar <- if (!is.null(t3)) t3[grepl("^Treatment from March", analysis)] else NULL
if (!is.null(t3)) t3 <- t3[!grepl("^Treatment from March", analysis)]
if (!is.null(t3)) J$sensitivity <- list(n_analyses = nrow(t3), agegroup_RRR = c(min(t3$RRR, na.rm = TRUE), max(t3$RRR, na.rm = TRUE)),
                                          rows = t3[, .(analysis, RRR, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p)])
expd <- tab("table_specifications_regional.csv"); if (!is.null(expd)) J$specifications_regional <- expd[, .(specification, RRR, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p)]
sc <- tab("table_serotype_completeness_year.csv"); if (!is.null(sc)) J$serotype_completeness <- sc[, .(year, cases, typed, pct_typed, ufs_with_30_typed, ufmes_with_min)]
sn <- tab("table_serotype_composition_national.csv"); if (!is.null(sn)) J$serotype_composition_national <- sn[, .(year, typed, f_denv1, f_denv2, f_denv3, f_denv4)]
pb <- tab("table_placebos_season.csv"); pp <- tab("table_placebos_permutation.csv")
if (!is.null(pb)) {
  jan <- pb[grepl("^3 seasons", window)]
  J$placebos_season <- list(window = jan$window[1], seasons = jan[, .(season, vaccine, RRT, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p)],
                               agegroup_placebos = c(min(jan[vaccine == FALSE, RRT]), max(jan[vaccine == FALSE, RRT])),
                               placebos_significant = jan[vaccine == FALSE & p < 0.05, season],
                               sensitivity_all_pre = pb[!grepl("^3 seasons", window), .(season, RRT, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p)])
  if (!is.null(pp)) { x <- pp[grepl("^3 seasons", window)]
    J$placebos_season$permutation <- list(n_placebos = x$n_placebos, p_unilateral = x$p_permutation_unilateral,
                                            p_twosided = x$p_permutation_twosided, p_minimum_possible = x$p_minimum_possible) }
}
dc <- tab("table_diagnostic_cohort.csv"); if (!is.null(dc)) J$diagnostic_cohort <- dc

# ---- 5b. coverage, doses, ICD codes, serotype and permutation by strata -------------------
cm <- tab("table_coverage_d1_d2_municipal.csv"); ci <- tab("table_coverage_single_age.csv")
cc <- tab("table_coverage_cohorts.csv")
if (!is.null(cm) || !is.null(cc)) J$coverage_single_age <- list(
  method = paste("SI-PNI doses cumulated through the end of Dec/2024 and Dec/2025. By municipality (10-14): / population aged",
                 "10-14 in 2024, capped at 100%. By single year of age: cohort = year of dose - age at dose; population",
                 "by single year of age by UNIFORM SPLIT of the age group (brpop has no single-year estimates)."),
  municipal_d1_d2_10a14 = cm, cohorts = cc, single_age = ci)
dmr <- tab("table_doses_summ.csv"); dmm <- tab("table_doses_monthly.csv")
if (!is.null(dmr)) J$doses_monthly <- list(
  pct_d1_after_may2024_first_phase_10a14 = dmr[grepl("^% of D1 given after May 2024 \\(Jun", item), val],
  summ = dmr, series = dmm)
d24 <- tab("table_crosssec_dose_cov_dec2024.csv"); d24p <- tab("table_crosssec_dose_cov_dec2024_permutation.csv")
if (!is.null(d24)) J$crosssec_dose_cov_dec2024 <- list(
  year_2025 = c(est(d24[year == 2025L], "RR_10pp", "CI95_lwr", "CI95_upr"),
               list(p_calibrated_unilateral = d24[year == 2025L, p_calibrated_unilateral],
                    p_permutation_unilateral = if (!is.null(d24p)) d24p[year == 2025L, p_unilateral] else NULL)),
  annual = d24[, .(year, RR_10pp, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p_calibrated_unilateral, municipalities,
                  cov_dec24_median, cov_final_median)],
  permutation = d24p)
cgs <- tab("table_cohort_group_sensitivity.csv")
if (!is.null(cgs)) J$birth_cohort_sensitivity <- cgs[, .(analysis, sensitivity, year, RR, cint_lwr = CI95_lwr, cint_upr = CI95_upr,
                                                                  p_calibrated_unilateral)]
pes <- tab("table_permutation_strata.csv"); pet <- tab("table_spatial_permutation_t.csv")
if (!is.null(pes)) J$spatial_permutation_strata <- c(as.list(pes), list(p_over_t = pet))
sca <- tab("table_sih_codes_year.csv"); if (!is.null(sca)) J$sih_codes_year <- sca
suf <- tab("table_serotype_uf.csv")
if (!is.null(suf)) J$serotype_uf <- list(highlights = suf[highlight == TRUE, .(year, abbr, typed, predominant, f_predominant,
                                                                            f_denv1, f_denv2, f_denv3, f_denv4)],
                                         all_ufs = suf[, .(year, abbr, typed, pct_typed, predominant, f_predominant,
                                                             f_denv1, f_denv2, f_denv3, f_denv4)])

# ---- 5c. birth cohort, plausibility, sensitivity analyses and negative controls ------------
if (!is.null(t3_mar) && nrow(t3_mar)) J$sensitivity_march_start <- t3_mar[, .(analysis, RRR, cint_lwr = CI95_lwr, cint_upr = CI95_upr, p, n_obs)]
cnx <- tab("table_birth_cohort_exact.csv"); cne <- tab("table_birth_cohort_exact_entry.csv")
if (!is.null(cnx)) J$birth_cohort_exact <- list(
  method_cohort = cnx$method_cohort[1],
  note = "NASC missing from the SIH cache: cohort expected from age and admission date (birthdays uniform over the year)",
  estimates = cnx[, .(cohortset, contrast, specification, gl, year, RR, cint_lwr = CI95_lwr, cint_upr = CI95_upr,
                        t_calibrated, p_calibrated_unilateral, municipalities, admissions_exposed_2025)],
  entry_cohorts_cohortset_ii = cne)
mde <- tab("table_cohort_mde.csv"); if (!is.null(mde)) J$cohort_mde <- mde
cvf <- tab("table_inforce_coverage_agegroup.csv"); if (!is.null(cvf)) J$inforce_coverage_agegroup <- cvf
pv2 <- tab("table_plausibility_v2.csv"); if (!is.null(pv2)) J$plausibility_v2 <- pv2
evt <- tab("table_estimate_vs_ceiling.csv"); if (!is.null(evt)) J$estimate_vs_ceiling <- evt
seas <- tab("table_event_monthly_seasonal.csv")
if (!is.null(seas)) J$event_study_seasonal <- list(means = seas[type != "monthly", .(type, RR, cint_lwr = CI95_lwr, cint_upr = CI95_upr)],
                                                    mar_jul_2024 = seas[mar_jul_2024 == TRUE, .(month, RR, cint_lwr = CI95_lwr, cint_upr = CI95_upr)],
                                                    series = seas[type == "monthly", .(k, month, RR, cint_lwr = CI95_lwr, cint_upr = CI95_upr)])
cej <- tab("table_composition_age_windows.csv"); if (!is.null(cej)) J$composition_age_windows <- cej
exu <- tab("table_exclusion_uf.csv"); if (!is.null(exu)) J$exclusion_uf_anomaly <- exu
loo <- tab("table_loo_health_region.csv"); lr <- tab("table_loo_health_region_summ.csv")
if (!is.null(lr)) J$loo_health_region <- list(summ = lr, by_region = loo[, .(year, health_reg, name_rs, municipalities_excluded,
                                                                              treated_excluded, RR, cint_lwr = CI95_lwr, cint_upr = CI95_upr)])
s1r <- tab("table_s1_reconciliation.csv"); if (!is.null(s1r)) J$s1_reconciliation <- s1r
cb <- tab("table_controls_sih_binary.csv"); cgr <- tab("table_controls_sih_gradient.csv")
if (!is.null(cb) || !is.null(cgr)) J$negative_controls_sih <- list(binary = cb, gradient_expanded = cgr)
prs <- tab("table_permutation_health_region.csv"); if (!is.null(prs)) J$permutation_health_region <- prs

# ---- 6. meta ------------------------------------------------------------------------------
manif <- file.path(ROOT, "outputs", "data_manifest.csv")
J$`_meta` <- list(generated_em = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), script = "R/43_manuscript_numbers.R",
                  r_version = R.version.string,
                  manifest = if (file.exists(manif)) { m <- fread(manif); list(fnames = nrow(m), gigabytes = round(sum(m$bytes) / 2^30, 2),
                                                                             generated_em = format(file.mtime(manif), "%Y-%m-%d %H:%M:%S")) } else NULL,
                  missing_tables = missing)
# extraction dates of each database (file write dates, read from the manifest) and versions
if (file.exists(manif)) {
  m <- fread(manif)
  m[, base := fifelse(grepl("^data/cache_datasus/", path), sub("^data/cache_datasus/([^/]+)/.*$", "\\1", path),
         fifelse(grepl("^data/raw/sipni/reduced/", path), "SI-PNI (reduced)", NA_character_))]
  J$`_meta`$extraction_sources <- m[!is.na(base), .(fnames = .N, first = min(modified_em), last = max(modified_em)), by = base]
}
J$`_meta`$versoes <- list(R = R.version.string, R_session = utils::sessionInfo()$R.version$version.string,
                          data.table = as.character(packageVersion("data.table")),
                          fixest = as.character(packageVersion("fixest")),
                          microdatasus = tryCatch(as.character(packageVersion("microdatasus")), error = function(e) "missing"))
target_dir <- file.path(ROOT, "outputs", "manuscript_numbers.json")
write_json(J, target_dir, auto_unbox = TRUE, digits = 8, pretty = TRUE, na = "null")
log_msg("manuscript_numbers.json: ", length(J) - 1L, " sections", if (length(missing)) paste0(" — MISSING: ", paste(missing, collapse = ", ")))
