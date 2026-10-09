# =============================================================================
# 40_tables.R — manuscript tables (flextable -> .docx)
#
# Reads the intermediate tables from outputs/models/ and builds:
#   Table 1  annual rates by municipality group and age group
#   Table 2  A: pooled effect by outcome; B: by period × coverage threshold
#   Table 3  A: periods with a common reference; B: dose-response; C: cohorts
#   Table S1 sensitivity; S2 placebos by season; S3 dose-response robustness
# Outputs: outputs/tables/tables.docx and one .csv per table.
# The full notes for each table are in the manuscript; here only
# the title and a short note.
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(flextable); library(officer)

read <- function(name) fread(file.path(DIR$mod, name))
t1   <- read("table1_rates_year_agegroup_group.csv");   t2 <- read("table2_main_effect.csv")
pl   <- read("table_plausibility_coverage.csv"); ea <- read("table_event_annual.csv")
et   <- read("table_event_season.csv");         cal <- read("table_calibration_placebo.csv")
pe   <- read("table_spatial_permutation.csv");      da <- read("table_crosssec_dose_annual.csv")
dd   <- read("table_crosssec_dose_2024_2025.csv"); dp <- read("table_crosssec_dose_permutation.csv")
cg   <- read("table_cohort_group.csv");             cna <- read("table_birth_cohort_annual.csv")
s3   <- read("table3_sensitivity.csv");           cl <- read("table_cluster_alternative.csv")
expd  <- read("table_specifications_regional.csv"); plc <- read("table_placebos_season.csv")
rob  <- read("table_robustness_crosssec_dose.csv")

{
  # ---- formatting and texts ------------------------------------------------------
  dm <- "."; bm <- ","
  f  <- function(x, d = 2) formatC(x, format = "f", digits = d, decimal.mark = dm, big.mark = bm)
  cint <- function(rr, lo, hi, d = 2) sprintf("%s (%s–%s)", f(rr, d), f(lo, d), f(hi, d))
  agegroup <- function(x, d = 2) sprintf("%s %s %s", f(min(x), d), "to", f(max(x), d))
  TXT <- list(
    t1 = "Table 1. Annual dengue hospitalisation rate per 100,000 inhabitants by municipality group, age group and year, Brazil, 2017–2025",
    t1h = c("Year", "First-phase\n10–14 years", "First-phase\ncomparison", "Ratio", "Never-included\n10–14 years", "Never-included\ncomparison", "Ratio"),
    t1n = "First-phase: 521 municipalities; never-included: 2,819. Comparison: 5–9, 15–19 and 20–29 years (the model uses all ages 5–59). Source: SIH/SUS and DATASUS.",
    t2 = "Table 2. Association estimates and the maximum reduction compatible with vaccine coverage, first-phase municipalities, Brazil, 2017–2025",
    t2a = "Panel A. Pooled estimate by outcome, February 2024 to December 2025",
    t2ah = c("Outcome (source)", "RRR", "95% CI", "p value", "Observations (n)"),
    t2b = "Panel B. Dengue hospitalisations by period against the plausibility bound",
    t2bh = c("Period", "Effective D1 coverage", "Observed RRR (95% CI)", "Smallest possible RRR: efficacy 100% / 67.5% / 84.1%", "Compatible with the vaccine?"),
    t2n = "RRR: ratio of rate ratios (Poisson with municipality × month, municipality × age and age × month fixed effects; SE clustered by municipality). Effective coverage: first-dose coverage in the target group weighted by each month's admissions. Bound = 1 − efficacy × coverage.",
    outc = c("Admissions (SIH)" = "Dengue hospitalisations (SIH/SUS)", "Probable cases (SINAN)" = "Probable cases (SINAN)",
             "Hospitalized (SINAN)" = "Hospitalised cases (SINAN)", "Dengue with warning signs/severe (SINAN)" = "Warning signs or severe dengue (SINAN)"),
    per = c("2024 (Feb-Dec)" = "2024 (Feb–Dec)", "2025" = "2025", "Whole period (Feb/2024-Dec/2025)" = "Pooled period"),
    sim = c("yes", "no"),
    t3 = "Table 3. Validity analyses: period estimates with a common reference, dose–response by final coverage and birth-cohort comparison, Brazil, 2017–2025",
    t3a = "Panel A. Treated × target-group contrast by period, re-anchored to the mean of the seven pre-vaccine periods",
    t3ah = c("Calendar year", "RRR (95% CI)", "Season (Feb–Jan)", "RRR (95% CI)"),
    vac = "(vaccine)", until_dec = "to Dec.",
    t3b = "Panel B. Dose–response: ratio per 10 percentage points of final D1 coverage in the target group, by year, re-anchored to the 2017–2023 mean",
    t3bh = c("Sample", "Pre-vaccine years (2017–2023)", "2024", "2025", "Calibrated p, 2025", "Permutation p, 2025"),
    samples = c("expanded (2,230)" = "Expansion (n = 2,230)", "all municipalities" = "All municipalities (n = 5,570)", "first phase (521)" = "First phase (n = 521)"),
    outc_b = c("Probable cases (SINAN)" = "Expansion: probable cases (SINAN)*", "Hospitalized (SINAN)" = "Expansion: hospitalised cases (SINAN)*",
               "Negative control: chikungunya (SINAN)" = "Expansion: chikungunya (negative control)*"),
    t3c = "Panel C. Birth-cohort comparison: eligible cohorts against non-eligible neighbouring cohorts, same municipality and month",
    t3ch = c("Contrast", "Neighbouring cohorts", "2017–2023 (no vaccine)", "2024", "2025"),
    c_f1 = "First phase × never included", c_expn = "Expansion: final coverage, per 10 p.p.",
    adj20 = "2007–2009 and 2016–2018", adjl = "2005–2007 and 2017–2019",
    t3n = "Calibrated p: t test with 6 df on the dispersion of pre-vaccine periods. Permutation: 400 spatial draws (panel A) and 300 permutations of coverage across municipalities (panel B). * Model with 2024 and 2025 indicators only.",
    s1 = "Table S1. Negative controls and sensitivity analyses of the pooled effect on dengue hospitalisations in 10–14-year-olds, Brazil, 2017–2025",
    s1h = c("Analysis", "RRR", "95% CI", "p value"),
    s1rot = c("Main model", "SE clustered by health region", "SE clustered by state", "Negative control: chikungunya cases",
              "Temporal placebo: adoption 24 months earlier", "Comparison with 20–59 years only", "Comparison with 5–9 and 15–19 years only",
              "Exposure by SI-PNI coverage (≥ 20% D1)", "Excluding February–April 2025", "Expansion municipalities as treated",
              "Expansion municipalities kept as controls", "Negative binomial model", "Municipalities with ≥ 100,000 inhabitants",
              "Aggregation by health region (449)", "+ age × state × year fixed effects", "+ age × region × month fixed effects", "+ both",
              "Northeast region", "Southeast region", "Central-West region", "South region", "North region"),
    s2 = "Table S2. Placebo seasons, each estimated against the three preceding seasons and with a common reference",
    s2h = c("Treated season", "Vaccination", "RRR (95% CI), pre = 3 seasons", "RRR (95% CI), common reference (Table 3A)"),
    s3 = "Table S3. Robustness of the 2025 dose–response in expansion municipalities: selection on prior epidemic, regional heterogeneity and expected 2024 gradient",
    s3h = c("Check", "Result"),
    s3rot = c("Spearman correlation of final coverage with the 10–14 share of admissions in 2023 / 2024",
              "Correlation of final coverage with the total admission rate in 2024 / with log population",
              "2025 gradient (reference)", "2025 gradient controlling for target group × 10–14 share in 2024",
              "2025 gradient controlling for target group × log(10–14 rate in 2024)",
              "2025 gradient in the Northeast / Southeast / South only", "2025 gradient in the Central-West / North only (few events)",
              "2025 gradient excluding each macro-region (min–max) / excluding São Paulo",
              "Effective 2024 coverage (admission-weighted) / final coverage in the expansion",
              "Expected 2024 gradient under 84.1% efficacy / observed"))
  fp <- function(p) ifelse(p < 0.001, paste0("< 0", dm, "001"), f(p, 3))

  # ---- Table 1 --------------------------------------------------------------------
  cmpr <- c("05-09", "15-19", "20-29")
  rt <- t1[, .(target = 1e5 * sum(admissions[agegroup == "10-14"]) / sum(pop_mean[agegroup == "10-14"]),
               cmp  = 1e5 * sum(admissions[agegroup %in% cmpr]) / sum(pop_mean[agegroup %in% cmpr])), by = .(group, year)]
  rt <- dcast(rt[, ratio := target / cmp], year ~ group, value.var = c("target", "cmp", "ratio"))
  tab1 <- data.table(year = rt$year, a1 = f(rt$`target_first phase`, 1), c1 = f(rt$`cmp_first phase`, 1), r1 = f(rt$`ratio_first phase`, 2),
                     a0 = f(rt$`target_never included`, 1), c0 = f(rt$`cmp_never included`, 1), r0 = f(rt$`ratio_never included`, 2))

  # ---- Table 2 --------------------------------------------------------------------
  t2 <- t2[match(names(TXT$outc), outcome)]
  tab2a <- data.table(TXT$outc[t2$outcome], f(t2$RRR), sprintf("%s–%s", f(t2$CI95_lwr), f(t2$CI95_upr)), fp(t2$p), f(t2$n_obs, 0))
  ag <- pl[level == "aggregate"][match(names(TXT$per), label_txt)]
  tab2b <- data.table(TXT$per[ag$label_txt], paste0(f(100 * ag$coverage_d1, 1), "%"), cint(ag$obs_RR, ag$obs_lwr, ag$obs_upr),
                      sprintf("%s / %s / %s", f(ag$bound_ef100), f(ag$bound_ef675), f(ag$bound_ef841)),
                      ifelse(ag$violates_ef100, TXT$sim[2], TXT$sim[1]))

  # ---- Table 3 --------------------------------------------------------------------
  lbl_year <- ifelse(ea$period >= "2024", paste(ea$period, TXT$vac), ea$period)
  lbl_tmp <- ifelse(et$period >= "2024", paste(et$period, TXT$vac), et$period)
  lbl_tmp[et$period == "2025/26"] <- paste("2025/26", TXT$until_dec, TXT$vac)
  tab3a <- data.table(lbl_year, cint(ea$RRR, ea$CI95_lwr, ea$CI95_upr), lbl_tmp, cint(et$RRR, et$CI95_lwr, et$CI95_upr))

  b <- rbindlist(lapply(names(TXT$samples), function(a) {
    x <- da[samp == a]
    data.table(TXT$samples[[a]], agegroup(x[year <= 2023, RR_10pp]), cint(x[year == 2024, RR_10pp], x[year == 2024, CI95_lwr], x[year == 2024, CI95_upr]),
               cint(x[year == 2025, RR_10pp], x[year == 2025, CI95_lwr], x[year == 2025, CI95_upr]), fp(x[year == 2025, p_calibrated_unilateral]),
               if (a == "expanded (2,230)") f(dp[year == 2025, p_unilateral], 3) else "—")
  }))
  b2 <- rbindlist(lapply(names(TXT$outc_b), function(d) {
    x <- dd[samp == "expanded (2,230)" & outcome == d]
    data.table(TXT$outc_b[[d]], "—", cint(x[year == 2024, RR_10pp], x[year == 2024, CI95_lwr], x[year == 2024, CI95_upr]),
               cint(x[year == 2025, RR_10pp], x[year == 2025, CI95_lwr], x[year == 2025, CI95_upr]), "—", "—")
  }))
  tab3b <- rbind(b, b2, use.names = FALSE)

  row_c <- function(x, lbl, adj, RRc, lo, hi)
    data.table(lbl, adj, agegroup(x[year <= 2023][[RRc]]), cint(x[year == 2024][[RRc]], x[year == 2024][[lo]], x[year == 2024][[hi]]),
               cint(x[year == 2025][[RRc]], x[year == 2025][[lo]], x[year == 2025][[hi]]))
  tab3c <- rbind(
    row_c(cna, TXT$c_f1, TXT$adj20, "RRR", "CI95_lwr", "CI95_upr"),
    row_c(cg[grepl("^first phase", analysis) & grepl("clean", analysis)], TXT$c_f1, TXT$adjl, "RR", "CI95_lwr", "CI95_upr"),
    row_c(cg[grepl("^Expanded", analysis) & grepl("script 20", analysis)], TXT$c_expn, TXT$adj20, "RR", "CI95_lwr", "CI95_upr"),
    row_c(cg[grepl("^Expanded", analysis) & grepl("clean", analysis)], TXT$c_expn, TXT$adjl, "RR", "CI95_lwr", "CI95_upr"),
    use.names = FALSE)

  # ---- Table S1 -------------------------------------------------------------------
  pick <- function(d, col, default) { x <- d[grepl(default, get(col))][1]; list(x$RRR, x$CI95_lwr, x$CI95_upr, x$p) }
  sources <- list(
    pick(t2, "outcome", "^Admissions"), pick(cl, "clustering", "^health region"), pick(cl, "clustering", "^UF$"),
    pick(s3, "analysis", "^Negative control"), pick(s3, "analysis", "^Placebo"), pick(s3, "analysis", "adults 20-59"),
    pick(s3, "analysis", "5-9 and 15-19"), pick(s3, "analysis", "^Adoption by D1 coverage"), pick(s3, "analysis", "^Excluding Feb"),
    pick(s3, "analysis", "^Including expanded municipalities"), pick(s3, "analysis", "^Expanded municipalities kept"),
    pick(s3, "analysis", "^Negative binomial"), pick(s3, "analysis", "100 thousand"), pick(s3, "analysis", "^Aggregated by health region"),
    pick(expd, "specification", "^\\+ agegroup\\^uf\\^year$"), pick(expd, "specification", "^\\+ agegroup\\^region\\^month$"),
    pick(expd, "specification", "^\\+ agegroup\\^uf\\^year and"), pick(s3, "analysis", "^Region Northeast$"), pick(s3, "analysis", "^Region Southeast$"),
    pick(s3, "analysis", "^Region Central-West$"), pick(s3, "analysis", "^Region South$"), pick(s3, "analysis", "^Region North$"))
  tabS1 <- rbindlist(lapply(seq_along(sources), function(i) {
    x <- sources[[i]]; data.table(TXT$s1rot[i], f(x[[1]]), sprintf("%s–%s", f(x[[2]]), f(x[[3]])), fp(x[[4]])) }))

  # ---- Table S2 -------------------------------------------------------------------
  p3 <- plc[grepl("^3 seasons", window)][season >= "2019/20"]
  p3 <- merge(p3, et[, .(season = period, RRc = RRR, loc = CI95_lwr, hic = CI95_upr)], by = "season")
  tabS2 <- data.table(p3$season, ifelse(p3$vaccine, TXT$sim[1], TXT$sim[2]), cint(p3$RRT, p3$CI95_lwr, p3$CI95_upr), cint(p3$RRc, p3$loc, p3$hic))

  # ---- Table S3 -------------------------------------------------------------------
  v <- function(default) rob[grepl(default, item)][1]
  vic <- function(default) { x <- v(default); if (is.na(x$val)) return("\u2014")
    m <- regmatches(x$observation, regexec("95% CI ([0-9.]+)-([0-9.]+)", x$observation))[[1]]
    cint(x$val, as.numeric(m[2]), as.numeric(m[3])) }
  excl <- rob[grepl("^2025 gradient, excluding (North|Northeast|Southeast|South|Central-West)$", item), val]
  tabS3 <- data.table(TXT$s3rot, c(
    sprintf("%s / %s", f(v("share 10-14 in 2023")$val), f(v("share 10-14 in 2024")$val)),
    sprintf("%s / %s", f(v("total rate in 2024")$val), f(v("log population")$val)),
    vic("expanded \\(reference\\)"), vic("controlling for target x share"), vic("controlling for target x log"),
    sprintf("%s / %s / %s", vic("Northeast only"), vic("Southeast only"), vic("South only")),
    sprintf("%s / %s", vic("Central-West only"), vic("North only")),
    sprintf("%s–%s / %s", f(min(excl)), f(max(excl)), vic("excluding Sao Paulo")),
    sprintf("%s%% / %s%%", f(100 * v("effective coverage 2024")$val, 1), f(100 * v("final coverage Dec/2025")$val, 1)),
    sprintf("%s / %s", f(v("EXPECTED")$val), vic("2024 gradient, expanded \\(observed\\)"))))

  # ---- Tables S4-S6, only if the intermediate tables exist -----------------
  TXT3 <- list(
    grp = c("first phase" = "First phase", "expanded" = "Expansion", "never included" = "Never included"),
    s4 = "Table S4. Dengue vaccine (TAK-003) coverage: first and second dose at ages 10–14 by municipality group, and first dose by birth cohort, December 2024 and December 2025, Brazil",
    s4a = "Panel A. Coverage at ages 10–14, population-weighted (doses through the end of the month)",
    s4ah = c("Municipality group", "D1 Dec 2024", "D1 Dec 2025", "D2 Dec 2024", "D2 Dec 2025"),
    s4b = "Panel B. D1 coverage by set of birth cohorts (cohort = year of dose − age; single-year population by uniform split of the 5-year age group)",
    s4bh = c("Cohorts", "Municipality group", "Dec 2024", "Dec 2025"),
    s5 = "Table S5. Dengue hospitalizations in SIH/SUS by principal diagnosis code and year of admission, Brazil, 2017–2025",
    s5h = c("Year", "A90", "A91", "A97", "Total", "% A97", "Ages 10–14"),
    s6 = "Table S6. Predominant serotype among probable dengue cases with an identified serotype, by state, 2023 and 2024, Brazil (* DF, RJ and BA highlighted)",
    s6h = c("State", "2023: typed", "2023: predominant (%)", "2024: typed", "2024: predominant (%)"))
  newer <- list(); HEADS_NEW <- list()
  fpath_m <- file.path(DIR$mod, c(cm = "table_coverage_d1_d2_municipal.csv", cc = "table_coverage_cohorts.csv",
                                sc = "table_sih_codes_year.csv", su = "table_serotype_uf.csv"))
  names(fpath_m) <- c("cm", "cc", "sc", "su")
  if (file.exists(fpath_m[["cm"]])) {
    x <- fread(fpath_m[["cm"]]); x <- dcast(x, group ~ month_ref, value.var = c("d1_weighted_pct", "d2_weighted_pct"))
    x <- x[match(names(TXT3$grp), group)]
    pc <- function(v) paste0(f(v, 1), "%")
    tabS4a <- data.table(TXT3$grp[x$group], pc(x$`d1_weighted_pct_2024-12`), pc(x$`d1_weighted_pct_2025-12`),
                         pc(x$`d2_weighted_pct_2024-12`), pc(x$`d2_weighted_pct_2025-12`))
    newer$tabS4a <- list(ttl = TXT3$s4, sub = TXT3$s4a); HEADS_NEW$tabS4a <- TXT3$s4ah
  }
  if (file.exists(fpath_m[["cc"]])) {
    x <- fread(fpath_m[["cc"]]); x <- dcast(x, cohorts + group ~ year_ref, value.var = "coverage_d1")
    x <- x[order(match(cohorts, unique(fread(fpath_m[["cc"]])$cohorts)), match(group, names(TXT3$grp)))]
    tabS4b <- data.table(x$cohorts, TXT3$grp[x$group], paste0(f(100 * x$`2024`, 1), "%"), paste0(f(100 * x$`2025`, 1), "%"))
    newer$tabS4b <- list(ttl = if (is.null(newer$tabS4a)) TXT3$s4 else NULL, sub = TXT3$s4b); HEADS_NEW$tabS4b <- TXT3$s4bh
  }
  if (file.exists(fpath_m[["sc"]])) {
    x <- fread(fpath_m[["sc"]])
    tabS5 <- data.table(x$year, f(x$A90, 0), f(x$A91, 0), f(x$A97, 0), f(x$total, 0), f(x$pct_A97, 1), f(x$total_10a14, 0))
    newer$tabS5 <- list(ttl = TXT3$s5, sub = NULL); HEADS_NEW$tabS5 <- TXT3$s5h
  }
  if (file.exists(fpath_m[["su"]])) {
    x <- fread(fpath_m[["su"]])
    pr <- function(a) { y <- x[year == a]; setNames(sprintf("%s (%s%%)", y$predominant, f(100 * y$f_predominant, 0)), y$abbr) }
    tp <- function(a) { y <- x[year == a]; setNames(f(y$typed, 0), y$abbr) }
    ufs <- sort(unique(x$abbr)); dest <- unique(x[highlight == TRUE, abbr])
    tabS6 <- data.table(ifelse(ufs %in% dest, paste0(ufs, "*"), ufs), tp(2023)[ufs], pr(2023)[ufs], tp(2024)[ufs], pr(2024)[ufs])
    for (j in names(tabS6)) set(tabS6, which(is.na(tabS6[[j]])), j, "—")
    newer$tabS6 <- list(ttl = TXT3$s6, sub = NULL); HEADS_NEW$tabS6 <- TXT3$s6h
  }

  # ---- reconciled Table 1; Tables S8-S10 --------------------------------
  TXT4 <- list(
    t1h_extra = c("First-phase\ncomparison 5–59", "Ratio\n(5–59)", "Never-included\ncomparison 5–59", "Ratio\n(5–59)",
                  "Crude DDD\n(5–9, 15–19, 20–29)", "Crude DDD\n(5–59)"),
    t1n_extra = " Comparison 5–59: all ages 5–59 except 10–14. Crude DDD = (first-phase ratio in the year / 2017–2023 mean) ÷ (the same in never-included municipalities).",
    s8 = "Table S8. Estimates against the smallest value compatible with the vaccine coverage achieved (plausibility bound), Brazil, 2024–2025",
    s8h = c("Analysis", "Period", "Estimate (95% CI)", "Coverage used", "Effective coverage", "Bound: efficacy 100% / 74.6% / 84.1% / 87.9%", "Below the 100% bound: point / whole CI"),
    s9 = "Table S9. Share of dengue hospitalizations at ages 5–59 that occurred at ages 10–14, by municipality group, seasonal window and year, Brazil, 2017–2025",
    s9h = c("Window", "Year", "First phase, % (95% CI)", "Never included, % (95% CI)", "Ratio (95% CI)"),
    s9n = "95% CI by municipality bootstrap (500 resamples). June–February 2025 runs only to December 2025.",
    s10 = "Table S10. Negative controls in SIH/SUS: binary model (first phase × never included) and 2025 gradient per 10 p.p. of coverage in expansion municipalities, Brazil, 2017–2025",
    s10h = c("Outcome", "Binary, pooled", "Binary, 2024", "Binary, 2025", "Gradient, 2024", "Gradient, 2025"),
    windows = c("jan-may" = "Jan–May", "jun-feb" = "Jun–Feb", "apr-jul" = "Apr–Jul"), yes_no = c("yes", "no"))
  T1H <- TXT$t1h
  fpath_rec <- file.path(DIR$mod, "table_s1_reconciliation.csv")
  if (file.exists(fpath_rec)) {
    rec <- fread(fpath_rec)[match(tab1$year, year)]
    tab1[, `:=`(ct1 = f(rec$rt_cmpr_all_first_phase, 1), rt1 = f(rec$ratio_cmpr_all_first_phase, 2),
                ct0 = f(rec$rt_cmpr_all_never, 1), rt0 = f(rec$ratio_cmpr_all_never, 2),
                d3 = f(rec$ddd_crude_cmpr3, 2), dt = f(rec$ddd_crude_cmpr_all, 2))]
    T1H <- c(TXT$t1h, TXT4$t1h_extra); TXT$t1n <- paste0(TXT$t1n, TXT4$t1n_extra)
  }
  sn <- function(x) ifelse(is.na(x), "—", ifelse(x, TXT4$yes_no[1], TXT4$yes_no[2]))
  fpath_r2 <- file.path(DIR$mod, c(tt = "table_estimate_vs_ceiling.csv", ce = "table_composition_age_windows.csv",
                                 cb = "table_controls_sih_binary.csv", cg = "table_controls_sih_gradient.csv"))
  names(fpath_r2) <- c("tt", "ce", "cb", "cg")
  if (file.exists(fpath_r2[["tt"]])) {
    x <- fread(fpath_r2[["tt"]])
    tabS8 <- data.table(x$analysis, x$period, cint(x$estimate, x$cint_lwr, x$cint_upr), x$coverage_used, paste0(f(100 * x$coverage_effective, 1), "%"),
                        sprintf("%s / %s / %s / %s", f(x$ceiling_ve100), f(x$ceiling_ve746), f(x$ceiling_ve841), f(x$ceiling_ve879)),
                        sprintf("%s / %s", sn(x$violates_point_ve100), sn(x$violates_cint_ve100)))
    newer$tabS8 <- list(ttl = TXT4$s8, sub = NULL); HEADS_NEW$tabS8 <- TXT4$s8h
  }
  if (file.exists(fpath_r2[["ce"]])) {
    x <- fread(fpath_r2[["ce"]])
    pct <- function(a, lo, hi) sprintf("%s (%s–%s)", f(100 * a, 1), f(100 * lo, 1), f(100 * hi, 1))
    tabS9 <- data.table(TXT4$windows[x$window], paste0(x$year, ifelse(x$complete, "", "*")), pct(x$frac_first_phase, x$frac_first_phase_lwr, x$frac_first_phase_upr),
                        pct(x$frac_never, x$frac_never_lwr, x$frac_never_upr), cint(x$ratio, x$ratio_lwr, x$ratio_upr))
    newer$tabS9 <- list(ttl = TXT4$s9, sub = TXT4$s9n); HEADS_NEW$tabS9 <- TXT4$s9h
  }
  if (file.exists(fpath_r2[["cb"]]) && file.exists(fpath_r2[["cg"]])) {
    xb <- fread(fpath_r2[["cb"]]); xg <- fread(fpath_r2[["cg"]])
    icb <- function(d) if (nrow(d) && !is.na(d$RRR)) cint(d$RRR, d$CI95_lwr, d$CI95_upr) else "—"
    icg <- function(d) if (nrow(d) && !is.na(d$RR_10pp)) cint(d$RR_10pp, d$CI95_lwr, d$CI95_upr) else "—"
    tabS10 <- rbindlist(lapply(unique(xb$outcome), function(dsf) data.table(
      dsf, icb(xb[outcome == dsf & period == "pooled Feb/2024-Dec/2025"]), icb(xb[outcome == dsf & period == "2024"]),
      icb(xb[outcome == dsf & period == "2025"]), icg(xg[outcome == dsf & year == 2024L]), icg(xg[outcome == dsf & year == 2025L]))))
    newer$tabS10 <- list(ttl = TXT4$s10, sub = NULL); HEADS_NEW$tabS10 <- TXT4$s10h
  }

  # ---- docx + csv -----------------------------------------------------------------
  HEADS <- list(tab1 = T1H, tab2a = TXT$t2ah, tab2b = TXT$t2bh, tab3a = TXT$t3ah, tab3b = TXT$t3bh, tab3c = TXT$t3ch,
              tabS1 = TXT$s1h, tabS2 = TXT$s2h, tabS3 = TXT$s3h)
  HEADS <- c(HEADS, HEADS_NEW)
  ft <- function(nm) {                       # headers may repeat ("Razão"), hence labels rather than names
    d <- get(nm); setnames(d, paste0("c", seq_len(ncol(d))))
    wide <- ncol(d) > 9                      # Table 1 with reconciliation: smaller font, narrow first column
    w1 <- if (wide) 0.45 else if (ncol(d) <= 4) 3 else 2
    flextable(d) |> set_header_labels(values = setNames(as.list(HEADS[[nm]]), names(d))) |>
      font(fontname = "Times New Roman", part = "all") |> fontsize(size = if (wide) 7 else 9, part = "all") |>
      bold(part = "header") |> align(j = 2:ncol(d), align = "center", part = "all") |>
      width(j = 1, width = w1) |>
      width(j = 2:ncol(d), width = (6.3 - w1) / (ncol(d) - 1)) |>
      set_table_properties(layout = "fixed")
  }
  par <- function(doc, x, bold = FALSE, size = 10)
    body_add_fpar(doc, fpar(ftext(x, fp_text(font.family = "Times New Roman", font.size = size, bold = bold))))
  doc <- read_docx() |>
    par(TXT$t1, TRUE) |> body_add_flextable(ft("tab1")) |> par(TXT$t1n, size = 8) |> body_add_break() |>
    par(TXT$t2, TRUE) |> par(TXT$t2a, size = 9) |> body_add_flextable(ft("tab2a")) |> par(TXT$t2b, size = 9) |>
      body_add_flextable(ft("tab2b")) |> par(TXT$t2n, size = 8) |> body_add_break() |>
    par(TXT$t3, TRUE) |> par(TXT$t3a, size = 9) |> body_add_flextable(ft("tab3a")) |> par(TXT$t3b, size = 9) |>
      body_add_flextable(ft("tab3b")) |> par(TXT$t3c, size = 9) |> body_add_flextable(ft("tab3c")) |> par(TXT$t3n, size = 8) |> body_add_break() |>
    par(TXT$s1, TRUE) |> body_add_flextable(ft("tabS1")) |> body_add_break() |>
    par(TXT$s2, TRUE) |> body_add_flextable(ft("tabS2")) |> body_add_break() |>
    par(TXT$s3, TRUE) |> body_add_flextable(ft("tabS3"))
  for (nm in names(newer)) {                  # S4-S10
    if (!is.null(newer[[nm]]$ttl)) doc <- doc |> body_add_break() |> par(newer[[nm]]$ttl, TRUE)
    if (!is.null(newer[[nm]]$sub)) doc <- doc |> par(newer[[nm]]$sub, size = 9)
    doc <- doc |> body_add_flextable(ft(nm))
  }
  print(doc, target = file.path(DIR$tab, "tables.docx"))
  for (nm in names(HEADS)) {
    d <- copy(get(nm)); setnames(d, HEADS[[nm]])
    fwrite(d, file.path(DIR$tab, paste0(sub("^tab", "table", nm), ".csv")), bom = TRUE)
  }
  log_msg("tables.docx written")
}
