# =============================================================================
# 33_dose_response.R — dose-response by vaccine coverage (SI-PNI)
#
# A. CROSS-SECTIONAL (main): FINAL D1 coverage in the target age group (Dec/2025),
#    fixed per municipality, × target age group × year, with the fixed effects of
#    the main model; each coefficient re-anchored on the 2017-2023 mean. In years
#    without the vaccine, future coverage should predict nothing; in 2024 a weak
#    gradient is expected (low coverage during the epidemic); in 2025, the full one.
#    Sample of interest: the 2,230 EXPANSION municipalities (highly heterogeneous
#    coverage, selection independent of the 2023 epidemic). Also: all
#    municipalities, official list, first phase, main sample; SINAN outcomes and
#    chikungunya (negative control); coverage categories; permutation of
#    coverage across municipalities (within population deciles). -> Table 3B, Fig. 2
# B. ROBUSTNESS of the 2025 gradient in the expansion: correlation of coverage with
#    the prior epidemic in the target age group, explicit control for it, exclusion
#    of each region/the state of São Paulo, and expected gradient in 2024. -> Table S3
# C. WITHIN each year, only the 521 first-phase municipalities (lagged monthly
#    coverage), with falsification (2025 coverage transposed to 2019 and 2023).
# D. LONGITUDINAL (lagged monthly coverage × target age group, main sample):
#    D1 alone, categories and D1 + D2 (exploratory).
# Fixed seed (20260923) for the permutation. Outputs in outputs/models/.
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(fixest)

N_PERM <- as.integer(Sys.getenv("N_PERM_DOSE", "300"))
NEVER  <- 10000L
SCALE <- 0.10                                   # effect per 10 percentage points
set.seed(20260923)

panel <- readRDS(file.path(DIR$processed, "panel.rds"))
panel[, cov_final := max(cov_d1_l1), by = code_muni]      # cumulative coverage in Dec/2025
panel[, `:=`(phase1 = as.integer(g < NEVER), lst = as.integer(g < NEVER | expn == 1L))]
years <- PARAM$year_start:PARAM$year_end

# keep only municipality×month cells with some outcome event (faster; same estimate)
prep <- function(d, y = "admissions") {
  d <- copy(d)[, tot_mt := sum(get(y)), by = .(code_muni, t)][tot_mt > 0]
  d[, dose := target * cov_final][]
}
SAMPLES <- list("expanded (2,230)" = quote(expn == 1L), "all municipalities" = quote(TRUE),
                 "official list (first phase + expanded)" = quote(lst == 1L), "first phase (521)" = quote(phase1 == 1L),
                 "main sample (521 + never)" = quote(expn == 0L))

# ---- A1. coefficient by year, re-anchored on 2017-2023, by sample ------------
annual <- list()
for (nm in names(SAMPLES)) {
  d <- prep(panel[eval(SAMPLES[[nm]])])
  m <- fepois(as.formula(sprintf("admissions ~ i(year, dose, ref = %d) | %s", years[1], FE_MAIN)),
              data = d, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
  tab <- reanchor_by_year(m, "dose", scale = SCALE)
  tab[, `:=`(samp = nm, municipalities = uniqueN(d$code_muni),
             admissions_target_2025 = d[target == 1L & year == 2025L, sum(admissions)])]
  annual[[nm]] <- tab
  cat("\n== ", nm, "\n", sep = ""); print(tab[, .(year, RR_10pp = round(RR, 3), CI95_lwr = round(CI95_lwr, 3),
                                                   CI95_upr = round(CI95_upr, 3), p_calibrated_unilateral = signif(p_calibrated_unilateral, 2))])
}
annual <- rbindlist(annual)
setnames(annual, "RR", "RR_10pp")
fwrite(annual[, .(samp, year, log_b, se, RR_10pp, CI95_lwr, CI95_upr, t_calibrated, p_calibrated_unilateral,
                 municipalities, admissions_target_2025)], file.path(DIR$mod, "table_crosssec_dose_annual.csv"))

# ---- A2. only 2024 and 2025 (reference: other years), by outcome -------------------
OUTC <- c(admissions = "Admissions (SIH)", den_cases = "Probable cases (SINAN)",
          den_hosp = "Hospitalized (SINAN)", den_severe = "Warning signs/severe (SINAN)",
          chik_cases = "Negative control: chikungunya (SINAN)")
two <- list()
for (nm in names(SAMPLES)[1:2]) for (y in names(OUTC)) {
  d <- prep(panel[eval(SAMPLES[[nm]])], y)
  d[, `:=`(dose24 = dose * (year == 2024L), dose25 = dose * (year == 2025L))]
  m <- fepois(as.formula(sprintf("%s ~ dose24 + dose25 | %s", y, FE_MAIN)), data = d,
              offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
  for (k in c("dose24", "dose25")) {
    r <- rr_row(m, k, nm, SCALE)
    two[[paste(nm, y, k)]] <- data.table(samp = nm, outcome = OUTC[[y]], year = if (k == "dose24") 2024L else 2025L,
                                          RR_10pp = r$RR, CI95_lwr = r$CI95_lwr, CI95_upr = r$CI95_upr, p = r$p, n_obs = r$n_obs)
  }
}
two <- rbindlist(two)
fwrite(two, file.path(DIR$mod, "table_crosssec_dose_2024_2025.csv"))
cat("\n== 2024 and 2025 by outcome\n"); print(two[, .(samp, outcome, year, RR_10pp = round(RR_10pp, 3), CI95 = sprintf("%.3f-%.3f", CI95_lwr, CI95_upr), p = signif(p, 2))])

# ---- A3. final coverage categories (expansion; ref. < 10%) --------------------
d <- prep(panel[expn == 1L])
d[, cat := cut(cov_final, c(-Inf, 0.10, 0.25, 0.40, 0.55, Inf), labels = c("<10%", "10-25%", "25-40%", "40-55%", ">=55%"), right = FALSE)]
d[, `:=`(cat24 = fifelse(year == 2024L & target == 1L, as.character(cat), "ref"),
         cat25 = fifelse(year == 2025L & target == 1L, as.character(cat), "ref"))]
m_cat <- fepois(as.formula(sprintf("admissions ~ i(cat24, ref = c('ref', '<10%%')) + i(cat25, ref = c('ref', '<10%%')) | %s", FE_MAIN)),
                data = d, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
ct <- as.data.table(coeftable(m_cat), keep.rownames = "term")
setnames(ct, c("Estimate", "Std. Error"), c("b", "se"))
ct <- ct[, .(year = fifelse(grepl("^cat24", term), 2024L, 2025L), category = sub("^cat2[45]::", "", term),
             RR = exp(b), CI95_lwr = exp(b - 1.96 * se), CI95_upr = exp(b + 1.96 * se), p = ct[[5]])]
ct <- merge(ct, d[target == 1L & year == 2025L, .(municipalities = uniqueN(code_muni)), by = .(category = as.character(cat))],
            by = "category", all.x = TRUE)[order(year, category)]
fwrite(ct, file.path(DIR$mod, "table_crosssec_dose_categories.csv"))
cat("\n== coverage categories (expanded)\n"); print(ct)

# ---- A4. permutation of coverage across expansion municipalities --------------------
d <- prep(panel[expn == 1L])
mun <- panel[expn == 1L & year == 2022L, .(pop = sum(pop) / 12, cov_final = cov_final[1]), by = code_muni]
mun <- mun[code_muni %in% unique(d$code_muni)]
mun[, dec := cut(pop, unique(quantile(pop, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
estfit_dose <- function(cf) {
  dd <- merge(d[, -"cov_final"], cf, by = "code_muni")
  dd[, `:=`(dose24 = target * cov_final * (year == 2024L), dose25 = target * cov_final * (year == 2025L))]
  m <- fepois(as.formula(sprintf("admissions ~ dose24 + dose25 | %s", FE_MAIN)), data = dd,
              offset = ~log(pop), vcov = "iid", notes = FALSE)
  coef(m)[c("dose24", "dose25")]
}
real <- estfit_dose(mun[, .(code_muni, cov_final)])
fpath <- file.path(DIR$mod, "permutation_crosssec_dose.csv")
done <- if (file.exists(fpath)) nrow(fread(fpath)) else 0L
if (done < N_PERM) for (i in (done + 1):N_PERM) {
  cf <- copy(mun)[, cov_final := sample(cov_final), by = dec][, .(code_muni, cov_final)]
  r <- tryCatch(estfit_dose(cf), error = function(e) c(NA_real_, NA_real_))
  fwrite(data.table(perm = i, b2024 = r[[1]], b2025 = r[[2]]), fpath, append = file.exists(fpath))
  if (i %% 25 == 0) log_msg("coverage permutation: ", i, "/", N_PERM)
}
pl <- fread(fpath)[!is.na(b2025)]
perm <- rbindlist(lapply(c(2024L, 2025L), function(a) {
  x <- pl[[if (a == 2024L) "b2024" else "b2025"]]; rv <- real[[if (a == 2024L) "dose24" else "dose25"]]
  data.table(year = a, n_draws = length(x), RR_10pp_real = exp(SCALE * rv),
             RR_10pp_placebo_p05 = exp(SCALE * quantile(x, .05)), RR_10pp_placebo_p95 = exp(SCALE * quantile(x, .95)),
             sd_log_placebos = sd(x), p_unilateral = (1 + sum(x <= rv)) / (1 + length(x)))
}))
fwrite(perm, file.path(DIR$mod, "table_crosssec_dose_permutation.csv"))
cat("\n== coverage permutation (expanded)\n"); print(perm)

# ---- A5. gradient with DEC/2024 coverage instead of the final one -------------
# Same convention as the final coverage (cov_d1_l1 in the month: cumulative up to the previous month),
# taken in Dec/2024 — the coverage with which the expansion municipalities enter the 2025
# season. Same model (A1), same re-anchoring and calibrated p; permutation with the SAME
# draws as A4: seed reset and the same `mun` table (same order and deciles), changing
# only the permuted value. Placed after A4 so as not to alter its random sequence.
cov24 <- panel[expn == 1L & month == as.Date("2024-12-01"), .(cov_dec24 = cov_d1_l1[1]), by = code_muni]
d5 <- prep(panel[expn == 1L])
d5[, cov_dec24 := cov24$cov_dec24[match(code_muni, cov24$code_muni)]][, dose := target * cov_dec24]
m5 <- fepois(as.formula(sprintf("admissions ~ i(year, dose, ref = %d) | %s", years[1], FE_MAIN)),
             data = d5, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
a5 <- reanchor_by_year(m5, "dose", scale = SCALE)
setnames(a5, "RR", "RR_10pp")
a5[, `:=`(exposure = "D1 coverage Dec/2024", municipalities = uniqueN(d5$code_muni),
          cov_dec24_median = median(cov24$cov_dec24), cov_final_median = median(mun$cov_final))]
fwrite(a5[, .(exposure, year, log_b, se, RR_10pp, CI95_lwr, CI95_upr, t_calibrated, p_calibrated_unilateral,
              p_calibrated_twosided, municipalities, cov_dec24_median, cov_final_median)],
       file.path(DIR$mod, "table_crosssec_dose_cov_dec2024.csv"))
cat("\n== gradient with Dec/2024 coverage (expanded)\n")
print(a5[, .(year, RR_10pp = round(RR_10pp, 3), CI95_lwr = round(CI95_lwr, 3), CI95_upr = round(CI95_upr, 3),
             p_calibrated_unilateral = signif(p_calibrated_unilateral, 2))])
mun24 <- copy(mun)[, cov_dec24 := cov24$cov_dec24[match(code_muni, cov24$code_muni)]]
real24 <- estfit_dose(mun24[, .(code_muni, cov_final = cov_dec24)])
fpath24 <- file.path(DIR$mod, "permutation_crosssec_dose_cov_dec2024.csv")
done24 <- if (file.exists(fpath24)) nrow(fread(fpath24)) else 0L
set.seed(20260923)
if (done24 > 0 && done24 < N_PERM) log_msg("WARNING: resuming the Dec/2024 permutation; the draws no longer match those of A4")
if (done24 < N_PERM) for (i in (done24 + 1):N_PERM) {
  cf <- copy(mun24)[, cov_final := sample(cov_dec24), by = dec][, .(code_muni, cov_final)]
  r <- tryCatch(estfit_dose(cf), error = function(e) c(NA_real_, NA_real_))
  fwrite(data.table(perm = i, b2024 = r[[1]], b2025 = r[[2]]), fpath24, append = file.exists(fpath24))
  if (i %% 25 == 0) log_msg("Dec/2024 coverage permutation: ", i, "/", N_PERM)
}
pl24 <- fread(fpath24)[!is.na(b2025)]
perm24 <- rbindlist(lapply(c(2024L, 2025L), function(a) {
  x <- pl24[[if (a == 2024L) "b2024" else "b2025"]]; rv <- real24[[if (a == 2024L) "dose24" else "dose25"]]
  data.table(exposure = "D1 coverage Dec/2024", year = a, n_draws = length(x), RR_10pp_real = exp(SCALE * rv),
             RR_10pp_placebo_p05 = exp(SCALE * quantile(x, .05)), RR_10pp_placebo_p95 = exp(SCALE * quantile(x, .95)),
             sd_log_placebos = sd(x), p_unilateral = (1 + sum(x <= rv)) / (1 + length(x)))
}))
fwrite(perm24, file.path(DIR$mod, "table_crosssec_dose_cov_dec2024_permutation.csv"))
cat("\n== Dec/2024 coverage permutation (expanded)\n"); print(perm24)

# ---- B. robustness of the 2025 gradient in the expansion ---------------------------------
expn <- panel[expn == 1L]
an <- expn[, .(i = sum(admissions), pop = sum(pop) / 12), by = .(code_muni, year, target)]
an <- dcast(an, code_muni + year ~ target, value.var = c("i", "pop"))
setnames(an, c("i_0", "i_1", "pop_0", "pop_1"), c("i0", "i1", "p0", "p1"))
an[, `:=`(part = i1 / (i0 + i1), rt_target = 1e5 * i1 / p1, rt_tot = 1e5 * (i0 + i1) / (p0 + p1))]
mun <- expn[, .(cov = cov_final[1], pop = max(pop)), by = code_muni]
rob <- list()
for (a in c(2023L, 2024L)) {
  x <- merge(an[year == a & (i0 + i1) >= 20], mun, by = "code_muni")
  rob[[paste("part", a)]] <- data.table(item = sprintf("Spearman final coverage x share 10-14 in %d", a),
                                        val = cor(x$cov, x$part, method = "spearman"),
                                        observation = sprintf("n = %d municipalities with >= 20 admissions", nrow(x)))
  rob[[paste("rt", a)]] <- data.table(item = sprintf("Spearman final coverage x total rate in %d", a),
                                      val = cor(x$cov, x$rt_tot, method = "spearman"), observation = "")
}
rob$pop <- data.table(item = "Spearman final coverage x log population", val = cor(mun$cov, log(mun$pop), method = "spearman"), observation = "")

d <- prep(expn)
d[, `:=`(dose24 = target * cov_final * (year == 2024L), dose25 = target * cov_final * (year == 2025L))]
d <- merge(d, an[year == 2024L, .(code_muni, part24 = part, rt24 = rt_target)], by = "code_muni", all.x = TRUE)
d[is.na(part24), part24 := median(an[year == 2024L, part], na.rm = TRUE)][is.na(rt24), rt24 := 0]
d[, `:=`(ctrl_part = target * part24 * (year == 2025L), ctrl_rt = target * log1p(rt24) * (year == 2025L))]
SPEC <- list(list(dd = quote(d), extra = "", lbl = "2025 gradient, expanded (reference)"),
              list(dd = quote(d), extra = " + ctrl_part", lbl = "2025 gradient controlling for target x share 10-14 in 2024"),
              list(dd = quote(d), extra = " + ctrl_rt", lbl = "2025 gradient controlling for target x log(rate 10-14 in 2024)"))
for (r in intersect(REGION_ORDER, unique(d$region))) SPEC <- c(SPEC,
  list(list(dd = bquote(d[region == .(r)]), extra = "", lbl = sprintf("2025 gradient, %s only", r)),
       list(dd = bquote(d[region != .(r)]), extra = "", lbl = sprintf("2025 gradient, excluding %s", r))))
SPEC <- c(SPEC, list(list(dd = quote(d[uf != "35"]), extra = "", lbl = "2025 gradient, excluding Sao Paulo (state)")))
for (e in SPEC) {
  dd <- eval(e$dd)
  m <- tryCatch(fepois(as.formula(sprintf("admissions ~ dose24 + dose25%s | %s", e$extra, FE_MAIN)), data = dd,
                       offset = ~log(pop), cluster = ~code_muni, notes = FALSE), error = function(err) NULL)
  if (is.null(m)) next
  r25 <- rr_row(m, "dose25", e$lbl, SCALE)
  rob[[e$lbl]] <- data.table(item = e$lbl, val = r25$RR,
                             observation = sprintf("95%% CI %.3f-%.3f; municipalities = %d", r25$CI95_lwr, r25$CI95_upr, uniqueN(dd$code_muni)))
  if (e$lbl == "2025 gradient, expanded (reference)") {
    r24 <- rr_row(m, "dose24", e$lbl, SCALE)
    rob[["obs24"]] <- data.table(item = "2024 gradient, expanded (observed)", val = r24$RR,
                                 observation = sprintf("95%% CI %.3f-%.3f", r24$CI95_lwr, r24$CI95_upr))
  }
}
# expected gradient in 2024: EFFECTIVE 2024 coverage (weighted by hospitalizations) / final coverage
ef <- expn[target == 1L & year == 2024L, .(cov = sum(cov_d1_l1 * pop) / sum(pop), i = sum(admissions)), by = month]
cov_ef24 <- sum(ef$cov * ef$i) / sum(ef$i)
end <- expn[target == 1L & month == max(month)]; cov_fin <- sum(end$cov_d1_l1 * end$pop) / sum(end$pop)
rob$ef24 <- data.table(item = "effective coverage 2024, expanded (admission-weighted)", val = cov_ef24, observation = "")
rob$fin  <- data.table(item = "final coverage Dec/2025, expanded (population-weighted)", val = cov_fin, observation = "")
rob$expd  <- data.table(item = "2024 gradient EXPECTED per 10 p.p. under 84.1% efficacy", val = 1 - 0.841 * SCALE * cov_ef24 / cov_fin,
                       observation = "1 - VE x 0.10 x (effective/final); compare with the observed")
rob <- rbindlist(rob)
fwrite(rob, file.path(DIR$mod, "table_robustness_crosssec_dose.csv"))
cat("\n== robustness (expanded)\n"); print(rob)

# ---- C. within each year, only the 521 first-phase municipalities ------------------
tr1 <- panel[expn == 0L & g < NEVER]                   # the 521 first-phase ones
intra <- list()
fit_intra <- function(dd, label_txt) {                    # dd has `cov` and `dose`
  m <- tryCatch(fepois(as.formula(sprintf("admissions ~ dose | %s", FE_MAIN)), data = dd,
                       offset = ~log(pop), cluster = ~code_muni, notes = FALSE), error = function(e) NULL)
  q <- quantile(dd[target == 1L, cov], c(.25, .5, .75), na.rm = TRUE)
  base <- data.table(analysis = label_txt, municipalities = uniqueN(dd$code_muni), admissions_target = dd[target == 1L, sum(admissions)],
                     cov_p25 = q[[1]], cov_p50 = q[[2]], cov_p75 = q[[3]])
  if (is.null(m)) return(cbind(base, data.table(RR_10pp = NA_real_, CI95_lwr = NA_real_, CI95_upr = NA_real_, p = NA_real_, n_obs = NA_integer_)))
  r <- rr_row(m, "dose", label_txt, SCALE)
  cbind(base, data.table(RR_10pp = r$RR, CI95_lwr = r$CI95_lwr, CI95_upr = r$CI95_upr, p = r$p, n_obs = r$n_obs))
}
for (a in c(2025L, 2024L)) {
  dd <- tr1[year == a][, `:=`(cov = cov_d1_l1, dose = target * cov_d1_l1)]
  intra[[as.character(a)]] <- fit_intra(dd, sprintf("Season %d (treated municipalities only)", a))
}
dd <- tr1[year %in% 2024:2025][, `:=`(cov = cov_d1_l1, dose = target * cov_d1_l1)]
intra$joint <- fit_intra(dd, "2024 and 2025 pooled (with temporal confounding)")
traj <- unique(tr1[year == 2025L, .(code_muni, m = as.integer(format(month, "%m")), cov_fake = cov_d1_l1)])
for (a in c(2019L, 2023L)) {                                # falsification: 2025 trajectory transposed
  dd <- tr1[year == a][, m := as.integer(format(month, "%m"))]
  dd <- merge(dd, traj, by = c("code_muni", "m"), all.x = TRUE)
  dd[is.na(cov_fake), cov_fake := 0][, `:=`(cov = cov_fake, dose = target * cov_fake)]
  intra[[paste0("fake_", a)]] <- fit_intra(dd, sprintf("Falsification: 2025 coverage applied to %d", a))
}
intra <- rbindlist(intra, fill = TRUE)
fwrite(intra, file.path(DIR$mod, "table_dose_response_withinseason.csv"))
cat("\n== within each year (first phase)\n"); print(intra[, .(analysis, municipalities, admissions_target, RR_10pp = round(RR_10pp, 3), CI95 = sprintf("%.3f-%.3f", CI95_lwr, CI95_upr))])

# ---- D. longitudinal: lagged monthly coverage × target age group (main sample)
p <- panel[expn == 0L]
p[, `:=`(cov1_target = target * cov_d1_l1, cov2_target = target * cov_d2_l1)]
m1 <- fepois(as.formula(sprintf("admissions ~ cov1_target | %s", FE_MAIN)), data = p, offset = ~log(pop), cluster = ~code_muni)
cov_med <- p[target == 1L & g < NEVER & month >= as.Date("2024-12-01"), mean(cov_d1_l1)]
p[, cov_cat := fifelse(cov1_target <= 0, "0 (no coverage)", fifelse(cov1_target < 0.20, "1 (<20%)",
                fifelse(cov1_target < 0.50, "2 (20-50%)", "3 (>=50%)")))]
m2 <- fepois(as.formula(sprintf("admissions ~ i(cov_cat, ref = '0 (no coverage)') | %s", FE_MAIN)),
             data = p, offset = ~log(pop), cluster = ~code_muni)
m3 <- fepois(as.formula(sprintf("admissions ~ cov1_target + cov2_target | %s", FE_MAIN)), data = p, offset = ~log(pop), cluster = ~code_muni)
ct2 <- as.data.table(coeftable(m2), keep.rownames = "term"); setnames(ct2, c("Estimate", "Std. Error"), c("b", "se"))
lon <- rbind(
  rr_row(m1, "cov1_target", "D1 alone — per 10 p.p. of coverage", 0.10)[, .(model = analysis, term = "cov1_target", scale = 0.10, RR, CI95_lwr, CI95_upr, p)],
  rr_row(m1, "cov1_target", sprintf("D1 alone — at the mean observed coverage (%.0f%%)", 100 * cov_med), cov_med)[, .(model = analysis, term = "cov1_target", scale = cov_med, RR, CI95_lwr, CI95_upr, p)],
  ct2[, .(model = "Coverage categories (ref. no coverage)", term = sub(".*::", "", term), scale = 1,
          RR = exp(b), CI95_lwr = exp(b - 1.96 * se), CI95_upr = exp(b + 1.96 * se), p = ct2[[5]])],
  rr_row(m3, "cov1_target", "D1+D2 (exploratory) — D1 per 10 p.p.", 0.10)[, .(model = analysis, term = "cov1_target", scale = 0.10, RR, CI95_lwr, CI95_upr, p)],
  rr_row(m3, "cov2_target", "D1+D2 (exploratory) — D2 per 10 p.p.", 0.10)[, .(model = analysis, term = "cov2_target", scale = 0.10, RR, CI95_lwr, CI95_upr, p)])
fwrite(lon, file.path(DIR$mod, "table7_dose_response.csv"))
saveRDS(list(d1 = m1, cat = m2, d1d2 = m3), file.path(DIR$mod, "models_dose_response.rds"))
cat("\n== longitudinal\n"); print(lon[, .(model, term, RR = round(RR, 3), CI95_lwr = round(CI95_lwr, 3), CI95_upr = round(CI95_upr, 3))])

# ---- SIH negative controls: cross-sectional gradient in the expansion -----------
# same specification as A2 (two indicators, 2024 and 2025, per 10 p.p. of final coverage)
FILE_CTRL <- file.path(DIR$processed, "sih_controls_aggregate.rds")
CTRL <- c(ctrl_all = "All causes except A90/A91 and chapters O, P, Z", ctrl_a00a09 = "A00-A09 (intestinal infectious diseases)",
          ctrl_j00j22 = "J00-J22 (acute respiratory infections)", ctrl_a920 = "A92.0 (chikungunya)")
if (file.exists(FILE_CTRL)) {
  pc <- merge(panel[expn == 1L], readRDS(FILE_CTRL), by = c("code_muni", "agegroup", "month"), all.x = TRUE)
  for (v in names(CTRL)) set(pc, which(is.na(pc[[v]])), v, 0L)
  grc <- rbindlist(lapply(names(CTRL), function(y) {
    d <- prep(pc, y)
    d[, `:=`(dose24 = dose * (year == 2024L), dose25 = dose * (year == 2025L))]
    m <- tryCatch(fepois(as.formula(sprintf("%s ~ dose24 + dose25 | %s", y, FE_MAIN)), data = d,
                         offset = ~log(pop), cluster = ~code_muni, notes = FALSE), error = function(e) NULL)
    if (is.null(m)) return(data.table(outcome = CTRL[[y]], year = c(2024L, 2025L)))
    rbindlist(lapply(c(dose24 = 2024L, dose25 = 2025L), function(a) {
      r <- rr_row(m, if (a == 2024L) "dose24" else "dose25", CTRL[[y]], SCALE)
      data.table(outcome = CTRL[[y]], var = y, year = a, RR_10pp = r$RR, CI95_lwr = r$CI95_lwr, CI95_upr = r$CI95_upr,
                 p = r$p, n_obs = r$n_obs, events_target_2025 = d[target == 1L & year == 2025L, sum(get(y))])
    }))
  }), fill = TRUE)
  fwrite(grc, file.path(DIR$mod, "table_controls_sih_gradient.csv"))
  cat("\n== SIH negative controls: gradient in the expanded municipalities\n")
  print(grc[, .(outcome, year, RR_10pp = round(RR_10pp, 3), CI95 = sprintf("%.3f-%.3f", CI95_lwr, CI95_upr))])
} else log_msg("WARNING: ", FILE_CTRL, " missing; negative-control gradient not estimated")
log_msg("33 done")
