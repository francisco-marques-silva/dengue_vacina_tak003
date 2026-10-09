# =============================================================================
# 31_validity.R — design validity: periods with a common reference,
# calibrated inference, alternative clusterings and coverage bound
#
# 1. One treated × target-age-group coefficient per YEAR (and per Feb-Jan SEASON),
#    all in the same model, each re-anchored on the mean of the seven periods without
#    vaccine (2017-2023). The dispersion of the placebo periods calibrates the inference
#    for 2024 and 2025 (t with 6 df). -> Table 3A, Figure S2
# 2. Standard error of the main model clustered by UF and by health region (the
#    unit at which the first phase was assigned). -> Table S1
# 3. Monthly event study without bins and three re-anchorings of the ATT.
# 4. Plausibility bound: under intention to treat, the maximum reduction in a
#    month is 1 - efficacy × current coverage. Coverage is aggregated by period
#    weighting each month by observed HOSPITALIZATIONS (effective coverage) and
#    compared with the observed RRR (by period: model from script 30; by month: the
#    monthly series from block 3). -> Table 2B
# Outputs (outputs/models/): table_event_annual.csv, table_event_season.csv,
#   table_calibration_placebo.csv, table_cluster_alternative.csv,
#   table_att_reanchorings.csv, table_event_monthly_all_pre.csv,
#   table_plausibility_coverage.csv
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(fixest)

NEVER <- 10000L
panel <- readRDS(file.path(DIR$processed, "panel.rds"))[expn == 0L]
# only municipality×month cells with some hospitalization identify the model (fepois
# would drop the rest); filtering beforehand only speeds things up
pa <- copy(panel)[, tot_mt := sum(admissions), by = .(code_muni, t)][tot_mt > 0]
pa[, ta := as.integer(g < NEVER & target == 1L)]                        # treated × target age group
years <- PARAM$year_start:PARAM$year_end

# ---- 1. by year and by season, common reference -----------------------------
m_year <- fepois(as.formula(sprintf("admissions ~ i(year, ta, ref = %d) | %s", years[1], FE_MAIN)),
                data = pa, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
ev_year <- reanchor_by_year(m_year, "ta")
ev_year[, period := as.character(year)]

pa[, season := fifelse(as.integer(format(month, "%m")) >= 2L, year, year - 1L)]
pa[season < PARAM$year_start, season := PARAM$year_start]        # Jan/2017 joins 2017/18
temps <- sort(unique(pa$season))
m_temp <- fepois(as.formula(sprintf("admissions ~ i(season, ta, ref = %d) | %s", temps[1], FE_MAIN)),
                 data = pa, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
# reanchor_by_year expects the name "year::"; the season uses the same calculation with another prefix
nms <- sprintf("season::%d:ta", temps[-1]); b <- coef(m_temp)[nms]; V <- vcov(m_temp)[nms, nms]
K <- length(b); n_pre <- sum(temps <= 2023L); wpre <- rep(0, K); wpre[temps[-1] <= 2023L] <- 1 / n_pre
ev_temp <- rbindlist(lapply(seq_along(temps), function(i) {
  w <- rep(0, K); if (i > 1) w[i - 1] <- 1; w <- w - wpre; r <- lincomb(w, b, V)
  data.table(year = temps[i], period = sprintf("%d/%02d", temps[i], (temps[i] + 1L) %% 100L),
             log_b = r[["est"]], se = r[["se"]], RR = exp(r[["est"]]),
             CI95_lwr = exp(r[["est"]] - 1.96 * r[["se"]]), CI95_upr = exp(r[["est"]] + 1.96 * r[["se"]]))
}))
sd_pre <- sd(ev_temp[year <= 2023L, log_b])
ev_temp[, t_calibrated := log_b / (sd_pre * sqrt(1 + 1 / n_pre))]
ev_temp[, p_calibrated_unilateral := fifelse(year > 2023L, pt(t_calibrated, df = n_pre - 1), NA_real_)]
ev_temp[, p_calibrated_twosided := fifelse(year > 2023L, 2 * pt(-abs(t_calibrated), df = n_pre - 1), NA_real_)]
ev_temp[, sd_log_placebos := sd_pre]

cols <- c("period", "log_RRR", "se", "RRR", "CI95_lwr", "CI95_upr", "t_calibrated",
          "p_calibrated_unilateral", "p_calibrated_twosided", "sd_log_placebos")
setnames(ev_year, c("log_b", "RR"), c("log_RRR", "RRR")); setnames(ev_temp, c("log_b", "RR"), c("log_RRR", "RRR"))
fwrite(ev_year[, ..cols], file.path(DIR$mod, "table_event_annual.csv"))
fwrite(ev_temp[, ..cols], file.path(DIR$mod, "table_event_season.csv"))
cal <- rbind(ev_year[year >= 2024L, .(subset = "calendar year", period, sd_log_placebos, t = t_calibrated, gl = 6L,
                                   p_unilateral = p_calibrated_unilateral, p_twosided = p_calibrated_twosided)],
             ev_temp[year >= 2024L, .(subset = "season Feb-Jan", period, sd_log_placebos, t = t_calibrated, gl = 6L,
                                    p_unilateral = p_calibrated_unilateral, p_twosided = p_calibrated_twosided)])
fwrite(cal, file.path(DIR$mod, "table_calibration_placebo.csv"))
cat("\n== by year (reference: 2017-2023 mean)\n"); print(ev_year[, .(period, RRR = round(RRR, 3), CI95_lwr = round(CI95_lwr, 3), CI95_upr = round(CI95_upr, 3), p_calibrated_unilateral = signif(p_calibrated_unilateral, 2))])
cat("\n== by season\n"); print(ev_temp[, .(period, RRR = round(RRR, 3), CI95_lwr = round(CI95_lwr, 3), CI95_upr = round(CI95_upr, 3), p_calibrated_unilateral = signif(p_calibrated_unilateral, 2))])

# ---- 2. standard error under alternative clusterings ------------------------
m_main <- fepois(as.formula(sprintf("admissions ~ treat | %s", FE_MAIN)), data = pa,
                  offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
pa <- merge(pa, health_region(), by = "code_muni", all.x = TRUE)
pa[is.na(health_reg), health_reg := paste0("uf", uf)]
m_rs <- fepois(as.formula(sprintf("admissions ~ treat | %s", FE_MAIN)), data = pa,
               offset = ~log(pop), cluster = ~health_reg, notes = FALSE)
cl <- rbind(rr_row(m_main, "treat", "municipality (main)"),
            rr_row(summary(m_main, cluster = ~uf), "treat", "UF"),
            rr_row(m_rs, "treat", "health region (unit of assignment)"))
setnames(cl, c("analysis", "RR"), c("clustering", "RRR"))
fwrite(cl, file.path(DIR$mod, "table_cluster_alternative.csv"))
cat("\n== standard error under alternative clusterings\n"); print(cl)

# ---- 3. monthly event study without bins; three re-anchorings of the ATT ----
t_ini <- (as.integer(format(PARAM$vaccine_start, "%Y")) - PARAM$year_start) * 12L + as.integer(format(PARAM$vaccine_start, "%m"))
pa[, rel := fifelse(ta == 1L, t - t_ini, -1000L)]
m_es <- fepois(as.formula(sprintf("admissions ~ i(rel, ref = c(-1, -1000)) | %s", FE_MAIN)),
               data = pa, offset = ~log(pop), cluster = ~code_muni, notes = FALSE)
nm <- names(coef(m_es)); ks <- as.integer(sub("^rel::(-?[0-9]+).*$", "\\1", nm)); ok <- !is.na(ks)
nm <- nm[ok]; ks <- ks[ok]; b <- coef(m_es)[nm]; V <- vcov(m_es)[nm, nm]
att <- list()
for (r in list(list(sel = ks < 0 & ks >= -24, lbl = "mean of the previous 24 months (Feb/2022-Dec/2023; ~07b)"),
               list(sel = ks < 0,             lbl = "mean of all pre months (Jan/2017-Dec/2023)"),
               list(sel = ks < -24,           lbl = "mean of Jan/2017-Jan/2022 (excludes 2022-2023)"))) {
  w <- rep(0, length(b)); w[ks >= 0] <- 1 / sum(ks >= 0); w[r$sel] <- w[r$sel] - 1 / sum(r$sel)
  e <- lincomb(w, b, V)
  att[[r$lbl]] <- data.table(reanchoring = r$lbl, months_pre = sum(r$sel), RRR = exp(e[["est"]]),
                             CI95_lwr = exp(e[["est"]] - 1.96 * e[["se"]]), CI95_upr = exp(e[["est"]] + 1.96 * e[["se"]]))
}
att <- rbindlist(att)
fwrite(att, file.path(DIR$mod, "table_att_reanchorings.csv"))
cat("\n== monthly ATT under three re-anchorings\n"); print(att)

w0 <- rep(0, length(b)); w0[ks < 0] <- 1 / sum(ks < 0)
series <- rbindlist(lapply(c(seq_along(b), 0L), function(i) {
  w <- -w0; if (i > 0) w[i] <- w[i] + 1; e <- lincomb(w, b, V)
  data.table(k = if (i > 0) ks[i] else -1L, RR_reanch = exp(e[["est"]]),
             lwr_reanch = exp(e[["est"]] - 1.96 * e[["se"]]), upr_reanch = exp(e[["est"]] + 1.96 * e[["se"]]))
}))[order(k)]
fwrite(series, file.path(DIR$mod, "table_event_monthly_all_pre.csv"))

# ---- 4. plausibility bound given coverage ----------------------------------
EFFICACY <- c(ef100 = 1.00, ef675 = 0.675, ef746 = 0.746, ef841 = 0.841, ef879 = 0.879)
# 100%: absolute bound; 67.5%: effectiveness against hospitalization at ages 10-14, SP,
# 2024 (Ranzani et al.); 74.6% / 87.9%: 1 and 2 doses, 1 year (Krug Mareto et al.);
# 84.1%: efficacy against hospitalization in the phase 3 trial, 54 months (Tricou et al.)
trs <- panel[g < NEVER, unique(code_muni)]
monthly <- panel[code_muni %in% trs & agegroup == TARGET_AGEGROUP & month >= PARAM$vaccine_start,
                 .(coverage_d1 = sum(cov_d1_l1 * pop) / sum(pop), admissions = sum(admissions),
                   pop_target = sum(pop)), by = month][order(month)]
# observed monthly = event study without bins re-anchored on the mean of ALL
# pre months (block 3; the series cited in the text), not the binned series from script 30
ev <- copy(series)
months_k <- seq(PARAM$vaccine_start, by = "month", length.out = max(ev$k) + 1L)
ev[k >= 0, month := months_k[k + 1L]]
monthly <- merge(monthly, ev[k >= 0, .(month, obs_RR = RR_reanch, obs_lwr = lwr_reanch, obs_upr = upr_reanch)],
                by = "month", all.x = TRUE)

ag <- rbindlist(lapply(list(list("2024 (Feb-Dec)", monthly[month <= as.Date("2024-12-01")]),
                            list("2025", monthly[month >= as.Date("2025-01-01")]),
                            list("Whole period (Feb/2024-Dec/2025)", monthly)), function(x)
  data.table(label_txt = x[[1]], months = nrow(x[[2]]), admissions = sum(x[[2]]$admissions),
             coverage_d1 = sum(x[[2]]$coverage_d1 * x[[2]]$admissions) / sum(x[[2]]$admissions),
             coverage_d1_wtd_time = mean(x[[2]]$coverage_d1))))
het <- fread(file.path(DIR$mod, "heterogeneity_season.csv"))
t2  <- fread(file.path(DIR$mod, "table2_main_effect.csv"))[outcome == "Admissions (SIH)"]
obs <- rbind(het[, .(label_txt = c(treat_2024 = "2024 (Feb-Dec)", treat_2025 = "2025")[term],
                     obs_RR = exp(Estimate), obs_lwr = exp(Estimate - 1.96 * `Std. Error`),
                     obs_upr = exp(Estimate + 1.96 * `Std. Error`))],
             data.table(label_txt = "Whole period (Feb/2024-Dec/2025)", obs_RR = t2$RRR, obs_lwr = t2$CI95_lwr, obs_upr = t2$CI95_upr))
ag <- merge(ag, obs, by = "label_txt", all.x = TRUE, sort = FALSE)
output <- rbind(monthly[, .(level = "monthly", label_txt = format(month, "%Y-%m"), months = 1L, admissions, coverage_d1,
                          coverage_d1_wtd_time = NA_real_, obs_RR, obs_lwr, obs_upr)],
               ag[, .(level = "aggregate", label_txt, months, admissions, coverage_d1, coverage_d1_wtd_time,
                      obs_RR, obs_lwr, obs_upr)])
for (e in names(EFFICACY)) output[, (paste0("bound_", e)) := 1 - EFFICACY[[e]] * coverage_d1]
output[, violates_ef100 := !is.na(obs_RR) & obs_RR < bound_ef100]
output[, violates_ef675 := !is.na(obs_RR) & obs_RR < bound_ef675]
fwrite(output, file.path(DIR$mod, "table_plausibility_coverage.csv"))
cat("\n== effective coverage and bound\n")
print(output[level == "aggregate", .(label_txt, admissions, coverage_pct = round(100 * coverage_d1, 1),
                                  bound_ef100 = round(bound_ef100, 2), bound_ef675 = round(bound_ef675, 2),
                                  bound_ef841 = round(bound_ef841, 2),
                                  observed = sprintf("%.2f (%.2f-%.2f)", obs_RR, obs_lwr, obs_upr), violates_ef100)])

# =====================================================================================
# Blocks 5-8: current coverage, effective coverage, seasonal re-anchoring and age composition.
# =====================================================================================
pf  <- readRDS(file.path(DIR$processed, "panel.rds"))                       # all municipalities
grp <- unique(pf[, .(code_muni, group = fifelse(g < NEVER, "first phase", fifelse(expn == 1L, "expanded", "never included")))])
month_of_k <- function(k) {                                                       # k = months since Feb/2024
  m <- as.integer(format(PARAM$vaccine_start, "%Y")) * 12L + as.integer(format(PARAM$vaccine_start, "%m")) - 1L + k
  as.Date(sprintf("%d-%02d-01", m %/% 12L, m %% 12L + 1L))
}

# ---- 5. CURRENT (in-force) coverage in the 10-14 age group -------------------
# D1 by single year of age (saved by script 15): cohort = year of dose - age at dose; in month M
# the 10-14 age group holds the cohorts with year(M) - cohort in 10..14; vaccinated in those cohorts
# up to the previous month (same 1-month lag as the panel) / 10-14 population in the year of M,
# truncated at 100%. The existing measure (cumulative doses at 10-14 / 2024 pop) counts in 2025
# those who have left the age group and misses those who entered it vaccinated at age 9.
vi <- readRDS(file.path(DIR$processed, "dengue_vaccination_age.rds"))[!is.na(age)]
vi[, cohort := as.integer(format(month, "%Y")) - age]
MV <- seq(as.Date("2024-01-01"), as.Date(sprintf("%d-12-01", PARAM$year_end)), by = "month")
vi <- vi[cohort %in% 2010:2015 & month %in% MV, .(d1 = sum(d1)), by = .(code_muni, cohort, month)]
gv <- CJ(code_muni = grp$code_muni, cohort = 2010:2015, month = MV)
gv <- merge(gv, vi, by = c("code_muni", "cohort", "month"), all.x = TRUE)[is.na(d1), d1 := 0L]
setorder(gv, code_muni, cohort, month)
gv[, cum_l1 := shift(cumsum(d1), 1, fill = 0), by = .(code_muni, cohort)]    # up to the previous month
gv <- gv[(as.integer(format(month, "%Y")) - cohort) %between% c(10L, 14L), .(vac = sum(cum_l1)), by = .(code_muni, month)]
pop_target <- readRDS(file.path(DIR$processed, "population.rds"))[agegroup == TARGET_AGEGROUP & year %in% 2024:2025, .(code_muni, year, pop_year = pop)]
gv[, year := as.integer(format(month, "%Y"))]
gv <- merge(gv, pop_target, by = c("code_muni", "year"))[pop_year > 0]
gv[, cov_inforce := pmin(vac / pop_year, 1)]
cm <- pf[agegroup == TARGET_AGEGROUP & month >= PARAM$vaccine_start, .(code_muni, month, cov_current = cov_d1_l1, pop, admissions)]
cm <- merge(cm, gv[, .(code_muni, month, cov_inforce)], by = c("code_muni", "month"), all.x = TRUE)[is.na(cov_inforce), cov_inforce := 0]
cm <- merge(cm, grp, by = "code_muni")
saveRDS(cm, file.path(DIR$processed, "inforce_coverage_municipality.rds"))       # used by script 35 (estimate × ceiling table)
cov_if <-cm[, .(municipalities = .N, cov_current = sum(cov_current * pop) / sum(pop), cov_inforce = sum(cov_inforce * pop) / sum(pop),
                  admissions_target = sum(admissions)), by = .(group, month)][order(match(group, GROUP_ORDER), month)]
fwrite(cov_if, file.path(DIR$mod, "table_inforce_coverage_agegroup.csv"))
cat("\n== in-force vs current coverage, Dec/2024 and Dec/2025\n")
print(cov_if[month %in% as.Date(c("2024-12-01", "2025-12-01")), .(group, month, cov_current = round(100 * cov_current, 1), cov_inforce = round(100 * cov_inforce, 1))])

# ---- 6. effective coverage with the current measure; observed and expected weighting ----
# expected in month m (group G) = target hospitalizations of never-included municipalities in m × ratio
# G / never-included of the target hospitalizations in 2017-2023
ia <- merge(pf[agegroup == TARGET_AGEGROUP, .(admissions = sum(admissions)), by = .(code_muni, month)], grp, by = "code_muni")
ia <- ia[, .(i = sum(admissions)), by = .(group, month)]
ratio_pre <- ia[month <= as.Date("2023-12-01"), .(i = sum(i)), by = group]
ratio_pre <- ratio_pre[, .(group, R = i / i[group == "never included"])]
i_never <- ia[group == "never included", .(month, i_never = i)]
expd <- ratio_pre[group != "never included", .(month = i_never$month, expected = i_never$i_never * R), by = group]
pv <- merge(cov_if[group != "never included"], expd, by = c("group", "month"))
PER <- list("2024 (Feb-Dec)" = quote(month <= as.Date("2024-12-01")), "2025" = quote(month >= as.Date("2025-01-01")),
            "Whole period (Feb/2024-Dec/2025)" = quote(TRUE))
plaus2 <- rbindlist(lapply(names(PER), function(p) {
  x <- pv[eval(PER[[p]])]
  x[, .(period = p, admissions = sum(admissions_target), expected = sum(expected),
        current_observed = sum(cov_current * admissions_target) / sum(admissions_target),
        current_expected = sum(cov_current * expected) / sum(expected),
        inforce_observed = sum(cov_inforce * admissions_target) / sum(admissions_target),
        inforce_expected = sum(cov_inforce * expected) / sum(expected)), by = group]
}))
plaus2 <- merge(plaus2, ratio_pre[, .(group, ratio_pre_2017_2023 = R)], by = "group")[order(match(group, GROUP_ORDER))]
fwrite(plaus2, file.path(DIR$mod, "table_plausibility_v2.csv"))
cat("\n== effective coverage: current vs in-force, observed vs expected\n"); print(plaus2)

# ---- 7. monthly event study without bins, seasonally matched re-anchoring ---------
# each month k (including k = -1, whose coefficient is 0) minus the mean of the coefficients of the
# same calendar month in 2017-2023, via a linear combination with the full covariance
nm_es <- names(coef(m_es)); k_es <- as.integer(sub("^rel::(-?[0-9]+).*$", "\\1", nm_es)); ok_es <- !is.na(k_es)
nm_es <- nm_es[ok_es]; k_es <- k_es[ok_es]; b_es <- coef(m_es)[nm_es]; V_es <- vcov(m_es)[nm_es, nm_es]
mk <- month_of_k(k_es); cal <- as.integer(format(mk, "%m")); year_k <- as.integer(format(mk, "%Y"))
weight_seas <- function(k) {                                                       # weight vector for re-anchored month k
  w <- rep(0, length(b_es)); if (k != -1L) w[match(k, k_es)] <- 1
  same <- which(year_k <= 2023L & cal == as.integer(format(month_of_k(k), "%m")))
  w[same] <- w[same] - 1 / length(same); w
}
row_seas <- function(w, type, k = NA_integer_) {
  e <- lincomb(w, b_es, V_es)
  data.table(type = type, k = k, month = if (is.na(k)) as.Date(NA) else month_of_k(k), log_b = e[["est"]], se = e[["se"]],
             RR = exp(e[["est"]]), CI95_lwr = exp(e[["est"]] - 1.96 * e[["se"]]), CI95_upr = exp(e[["est"]] + 1.96 * e[["se"]]))
}
all_k <- sort(c(k_es, -1L))
W_seas <- sapply(all_k, weight_seas)
seas <- rbindlist(lapply(seq_along(all_k), function(j) row_seas(W_seas[, j], "monthly", all_k[j])))
mean_w <- function(sel) rowMeans(W_seas[, sel, drop = FALSE])
seas <- rbind(seas,
             row_seas(mean_w(all_k >= 0), "post mean (Feb/2024-Dec/2025)"),
             row_seas(mean_w(all_k >= 0 & all_k <= 10), "2024 mean (Feb-Dec)"),
             row_seas(mean_w(all_k >= 11), "2025 mean"),
             row_seas(mean_w(all_k %in% 1:5), "Mar-Jul/2024 mean"))
seas[, mar_jul_2024 := type == "monthly" & k %in% 1:5]
fwrite(seas, file.path(DIR$mod, "table_event_monthly_seasonal.csv"))
cat("\n== seasonally matched re-anchoring\n")
print(seas[type != "monthly" | mar_jul_2024, .(type, month, RR = round(RR, 3), CI95_lwr = round(CI95_lwr, 3), CI95_upr = round(CI95_upr, 3))])

# ---- 8. age composition by window: share of 5-59 hospitalizations at ages 10-14 ----------
# windows from the text, by year Y: Jan-May of Y; Jun of Y to Feb of Y+1; Apr-Jul of Y. 95% CI by
# municipality bootstrap (500 resamples, within each group; fixed seed 20260925).
cp <- panel[, .(target = sum(admissions * (agegroup == TARGET_AGEGROUP)), total = sum(admissions)),
             by = .(code_muni, group = fifelse(g < NEVER, "first phase", "never included"), month)]
cp[, `:=`(a = as.integer(format(month, "%Y")), m = as.integer(format(month, "%m")))]
jan <- rbind(cp[m <= 5L, .(code_muni, group, window = "jan-may", year = a, target, total)],
             cp[m >= 6L, .(code_muni, group, window = "jun-feb", year = a, target, total)],
             cp[m <= 2L, .(code_muni, group, window = "jun-feb", year = a - 1L, target, total)],
             cp[m %between% c(4L, 7L), .(code_muni, group, window = "apr-jul", year = a, target, total)])
jan <- jan[year >= PARAM$year_start, .(target = sum(target), total = sum(total)), by = .(code_muni, group, window, year)]
jan[, cell := paste(window, year)]
set.seed(20260925)
B_BOOT <- 500L
boot_frac <- function(g) {                                     # matrix of resamples × cells
  x <- jan[group == g]; A <- dcast(x, code_muni ~ cell, value.var = "target", fill = 0)
  Tt <- dcast(x, code_muni ~ cell, value.var = "total", fill = 0)
  A <- as.matrix(A[, -1]); Tt <- as.matrix(Tt[, -1]); n <- nrow(A)
  cnt <- rmultinom(B_BOOT, n, rep(1 / n, n))                   # counts of each municipality in each resample
  list(est = colSums(A) / colSums(Tt), boot = t(cnt) %*% A / (t(cnt) %*% Tt))
}
bf <- boot_frac("first phase"); bn <- boot_frac("never included")
cel <- intersect(names(bf$est), names(bn$est))
q <- function(m) apply(m, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
qf <- q(bf$boot[, cel, drop = FALSE]); qn <- q(bn$boot[, cel, drop = FALSE]); qr <- q(bf$boot[, cel] / bn$boot[, cel])
cmpr <- data.table(cell = cel, frac_first_phase = bf$est[cel], frac_first_phase_lwr = qf[1, ], frac_first_phase_upr = qf[2, ],
                   frac_never = bn$est[cel], frac_never_lwr = qn[1, ], frac_never_upr = qn[2, ],
                   ratio = bf$est[cel] / bn$est[cel], ratio_lwr = qr[1, ], ratio_upr = qr[2, ])
cmpr[, `:=`(window = sub(" [0-9]{4}$", "", cell), year = as.integer(sub("^.* ", "", cell)))]
cmpr[, complete := !(window == "jun-feb" & year == PARAM$year_end)]            # Jun-Feb of 2025 only runs through Dec/2025
cmpr <- cmpr[order(match(window, c("jan-may", "jun-feb", "apr-jul")), year),
             .(window, year, complete, frac_first_phase, frac_first_phase_lwr, frac_first_phase_upr, frac_never, frac_never_lwr,
               frac_never_upr, ratio, ratio_lwr, ratio_upr)]
fwrite(cmpr, file.path(DIR$mod, "table_composition_age_windows.csv"))
cat("\n== age composition by window\n")
print(cmpr[, .(window, year, f1 = round(100 * frac_first_phase, 1), never = round(100 * frac_never, 1), ratio = round(ratio, 2))])
log_msg("31 done")
