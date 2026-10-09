# =============================================================================
# 30_main_model.R — description, triple difference and event study
#
# Model (Poisson with high-dimensional fixed effects, fixest::fepois):
#   log E[Y_mat] = log(pop) + β·Treat_mat + α_{m×t} + γ_{m×a} + δ_{a×t}
#   m municipality, a age group, t month; Trat = 1 for ages 10-14 in the 521
#   first-phase municipalities from Feb/2024 onward. exp(β) = ratio of rate
#   ratios (RRR). Standard errors clustered by municipality.
# Sample: 521 treated + 2,819 never-included municipalities (the expansion is excluded).
#
# Blocks: 1 description (Table 1); 2 main model by outcome (Table 2A) and
# averted hospitalizations; 3 effect by vaccination year (2024 × 2025, Table 2B);
# 4 monthly event study (Sun & Abraham; reference k = -1) with pre-trend tests
# and re-anchoring on the mean of the whole pre-period (Figure S1).
# Outputs in outputs/models/: table1_rates_year_agegroup_group.csv,
#   table2_main_effect.csv, averted_admissions.csv,
#   heterogeneity_season.csv, table_event_study_reanchored.csv,
#   table_pretrend.csv, main_models.rds, event_study.rds
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(fixest)

panel <- readRDS(file.path(DIR$processed, "panel.rds"))[expn == 0L]
NEVER <- 10000L
log_msg("Main sample: ", uniqueN(panel$code_muni), " municipalities (", panel[g < NEVER, uniqueN(code_muni)], " treated)")

# ---- 1. description: annual rates by group, age group and year (Table 1) ----
panel[, group := fifelse(g < NEVER, "first phase", "never included")]
tab1 <- panel[, .(pop_mean = sum(pop) / uniqueN(month), admissions = sum(admissions),   # mean population of the year
                   cases = sum(den_cases), deaths = sum(deaths_sim)), by = .(group, agegroup, year)]
tab1[, `:=`(rt_adm_100k = 1e5 * admissions / pop_mean, rt_cases_100k = 1e5 * cases / pop_mean,
            rt_deaths_100k = 1e5 * deaths / pop_mean)]
fwrite(tab1, file.path(DIR$mod, "table1_rates_year_agegroup_group.csv"))

# ---- 2. main model, by outcome -----------------------------------------------
OUTCOMES <- c(admissions = "Admissions (SIH)", den_cases = "Probable cases (SINAN)",
               den_severe = "Dengue with warning signs/severe (SINAN)", den_hosp = "Hospitalized (SINAN)")
mods <- list(); res <- list()
for (y in names(OUTCOMES)) {
  mods[[y]] <- fepois(as.formula(sprintf("%s ~ treat | %s", y, FE_MAIN)), data = panel,
                      offset = ~log(pop), cluster = ~code_muni)
  res[[y]] <- cbind(outcome = OUTCOMES[[y]], rr_row(mods[[y]], "treat", "main model"))
}
res <- rbindlist(res)[, .(outcome, RRR = RR, CI95_lwr, CI95_upr, p, reduction_pct = 100 * (1 - RR), n_obs)]
fwrite(res, file.path(DIR$mod, "table2_main_effect.csv"))
saveRDS(mods, file.path(DIR$mod, "main_models.rds"))
print(res)

# averted hospitalizations (counterfactual = observed / RRR)
b <- coef(mods$admissions)[["treat"]]
avert <- panel[treat == 1L, .(observed = sum(admissions))][, averted := observed * (exp(-b) - 1)][]
fwrite(avert, file.path(DIR$mod, "averted_admissions.csv"))

# ---- 3. effect by vaccination year (2024 and 2025 in the same model) ---------
panel[, `:=`(treat_2024 = treat * (year == 2024L), treat_2025 = treat * (year == 2025L))]
m_temp <- fepois(as.formula(sprintf("admissions ~ treat_2024 + treat_2025 | %s", FE_MAIN)),
                 data = panel, offset = ~log(pop), cluster = ~code_muni)
fwrite(as.data.table(coeftable(m_temp), keep.rownames = "term"), file.path(DIR$mod, "heterogeneity_season.csv"))

# ---- 4. monthly event study (Sun & Abraham; ref. k = -1; bins <= -24 and >= 22)
panel[, rel := fifelse(cohort < NEVER, t - cohort, NA_integer_)]
rels <- sort(unique(na.omit(panel$rel)))
BIN <- list("-24" = rels[rels <= -24], "22" = rels[rels >= 22]); BIN <- BIN[lengths(BIN) > 0]
m_es <- fepois(as.formula(sprintf("admissions ~ sunab(cohort, t, ref.p = -1, bin.rel = %s) | %s",
                                  paste(deparse(BIN, width.cutoff = 500L), collapse = ""), FE_MAIN)),
               data = panel, offset = ~log(pop), cluster = ~code_muni)
saveRDS(m_es, file.path(DIR$mod, "event_study.rds"))

ct <- as.data.table(coeftable(m_es), keep.rownames = "term")
ct[, k := as.integer(sub(".*::(-?[0-9]+).*", "\\1", term))]
ct <- ct[!is.na(k)][order(k)]
setnames(ct, c("Estimate", "Std. Error"), c("b", "se"))
# with a single adoption cohort, the sunab aggregation is the identity: it is enough
# to reorder the raw covariance by the same relative periods
Vb <- vcov(m_es); k_v <- as.integer(sub("^t::(-?[0-9]+).*", "\\1", colnames(Vb)))
stopifnot(!anyNA(k_v), setequal(k_v, ct$k))
V <- Vb[match(ct$k, k_v), match(ct$k, k_v), drop = FALSE]
b <- ct$b; pre <- which(ct$k < 0); post <- which(ct$k >= 0); K <- length(b)

# 4a. joint Wald: pre coefficients = 0 (reference k = -1)
R <- diag(K)[pre, , drop = FALSE]
W <- as.numeric(t(R %*% b) %*% solve(R %*% V %*% t(R)) %*% (R %*% b))
p_wald <- pchisq(W, length(pre), lower.tail = FALSE)
# 4b. anchor-free flatness of the pre-period: deviations of k <= -3 around their own mean
dist <- which(ct$k <= -3)
wd <- rep(0, K); wd[dist] <- 1 / length(dist)
Md <- (diag(K) - matrix(rep(wd, each = K), nrow = K))[dist, , drop = FALSE][-1, , drop = FALSE]
bd <- as.numeric(Md %*% b); Wd <- as.numeric(t(bd) %*% solve(Md %*% V %*% t(Md)) %*% bd)
p_flat <- pchisq(Wd, nrow(Md), lower.tail = FALSE)
# 4c. linear trend in the pre-period
k_pre <- ct$k[pre]; cl <- rep(0, K); cl[pre] <- (k_pre - mean(k_pre)) / sum((k_pre - mean(k_pre))^2)
lin <- lincomb(cl, b, V); p_lin <- 2 * pnorm(abs(lin[["est"]] / lin[["se"]]), lower.tail = FALSE)
log_msg(sprintf("Wald pre: X2 = %.1f (df %d), p = %.3g | flatness (k <= -3): p = %.3g | linear trend: p = %.3g",
                W, length(pre), p_wald, p_flat, p_lin))

# 4d. re-anchoring on the mean of the WHOLE pre-period: the "-24" bin pools 62 months
# (Jan/2017-Feb/2022) and gets a weight proportional to that number of months
n_pre_months <- (as.integer(format(PARAM$vaccine_start, "%Y")) - PARAM$year_start) * 12L +
  as.integer(format(PARAM$vaccine_start, "%m")) - 1L                   # 85 months before Feb/2024
weight <- ifelse(ct$k == min(ct$k), n_pre_months - length(pre), 1)
w <- rep(0, K); w[pre] <- weight[pre] / sum(weight[pre])
M <- diag(K) - matrix(rep(w, each = K), nrow = K)
b_r <- as.numeric(M %*% b); se_r <- sqrt(pmax(diag(M %*% V %*% t(M)), 0))
w2 <- rep(0, K); w2[post] <- 1 / length(post)
att <- lincomb(w2 - w, b, V)

fwrite(ct[, .(k, RR_ref_m1 = exp(b), lwr_ref_m1 = exp(b - 1.96 * se), upr_ref_m1 = exp(b + 1.96 * se),
              RR_reanch = exp(b_r), lwr_reanch = exp(b_r - 1.96 * se_r), upr_reanch = exp(b_r + 1.96 * se_r))],
       file.path(DIR$mod, "table_event_study_reanchored.csv"))
fwrite(data.table(
  measure = c("Re-anchored ATT (post mean - pre mean)", "CI95_lwr", "CI95_upr", "p ATT",
             "Wald pre (X2)", "gl", "p Wald pre", "Linear pre-trend (log-RR/month)", "p linear trend",
             "Wald pre flatness (k<=-3, no anchor)", "df flatness", "p flatness"),
  val  = c(exp(att[["est"]]), exp(att[["est"]] - 1.96 * att[["se"]]), exp(att[["est"]] + 1.96 * att[["se"]]),
             2 * pnorm(abs(att[["est"]] / att[["se"]]), lower.tail = FALSE),
             W, length(pre), p_wald, lin[["est"]], p_lin, Wd, nrow(Md), p_flat)),
  file.path(DIR$mod, "table_pretrend.csv"))
log_msg(sprintf("Re-anchored ATT: %.3f (%.3f-%.3f)", exp(att[["est"]]),
                exp(att[["est"]] - 1.96 * att[["se"]]), exp(att[["est"]] + 1.96 * att[["se"]])))

# ---- reconciliation of Table 1 (Table S1 in the manuscript) --------------
# Two definitions of the comparison: that of Table 1 (5-9, 15-19 and 20-29) and all age groups
# 5-59 except 10-14 (that of the model). Crude DDD per year = (first-phase ratio in the year / mean
# of the annual ratios 2017-2023) ÷ (the same in the never-included municipalities).
COMP3 <- c("05-09", "15-19", "20-29")
rec <- tab1[, .(rt_target = 1e5 * sum(admissions[agegroup == TARGET_AGEGROUP]) / sum(pop_mean[agegroup == TARGET_AGEGROUP]),
                rt_cmpr3 = 1e5 * sum(admissions[agegroup %in% COMP3]) / sum(pop_mean[agegroup %in% COMP3]),
                rt_cmpr_all = 1e5 * sum(admissions[agegroup != TARGET_AGEGROUP]) / sum(pop_mean[agegroup != TARGET_AGEGROUP])),
            by = .(group, year)]
rec[, `:=`(ratio_cmpr3 = rt_target / rt_cmpr3, ratio_cmpr_all = rt_target / rt_cmpr_all)]
rec[, `:=`(rel_cmpr3 = ratio_cmpr3 / mean(ratio_cmpr3[year <= 2023L]),
           rel_cmpr_all = ratio_cmpr_all / mean(ratio_cmpr_all[year <= 2023L])), by = group]
rw <- dcast(rec, year ~ group, value.var = c("rt_target", "rt_cmpr3", "rt_cmpr_all", "ratio_cmpr3", "ratio_cmpr_all",
                                           "rel_cmpr3", "rel_cmpr_all"))
setnames(rw, gsub("_first phase$", "_first_phase", gsub("_never included$", "_never", names(rw))))
rw[, `:=`(ddd_crude_cmpr3 = rel_cmpr3_first_phase / rel_cmpr3_never, ddd_crude_cmpr_all = rel_cmpr_all_first_phase / rel_cmpr_all_never)]
fwrite(rw, file.path(DIR$mod, "table_s1_reconciliation.csv"))
cat("\n== reconciliation of Table 1 (crude DDD by year)\n")
print(rw[, .(year, ratio_cmpr3_first_phase = round(ratio_cmpr3_first_phase, 2), ratio_cmpr_all_first_phase = round(ratio_cmpr_all_first_phase, 2),
             ddd_crude_cmpr3 = round(ddd_crude_cmpr3, 3), ddd_crude_cmpr_all = round(ddd_crude_cmpr_all, 3))])
