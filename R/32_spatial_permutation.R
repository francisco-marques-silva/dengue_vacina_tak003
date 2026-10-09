# =============================================================================
# 32_spatial_permutation.R — randomization inference over the set of
# treated municipalities
#
# N_PERM fictitious sets of "treated" municipalities are drawn with the same distribution
# of size and prior incidence as the 521 real ones (strata: 2022 population decile
# × 2017-2023 incidence above/below the decile median; in each stratum the same
# number of treated is drawn) and the model is re-estimated for each one.
# The position of the real coefficient in the fictitious distribution gives the Fisher p for
# 2019 (placebo), 2024, 2025 and the pooled effect; the standard deviation of the
# fictitious coefficients measures the real uncertainty of the design. Since each fictitious
# set repeats ~20% of the 521 real ones, the test is conservative.
# Fixed seed (20260923). Time: ~4-8 s per draw; ~1 h with N_PERM = 400.
# Outputs: outputs/models/spatial_permutation_municipality.csv (all draws;
#         resumable) and table_spatial_permutation.csv (summary)
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(fixest)

N_PERM <- as.integer(Sys.getenv("N_PERM", "400"))
NEVER  <- 10000L
set.seed(20260923)

panel <- readRDS(file.path(DIR$processed, "panel.rds"))[expn == 0L]
pa <- copy(panel)[, tot_mt := sum(admissions), by = .(code_muni, t)][tot_mt > 0]
t_ini <- (as.integer(format(PARAM$vaccine_start, "%Y")) - PARAM$year_start) * 12L + as.integer(format(PARAM$vaccine_start, "%m"))

# ---- draw strata (from the full panel: municipalities without hospitalizations are also eligible)
mun <- panel[year == 2022L, .(pop = sum(pop) / 12, treated = as.integer(g[1] < NEVER), uf = uf[1]), by = code_muni]
pre <- panel[year <= 2023L, .(rt = 1e5 * sum(admissions) / (sum(pop) / 12)), by = code_muni]
mun <- merge(mun, pre, by = "code_muni", all.x = TRUE)
mun[is.na(rt), rt := 0]
mun[, dec_pop := cut(pop, unique(quantile(pop, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
mun[, rt_high := as.integer(rt > median(rt)), by = dec_pop]
mun[, stratum := paste(dec_pop, rt_high, sep = "_")]
real_ids <- mun[treated == 1L, code_muni]

# ---- coefficients for a set of treated municipalities: 2019 (placebo), 2024, 2025 and pooled
estfit <- function(treated_ids) {
  pa[, tm := as.integer(code_muni %in% treated_ids & target == 1L)]
  pa[, `:=`(p2019 = tm * (year == 2019L), p2024 = tm * (year == 2024L), p2025 = tm * (year == 2025L),
            ppool = tm * (t >= t_ini))]
  m1 <- fepois(as.formula(sprintf("admissions ~ p2019 + p2024 + p2025 | %s", FE_MAIN)),
               data = pa, offset = ~log(pop), vcov = "iid", notes = FALSE)
  m2 <- fepois(as.formula(sprintf("admissions ~ ppool | %s", FE_MAIN)),
               data = pa, offset = ~log(pop), vcov = "iid", notes = FALSE)
  # the standard errors (iid, as in the fit) feed the p based on the t statistic
  c(coef(m1)[c("p2019", "p2024", "p2025")], ppool = unname(coef(m2)["ppool"]),
    setNames(se(m1)[c("p2019", "p2024", "p2025")], c("se_p2019", "se_p2024", "se_p2025")),
    se_ppool = unname(se(m2)["ppool"]))
}
real_set <- estfit(real_ids)
log_msg("Actual coefficients (log): ", paste(sprintf("%s=%.3f", names(real_set)[1:4], real_set[1:4]), collapse = ", "))

# ---- draws (resumable; the seed reproduces the sequence only within a single run)
fpath <- file.path(DIR$mod, "spatial_permutation_municipality.csv")
done <- if (file.exists(fpath)) nrow(fread(fpath)) else 0L
if (done > 0 && done < N_PERM)
  log_msg("Resuming from draw ", done + 1, ": the draw sequence is no longer that of a single run")
if (done < N_PERM) for (i in (done + 1):N_PERM) {
  S <- mun[, { k <- sum(treated); if (k == 0) character() else sample(code_muni, k) }, by = stratum]$V1
  r <- tryCatch(estfit(S), error = function(e) rep(NA_real_, 8))
  fwrite(data.table(perm = i, p2019 = r[[1]], p2024 = r[[2]], p2025 = r[[3]], ppool = r[[4]],
                    overlap = mean(real_ids %in% S),
                    se_p2019 = r[[5]], se_p2024 = r[[6]], se_p2025 = r[[7]], se_ppool = r[[8]]),
         fpath, append = file.exists(fpath))
  if (i %% 25 == 0) log_msg("draw ", i, "/", N_PERM)
}

# ---- summary -------------------------------------------------------------------
pl <- fread(fpath)[!is.na(ppool)]
summ <- rbindlist(lapply(c("p2019", "p2024", "p2025", "ppool"), function(v) {
  x <- pl[[v]]; real <- real_set[[v]]
  data.table(scheme = "municipality", coefficient = v, n_draws = length(x),
             overlap_mean = mean(pl$overlap), RRR_real = exp(real),
             RRR_placebo_median = exp(median(x)), RRR_placebo_p05 = exp(quantile(x, 0.05)),
             RRR_placebo_p95 = exp(quantile(x, 0.95)), sd_log_placebos = sd(x),
             p_unilateral = (1 + sum(x <= real)) / (1 + length(x)),
             p_twosided = (1 + sum(abs(x - median(x)) >= abs(real - median(x)))) / (1 + length(x)))
}))
fwrite(summ, file.path(DIR$mod, "table_spatial_permutation.csv"))
print(summ)

# ---- draw strata and p based on the t statistic (same draws) ----
strata <- data.table(deciles_population = uniqueN(mun$dec_pop), categories_incidence_pre = uniqueN(mun$rt_high),
                       strata_total = uniqueN(mun$stratum), strata_with_treated = mun[treated == 1L, uniqueN(stratum)],
                       municipalities_drawable = nrow(mun), treated_by_draw = sum(mun$treated))
fwrite(strata, file.path(DIR$mod, "table_permutation_strata.csv"))
if (all(paste0("se_", c("p2019", "p2024", "p2025", "ppool")) %in% names(pl))) {
  res_t <- rbindlist(lapply(c("p2019", "p2024", "p2025", "ppool"), function(v) {
    x <- pl[[v]] / pl[[paste0("se_", v)]]; x <- x[is.finite(x)]
    real <- real_set[[v]] / real_set[[paste0("se_", v)]]
    data.table(coefficient = v, n_draws = length(x), t_real = real, t_placebo_median = median(x),
               p_unilateral_t = (1 + sum(x <= real)) / (1 + length(x)),
               p_twosided_t = (1 + sum(abs(x - median(x)) >= abs(real - median(x)))) / (1 + length(x)))
  }))
  res_t <- merge(res_t, summ[, .(coefficient, p_unilateral_coef = p_unilateral, p_twosided_coef = p_twosided)],
                 by = "coefficient", sort = FALSE)
  fwrite(res_t, file.path(DIR$mod, "table_spatial_permutation_t.csv"))
  print(strata); print(res_t)
} else log_msg("WARNING: draws without standard error (file from an earlier version); p on t not computed")

# ---- randomization by HEALTH REGION ---------------------------------------------------
# Draw unit = health region; "treated" = regions with a first-phase municipality. In
# each draw, within each stratum (2022 population tercile of the region × 2017-2023
# incidence above/below the tercile median), the same number of treated regions is drawn
# and ALL eligible municipalities in them are treated (binary sample). The real set
# follows the same rule: all eligible municipalities of the treated regions (the 521 first-phase
# ones and the never-included municipalities located in those regions). Own seed, resumable.
rsm <- merge(unique(panel[, .(code_muni, g)]), health_region(), by = "code_muni")
reg_treat <- rsm[g < NEVER, unique(health_reg)]
pr <- merge(panel[, .(code_muni, year, pop, admissions)], rsm[, .(code_muni, health_reg)], by = "code_muni")
rg <- pr[, .(pop22 = sum(pop[year == 2022L]) / 12,
             rt = 1e5 * sum(admissions[year <= 2023L]) / (sum(pop[year <= 2023L]) / 12)), by = health_reg]
rg[, is_treated := as.integer(health_reg %in% reg_treat)]
rg[, terc := cut(pop22, unique(quantile(pop22, 0:3 / 3)), include.lowest = TRUE, labels = FALSE)]
rg[, rt_high := as.integer(rt > median(rt)), by = terc]
rg[, stratum := paste(terc, rt_high, sep = "_")]
munis_de <- function(regs) rsm[health_reg %in% regs, code_muni]
real_r <- estfit(munis_de(reg_treat))
log_msg(sprintf("regions: %d treated of %d eligible for drawing; %d strata; %d treated municipalities in the actual set",
                length(reg_treat), nrow(rg), uniqueN(rg$stratum), length(munis_de(reg_treat))))
fpath_r <- file.path(DIR$mod, "spatial_permutation_region.csv")
done_r <- if (file.exists(fpath_r)) nrow(fread(fpath_r)) else 0L
set.seed(20260925)
if (done_r > 0 && done_r < N_PERM) log_msg("WARNING: resuming the region draws; the sequence is no longer that of a single run")
if (done_r < N_PERM) for (i in (done_r + 1):N_PERM) {
  R <- rg[, { k <- sum(is_treated); if (k == 0) character() else sample(health_reg, k) }, by = stratum]$V1
  r <- tryCatch(estfit(munis_de(R)), error = function(e) rep(NA_real_, 8))
  fwrite(data.table(perm = i, p2019 = r[[1]], p2024 = r[[2]], p2025 = r[[3]], ppool = r[[4]],
                    se_p2019 = r[[5]], se_p2024 = r[[6]], se_p2025 = r[[7]], se_ppool = r[[8]],
                    regions_real_drawn = sum(R %in% reg_treat)), fpath_r, append = file.exists(fpath_r))
  if (i %% 50 == 0) log_msg("region draw ", i, "/", N_PERM)
}
plr <- fread(fpath_r)[!is.na(ppool)]
res_r <- rbindlist(lapply(c("p2019", "p2024", "p2025", "ppool"), function(v) {
  x <- plr[[v]]; rt_ <- x / plr[[paste0("se_", v)]]; real <- real_r[[v]]; real_t <- real / real_r[[paste0("se_", v)]]
  data.table(coefficient = v, n_draws = length(x), RRR_real = exp(real), RRR_main_model = exp(real_set[[v]]),
             RRR_placebo_p05 = exp(quantile(x, 0.05)), RRR_placebo_p95 = exp(quantile(x, 0.95)),
             p_unilateral_coef = (1 + sum(x <= real)) / (1 + length(x)),
             t_real = real_t, p_unilateral_t = (1 + sum(rt_ <= real_t, na.rm = TRUE)) / (1 + sum(is.finite(rt_))),
             regions_treated = length(reg_treat), regions_drawable = nrow(rg), strata = uniqueN(rg$stratum),
             municipalities_treated_real = length(munis_de(reg_treat)))
}))
fwrite(res_r, file.path(DIR$mod, "table_permutation_health_region.csv"))
cat("\n== permutation by health region\n"); print(res_r)
if (!requireNamespace("fwildclusterboot", quietly = TRUE))
  log_msg("fwildclusterboot missing (and it only supports feols, not fepois): wild bootstrap not run")
log_msg("32 done")
