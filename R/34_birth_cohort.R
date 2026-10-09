# =============================================================================
# 34_birth_cohort.R — comparison by BIRTH COHORT (no denominator)
#
# The fixed age groups of the main model do not separate direct protection from
# age dynamics: someone vaccinated at age 14 in 2024 is in the comparison group
# (15-19) in 2025, and someone aged 9 enters the target age group without having
# been eligible. Here the unit is the cohort (hospitalization year − age), which
# carries eligibility. With municipality × cohort and cohort × month fixed effects
# the size of each cohort is absorbed, and the model runs without an offset.
#
# A. Municipality × cohort × month panel (2000-2018), main sample: triple
#    difference pooled and by year, with three sets of comparison cohorts;
#    one coefficient per year re-anchored on 2017-2023 (adjacent cohorts).
# B. TWO GROUPS (eligible cohorts × non-eligible neighbouring cohorts), as in the
#    main model (municipality×month, municipality×group, group×month), with one
#    coefficient per year: (i) in the expansion, exposed group × final coverage (per 10 p.p.) —
#    the decisive test of the paper: if the 2025 gradient from script 33 is due to the
#    vaccine, it should appear against neighbouring cohorts; (ii) first phase × never-included municipalities.
#    Sets: exposed 2009-2015 vs clean neighbouring cohorts 2005-2007 and 2017-2019
#    (eligibility is by age AT VACCINATION: the 2009 cohort is mostly
#    eligible; and cohort = year − age shifts half of each cohort by one year, so
#    the immediately neighbouring cohorts are dropped); and, for comparison, exposed 2010-2015 vs
#    2007-2009 and 2016-2018 (the block A set). -> Table 3C
# Outputs (outputs/models/): table_birth_cohort.csv,
#   table_birth_cohort_annual.csv, table_cohort_group.csv
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(fixest)

NEVER <- 10000L
years  <- PARAM$year_start:PARAM$year_end
pc    <- readRDS(file.path(DIR$processed, "panel_cohort.rds"))
FE_C  <- "code_muni^t + code_muni^cohort + cohort^t"

# ---- A. cohort panel, main sample --------------------------------------
SET_A <- list("adjacent cohorts (2007-2009 and 2016-2018)"            = 2007:2018,
               "adjacent without boundary cohorts (2007-2008, 2017-2018)" = c(2007:2008, 2010:2015, 2017:2018),
               "all cohorts 2000-2018"                              = 2000:2018)
pc[, `:=`(treat_2024 = treat * (year == 2024L), treat_2025 = treat * (year == 2025L))]
pooled <- list(); by_year <- list()
for (nm in names(SET_A)) {
  d <- pc[cohort %in% SET_A[[nm]]][, tot_mt := sum(admissions), by = .(code_muni, t)][tot_mt > 0]
  for (expd in list(c("treat", "treat", nm), c("treat_2024 + treat_2025", "treat_2024", paste(nm, "— 2024")),
                   c("treat_2024 + treat_2025", "treat_2025", paste(nm, "— 2025")))) {
    m <- tryCatch(fepois(as.formula(sprintf("admissions ~ %s | %s", expd[1], FE_C)), data = d,
                         cluster = ~code_muni, notes = FALSE),
                  error = function(e) { warning(expd[3], ": ", conditionMessage(e)); NULL })
    if (is.null(m)) next
    r <- rr_row(m, expd[2], expd[3])
    row <- data.table(analysis = expd[3], RRR = r$RR, CI95_lwr = r$CI95_lwr, CI95_upr = r$CI95_upr, p = r$p,
                        n_obs = r$n_obs, admissions_exposed_post = d[treat == 1L, sum(admissions)])
    if (expd[1] == "treat") pooled[[expd[3]]] <- row else by_year[[expd[3]]] <- row
  }
}
res <- rbind(rbindlist(pooled), rbindlist(by_year))     # pooled first, then by year (reference order)
fwrite(res, file.path(DIR$mod, "table_birth_cohort.csv"))
cat("\n== triple difference by birth cohort\n"); print(res)

# by year, adjacent cohorts, common reference 2017-2023
d <- pc[cohort %in% SET_A[[1]]][, tot_mt := sum(admissions), by = .(code_muni, t)][tot_mt > 0]
d[, tc := treated * cohort_exposed]
m_year <- fepois(as.formula(sprintf("admissions ~ i(year, tc, ref = %d) | %s", years[1], FE_C)), data = d,
                cluster = ~code_muni, notes = FALSE)
ev <- reanchor_by_year(m_year, "tc")
setnames(ev, c("log_b", "RR"), c("log_RRR", "RRR"))
fwrite(ev[, .(year, log_RRR, se, RRR, CI95_lwr, CI95_upr, p_calibrated_unilateral)],
       file.path(DIR$mod, "table_birth_cohort_annual.csv"))
cat("\n== by year, adjacent cohorts\n"); print(ev[, .(year, RRR = round(RRR, 3), CI95_lwr = round(CI95_lwr, 3), CI95_upr = round(CI95_upr, 3))])

# ---- B. two groups: exposed × neighbouring cohorts ------------------------------------
sih <- readRDS(file.path(DIR$processed, "sih_age_aggregate.rds"))
sih[, year := as.integer(format(month, "%Y"))][, cohort := year - age]
panel <- readRDS(file.path(DIR$processed, "panel.rds"))
panel[, cov_final := max(cov_d1_l1), by = code_muni]
mun <- panel[, .(expn = expn[1], treated = as.integer(g[1] < NEVER), cov_final = cov_final[1]), by = code_muni]
months <- seq(as.Date(sprintf("%d-01-01", PARAM$year_start)), as.Date(sprintf("%d-12-01", PARAM$year_end)), by = "month")
FE_G <- "code_muni^t + code_muni^group + group^t"
SET_B <- list("clean adjacent (2005-2007 and 2017-2019)"       = list(exp = 2009:2015, adj = c(2005:2007, 2017:2019)),
               "adjacent from script 20 (2007-2009 and 2016-2018)" = list(exp = 2010:2015, adj = c(2007:2009, 2016:2018)))
res_g <- list(); sens_g <- list()
for (nm in names(SET_B)) for (samp in c("expn", "phase1")) {
  cj <- SET_B[[nm]]; stopifnot(!any(cj$exp %in% cj$adj))
  munis <- if (samp == "expn") mun[expn == 1L, code_muni] else mun[expn == 0L, code_muni]
  d <- sih[code_muni %in% munis & (cohort %in% cj$exp | cohort %in% cj$adj)]
  d[, group := as.integer(cohort %in% cj$exp)]
  g <- CJ(code_muni = munis, group = 0:1, month = months)
  g <- merge(g, d[, .(admissions = sum(admissions)), by = .(code_muni, group, month)],
             by = c("code_muni", "group", "month"), all.x = TRUE)
  g[is.na(admissions), admissions := 0L]
  g[, `:=`(year = as.integer(format(month, "%Y")),
           t = (as.integer(format(month, "%Y")) - PARAM$year_start) * 12L + as.integer(format(month, "%m")))]
  g <- merge(g, mun, by = "code_muni")
  g[, tot_mt := sum(admissions), by = .(code_muni, t)]; g <- g[tot_mt > 0]
  if (samp == "expn") {                       # exposed group × final coverage × year, per 10 p.p.
    g[, v := group * cov_final]; scale <- 0.10
    lbl <- sprintf("Expanded: exposed cohort x final coverage, per 10 p.p. — %s", nm)
  } else {                                      # exposed group × treated × year
    g[, v := group * treated]; scale <- 1
    lbl <- sprintf("first phase vs never: exposed cohort x treated — %s", nm)
  }
  m <- fepois(as.formula(sprintf("admissions ~ i(year, v, ref = %d) | %s", years[1], FE_G)), data = g,
              cluster = ~code_muni, notes = FALSE)
  tab <- reanchor_by_year(m, "v", scale = scale)
  tab[, `:=`(analysis = lbl, municipalities = uniqueN(g$code_muni),
             admissions_exposed_2025 = g[group == 1L & year == 2025L, sum(admissions)])]
  res_g[[paste(samp, nm)]] <- tab
  # sensitivity analyses of the same contrast: (i) same model, 2024 and 2025 re-anchored on the
  # 2017-2022 mean (2023 out of the reference; calibrated p with 5 df); (ii) same model without
  # the months Feb-Apr/2025 (near-expiry doses given outside the age group)
  s1 <- reanchor_by_year(m, "v", scale = scale, last_without_vaccine = 2022L)
  m_without <- fepois(as.formula(sprintf("admissions ~ i(year, v, ref = %d) | %s", years[1], FE_G)),
                  data = g[!(month >= as.Date("2025-02-01") & month <= as.Date("2025-04-01"))],
                  cluster = ~code_muni, notes = FALSE)
  s2 <- reanchor_by_year(m_without, "v", scale = scale)
  sens_g[[paste(samp, nm)]] <- rbind(
    s1[year %in% 2023:2025][, sensitivity := "reference 2017-2022 (without 2023)"],
    s2[year %in% 2024:2025][, sensitivity := "excluding Feb-Apr/2025"], fill = TRUE)[, analysis := lbl]
}
res_g <- rbindlist(res_g)
fwrite(res_g[, .(analysis, year, log_b, se, RR, CI95_lwr, CI95_upr, p_calibrated_unilateral, municipalities, admissions_exposed_2025)],
       file.path(DIR$mod, "table_cohort_group.csv"))
cat("\n== two groups (exposed vs neighbors)\n")
print(res_g[, .(analysis, year, RR = round(RR, 3), CI95_lwr = round(CI95_lwr, 3), CI95_upr = round(CI95_upr, 3), p_cal = signif(p_calibrated_unilateral, 2))])

sens_g <- rbindlist(sens_g)
fwrite(sens_g[, .(analysis, sensitivity, year, log_b, se, RR, CI95_lwr, CI95_upr, p_calibrated_unilateral, p_calibrated_twosided)],
       file.path(DIR$mod, "table_cohort_group_sensitivity.csv"))
cat("\n== sensitivity analyses of the two groups\n")
print(sens_g[, .(analysis, sensitivity, year, RR = round(RR, 3), CI95_lwr = round(CI95_lwr, 3), CI95_upr = round(CI95_upr, 3),
                 p_cal = signif(p_calibrated_unilateral, 2))])

# ---- C. cohort by YEAR OF BIRTH (sih_birth_cohort_aggregate.rds, from script 12) -----
# The cache has no NASC: the cohort is the one EXPECTED from age and hospitalization date (see 12);
# the aggregate's `method` column says which was used. Eligible = born 2010-2015 (aged
# 10-14 in completed years at some point in Feb/2024-Dec/2025; the 2009 cohort is left out).
#   (i)  older neighbouring cohorts 2005-2008, full panel 2017-2025 (main analysis of the block)
#   (ii) neighbouring cohorts 2005-2008 and 2016-2019, without cohort × month cells with age < 5
# Same two contrasts as block B; by year re-anchored on 2017-2023 (6 df) and 2017-2022 (5 df);
# pooled Feb/2024-Dec/2025 (model with one post indicator, and re-anchored 2024-2025 mean
# with calibrated p); sensitivity without Feb-Apr/2025. + minimum detectable effect (MDE).
FILE_BC <- file.path(DIR$processed, "sih_birth_cohort_aggregate.rds")
if (file.exists(FILE_BC)) {
  cn0 <- readRDS(FILE_BC)
  method_cohort <- paste(unique(cn0$method), collapse = "; ")
  cn <- cn0[, .(admissions = sum(admissions)), by = .(code_muni, cohort, month)]; rm(cn0)
  EXP <- 2010:2015; OLDER <- 2005:2008; YOUNGER <- 2016:2019
  SET_C <- list("(i) older neighbors 2005-2008, full panel" = list(adj = OLDER, min_age = -Inf),
                 "(ii) neighbors 2005-2008 and 2016-2019, without cells with age < 5" = list(adj = c(OLDER, YOUNGER), min_age = 5L))
  mean_reanch <- function(m, tab, scale, latest) {        # re-anchored mean of 2024 and 2025, with calibrated p
    nms <- sprintf("year::%d:v", years[-1]); b <- coef(m)[nms]; V <- vcov(m)[nms, nms]
    n_pre <- sum(years <= latest); w <- rep(0, length(b)); w[years[-1] <= latest] <- -1 / n_pre
    w[years[-1] %in% 2024:2025] <- w[years[-1] %in% 2024:2025] + 0.5
    e <- lincomb(w, b, V); sdp <- tab$sd_log_placebos[1]; tt <- e[["est"]] / (sdp * sqrt(1 + 1 / n_pre))
    data.table(year = NA_integer_, log_b = e[["est"]], se = e[["se"]], RR = exp(scale * e[["est"]]),
               CI95_lwr = exp(scale * (e[["est"]] - 1.96 * e[["se"]])), CI95_upr = exp(scale * (e[["est"]] + 1.96 * e[["se"]])),
               t_calibrated = tt, p_calibrated_unilateral = pt(tt, df = n_pre - 1))
  }
  resC <- list()
  for (nm in names(SET_C)) for (samp in c("expn", "phase1")) {
    cj <- SET_C[[nm]]
    munis <- if (samp == "expn") mun[expn == 1L, code_muni] else mun[expn == 0L, code_muni]
    d <- cn[code_muni %in% munis & cohort %in% c(EXP, cj$adj)]
    d <- d[(as.integer(format(month, "%Y")) - cohort) >= cj$min_age]
    d[, group := as.integer(cohort %in% EXP)]
    g <- CJ(code_muni = munis, group = 0:1, month = months)
    g <- merge(g, d[, .(admissions = sum(admissions)), by = .(code_muni, group, month)], by = c("code_muni", "group", "month"), all.x = TRUE)
    g[is.na(admissions), admissions := 0]
    g[, `:=`(year = as.integer(format(month, "%Y")),
             t = (as.integer(format(month, "%Y")) - PARAM$year_start) * 12L + as.integer(format(month, "%m")))]
    g <- merge(g, mun, by = "code_muni")
    g[, tot_mt := sum(admissions), by = .(code_muni, t)]; g <- g[tot_mt > 0]
    if (samp == "expn") { g[, v := group * cov_final]; scale <- 0.10
      lbl <- "Expanded: exposed cohort x final coverage, per 10 p.p."
    } else { g[, v := group * treated]; scale <- 1; lbl <- "first phase vs never: exposed cohort x treated" }
    fml_year <- as.formula(sprintf("admissions ~ i(year, v, ref = %d) | %s", years[1], FE_G))
    m <- fepois(fml_year, data = g, cluster = ~code_muni, notes = FALSE)
    r6 <- reanchor_by_year(m, "v", scale = scale); r5 <- reanchor_by_year(m, "v", scale = scale, last_without_vaccine = 2022L)
    m_without <- fepois(fml_year, data = g[!(month >= as.Date("2025-02-01") & month <= as.Date("2025-04-01"))], cluster = ~code_muni, notes = FALSE)
    rs <- reanchor_by_year(m_without, "v", scale = scale)
    g[, vpos := v * (month >= PARAM$vaccine_start)]
    mp <- fepois(as.formula(sprintf("admissions ~ vpos | %s", FE_G)), data = g, cluster = ~code_muni, notes = FALSE)
    lp <- rr_row(mp, "vpos", lbl, scale)
    base <- list(cohortset = nm, contrast = lbl, municipalities = uniqueN(g$code_muni),
                 admissions_exposed_2025 = g[group == 1L & year == 2025L, sum(admissions)])
    resC[[paste(nm, samp)]] <- rbind(
      r6[, .(specification = "by year, reference 2017-2023", gl = 6L, year, log_b, se, RR, CI95_lwr, CI95_upr, t_calibrated, p_calibrated_unilateral)],
      r5[, .(specification = "by year, reference 2017-2022", gl = 5L, year, log_b, se, RR, CI95_lwr, CI95_upr, t_calibrated, p_calibrated_unilateral)],
      rs[, .(specification = "by year, without Feb-Apr/2025 (ref. 2017-2023)", gl = 6L, year, log_b, se, RR, CI95_lwr, CI95_upr, t_calibrated, p_calibrated_unilateral)],
      data.table(specification = "pooled 2024-2025: re-anchored mean (ref. 2017-2023)", gl = 6L, mean_reanch(m, r6, scale, 2023L)),
      data.table(specification = "pooled 2024-2025: re-anchored mean (ref. 2017-2022)", gl = 5L, mean_reanch(m, r5, scale, 2022L)),
      data.table(specification = "pooled Feb/2024-Dec/2025: post indicator (model)", gl = NA_integer_, year = NA_integer_,
                 log_b = log(lp$RR) / scale, se = lp$se_log, RR = lp$RR, CI95_lwr = lp$CI95_lwr, CI95_upr = lp$CI95_upr,
                 t_calibrated = NA_real_, p_calibrated_unilateral = NA_real_), fill = TRUE)[, (names(base)) := base]
  }
  resC <- rbindlist(resC, fill = TRUE)[, method_cohort := method_cohort]
  fwrite(resC, file.path(DIR$mod, "table_birth_cohort_exact.csv"))
  entry <- data.table(cohort = c(OLDER, EXP, YOUNGER),
                        role = rep(c("older neighbour", "eligible", "younger neighbour"), c(length(OLDER), length(EXP), length(YOUNGER))))
  entry[, year_entry_cohortset_ii := pmax(PARAM$year_start, cohort + 5L)]
  fwrite(entry, file.path(DIR$mod, "table_birth_cohort_exact_entry.csv"))
  cat("\n== cohort by year of birth (", method_cohort, ")\n", sep = "")
  print(resC[is.na(year) | year >= 2024L, .(cohortset = substr(cohortset, 1, 4), contrast = substr(contrast, 1, 10), specification, year,
                                          RR = round(RR, 3), CI95 = sprintf("%.3f-%.3f", CI95_lwr, CI95_upr), p_cal = signif(p_calibrated_unilateral, 2))])

  # ---- minimum detectable effect of the gradient per 10 p.p. in the expansion (2025) ---------
  # MDE = 2.80 × SE(log per 10 p.p.): 80% power, two-sided alpha 5%
  mde_de <- function(cset, se_log_crude) data.table(cohortset = cset, se_log_10pp = 0.10 * se_log_crude,
                                                    mde_log = 2.80 * 0.10 * se_log_crude, mde_RR_reduction = exp(-2.80 * 0.10 * se_log_crude),
                                                    mde_RR_increase = exp(2.80 * 0.10 * se_log_crude))
  mde <- rbind(
    rbindlist(lapply(names(SET_C), function(nm) mde_de(paste("C", nm),
      resC[cohortset == nm & grepl("^Expanded", contrast) & specification == "by year, reference 2017-2023" & year == 2025L, se]))),
    rbindlist(lapply(names(SET_B), function(nm) mde_de(paste("B (current)", nm),
      res_g[grepl("^Expanded", analysis) & grepl(nm, analysis, fixed = TRUE) & year == 2025L, se]))))
  fwrite(mde, file.path(DIR$mod, "table_cohort_mde.csv"))
  cat("\n== minimum detectable effect (2025 gradient in the expanded municipalities)\n"); print(mde)
} else log_msg("WARNING: ", FILE_BC, " missing; block C not estimated")
log_msg("34 done")
