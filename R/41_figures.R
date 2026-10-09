# =============================================================================
# 41_figures.R — manuscript figures (ggplot2)
#
#   Figure 1  monthly hospitalization rates (10-14 vs comparison) by municipality
#             group and the ratio between them             <- panel.rds
#   Figure 2  dose-response by year (expansion, all, first phase) <- 33
#   Figure S1 monthly event study, two anchorings              <- 30
#   Figure S2 estimates by year and by season, common reference <- 31
#   Figure S3 monthly doses (first phase and expansion) and target-age hospitalizations <- 35
#   Figure S4 selection flow chart of municipalities and hospitalizations <- panel.rds
# No title, subtitle or descriptive caption (those go in the manuscript).
#
# Visual conventions (the same in all figures and in the maps of script 44):
#   - color follows the ENTITY, never a local meaning: first phase = blue #2a78d6,
#     expansion = orange #eb6834, never included = gray; the blue/orange pair
#     was validated (color blindness and normal vision) with the dataviz skill validator;
#   - estimates do not use group color: years with vaccine in dark ink, years without
#     vaccine in gray (emphasis, not identity);
#   - one y axis per panel (never a dual axis); the vaccination period is a neutral
#     background band, not a dashed line;
#   - direct labels at the end of the lines when there are few series; grid and axes in
#     thin gray strokes; text always in ink, never in the series color.
# Outputs: outputs/figures/Fig1.png|.tif|.pdf, Fig2, FigS1-S5 (300 dpi;
# FigS3 png and tif only)
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(ggplot2); library(patchwork)

# ---- palette and ink (reference instance of the dataviz skill) ------------------
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"                                  # entities (validated)
BLUE_LIGHT <- "#9ec5f4"; ORANGE_LIGHT <- "#f5b89e"; GRAY_LIGHT <- "#c3c2b7"
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRAY <- "#898781"; GRID <- "#e1e0d9"; BASE <- "#c3c2b7"
BACKGROUND_VAC <- "#f0efec"                                                   # vaccination period band
NEVER <- 10000L
INI <- PARAM$vaccine_start
thm <- theme_minimal(base_size = 9) +
  theme(text = element_text(colour = INK), axis.text = element_text(colour = INK2, size = 7.5),
        axis.title = element_text(colour = INK2, size = 8.5),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = GRID, linewidth = 0.25),
        axis.line.x = element_line(colour = BASE, linewidth = 0.3), axis.ticks.x = element_line(colour = BASE, linewidth = 0.3),
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        legend.text = element_text(colour = INK2, size = 8), legend.key.width = unit(0.6, "cm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -4, 0),
        strip.text = element_text(face = "bold", hjust = 0, size = 8.5, colour = INK),
        plot.tag = element_text(face = "bold", size = 11), panel.spacing = unit(0.45, "cm"),
        plot.caption = element_text(colour = INK2, size = 7, hjust = 0),
        plot.background = element_rect(fill = "white", colour = NA))
save <- function(p, name, wd, alt, exts = c("png", "tif", "pdf")) {
  for (ext in exts) {
    fpath <- file.path(DIR$fig, paste0(name, ".", ext))
    if (ext == "pdf") ggsave(fpath, p, width = wd, height = alt, units = "cm", device = cairo_pdf, bg = "white")
    else if (ext == "tif") ggsave(fpath, p, width = wd, height = alt, units = "cm", dpi = 300, bg = "white", device = "tiff", compression = "lzw")
    else ggsave(fpath, p, width = wd, height = alt, units = "cm", dpi = 300, bg = "white", device = "png")
  }
}
# ymin = 0 (not -Inf): on a log axis, log10(0) = -Inf reaches the edge; -Inf would become NaN and the band would vanish
background_vac <- function(x0, x1) annotate("rect", xmin = x0, xmax = x1, ymin = 0, ymax = Inf, fill = BACKGROUND_VAC)

panel <- readRDS(file.path(DIR$processed, "panel.rds"))[expn == 0L]
da <- fread(file.path(DIR$mod, "table_crosssec_dose_annual.csv"))
ev <- fread(file.path(DIR$mod, "table_event_study_reanchored.csv"))
ea <- fread(file.path(DIR$mod, "table_event_annual.csv"))
et <- fread(file.path(DIR$mod, "table_event_season.csv"))
# monthly series without the bin from 31 (whole pre-period) and the seasonally matched one (Figures S1 and S5)
le <- function(nm) { f <- file.path(DIR$mod, nm); if (file.exists(f)) fread(f) else NULL }
ev31 <- le("table_event_monthly_all_pre.csv"); seas31 <- le("table_event_monthly_seasonal.csv")
# Figure S3 (monthly doses, from 35) and Figure S4 (selection flow chart)
fpath_dm <- file.path(DIR$mod, "table_doses_monthly.csv")
dm <- if (file.exists(fpath_dm)) fread(fpath_dm)[, month := as.Date(month)] else NULL
pfull <- readRDS(file.path(DIR$processed, "panel.rds"))
flow_n <- list(total = uniqueN(pfull$code_muni), f1 = pfull[g < NEVER, uniqueN(code_muni)],
                expn = pfull[expn == 1L, uniqueN(code_muni)], never = pfull[g >= NEVER & expn == 0L, uniqueN(code_muni)],
                int_5a59 = panel[, sum(admissions)], int_target = panel[agegroup == "10-14", sum(admissions)])
rm(pfull)
END <- max(panel$month)

{
  big <- ","; dec <- "."
  fmtn <- function(x) formatC(as.numeric(x), format = "f", digits = 0, big.mark = big)
  fmtr <- function(x) sub("[.,]$", "", sub("0+$", "", formatC(as.numeric(x), format = "f", digits = 2, decimal.mark = dec)))
  fmt1 <- function(x) formatC(as.numeric(x), format = "f", digits = 1, decimal.mark = dec)
  TXT <- list(
    g1 = "First-phase municipalities (n = 521)", g0 = "Never-included municipalities (n = 2,819)",
    g1c = "First phase", g0c = "Never\nincluded",
    target = "Ages 10–14", cmpr = "Comparison",
    y_rt = "Dengue hospitalizations\nper 100,000 per month", y_rto = "Rate ratio\nages 10–14 / comparison",
    cap_A = "Comparison: ages 5–9, 15–19 and 20–29. Shaded area: vaccination period (from February 2024).",
    cap_B = "Thick line: 3-month moving average; thin line: monthly value.",
    samples = c("expanded (2,230)" = "Expansion municipalities (n = 2,230)", "all municipalities" = "All municipalities (n = 5,570)",
                 "first phase (521)" = "First-phase municipalities (n = 521)"),
    y_dose = "RR per 10 p.p. of final D1 coverage (95% CI)", with = "Vaccine years", without = "Pre-vaccine years",
    anc1 = "Reference: month −1 (January 2024)", anc2 = "Reference: pre-period mean (weighted by months)",
    x_k = "Months relative to the start of vaccination (February 2024 = 0); month −24 pools all earlier months",
    with_k = "Vaccine months", without_k = "Pre-vaccine months",
    y_rrt = "Ratio of rate ratios (95% CI)", year_cal = "Calendar year", temp = "Season (February–January)")
  x_dates <- scale_x_date(breaks = as.Date(sprintf("%d-01-01", PARAM$year_start:PARAM$year_end)), date_labels = "%Y",
                          expand = expansion(mult = c(0.01, 0.12)))
  colors_vac <- function(without, with) scale_colour_manual(values = setNames(c(GRAY, INK), c(without, with)))

  # ---- Figure 1 ----------------------------------------------------------------
  # A: color = municipality group (entity); intensity = age group (strong = 10-14,
  #    light = comparison), identified by a direct label at the end of each line
  d <- panel[agegroup %in% c("05-09", "10-14", "15-19", "20-29")]
  d[, group := factor(fifelse(g < NEVER, TXT$g1, TXT$g0), levels = c(TXT$g1, TXT$g0))]
  d[, cat := fifelse(agegroup == "10-14", "target", "cmpr")]
  s <- d[, .(rt = 1e5 * sum(admissions) / sum(pop)), by = .(group, cat, month)]
  s[, series := paste(as.integer(group), cat)]
  COLORS_A <- c("1 target" = BLUE, "1 cmpr" = BLUE_LIGHT, "2 target" = INK2, "2 cmpr" = GRAY_LIGHT)
  # end labels: the largest value on top, the smallest below (no overlap)
  fimA <- s[month == END][order(group, -rt)][, `:=`(lab = fifelse(cat == "target", TXT$target, TXT$cmpr),
                                                   vj = c(-0.3, 1.3)[seq_len(.N)]), by = group]
  pA <- ggplot(s, aes(month, rt, colour = series)) + background_vac(INI, END + 15) +
    geom_line(aes(linewidth = cat)) + facet_wrap(~group, ncol = 1, scales = "free_y") +
    geom_text(data = fimA, aes(x = END + 25, y = rt, label = lab, vjust = vj), hjust = 0, size = 2.4, colour = INK2) +
    scale_colour_manual(values = COLORS_A, guide = "none") +
    scale_linewidth_manual(values = c(target = 0.7, cmpr = 0.45), guide = "none") +
    x_dates + scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    coord_cartesian(clip = "off") + labs(x = NULL, y = TXT$y_rt, tag = "A", caption = TXT$cap_A) + thm
  # B: ratio by group; thin line = monthly, thick = 3-month moving average
  r <- dcast(s, group + month ~ cat, value.var = "rt")
  r[, rto := target / cmpr]; setorder(r, group, month)
  r[, rto_s := frollmean(rto, 3, align = "center"), by = group]
  fimB <- r[!is.na(rto_s), .SD[.N], by = group][, lab := fifelse(group == TXT$g1, TXT$g1c, TXT$g0c)]
  fimB[, vj := fifelse(rto_s == max(rto_s), -0.2, 1.1)]
  pB <- ggplot(r, aes(month, rto, colour = group)) + background_vac(INI, END + 15) +
    geom_hline(yintercept = 1, colour = BASE, linewidth = 0.35) +
    geom_line(linewidth = 0.25, alpha = 0.45) + geom_line(aes(y = rto_s), linewidth = 0.75) +
    geom_text(data = fimB, aes(x = END + 25, y = rto_s, label = lab, vjust = vj), hjust = 0, size = 2.4,
              colour = INK2, lineheight = 0.85) +
    scale_colour_manual(values = setNames(c(BLUE, INK2), c(TXT$g1, TXT$g0)), guide = "none") +
    x_dates + scale_y_continuous(limits = c(0, 3.2), breaks = seq(0, 3, 0.5), labels = fmt1) +
    coord_cartesian(clip = "off") + labs(x = NULL, y = TXT$y_rto, tag = "B", caption = TXT$cap_B) + thm
  save(pA / pB + plot_layout(heights = c(2, 1.15)), "Fig1", 18, 19.5)

  # ---- Figure 2: estimates by year, emphasis on years with vaccine --------------------
  pf <- da[samp %in% names(TXT$samples)]
  pf[, `:=`(panel = factor(TXT$samples[samp], levels = TXT$samples),
            vaccine = factor(fifelse(year >= 2024, TXT$with, TXT$without), levels = c(TXT$without, TXT$with)))]
  p2 <- ggplot(pf, aes(year, RR_10pp, colour = vaccine)) + background_vac(2023.5, 2025.5) +
    geom_hline(yintercept = 1, colour = BASE, linewidth = 0.35) +
    geom_linerange(aes(ymin = CI95_lwr, ymax = CI95_upr), linewidth = 0.5) + geom_point(size = 1.9) +
    facet_wrap(~panel, ncol = 1) + colors_vac(TXT$without, TXT$with) +
    scale_y_log10(breaks = c(0.85, 0.9, 0.95, 1, 1.05, 1.1), labels = fmtr) +
    scale_x_continuous(breaks = 2017:2025, expand = expansion(add = 0.4)) + labs(x = NULL, y = TXT$y_dose) + thm
  save(p2, "Fig2", 16, 18.5)

  # ---- Figure S1: one panel per anchoring, rather than two overlaid series ----------
  e <- rbind(ev[, .(k, RR = RR_ref_m1, lwr = lwr_ref_m1, upr = upr_ref_m1, anchor = TXT$anc1)],
             data.table(k = -1L, RR = 1, lwr = 1, upr = 1, anchor = TXT$anc1),
             ev[, .(k, RR = RR_reanch, lwr = lwr_reanch, upr = upr_reanch, anchor = TXT$anc2)])
  # bottom panel: month -1 comes from the series without the bin from 31 (the one from 30 lacks it)
  if (!is.null(ev31)) e <- rbind(e, ev31[k == -1L, .(k, RR = RR_reanch, lwr = lwr_reanch, upr = upr_reanch, anchor = TXT$anc2)])
  e[, `:=`(anchor = factor(anchor, levels = c(TXT$anc1, TXT$anc2)),
           vaccine = factor(fifelse(k >= 0, TXT$with_k, TXT$without_k), levels = c(TXT$without_k, TXT$with_k)))]
  pS1 <- ggplot(e, aes(k, RR, colour = vaccine)) + background_vac(-0.5, Inf) +
    geom_hline(yintercept = 1, colour = BASE, linewidth = 0.35) +
    geom_linerange(aes(ymin = lwr, ymax = upr), linewidth = 0.4) + geom_point(size = 1.3) +
    facet_wrap(~anchor, ncol = 1) + colors_vac(TXT$without_k, TXT$with_k) +
    scale_y_log10(breaks = c(0.25, 0.5, 1, 2), labels = fmtr) +
    scale_x_continuous(breaks = seq(-24, 20, 4), expand = expansion(add = 0.8)) +
    labs(x = TXT$x_k, y = TXT$y_rrt) + thm
  save(pS1, "FigS1", 18.5, 13)

  # ---- Figure S2: calendar year and season, each with its own x axis --------------
  panel_s2 <- function(tab, title) {
    x <- copy(tab)[, `:=`(post = seq_len(.N),
                          vaccine = factor(fifelse(grepl("^202[45]", period), TXT$with, TXT$without), levels = c(TXT$without, TXT$with)))]
    ggplot(x, aes(post, RRR, colour = vaccine)) + background_vac(min(x[vaccine == TXT$with, post]) - 0.5, max(x$post) + 0.5) +
      geom_hline(yintercept = 1, colour = BASE, linewidth = 0.35) +
      geom_linerange(aes(ymin = CI95_lwr, ymax = CI95_upr), linewidth = 0.5) + geom_point(size = 1.9) +
      colors_vac(TXT$without, TXT$with) +
      scale_x_continuous(breaks = x$post, labels = x$period, expand = expansion(add = 0.4)) +
      scale_y_log10(breaks = c(0.6, 0.8, 1, 1.2, 1.4), labels = fmtr, limits = range(c(pl_lim, 1))) +
      labs(x = NULL, y = TXT$y_rrt, subtitle = title) + thm +
      theme(plot.subtitle = element_text(face = "bold", size = 8.5, colour = INK))
  }
  pl_lim <- range(c(ea$CI95_lwr, ea$CI95_upr, et$CI95_lwr, et$CI95_upr))
  pS2 <- panel_s2(ea, TXT$year_cal) / panel_s2(et, TXT$temp) + plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = "top", legend.justification = "left")
  save(pS2, "FigS2", 16, 14)

  # ---- texts for figures S3 and S4 (American English in the English version) ----
  TXT2 <- list(
    f3_y1 = "Dengue hospitalizations\nages 10–14, first phase", f3_y2 = "Doses administered\nages 10–14",
    f3_lab = c(d1_first_phase = "D1 · first phase", d2_first_phase = "D2 · first phase",
               d1_expanded = "D1 · expansion", d2_expanded = "D2 · expansion"),
    f3_cap = sprintf("First phase: %s municipalities; expansion: %s. Shaded area: vaccination period.",
                     fmtn(flow_n$f1), fmtn(flow_n$expn)),
    f4 = c(sprintf("%s municipalities", fmtn(flow_n$total)),
           sprintf("%s first-phase\nmunicipalities", fmtn(flow_n$f1)),
           sprintf("%s never-included\nmunicipalities", fmtn(flow_n$never)),
           sprintf("%s expansion\nmunicipalities", fmtn(flow_n$expn)),
           sprintf("Main model (triple difference)\n%s municipalities", fmtn(flow_n$f1 + flow_n$never)),
           "Dose–response by coverage\n(outside the main model)",
           sprintf("Dengue hospitalizations, ages 5–59\n(main sample): %s", fmtn(flow_n$int_5a59)),
           sprintf("Target group, ages 10–14: %s", fmtn(flow_n$int_target))))

  # ---- Figure S3: two panels aligned in time, each with its own axis -------
  # (a dual y axis would suggest an arbitrary correlation between the two scales)
  if (!is.null(dm)) {
    lim_x <- c(min(dm$month) - 15, END + 15)
    x_s3 <- scale_x_date(limits = c(lim_x[1], lim_x[2] + 150), date_breaks = "6 months",
                         labels = function(x) format(x, "%m/%Y"), expand = expansion(mult = 0))
    top <- ggplot(dm, aes(month, admissions_target_first_phase)) + background_vac(INI - 15, END + 15) +
      geom_col(fill = BLUE, width = 22) +
      scale_y_continuous(labels = fmtn, expand = expansion(mult = c(0, 0.05))) + x_s3 +
      labs(x = NULL, y = TXT2$f3_y1) + thm + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    dl <- melt(dm[month >= INI, .(month, d1_first_phase, d2_first_phase, d1_expanded, d2_expanded)], id.vars = "month",
               variable.factor = FALSE)
    dl[, `:=`(lab = TXT2$f3_lab[variable], dose = substr(variable, 1, 2))]
    COLORS_S3 <- c(d1_first_phase = BLUE, d2_first_phase = BLUE_LIGHT, d1_expanded = ORANGE, d2_expanded = ORANGE_LIGHT)
    end3 <- dl[month == max(month)][order(-value)]
    # spaces out the end labels so they do not overlap (minimum of 6% of the scale between them)
    step <- 0.06 * max(dl$value); y <- end3$value
    for (i in seq_along(y)[-1]) y[i] <- min(y[i], y[i - 1] - step)
    end3[, y_lab := y]
    bot <- ggplot(dl, aes(month, value, colour = variable)) + background_vac(INI - 15, END + 15) +
      geom_line(aes(linewidth = dose)) +
      geom_segment(data = end3, aes(x = month, xend = month + 20, y = value, yend = y_lab), linewidth = 0.25, colour = BASE) +
      geom_text(data = end3, aes(x = month + 25, y = y_lab, label = lab), hjust = 0, size = 2.4, colour = INK2) +
      scale_colour_manual(values = COLORS_S3, labels = TXT2$f3_lab, breaks = names(COLORS_S3)) +
      scale_linewidth_manual(values = c(d1 = 0.75, d2 = 0.5), guide = "none") +
      scale_y_continuous(labels = fmtn, expand = expansion(mult = c(0, 0.05))) + x_s3 +
      guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
      labs(x = NULL, y = TXT2$f3_y2, caption = TXT2$f3_cap) + thm + theme(legend.position = "bottom")
    save(top / bot + plot_layout(heights = c(1, 1.3)), "FigS3", 18, 14, exts = c("png", "tif"))
  }

  # ---- Figure S4: selection flow chart ---------------------------------
  nb <- data.table(x = c(0, -2.3, 0, 2.3, -1.15, 2.3, -1.15, -1.15), y = c(4, 2.6, 2.6, 2.6, 1.2, 1.2, -0.1, -1.3),
                   lab = TXT2$f4, border = c(INK2, BLUE, GRAY, ORANGE, INK2, INK2, INK2, INK2))
  seg <- data.table(de = c(1, 1, 1, 2, 3, 4, 5, 7), to = c(2, 3, 4, 5, 5, 6, 7, 8))
  seg[, `:=`(x = nb$x[de], y = nb$y[de] - 0.42, xend = nb$x[to], yend = nb$y[to] + 0.42)]
  pS4 <- ggplot() +
    geom_segment(data = seg, aes(x, y, xend = xend, yend = yend), colour = GRAY, linewidth = 0.35,
                 arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
    geom_label(data = nb, aes(x, y, label = lab, colour = I(border)), size = 2.9, lineheight = 0.95, fill = "white",
               label.padding = unit(2.5, "mm"), label.r = unit(1, "mm"), text.colour = INK) +
    coord_cartesian(xlim = c(-3.6, 3.6), ylim = c(-1.8, 4.5)) + theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  save(pS4, "FigS4", 18, 12)

  # ---- Figure S5: monthly coefficients Mar/2023-Dec/2025, two re-anchorings ----
  if (!is.null(ev31) && !is.null(seas31)) {
    TXT5 <- list(r1 = "Reference: mean of all pre-period months (Jan 2017–Dec 2023)",
                      r2 = "Reference: mean of the same calendar months in 2017–2023 (seasonally matched)")
    month_k <- function(k) { m <- as.integer(format(INI, "%Y")) * 12L + as.integer(format(INI, "%m")) - 1L + k
                           as.Date(sprintf("%d-%02d-01", m %/% 12L, m %% 12L + 1L)) }
    s5 <- rbind(ev31[, .(k, RR = RR_reanch, lwr = lwr_reanch, upr = upr_reanch, ref = TXT5$r1)],
                seas31[type == "monthly", .(k, RR, lwr = CI95_lwr, upr = CI95_upr, ref = TXT5$r2)])
    s5 <- s5[k >= -11L]                                                      # from Mar/2023 onwards
    s5[, `:=`(month = month_k(k), ref = factor(ref, levels = c(TXT5$r1, TXT5$r2)),
              vaccine = factor(fifelse(k >= 0, TXT$with_k, TXT$without_k), levels = c(TXT$without_k, TXT$with_k)))]
    pS5 <- ggplot(s5, aes(month, RR, colour = vaccine)) + background_vac(INI - 15, END + 15) +
      geom_hline(yintercept = 1, colour = BASE, linewidth = 0.35) +
      geom_linerange(aes(ymin = lwr, ymax = upr), linewidth = 0.4) + geom_point(size = 1.4) +
      facet_wrap(~ref, ncol = 1) + colors_vac(TXT$without_k, TXT$with_k) +
      scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = fmtr) +
      scale_x_date(date_breaks = "3 months", labels = function(x) format(x, "%m/%Y"), expand = expansion(add = 20)) +
      labs(x = NULL, y = TXT$y_rrt) + thm
    save(pS5, "FigS5", 18.5, 13)
  }
  log_msg("figures written to ", DIR$fig)
}
