# =============================================================================
# 44_maps.R — maps of Brazil (sf + ggplot2)
#
#   Map 1   study design: first-phase, expansion and never-included
#           municipalities                                      <- panel.rds
#   Map 2   final D1 coverage at ages 10-14 (Dec 2025), by included
#           municipality                                        <- panel.rds
#   Map 3   2017-2025 grid: dengue hospitalizations at ages 10-14 per 100,000
#           person-years, by health region; outline on regions with a
#           first-phase municipality                            <- panel.rds
#   Map 4   2023 -> 2024 change in the hospitalization rate, by health region,
#           at ages 10-14 and in neighboring cohorts (05-09, 15-19) <- panel.rds
# Boundaries: geobr (IBGE 2022, simplified), downloaded once and stored in
# data/raw/geo/ (municipalities ~15 MB, states < 1 MB); health regions are
# dissolved from municipalities with the same brpop table used in the analyses
# (health_region() in 01_functions.R), so that map and model use the same unit.
# No title, subtitle or descriptive caption (they go in the manuscript).
#
# Visual conventions (the same as in 41; dataviz skill):
#   - entity colour: first phase blue #2a78d6, expansion orange #eb6834, never
#     included light grey; blue/orange pair validated for colour blindness;
#   - magnitude = a single hue (blue), light -> dark; the ordinal ramp of Map 2
#     was validated (monotonicity, steps, contrast of the light end >= 2:1);
#   - polarity (Map 4) = blue <-> orange with neutral grey in the middle, the same
#     number of classes in each arm and symmetric cut points on a log scale;
#   - "no data" in white with a thin outline, never a grey resembling a class;
#   - mainland territory only (oceanic islands outside the frame), legend
#     inside the frame, in the empty southwest corner; state borders in thin white
#     (separate without weighing); the first-phase outline is single (dissolved regions).
# Outputs: outputs/figures/Map1..Map4 .png|.tif|.pdf (300 dpi)
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(sf); library(ggplot2); library(patchwork)
sf_use_s2(FALSE)

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; NEUTRAL <- "#e1e0d9"; INK <- "#0b0b0b"; INK2 <- "#52514e"; GRAY <- "#898781"
NEVER <- 10000L
COV   <- c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b")               # ordinal, validated
SEQ   <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#1c5cab", "#0d366b")    # sequential (0 recedes into the background)
DIVERG <- c("#1c5cab", "#5598e7", "#b7d3f6", NEUTRAL, "#f8c9b3", "#eb6834", "#b84a1f")
SUBSET <- coord_sf(xlim = c(-74.2, -34.6), ylim = c(-33.9, 5.4), expand = FALSE, datum = NA)
thm_map <- theme_void(base_size = 9) +
  theme(text = element_text(colour = INK), legend.title = element_text(size = 8, colour = INK, face = "bold"),
        legend.text = element_text(size = 7.5, colour = INK2), legend.key.size = unit(0.4, "cm"),
        legend.key.spacing.y = unit(0.05, "cm"),
        strip.text = element_text(face = "bold", hjust = 0, size = 8.5, colour = INK, margin = margin(0, 0, 2, 0)),
        panel.spacing = unit(0.3, "cm"),
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(4, 4, 4, 4))
leg_within <- theme(legend.position = "inside", legend.position.inside = c(0.02, 0.03), legend.justification = c(0, 0))
leg_bottom  <- theme(legend.position = "bottom", legend.justification = "left", legend.box = "vertical",
                    legend.box.just = "left", legend.title.position = "top")
save <- function(p, name, wd, alt) {
  for (ext in c("png", "tif", "pdf")) {
    fpath <- file.path(DIR$fig, paste0(name, ".", ext))
    if (ext == "pdf") ggsave(fpath, p, width = wd, height = alt, units = "cm", device = cairo_pdf, bg = "white")
    else if (ext == "tif") ggsave(fpath, p, width = wd, height = alt, units = "cm", dpi = 300, bg = "white", device = "tiff", compression = "lzw")
    else ggsave(fpath, p, width = wd, height = alt, units = "cm", dpi = 300, bg = "white", device = "png")
  }
}

# ---- 1. boundaries (downloaded once; read thereafter) ------------------------
DIR_GEO <- file.path(DIR$raw, "geo"); dir.create(DIR_GEO, showWarnings = FALSE, recursive = TRUE)
f_mun <- file.path(DIR_GEO, "municipalities_2022.rds"); f_uf <- file.path(DIR_GEO, "uf_2022.rds"); f_rs <- file.path(DIR_GEO, "health_regions.rds")
if (!file.exists(f_mun) || !file.exists(f_uf)) {
  log_msg("Downloading IBGE boundaries via geobr (once)")
  m <- geobr::read_municipality(year = 2022, simplified = TRUE, showProgress = FALSE)
  m <- st_transform(m, 4674)[, "code_muni"]; m$code_muni <- mun6(m$code_muni)
  saveRDS(m, f_mun)
  u <- st_transform(geobr::read_state(year = 2022, simplified = TRUE, showProgress = FALSE), 4674)[, "abbrev_state"]
  saveRDS(u, f_uf)
}
mun_sf <- readRDS(f_mun); uf_sf <- readRDS(f_uf)
if (!file.exists(f_rs)) {                       # health regions = union of municipalities (brpop table)
  rs <- merge(mun_sf, health_region(), by = "code_muni")
  rs <- aggregate(rs[, "health_reg"], by = list(health_reg = rs$health_reg), FUN = function(x) x[1])
  rs <- st_make_valid(rs[, "health_reg"]); rs$health_reg <- as.character(rs$health_reg)
  saveRDS(rs, f_rs)
}
rs_sf <- readRDS(f_rs)
log_msg(nrow(mun_sf), " municipalities, ", nrow(rs_sf), " health regions, ", nrow(uf_sf), " states")

# ---- 2. data by municipality and by health region ---------------------------
panel <- readRDS(file.path(DIR$processed, "panel.rds"))
mun <- unique(panel[, .(code_muni, g, expn)])
mun[, group := fifelse(g < NEVER, "phase1", fifelse(expn == 1L, "expanded", "never"))]
mun <- merge(mun, panel[agegroup == TARGET_AGEGROUP, .(cov_final = max(cov_d1_l1)), by = code_muni], by = "code_muni")
mun[group == "never", cov_final := NA_real_]
# categories with an explicit "na" level (municipality not included / no data)
with_na <- function(x) factor(fifelse(is.na(x), "na", as.character(x)), levels = c(levels(x), "na"))
LEV_COV <- c("<10%", "10–25%", "25–40%", "40–55%", "≥55%")
mun[, cov_cat := with_na(cut(cov_final, c(-Inf, 0.10, 0.25, 0.40, 0.55, Inf), labels = LEV_COV, right = FALSE))]

reg <- merge(panel[agegroup %in% c("05-09", "10-14", "15-19")], health_region(), by = "code_muni")
reg <- reg[, .(admissions = sum(admissions), pop = sum(pop) / 12), by = .(health_reg, year, target = agegroup == TARGET_AGEGROUP)]
reg[, rate := admissions / pop * 1e5]                                 # per 100,000 person-years
LEV_RT <- c("0–5", "5–10", "10–25", "25–50", "50–100", "≥100")
reg[, rate_cat := with_na(cut(rate, c(0, 5, 10, 25, 50, 100, Inf), labels = LEV_RT, right = FALSE, include.lowest = TRUE))]
var <- dcast(reg[year %in% 2023:2024], health_reg + target ~ year, value.var = "rate")
setnames(var, c("2023", "2024"), c("t23", "t24"))
var[, ratio := fifelse(t23 > 0, t24 / t23, NA_real_)]
# 7 symmetric classes on a log scale (0.8 <-> 1.25; 0.5 <-> 2; 0.25 <-> 4), 3 per arm
LEV_RTO <- paste0("r", 1:7)
var[, ratio_cat := with_na(cut(ratio, c(0, 0.25, 0.5, 0.8, 1.25, 2, 4, Inf), labels = LEV_RTO, right = FALSE))]
# single outline of regions with a first-phase municipality (dissolved: no inner borders)
rs_phase1 <- rs_sf[rs_sf$health_reg %in% health_region()[code_muni %in% mun[group == "phase1", code_muni], health_reg], ]
# the simplified boundaries have gaps between neighboring municipalities; a closing (expand and
# shrink by 0.02°, ~2 km) removes them, so the outline is only the outer perimeter
outline_f1 <- suppressWarnings(suppressMessages(
  st_sf(geometry = st_buffer(st_union(st_buffer(st_make_valid(rs_phase1), 0.02)), -0.02))))

mun_sf  <- merge(mun_sf, mun[, .(code_muni, group, cov_cat)], by = "code_muni", all.x = TRUE)
mun_sf$group[is.na(mun_sf$group)] <- "never"                        # municipality with no population in the panel
mun_sf$cov_cat[is.na(mun_sf$cov_cat)] <- "na"
N_GROUP <- mun[, .N, by = group][match(c("phase1", "expanded", "never"), group), N]
layer_uf <- geom_sf(data = uf_sf, fill = NA, colour = "white", linewidth = 0.3)
layer_uf_thin <- geom_sf(data = uf_sf, fill = NA, colour = "white", linewidth = 0.15)
# "no data" = white with a thin outline; the level enters the legend only if present
class_levels <- function(x, base) c(base, if (any(x == "na")) "na")

{
  fmt_n <- function(x) formatC(x, big.mark = ",", decimal.mark = ".", format = "d")
  TXT <- list(
    t_group = "Vaccination strategy",
    groups = sprintf(c(phase1 = "First phase (n = %s)", expanded = "Expansion (n = %s)", never = "Never included (n = %s)"),
                     fmt_n(N_GROUP)),
    t_cov = "D1 coverage, ages 10–14\n(December 2025)", cov_na = "Not included",
    t_rt = "Dengue hospitalizations, ages 10–14, per 100,000 per year", rt_na = "No population in the age group",
    t_rto = "Ratio of hospitalization rates, 2024 / 2023", rto_na = "No hospitalizations in 2023",
    outline = "Health regions with a first-phase municipality",
    target = "Ages 10–14", viz = "Ages 5–9 and 15–19 (neighboring cohorts)",
    ratio_lab = c("<0.25", "0.25–0.5", "0.5–0.8", "0.8–1.25", "1.25–2", "2–4", "≥4")
  )
  leg_outline <- list(
    geom_sf(data = outline_f1, aes(colour = "c"), fill = NA, linewidth = 0.25, key_glyph = "path"),
    scale_colour_manual(values = c(c = INK), labels = TXT$outline, name = NULL))

  # Map 1 — study design
  p1 <- ggplot(mun_sf) +
    geom_sf(aes(fill = group), colour = NA) + layer_uf + SUBSET +
    scale_fill_manual(values = c(phase1 = BLUE, expanded = ORANGE, never = NEUTRAL), labels = TXT$groups,
                      breaks = c("phase1", "expanded", "never"), name = TXT$t_group) +
    thm_map + leg_within
  save(p1, "Map1", 16, 16)

  # Map 2 — final D1 coverage (included municipalities)
  lv2 <- class_levels(mun_sf$cov_cat, LEV_COV)
  p2 <- ggplot(mun_sf) +
    geom_sf(aes(fill = cov_cat), colour = NA) + layer_uf + SUBSET +
    scale_fill_manual(values = setNames(c(COV, NEUTRAL), c(LEV_COV, "na")), limits = lv2,
                      labels = setNames(c(LEV_COV, TXT$cov_na), c(LEV_COV, "na")), name = TXT$t_cov) +
    thm_map + leg_within
  save(p2, "Map2", 16, 16)

  # Map 3 — annual grid of the rate at ages 10-14 by health region
  g3 <- merge(rs_sf, reg[target == TRUE, .(health_reg, year, rate_cat)], by = "health_reg")
  lv3 <- class_levels(g3$rate_cat, LEV_RT)
  p3 <- ggplot(g3) +
    geom_sf(aes(fill = rate_cat), colour = "white", linewidth = 0.04) + layer_uf_thin +
    leg_outline + facet_wrap(~ year, ncol = 3) + SUBSET +
    scale_fill_manual(values = setNames(c(SEQ, "white"), c(LEV_RT, "na")), limits = lv3,
                      labels = setNames(c(LEV_RT, TXT$rt_na), c(LEV_RT, "na")), name = TXT$t_rt) +
    guides(fill = guide_legend(nrow = 1, order = 1, override.aes = list(colour = ifelse(lv3 == "na", GRAY, NA))),
           colour = guide_legend(order = 2)) +
    thm_map + leg_bottom
  save(p3, "Map3", 18, 21)

  # Map 4 — change 2023 -> 2024, ages 10-14 vs neighboring cohorts
  g4 <- merge(rs_sf, var[, .(health_reg, target, ratio_cat)], by = "health_reg")
  g4$panel <- factor(ifelse(g4$target, TXT$target, TXT$viz), levels = c(TXT$target, TXT$viz))
  lv4 <- class_levels(g4$ratio_cat, LEV_RTO)
  p4 <- ggplot(g4) +
    geom_sf(aes(fill = ratio_cat), colour = "white", linewidth = 0.04) +
    geom_sf(data = g4[g4$ratio_cat == "na", ], fill = "white", colour = GRAY, linewidth = 0.1) +
    layer_uf_thin + leg_outline + facet_wrap(~ panel, ncol = 2) + SUBSET +
    scale_fill_manual(values = setNames(c(DIVERG, "white"), c(LEV_RTO, "na")), limits = lv4,
                      labels = setNames(c(TXT$ratio_lab, TXT$rto_na), c(LEV_RTO, "na")), name = TXT$t_rto) +
    guides(fill = guide_legend(nrow = 1, order = 1, override.aes = list(colour = ifelse(lv4 == "na", GRAY, NA))),
           colour = guide_legend(order = 2)) +
    thm_map + leg_bottom
  save(p4, "Map4", 18, 11.5)
}
log_msg("Maps 1-4 written to ", DIR$fig)
