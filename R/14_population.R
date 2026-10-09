# =============================================================================
# 14_population.R — population by municipality × age group × year (denominators)
#
# Source: revised DATASUS estimates (brpop, source = "datasus2024", up to 2024).
# For 2025, the 2024 estimate is multiplied by the 2024 -> 2025 change in the
# UFRN-PPGDem-LEPP projections (brpop, source = "ufrn"), municipality by municipality and
# age group by age group. The data come bundled in the brpop package (no download).
# Output: data/processed/population.rds
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })

pops <- list()
for (src in c("datasus2024", "ufrn")) {
  p <- as.data.table(brpop::mun_pop_age(source = src))
  p <- p[grepl("^From \\d+", age_group)]                     # drops the "Total" row
  p[, agegroup := to_agegroup(as.integer(sub("^From (\\d+).*", "\\1", age_group)))]
  pops[[src]] <- p[!is.na(agegroup), .(pop = sum(pop, na.rm = TRUE)),
                     by = .(code_muni = mun6(code_muni), year = as.integer(year), agegroup)]
}
years <- PARAM$year_start:PARAM$year_end
pop  <- pops$datasus2024[year %in% years]

ult <- max(pop$year)                                           # last year with a DATASUS estimate (2024)
for (a in setdiff(years, unique(pop$year))) {                  # years without an estimate: 2025
  fac <- merge(pops$ufrn[year == a,   .(code_muni, agegroup, p_a = pop)],
                 pops$ufrn[year == ult, .(code_muni, agegroup, p_u = pop)],
                 by = c("code_muni", "agegroup"))[, .(code_muni, agegroup, f = p_a / p_u)]
  proj <- merge(pop[year == ult], fac, by = c("code_muni", "agegroup"), all.x = TRUE)
  proj[is.na(f) | !is.finite(f), f := 1]
  pop <- rbind(pop, proj[, .(code_muni, year = a, agegroup, pop = round(pop * f))])
  log_msg("Population for ", a, " projected from ", ult, " with the UFRN change")
}

saveRDS(pop, file.path(DIR$processed, "population.rds"))
log_msg("Population: ", uniqueN(pop$code_muni), " municipalities, ", paste(range(pop$year), collapse = "-"))
