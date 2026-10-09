# =============================================================================
# 16_municipality_list.R — official list of municipalities included in the
# dengue vaccination (intention-to-treat exposure)
#
# The list is the versioned INPUT in data/raw/lists/:
#   dengue_vaccine_municipalities_full.csv  2,751 municipalities: code_muni (IBGE, 6 digits),
#                                          uf, name, health_region, phase ("1" | "expanded"),
#                                          start_date ("2024-02-01" in the first phase; "2024-07-01" in the expansion)
#   dengue_vaccine_municipalities.csv           the 521 of the first phase (subset of the above)
# This script only checks that it is intact; it is what 20_panel.R reads.
#
# Origin (documents kept in the same folder, not read by the code):
#   municipalities_included.xlsx     Ministry of Health spreadsheet with the 2,751
#                                    included municipalities (state, municipality, health region)
#   informe_tecnico_dengue_2024.pdf  Operational Technical Report of the strategy; Annex II
#                                    lists the 521 first-phase municipalities (vaccination from
#                                    Feb/2024)
# Phase 1 = listed in Annex II; expansion = in the spreadsheet but not in Annex II. Names
# were matched to the IBGE code using the DATASUS municipality table; two renamed
# municipalities were resolved by hand (Campo Grande/RN = formerly Augusto Severo, 240130;
# Tabocão/TO = formerly Fortaleza do Tabocão, 170825).
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })

DIR_LIST  <- file.path(DIR$raw, "lists")
F_FULL <- file.path(DIR_LIST, "dengue_vaccine_municipalities_full.csv")
F_PHASE1    <- file.path(DIR_LIST, "dengue_vaccine_municipalities.csv")
COLS <- c("code_muni", "uf", "name", "health_region", "phase", "start_date")

if (!file.exists(F_FULL)) stop("Missing ", F_FULL, " — it is part of the repository (data/raw/lists/)")
lst <- fread(F_FULL, colClasses = "character")
stopifnot(all(COLS %in% names(lst)),
          nrow(lst) == 2751L, sum(lst$phase == "1") == 521L, sum(lst$phase == "expanded") == 2230L,
          !anyDuplicated(lst$code_muni), all(nchar(lst$code_muni) == 6L),
          all(lst[phase == "1", start_date] == "2024-02-01"),
          all(lst[phase == "expanded", start_date] == "2024-07-01"))
if (!file.exists(F_PHASE1)) fwrite(lst[phase == "1"], F_PHASE1)      # derived; recreated if missing
log_msg("Official list checked: ", nrow(lst), " municipalities (", sum(lst$phase == "1"),
        " in the first phase, ", sum(lst$phase == "expanded"), " in the expansion)")
