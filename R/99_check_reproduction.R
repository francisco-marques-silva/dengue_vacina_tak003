# =============================================================================
# 99_check_reproduction.R — did the run reproduce the published results?
#
# Compares outputs/manuscript_numbers.json (freshly generated) with the frozen version
# in docs/reference/manuscript_numbers_reference.json (generated on 2026-09-23,
# the one in the manuscripts), number by number:
#   - equal: relative difference <= TOL (default 1e-6) or identical text;
#   - DIFFERENT: above the tolerance (this is what matters);
#   - only in reference / only in new: keys that one of the versions lacks
#     (expected for what was reorganized; listed for checking).
# Output: outputs/reproduction_check.csv and a summary on the console. Ends with
# an error if any number is DIFFERENT, so that run_all.R flags it.
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })
library(jsonlite)

TOL <- as.numeric(Sys.getenv("TOL_REPRODUCTION", "1e-6"))

# ---- manuscript numbers vs frozen reference -------------------------------
ref  <- fromJSON(file.path(ROOT, "docs", "reference", "manuscript_numbers_reference.json"), simplifyVector = FALSE)
new <- fromJSON(file.path(ROOT, "outputs", "manuscript_numbers.json"), simplifyVector = FALSE)

# flattens the nested list into path -> scalar value pairs
flatten <- function(x, prefix = "") {
  if (is.list(x)) {
    if (!length(x)) return(list())
    nms <- if (is.null(names(x))) sprintf("[%d]", seq_along(x)) else names(x)
    out <- list()
    for (i in seq_along(x)) out <- c(out, flatten(x[[i]], paste0(prefix, if (nzchar(prefix) && !startsWith(nms[i], "[")) "." else "", nms[i])))
    return(out)
  }
  if (is.null(x)) return(setNames(list(NA), prefix))
  setNames(list(x), prefix)
}
a <- flatten(ref); b <- flatten(new)
# keys excluded from the comparison: `_meta` (date/version), `diagnostic_cohort` (table
# reformatted) and `plausibilidade.mensal[*].observed_*` (in the reference they came from
# an earlier version of the monthly series; the valid coefficients are in
# study_events.coefficients_reanchored, which is compared)
EXCLUDE <- "^(_meta|diagnostic_cohort|plausibility\\.monthly\\[[0-9]+\\]\\.observed)"
a <- a[!grepl(EXCLUDE, names(a))]; b <- b[!grepl(EXCLUDE, names(b))]

common <- intersect(names(a), names(b))
cmpr <- rbindlist(lapply(common, function(k) {
  x <- a[[k]]; y <- b[[k]]
  if (is.numeric(x) && is.numeric(y)) {
    d <- if (is.na(x) && is.na(y)) 0 else if (is.na(x) || is.na(y)) Inf else abs(x - y) / max(abs(x), 1e-12)
    data.table(json_key = k, reference = as.character(x), new = as.character(y), dif_relative = d, status = if (d <= TOL) "equal" else "DIFFERENT")
  } else data.table(json_key = k, reference = as.character(x), new = as.character(y), dif_relative = NA_real_,
                    status = if (identical(as.character(x), as.character(y))) "equal" else "DIFFERENT")
}))
only_ref  <- data.table(json_key = setdiff(names(a), names(b)))[, status := rep("only in reference", .N)]
only_new <- data.table(json_key = setdiff(names(b), names(a)))[, status := rep("only in new", .N)]
everything <- rbind(cmpr, only_ref, only_new, fill = TRUE)
fwrite(everything, file.path(ROOT, "outputs", "reproduction_check.csv"))

cat("\n== reproduction check (relative tolerance ", TOL, ")\n", sep = "")
print(everything[, .N, by = status])
dif <- cmpr[status == "DIFFERENT"]
if (nrow(dif)) { cat("\nDIFFERING NUMBERS:\n"); print(dif[order(-dif_relative)], nrows = 200) }
if (nrow(only_ref)) { cat("\nOnly in the reference (sections reorganized or removed):\n"); print(unique(sub("\\[.*$", "", sub("\\..*$", "", only_ref$json_key)))) }
if (nrow(only_new)) { cat("\nOnly in the new run (sections added):\n"); print(unique(sub("\\[.*$", "", sub("\\..*$", "", only_new$json_key)))) }
if (nrow(dif)) stop(nrow(dif), " number(s) differ from the reference — see outputs/reproduction_check.csv")
log_msg("Reproduction confirmed: ", nrow(cmpr), " numbers identical to the reference")
