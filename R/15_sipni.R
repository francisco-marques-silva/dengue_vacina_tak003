# =============================================================================
# 15_sipni.R — SI-PNI: dengue vaccine doses administered, by municipality of
# residence × age group × month (Feb/2024 to Dec/2025)
#
# The national monthly SI-PNI microdata (SUS Open Data Portal,
# dataset "doses aplicadas pelo PNI <year>") are ~1.5-2 GB compressed and cover all
# vaccines. One month at a time: download -> keep only dengue vaccine rows ->
# aggregate -> save the (small) aggregate in data/raw/sipni/reduced/ -> delete the
# raw file. Months already reduced are skipped, so with the reduced/ folder filled
# the script downloads nothing and only consolidates.
#
# Row filtering: if tar+findstr (Windows) or unzip+grep (Linux/Mac) are available,
# the zip is filtered with them (~1 min per month); otherwise the zip is read in chunks
# inside R (~20 min per month). The result is the same.
# Output: data/processed/dengue_vaccination.rds (code_muni, agegroup, month, d1, d2)
# =============================================================================
local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
        if (!length(f)) for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }
        source(if (length(f)) file.path(dirname(normalizePath(f[1], mustWork = FALSE)), "00_setup.R")
               else "R/00_setup.R") })

DIR_RED <- file.path(DIR$raw, "sipni", "reduced")
DIR_TMP <- file.path(DIR$raw, "sipni", "tmp")
for (d in c(DIR_RED, DIR_TMP)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
# the same month also aggregated by SINGLE AGE, in its own folder (the consolidation
# by age group below reads only DIR_RED and stays identical) -> dengue_vaccination_age.rds (script 35)
DIR_RED_AGE <- file.path(DIR$raw, "sipni", "reduced_age")
dir.create(DIR_RED_AGE, recursive = TRUE, showWarnings = FALSE)
MONTHS <- c("janeiro","fevereiro","marco","abril","maio","junho","julho","agosto","setembro",
           "outubro","novembro","dezembro")
TARGET  <- rbind(data.table(year = 2024L, m = 2:12), data.table(year = 2025L, m = 1:12))
BASE  <- "https://dadosabertos.saude.gov.br"
# possible column names (they changed between portal versions)
CAND <- list(
  mun   = c("co_municipio_paciente", "paciente_endereco_coibgemunicipio",
            "co_municipio_residencia_paciente", "co_ibge_municipio_paciente", "codigo_municipio_paciente"),
  age = c("nu_idade_paciente", "paciente_idade", "nu_idade", "idade_paciente"),
  vac   = c("ds_vacina", "ds_nome", "no_vacina", "vacina_nome", "sg_vacina", "ds_vacina_nome",
            "descricao_vacina", "imunobiologico", "sg_imunobiologico"),
  dose  = c("ds_dose_vacina", "ds_tipo_dose", "vacina_descricao_dose", "no_dose", "co_dose_vacina",
            "descricao_dose"),
  data  = c("dt_vacina", "vacina_dataaplicacao", "dt_aplicacao", "data_aplicacao"))
key_txt <- function(x) gsub("[^a-z0-9]", "", tolower(stringi::stri_trans_general(enc2utf8(x), "Latin-ASCII")))

pending <- TARGET[!file.exists(file.path(DIR_RED, sprintf("sipni_dengue_%d_%02d.rds", year, m)))]
log_msg(nrow(TARGET) - nrow(pending), " months already reduced; ", nrow(pending), " to download")

# ---- 1. resources (files) for each year on the portal -------------------------
portal_resources <- list()
for (a in unique(pending$year)) {
  slug <- sprintf("doses-aplicadas-pelo-programa-de-nacional-de-imunizacoes-pni-%d", a)
  r <- tryCatch({                                            # CKAN API
    j <- jsonlite::fromJSON(sprintf("%s/api/3/action/package_show?id=%s", BASE, slug), simplifyVector = TRUE)
    as.data.table(j$result$resources)[, .(name = name, url = url, fmt = toupper(format))]
  }, error = function(e) NULL)
  if (is.null(r) || !nrow(r)) {                              # Next.js page: embedded JSON
    html <- tryCatch(paste(readLines(sprintf("%s/dataset/%s", BASE, slug), warn = FALSE, encoding = "UTF-8"),
                           collapse = "\n"), error = function(e) "")
    j <- sub("</script>.*", "", sub('.*id="__NEXT_DATA__"[^>]*>', "", html))
    res <- NULL
    search <- function(x) {                                 # descends until it finds "resources"
      if (is.list(x)) {
        if (!is.null(x$resources) && length(x$resources)) return(x$resources)
        for (e in x) { z <- search(e); if (!is.null(z)) return(z) }
      }
      NULL
    }
    res <- tryCatch(search(jsonlite::fromJSON(j, simplifyVector = FALSE)), error = function(e) NULL)
    if (!is.null(res)) {
      field <- function(k) vapply(res, function(x) if (is.null(x[[k]])) "" else as.character(x[[k]])[1], "")
      r <- data.table(name = field("name"), url = field("url"), fmt = toupper(field("format")))
    } else {                                                 # last resort: links on the page
      hrefs <- unique(regmatches(html, gregexpr('https?://[^"\']+_csv\\.zip', html))[[1]])
      r <- data.table(name = basename(hrefs), url = hrefs, fmt = "CSV")
    }
  }
  portal_resources[[as.character(a)]] <- r
  log_msg("Resources found for ", a, ": ", nrow(r))
}

# ---- 2. one month at a time --------------------------------------------------
has_tar     <- .Platform$OS.type == "windows" && file.exists(file.path(Sys.getenv("SystemRoot", "C:/Windows"), "System32", "tar.exe"))
has_unzip   <- .Platform$OS.type != "windows" && nzchar(Sys.which("unzip")) && nzchar(Sys.which("grep"))
for (i in seq_len(nrow(pending))) {
  year <- pending$year[i]; m <- pending$m[i]; tag <- sprintf("%d_%02d", year, m)
  fout <- file.path(DIR_RED, sprintf("sipni_dengue_%s.rds", tag))
  r <- portal_resources[[as.character(year)]]
  if (is.null(r) || !nrow(r)) { log_msg("[", tag, "] no resource list"); next }

  # finds the month's resource by name ("vacinacaofevereiro2024", "..._2024_02", "vaccination_fev_2024_csv.zip")
  k <- key_txt(r$name)
  hit <- which(grepl(key_txt(MONTHS[m]), k, fixed = TRUE) & grepl(as.character(year), k, fixed = TRUE))
  if (!length(hit)) hit <- which(grepl(sprintf("%d%02d", year, m), k, fixed = TRUE))
  if (!length(hit)) hit <- which(grepl(sprintf("%s%d", substr(MONTHS[m], 1, 3), year), k, fixed = TRUE))
  if (!length(hit)) { log_msg("[", tag, "] resource not found on the portal"); next }
  hit <- hit[order(match(r$fmt[hit], c("CSV", "ZIP", "CSV.ZIP")), na.last = TRUE)]
  url <- r$url[hit[1]]

  # download with up to 3 attempts
  dest <- file.path(DIR_TMP, sprintf("sipni_%s%s", tag, if (grepl("\\.zip", url, ignore.case = TRUE)) ".zip" else ".csv"))
  ok <- FALSE
  for (t in 1:3) {
    log_msg("[", tag, "] downloading", if (t > 1) paste0(" (attempt ", t, ")"), " ", url)
    ok <- tryCatch({
      curl::curl_download(url, dest, quiet = TRUE,
        handle = curl::new_handle(timeout = 0L, connecttimeout = 60L, low_speed_limit = 1024L,
                                  low_speed_time = 300L, tcp_keepalive = TRUE)); TRUE },
      error = function(e) { log_msg("[", tag, "] download error: ", conditionMessage(e)); FALSE })
    if (ok) break
    unlink(dest); if (t < 3) Sys.sleep(c(120, 300)[t])
  }
  if (!ok) { log_msg("[", tag, "] failed 3 times; skipping"); next }
  log_msg("[", tag, "] ", round(file.size(dest) / 1e6), " MB; filtering the dengue vaccine rows")

  # rows with "dengue" (plus the header, which contains "_paciente")
  filt <- file.path(DIR_TMP, sprintf("filtered_%s.csv", tag))
  zip  <- grepl("\\.zip$", dest, ignore.case = TRUE)
  if (zip && has_tar) {
    shell(sprintf('%s -xOf "%s" | findstr /i /c:dengue /c:_paciente > "%s"',
                  gsub("/", "\\\\", file.path(Sys.getenv("SystemRoot", "C:/Windows"), "System32", "tar.exe")),
                  normalizePath(dest, winslash = "\\"), normalizePath(filt, winslash = "\\", mustWork = FALSE)))
  } else if (zip && has_unzip) {
    system(sprintf('unzip -p "%s" | grep -i -e dengue -e _paciente > "%s"', dest, filt))
  }
  if (file.exists(filt) && file.size(filt) > 0) {
    d <- fread(filt, colClasses = "character", encoding = "Latin-1", showProgress = FALSE, fill = TRUE)
  } else {                                                   # portable path: reads in chunks
    con <- if (zip) unz(dest, utils::unzip(dest, list = TRUE)$Name[1]) else file(dest)
    open(con, "r")
    heads <- readLines(con, n = 1L, warn = FALSE); parts <- list(); n <- 0L
    repeat {
      x <- readLines(con, n = 500000L, warn = FALSE)
      if (!length(x)) break
      x <- x[grepl("dengue", x, ignore.case = TRUE)]
      if (length(x)) { n <- n + 1L; parts[[n]] <- x }
    }
    close(con)
    d <- fread(text = c(heads, unlist(parts, use.names = FALSE)), colClasses = "character",
               encoding = "Latin-1", showProgress = FALSE, fill = TRUE)
  }
  unlink(c(dest, filt))

  # columns: picks, among the candidates, the one present in the file
  heads <- names(d)
  cols <- vapply(names(CAND), function(p) {
    h <- CAND[[p]][tolower(CAND[[p]]) %in% tolower(heads)]
    if (!length(h)) NA_character_ else heads[match(tolower(h[1]), tolower(heads))] }, "")
  if (anyNA(cols)) { log_msg("[", tag, "] unrecognized columns: ", paste(heads, collapse = ", ")); next }
  d <- d[, unname(cols), with = FALSE]; setnames(d, unname(cols), names(cols))
  d <- d[grepl("dengue", vac, ignore.case = TRUE)]
  d[, dose_n := fifelse(grepl("^\\s*(1|D1|1ª|Primeira)", dose, ignore.case = TRUE), 1L,
              fifelse(grepl("^\\s*(2|D2|2ª|Segunda)",  dose, ignore.case = TRUE), 2L, NA_integer_))]
  d[, `:=`(code_muni = mun6(mun), month = month_of(as_date(data)),
           agegroup = to_agegroup(suppressWarnings(as.integer(age))))]
  ag <- d[!is.na(month) & !is.na(dose_n), .(d1 = sum(dose_n == 1L), d2 = sum(dose_n == 2L)),
          by = .(code_muni, agegroup, month)]
  saveRDS(d[!is.na(month) & !is.na(dose_n), .(d1 = sum(dose_n == 1L), d2 = sum(dose_n == 2L)),
            by = .(code_muni, age = suppressWarnings(as.integer(age)), month)],
          file.path(DIR_RED_AGE, sprintf("sipni_dengue_age_%s.rds", tag)))   # before fout: fout marks the month as done
  saveRDS(ag, fout)
  log_msg("[", tag, "] ok: ", nrow(ag), " rows; D1 = ", sum(ag$d1), ", D2 = ", sum(ag$d2))
  rm(d, ag); gc()
}

# ---- 3. consolidate ---------------------------------------------------------
fs <- list.files(DIR_RED, "\\.rds$", full.names = TRUE)
if (!length(fs)) stop("No reduced SI-PNI month in ", DIR_RED)
vac <- rbindlist(lapply(fs, readRDS), fill = TRUE)
vac <- vac[, .(d1 = sum(d1), d2 = sum(d2)), by = .(code_muni, agegroup, month)]
saveRDS(vac, file.path(DIR$processed, "dengue_vaccination.rds"))
log_msg("dengue_vaccination.rds: ", length(fs), " months, ", uniqueN(vac$code_muni), " municipalities, D1 = ",
        sum(vac$d1), ", D2 = ", sum(vac$d2))

# consolidates the single-age aggregate
fi <- list.files(DIR_RED_AGE, "\\.rds$", full.names = TRUE)
if (length(fi)) {
  vi <- rbindlist(lapply(fi, readRDS), fill = TRUE)[, .(d1 = sum(d1), d2 = sum(d2)), by = .(code_muni, age, month)]
  saveRDS(vi, file.path(DIR$processed, "dengue_vaccination_age.rds"))
  log_msg("dengue_vaccination_age.rds: ", length(fi), " months, D1 = ", sum(vi$d1), ", D2 = ", sum(vi$d2))
}
if (length(fi) != length(fs)) log_msg("WARNING: ", length(fs), " months by age group and ", length(fi), " by single age")
