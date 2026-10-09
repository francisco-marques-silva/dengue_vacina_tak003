# =============================================================================
# 01_functions.R — the project's only functions (loaded by 00_setup.R)
#
# Only what does not fit as linear code inside a script:
#   A. study constants (states, age groups, dengue ICD codes)
#   B. conversions used by several scripts (age, date, IBGE code, log)
#   C. local cache of DATASUS microdata (download once, read from disk later)
#   D. two statistical calculations repeated in four analyses: linear
#      combination of coefficients with CI (lincomb) and re-anchoring of per-year
#      coefficients on the mean of the pre-vaccine years (reanchor_by_year)
#   E. reading a coefficient as a rate ratio (rr_row) and the health region
#      table (health_region)
# =============================================================================

# ---- A. constants -----------------------------------------------------------
UFS <- c("AC","AL","AM","AP","BA","CE","DF","ES","GO","MA","MG","MS","MT","PA",
         "PB","PE","PI","PR","RJ","RN","RO","RR","RS","SC","SE","SP","TO")

# panel age groups: 10-14 is the target; the others are comparisons
AGEGROUP_CUTS  <- c(5, 10, 15, 20, 30, 40, 50, 60)
AGEGROUP_LABELS <- c("05-09", "10-14", "15-19", "20-29", "30-39", "40-49", "50-59")
TARGET_AGEGROUP     <- "10-14"

ICD_DENGUE <- "^A9(0|1|7)"     # A90, A91 (old ICD-10) and A97.x (current ICD-10)

# fixed effects of the main model: municipality×month, municipality×age group, age group×month
FE_MAIN <- "code_muni^t + code_muni^agegroup + agegroup^t"

# ---- B. conversions ---------------------------------------------------------
to_agegroup <- function(age_years)
  as.character(cut(age_years, breaks = AGEGROUP_CUTS, labels = AGEGROUP_LABELS, right = FALSE))

# SINAN — NU_IDADE_N: 4 digits; 1st = unit (1 hour, 2 day, 3 month, 4 year)
age_sinan <- function(x) {
  x <- suppressWarnings(as.integer(as.character(x)))
  fifelse(x %/% 1000L == 4L, x %% 1000L, fifelse(x %/% 1000L %in% 1:3, 0L, NA_integer_))
}
# SIM — IDADE: 3 digits; 1st = unit (0 min, 1 h, 2 day, 3 month, 4 year, 5 = 100+ years)
age_sim <- function(x) {
  x <- suppressWarnings(as.integer(as.character(x)))
  un <- x %/% 100L; val <- x %% 100L
  fifelse(un == 4L, val, fifelse(un == 5L, 100L + val, fifelse(un %in% 0:3, 0L, NA_integer_)))
}
# SIH — IDADE + COD_IDADE (2 day, 3 month, 4 year, 5 hundreds of years)
age_sih <- function(age, cod) {
  age <- suppressWarnings(as.integer(as.character(age))); cod <- as.character(cod)
  fifelse(cod == "4", age, fifelse(cod == "5", 100L + age,
          fifelse(cod %in% c("2", "3"), 0L, NA_integer_)))
}

# dates: accepts Date, "YYYY-MM-DD", "YYYYMMDD" or "DDMMYYYY"
as_date <- function(x, format_short = c("ymd", "dmy")) {
  format_short <- match.arg(format_short)
  if (inherits(x, "Date")) return(x)
  x <- trimws(as.character(x))
  out <- as.Date(rep(NA_character_, length(x)))
  iso <- grepl("^\\d{4}-\\d{2}-\\d{2}", x)
  out[iso] <- as.Date(substr(x[iso], 1, 10))
  eight <- grepl("^\\d{8}$", x)
  out[eight] <- as.Date(x[eight], format = if (format_short == "ymd") "%Y%m%d" else "%d%m%Y")
  out
}
month_of  <- function(d) as.Date(format(d, "%Y-%m-01"))               # 1st day of the month
mun6    <- function(x) substr(gsub("\\D", "", as.character(x)), 1, 6) # IBGE code, 6 digits
log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# ---- C. DATASUS microdata cache --------------------------------------------
# Each DATASUS file (SINAN: national year; SIM: state-year; SIH: state-month) is
# downloaded ONCE, with a broad set of columns, and saved as Parquet in
# data/cache_datasus/<system>/. After that it is read from disk.
#
# FROZEN cache (default): a file that already exists is never replaced, even
# if it is from a recent year DATASUS still revises. Only with DATASUS_UPDATE=1 are
# the last two years re-downloaded after 120 days. Missing files are
# always downloaded.
CACHE_COLUMNS <- list(
  "SIM-DO" = c("CODMUNRES", "CODMUNOCOR", "DTOBITO", "DTNASC", "IDADE", "SEXO",
               "RACACOR", "ESC", "ESTCIV", "OCUP", "LOCOCOR", "CODESTAB", "CAUSABAS",
               "CIRCOBITO", "ACIDTRAB", "LINHAA", "LINHAB", "LINHAC", "LINHAD", "LINHAII"),
  "SIH-RD" = c("UF_ZI", "ANO_CMPT", "MES_CMPT", "N_AIH", "IDENT", "CNES", "MUNIC_RES",
               "MUNIC_MOV", "DT_INTER", "DT_SAIDA", "IDADE", "COD_IDADE", "SEXO",
               "RACA_COR", "DIAG_PRINC", "DIAG_SECUN", "CID_ASSO", "CID_MORTE",
               paste0("DIAGSEC", 1:9), "PROC_REA", "CAR_INT", "COMPLEX", "MORTE",
               "DIAS_PERM", "UTI_MES_TO", "VAL_TOT"),
  "SINAN" = c("ID_AGRAVO", "ID_MUNICIP", "ID_MN_RESI", "DT_NOTIFIC", "DT_SIN_PRI",
              "SEM_PRI", "NU_IDADE_N", "CS_SEXO", "CS_RACA", "CS_GESTANT", "CLASSI_FIN",
              "CRITERIO", "HOSPITALIZ", "DT_INTERNA", "EVOLUCAO", "DT_OBITO", "SOROTIPO")
)
SINAN_PREFIX <- c("SINAN-DENGUE" = "DENG", "SINAN-CHIKUNGUNYA" = "CHIK")
FTP_SINAN <- c(final  = "ftp://ftp.datasus.gov.br/dissemin/publicos/SINAN/DADOS/FINAIS/",
               prelim = "ftp://ftp.datasus.gov.br/dissemin/publicos/SINAN/DADOS/PRELIM/")

# path of a file in the cache (the SIH "_full" suffix is inherited from the
# original cache: all hospitalizations, no ICD filter)
cache_file <- function(system, year, uf = NULL, month = NULL) {
  name <- paste(c(system, uf, year, if (!is.null(month)) sprintf("%02d", month)), collapse = "_")
  if (system == "SIH-RD") name <- paste0(name, "_completo")
  file.path(DIR$cache, system, paste0(name, ".parquet"))
}

cache_valid <- function(f, year) {
  if (!file.exists(f)) return(FALSE)
  if (Sys.getenv("DATASUS_UPDATE") != "1") return(TRUE)        # frozen
  recent <- year >= as.integer(format(Sys.Date(), "%Y")) - 2
  !recent || as.numeric(difftime(Sys.time(), file.mtime(f), units = "days")) <= 120
}

ftp_lst <- function(url_dir, timeout = 120) {
  h <- curl::new_handle(dirlistonly = TRUE, connecttimeout = 60, timeout = timeout)
  r <- curl::curl_fetch_memory(url_dir, handle = h)
  trimws(strsplit(rawToChar(r$content), "\r?\n")[[1]])
}

# Reads a DBF in chunks, keeping only the requested columns. Needed for national
# SINAN (DENGBR24 has > 6 million rows and ~120 columns; foreign::read.dbf
# would load everything into memory).
read_dbf_columns <- function(dbf, cols = NULL, block = 200000L) {
  con <- file(dbf, "rb"); on.exit(close(con), add = TRUE)
  h <- readBin(con, "raw", 32L)
  nrec <- readBin(h[5:8], "integer", size = 4L, endian = "little")
  hlen <- readBin(h[9:10], "integer", size = 2L, signed = FALSE, endian = "little")
  rlen <- readBin(h[11:12], "integer", size = 2L, signed = FALSE, endian = "little")
  seek(con, 0L); hdr <- readBin(con, "raw", hlen)
  nms <- character(); tam <- integer(); post <- 33L
  while (post + 31L <= hlen && hdr[post] != as.raw(0x0D)) {
    b <- hdr[post:(post + 31L)]
    nm <- b[1:11]; nms <- c(nms, rawToChar(nm[nm != as.raw(0)]))
    tam <- c(tam, as.integer(b[17])); post <- post + 32L
  }
  ini <- cumsum(c(2L, tam))[seq_along(tam)]          # byte 1 = deletion flag
  sel <- if (is.null(cols)) nms else intersect(cols, nms)
  if (!length(sel)) stop("None of the requested columns exists in the DBF.")
  j <- match(sel, nms)
  SEP <- rlen + 1L; NL <- rlen + 2L
  idx <- unlist(lapply(seq_along(j), function(k)
    c(ini[j[k]]:(ini[j[k]] + tam[j[k]] - 1L), if (k < length(j)) SEP)))
  idx <- c(idx, NL)
  tmp <- tempfile(fileext = ".txt"); out <- file(tmp, "wb")
  on.exit(unlink(tmp), add = TRUE)
  seek(con, hlen); loaded <- 0
  while (loaded < nrec) {
    k <- min(block, nrec - loaded)
    r <- readBin(con, "raw", rlen * k)
    k <- length(r) %/% rlen; if (k == 0) break
    m <- matrix(r[seq_len(rlen * k)], nrow = rlen)
    live <- m[1, ] != as.raw(0x2A)
    m <- rbind(m, as.raw(0x7C), as.raw(0x0A))
    o <- m[idx, live, drop = FALSE]
    o[o == as.raw(0) | o == as.raw(0x0D)] <- as.raw(0x20)
    writeBin(as.vector(o), out)
    loaded <- loaded + k
  }
  close(out)
  fread(tmp, sep = "|", header = FALSE, colClasses = "character", quote = "",
        strip.white = TRUE, encoding = "Latin-1", na.strings = "", fill = TRUE,
        col.names = sel, showProgress = FALSE)
}

# national SINAN straight from FTP (FINAIS, then PRELIM), reading the DBF in chunks
download_sinan <- function(system, year, timeout) {
  default <- sprintf("^%sBR%02d.*\\.dbc$", SINAN_PREFIX[[system]], year %% 100)
  for (rep in names(FTP_SINAN)) {
    lst <- tryCatch(ftp_lst(FTP_SINAN[[rep]]), error = function(e) character())
    fpaths <- sort(grep(default, lst, value = TRUE, ignore.case = TRUE))
    if (!length(fpaths)) next
    dbc2dbf <- get(".dbc2dbf", envir = asNamespace("microdatasus"))
    parts <- lapply(fpaths, function(a) {
      dbc <- tempfile(fileext = ".dbc"); dbf <- tempfile(fileext = ".dbf")
      on.exit(unlink(c(dbc, dbf)), add = TRUE)
      ok <- FALSE
      for (t in 1:3) {
        # keepalive: on long downloads the FTP control connection drops without it
        ok <- tryCatch({ curl::curl_download(paste0(FTP_SINAN[[rep]], a), dbc, quiet = TRUE,
                           handle = curl::new_handle(timeout = timeout * 3, connecttimeout = 60,
                                                     tcp_keepalive = 1L, tcp_keepidle = 30L,
                                                     tcp_keepintvl = 15L)); TRUE },
                       error = function(e) FALSE)
        if (ok) break; Sys.sleep(30)
      }
      if (!ok) stop("Download failed: ", a)
      dbc2dbf(dbc, dbf)
      read_dbf_columns(dbf, CACHE_COLUMNS$SINAN)
    })
    return(rbindlist(parts, fill = TRUE))
  }
  NULL
}

# SIM and SIH through microdatasus, dropping columns that do not exist in that year
download_microdatasus <- function(system, year, uf, month, timeout) {
  cols <- CACHE_COLUMNS[[system]]
  for (attempt in 1:4) {
    args <- list(year_start = year, year_end = year, information_system = system,
                 vars = cols, timeout = timeout, quiet = TRUE)
    if (!is.null(uf))  args$uf <- uf
    if (!is.null(month)) args <- c(args, list(month_start = month, month_end = month))
    res <- tryCatch(do.call(microdatasus::fetch_datasus, args),
                    microdatasus_unknown_vars = function(e) e,
                    error = function(e) { warning(system, " ", uf, " ", year, " ", month, ": ",
                                                  conditionMessage(e)); NULL })
    if (!inherits(res, "condition")) return(res)
    missing <- intersect(cols, regmatches(conditionMessage(res),
                                         gregexpr("[A-Z][A-Z0-9_]+", conditionMessage(res)))[[1]])
    if (!length(missing)) return(NULL)
    cols <- setdiff(cols, missing)
  }
  NULL
}

# Returns a data.table with the requested columns (those missing in that
# year come as NA). Downloads only if the file is not in the cache.
#   system: "SIM-DO" (uf), "SIH-RD" (uf and month), "SINAN-DENGUE"/"SINAN-CHIKUNGUNYA" (national)
get_datasus <- function(system, year, uf = NULL, month = NULL, vars = NULL, timeout = 600) {
  if (startsWith(system, "SINAN")) uf <- NULL
  f <- cache_file(system, year, uf, month)
  if (!cache_valid(f, year)) {
    d <- if (startsWith(system, "SINAN")) download_sinan(system, year, timeout)
         else download_microdatasus(system, year, uf, month, timeout)
    if (is.null(d)) return(NULL)
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(f, ".", Sys.getpid(), ".tmp")
    arrow::write_parquet(as.data.table(d), tmp, compression = "zstd")
    if (!suppressWarnings(file.rename(tmp, f)) && file.exists(f)) unlink(tmp)
  }
  d <- as.data.table(arrow::read_parquet(f, col_select = if (is.null(vars)) NULL else
    tidyselect::any_of(vars)))
  for (v in setdiff(vars, names(d))) set(d, j = v, value = NA_character_)
  d
}

# Joins the parts (SIH months or SIM states) of a year. Failed parts
# are retried up to 2 times; if any is still missing it returns NULL,
# so the year is not saved incomplete (the good parts stay in the cache).
get_parts <- function(f, keys, label_txts, o_que, pause = 30) {
  parts <- lapply(keys, f)
  for (t in 1:2) {
    failed <- vapply(parts, function(p) is.null(p) || !nrow(p), logical(1))
    if (!any(failed) || all(failed)) break
    Sys.sleep(pause); parts[failed] <- lapply(keys[failed], f)
  }
  failed <- vapply(parts, function(p) is.null(p) || !nrow(p), logical(1))
  if (any(failed) && !all(failed))
    message(o_que, ": failed for ", paste(label_txts[failed], collapse = ", "),
            " — not saved now (try again later)")
  if (any(failed)) return(NULL)
  rbindlist(parts, fill = TRUE)
}
get_sih <- function(uf, year, months = 1:12, vars = NULL, timeout = 600)
  get_parts(function(m) get_datasus("SIH-RD", year, uf = uf, month = m, vars = vars, timeout = timeout),
               months, sprintf("month %02d", months), paste("SIH-RD", uf, year))
get_sim <- function(year, vars = NULL, timeout = 600)
  get_parts(function(u) get_datasus("SIM-DO", year, uf = u, vars = vars, timeout = timeout),
               UFS, UFS, paste("SIM-DO", year))

# ---- D. linear combination and per-year re-anchoring -------------------------
# estimate and standard error of w'b, with the full covariance matrix
lincomb <- function(w, b, V) {
  est <- sum(w * b); c(est = est, se = sqrt(as.numeric(t(w) %*% V %*% w)))
}

# Model with one coefficient per year (i(year, x, ref = first year)): each
# coefficient is re-anchored on the mean of the pre-vaccine years (2017-2023, including the
# reference year, whose coefficient is 0), and the dispersion of those years calibrates the
# inference for 2024 and 2025 (t test with 6 degrees of freedom, in the manner of
# Conley and Taber). `scale` = 0.10 returns the effect per 10 percentage points.
#   m: fixest model; var: name of the variable interacted with i(year, ...)
reanchor_by_year <- function(m, var, years = PARAM$year_start:PARAM$year_end,
                             last_without_vaccine = 2023L, scale = 1) {
  nms <- sprintf("year::%d:%s", years[-1], var)
  b <- coef(m)[nms]; V <- vcov(m)[nms, nms]; K <- length(b)
  n_pre <- sum(years <= last_without_vaccine)
  wpre <- rep(0, K); wpre[years[-1] <= last_without_vaccine] <- 1 / n_pre
  tab <- rbindlist(lapply(seq_along(years), function(i) {
    w <- rep(0, K); if (i > 1) w[i - 1] <- 1; w <- w - wpre; r <- lincomb(w, b, V)
    data.table(year = years[i], log_b = r[["est"]], se = r[["se"]],
               RR = exp(scale * r[["est"]]),
               CI95_lwr = exp(scale * (r[["est"]] - 1.96 * r[["se"]])),
               CI95_upr = exp(scale * (r[["est"]] + 1.96 * r[["se"]])))
  }))
  sd_pre <- sd(tab[year <= last_without_vaccine, log_b])
  tab[, t_calibrated := log_b / (sd_pre * sqrt(1 + 1 / n_pre))]
  tab[, p_calibrated_unilateral := fifelse(year > last_without_vaccine, pt(t_calibrated, df = n_pre - 1), NA_real_)]
  tab[, p_calibrated_twosided := fifelse(year > last_without_vaccine, 2 * pt(-abs(t_calibrated), df = n_pre - 1), NA_real_)]
  tab[, sd_log_placebos := sd_pre]
  tab[]
}

# ---- E. coefficient as rate ratio; health regions ----------------------------
# one row with RR, 95% CI and p of a model term (scale 0.10 = per 10 p.p.)
rr_row <- function(m, term, analysis, scale = 1) {
  ct <- coeftable(m); b <- ct[term, 1]; se <- ct[term, 2]
  data.table(analysis = analysis, RR = exp(scale * b),
             CI95_lwr = exp(scale * (b - 1.96 * se)), CI95_upr = exp(scale * (b + 1.96 * se)),
             p = ct[term, 4], se_log = se, n_obs = nobs(m))
}

# municipality -> health region (449 regions, brpop table)
health_region <- function() {
  rs <- as.data.table(brpop::mun_reg_saude_449)
  cm <- intersect(c("code_muni", "cod_mun"), names(rs))[1]
  cr <- intersect(c("codi_reg_saude", "code_health_region", "cod_reg_saude", "reg_saude"), names(rs))[1]
  if (is.na(cm) || is.na(cr)) stop("brpop::mun_reg_saude_449 lacks the expected columns")
  unique(rs[, .(code_muni = mun6(get(cm)), health_reg = as.character(get(cr)))])
}
