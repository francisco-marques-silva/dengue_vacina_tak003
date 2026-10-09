# Dengue vaccination (TAK-003) at ages 10–14 and dengue hospitalizations — Brazil, 2017–2025

Code and data to reproduce the analyses in the article, from DATASUS microdata to the tables, figures and
numbers in the text.

## How to reproduce

Requires R ≥ 4.2 (results were produced with R 4.6.1).

```r
renv::restore()            # same package versions (renv.lock)
source("R/run_all.R")      # runs every script in order; logs in logs/
```

The repository already includes the aggregated data, so nothing is downloaded from DATASUS and the results
are identical. The last script (`99_check_reproduction.R`) compares the generated numbers with
`docs/reference/manuscript_numbers_reference.json` and flags any difference.

To rebuild everything from the microdata (~2 GB download, several hours), delete `data/` (except
`data/raw/lists/`) and `outputs/` and run with `Sys.setenv(CACHE_SOURCE = "none")`. Because DATASUS revises
its preliminary files, results may differ slightly; script 99 shows where.

## Design

Triple difference (municipality × age group × month), Poisson with fixed effects
`municipality×month + municipality×age group + age group×month`, log(population) offset and standard errors
clustered by municipality. Treated: ages 10–14 in the 521 first-phase municipalities, from February 2024.
Comparison: 2,819 municipalities never included in vaccination. Dengue hospitalizations: ICD-10 A90, A91 and
A97 as principal diagnosis. Parameters in `R/00_setup.R`.

## Scripts

| Script | Purpose |
|---|---|
| `00_setup.R`, `01_functions.R` | packages, parameters and functions |
| `10`–`16` | data: SINAN, SIH, SIM, population, SI-PNI and municipality list |
| `20_panel.R` | municipality × age group × month and birth-cohort panels |
| `30_main_model.R` | main model and event study |
| `31_validity.R` | calibrated inference, coverage bound and age composition |
| `32_spatial_permutation.R` | permutation of municipalities and of health regions |
| `33_dose_response.R` | dose-response by vaccine coverage |
| `34_birth_cohort.R` | eligible vs. neighbouring birth cohorts |
| `35_sensitivity.R` | sensitivity analyses and negative-control outcomes |
| `40`–`44` | tables, figures, maps and numbers in the text |
| `99_check_reproduction.R` | check against the reference |

## Included data

| Folder | Contents |
|---|---|
| `data/raw/lists/` | official list of municipalities included in vaccination |
| `data/raw/sipni/` | vaccine doses by municipality and month (SI-PNI) |
| `data/interim/` | SINAN, SIH and SIM aggregates by state and year |
| `data/processed/` | SIH aggregates by birth cohort and negative-control outcomes |
| `outputs/models/` | permutation draws |
| `docs/reference/` | reference numbers |

Tables and figures are written to `outputs/tables/` and `outputs/figures/`.
