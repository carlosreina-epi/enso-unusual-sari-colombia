# ENSO and unusual severe acute respiratory infection in Colombia, 2011–2019

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23924392.svg)](https://doi.org/10.5281/zenodo.23924392)

[![Preprint](https://img.shields.io/badge/Preprint-Research%20Square-blue)](https://doi.org/10.21203/rs.3.rs-11612511/v1)

Code and aggregated data for the manuscript:

> Reina-Bolaños CA, Guzmán S, Blanco C. El Niño–Southern Oscillation and unusual severe acute respiratory infection in Colombia, 2011–2019: a national time-series analysis using distributed lag non-linear models. *Research Square* [preprint]; 2026. https://doi.org/10.21203/rs.3.rs-11612511/v1 (under review at the *International Journal of Biometeorology*).

## What this repository does

Estimates the association between ENSO intensity (Oceanic Niño Index, ONI) and monthly cases of unusual severe acute respiratory infection (SARI, SIVIGILA event 348) in Colombia, using quasi-Poisson distributed lag non-linear models (lags 0–12 months) with control for seasonality and long-term trend and a population offset. Models are fitted nationally, by region and by departmental Multidimensional Poverty Index tertile, with prespecified sensitivity analyses.

## Data

| File | Content | Source |
|---|---|---|
| `irag348_semanal_estratos.csv` | Weekly event-348 case counts 2011–2019, national and by stratum | SIVIGILA, Instituto Nacional de Salud, via [datos.gov.co](https://www.datos.gov.co/Salud-y-Protecci-n-Social/Datos-de-Vigilancia-en-Salud-P-blica-de-Colombia/4hyg-wa9d) (dataset 4hyg-wa9d) |
| `raw/oni.csv` | Monthly Oceanic Niño Index, 2009–2020 | [NOAA Climate Prediction Center](https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt) (ERSSTv5) |
| `raw/pop_dpto_ano.csv` | Population by department and year | DANE projections, base CNPV 2018 |
| `raw/dept_strata.csv` | Region and Multidimensional Poverty Index tertile by department | IDEAM 2017; DANE 2023 |

All data are public and aggregated; no individual records are included. Annual totals of the case series were checked against department-year totals from the same source. Strata `Orinoquia_Amazonia` and `IPM_T3` were obtained by subtraction and include cases without an assigned department.

## How to run

Requirements: R ≥ 4.1 and the packages `ggplot2` and `dlnm`.

```r
install.packages(c("ggplot2", "dlnm"))
```

From this folder:

```bash
Rscript analisis_enso_irag_dlnm.R   # full analysis: tables and figures in resultados/
Rscript verificar_dlnm.R            # re-fits the main model with dlnm::crossbasis and compares
```

`analisis_enso_irag_dlnm.R` builds the cross-basis with base-R natural splines (equivalent to `dlnm::crossbasis`); `verificar_dlnm.R` reproduces the main model with the `dlnm` package (v2.4.10), giving the same cumulative relative risk (0.95, 95% CI 0.78–1.15 per 1 °C ONI).

## Outputs (`resultados/`)

- `tabla1_descriptivos.csv`: cases, monthly means and incidence by stratum
- `tabla2_rr_acumulados.csv`: cumulative relative risks for the main model, strata and sensitivity analyses
- `fig1`–`fig5`, `figS1_dag.png`: manuscript figures
- `acf_residuos_principal.csv`: residual autocorrelation of the main model
- `verificacion_dlnm_rezago.png`: lag-response curve from the `dlnm` package

## Citation
[![Preprint](https://img.shields.io/badge/Preprint-Research%20Square-blue)](https://doi.org/10.21203/rs.3.rs-11612511/v1)

**Preprint:** Reina-Bolaños CA, Guzmán S, Blanco C. El Niño–Southern Oscillation and unusual severe acute respiratory infection in Colombia, 2011–2019: a national time-series analysis using distributed lag non-linear models. *Research Square* [preprint]; 2026. https://doi.org/10.21203/rs.3.rs-11612511/v1 (under review at the *International Journal of Biometeorology*).

Archived version: Reina-Bolaños CA, Guzmán S, Blanco C. ENSO and unusual severe acute respiratory infection in Colombia, 2011-2019: code and aggregated data (v1.0.0). Zenodo; 2026. doi:10.5281/zenodo.23924392

See also `CITATION.cff`. Please cite the article once published.

## License

Code: MIT License (see `LICENSE`). Data remain subject to the terms of their original public sources.
