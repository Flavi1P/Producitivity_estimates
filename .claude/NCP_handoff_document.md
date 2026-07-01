# Handoff Document: NCP Derivation from BGC-Argo Floats (North Atlantic)

## Purpose

This document describes how Net Community Production (NCP) is estimated from BGC-Argo float data in the North Atlantic, using nitrate profiles predicted by the CANYON-B neural network. It is intended to provide a scientific writing agent with the methodological depth needed to draft or review manuscript sections describing this approach.

---

## Overview

The toolbox derives depth-integrated NCP from changes in CANYON-B-predicted nitrate profiles between successive time bins. The core principle is the nitrate budget: a net decrease in vertically integrated nitrate within the mixed layer between two time periods reflects biological nitrogen uptake, which is converted to carbon units via the Redfield ratio. The pipeline is orchestrated by Snakemake and is described below step by step.

---

## Study Region: North Atlantic

The analysis is spatially restricted to named basins defined by user-supplied polygons in a YAML configuration file. For the North Atlantic paper, the regions used include:

- **Iceland Basin**: roughly bounded by the Reykjanes Ridge to the west (~32°W) and extending to ~12°W, between ~55–63°N. One tested polygon spans from the Reykjanes Ridge east to approximately (−12°E, 63°N) and south to (−23°W, 55.7°N).
- **Irminger Sea**: west of the Reykjanes Ridge, bounded by Greenland to the northwest (~26–42°W, 57–65°N).
- **Labrador Sea**: roughly 50–60°W, 53–62°N.
- **Norwegian Sea**: roughly 3–8.5°E, 62–70°N.

The study period nominally spans 2020–2026, filtered from the float database.

---

## Pipeline Steps

### Step 1 — Float Discovery and Download (`download_sprof.py`)

BGC-Argo Sprof NetCDF files (one per float, containing all profiles) are retrieved from the GDAC via the `argo_gdac` utility. The bounding box for download is taken as the union of all basin polygon vertices defined in the config. Only floats carrying a dissolved oxygen (DOXY) sensor are selected. Files are saved to a shared directory to avoid re-downloading across runs.

### Step 2 — Profile Processing (`process_sprof.R`)

For each Sprof file:
- Variables extracted: longitude, latitude, date (from JULD in days since 1950-01-01), pressure, temperature (TEMP), salinity (PSAL), and adjusted dissolved oxygen (DOXY_ADJUSTED).
- DOXY quality control: only QC flags 1 and 2 (good and probably good data) are retained; all other values are set to NA.
- Profiles are interpolated to a regular 1-m depth grid from 0–2000 m using linear interpolation (`zoo::na.approx`). Depths with fewer than 10 valid measurements are left as NA rather than extrapolated.
- Output: per-float CSV with columns `lon`, `lat`, `date`, `depth` (1 m steps), `temp`, `sal`, `oxygen`, `prof_number` (profile rank by date).

### Step 3 — Nitrate Prediction with CANYON-B (`predict_nitrate.R`)

Nitrate profiles are not measured directly by the floats used here. Instead, they are estimated using **CANYON-B** (Carbonate system and Nutrients from hydrographic properties, Bayesian Neural Network; Bittig et al. 2018), a committee of neural networks trained on climatological biogeochemical data. The implementation used (`scripts/utils/fastr_canyon.R`) is a vectorized, in-house R port that pre-loads the published network weights (`scripts/utils/canyon_weights/wgts_NO3.txt` and companion files) and processes all samples in a single batched matrix operation, avoiding profile-by-profile loops.

**Inputs to CANYON-B per depth level:**
- Date (decimal year)
- Latitude, longitude
- Pressure (db)
- In-situ temperature (°C)
- Practical salinity (PSU)
- Dissolved oxygen (µmol kg⁻¹) — DOXY_ADJUSTED

**Outputs:**
- `canyon_nitrate`: predicted nitrate concentration (mmol m⁻³)
- `canyon_nitrate_ci`: total 1-sigma uncertainty (mmol m⁻³), combining measurement uncertainty, neural network committee uncertainty, and parameter uncertainty

Profiles with fewer than 10 valid oxygen measurements or lacking valid geographic coordinates are skipped. Processing is parallelised across floats using `furrr::future_map` (multisession workers).

### Step 4 — Merge and MLD (`merge_nitrate_files.R`)

All per-float nitrate CSVs are concatenated into a single merged file. Additional derived variables are computed:
- **Absolute Salinity (SA)** and **Conservative Temperature (CT)** via TEOS-10 equations (R package `gsw`)
- **Potential density anomaly (σ₀)** from SA and CT
- **Mixed Layer Depth (MLD)** per profile: computed using a threshold criterion Δσ₀ = 0.03 kg m⁻³ referenced to the near-surface layer (0–5 m average), implemented via `castr::mld()`. The MLD is the shallowest depth at which σ₀ exceeds the reference value by 0.03 kg m⁻³.

The merged CSV contains columns: `float_wmo`, `prof_number`, `date`, `lon`, `lat`, `canyon_nitrate`, `canyon_nitrate_ci`, `sa`, `ct`, `sigma0`, `MLD`, `sal`, `oxygen`, `depth`.

### Step 5 — NCP Computation (`compute_ncp.R`)

This is the core computation. The merged CSV is read and filtered first by bounding box then by exact polygon membership (using `sf::st_within`). The following steps are then applied for each configured time-step window (e.g. 10, 15, 20, 30 days):

#### 5a. MLD outlier removal
Profile-level MLD values are screened with a 1.5 × IQR rule. Profiles with MLD below `Q1 – 1.5×IQR` or above `Q3 + 1.5×IQR` are flagged and excluded from downstream computation. This prevents a single anomalously deep or shallow mixed layer (e.g. from a data artefact) from distorting the basin-mean MLD.

#### 5b. Time binning
Both MLD/Zeu and nitrate profiles are assigned to fixed-duration time bins starting from the earliest clean profile date. Within each bin, all available profiles from all floats within the basin polygon are averaged:
- For MLD/Zeu: arithmetic mean (and, when `mld_method = "winter_max"`, also bin maximum)
- For nitrate: depth-wise arithmetic mean across all profiles and floats within the bin

#### 5c. MLD method
Two options are available via `mld_method` in the config:
- **`"average"`**: bin-mean MLD is always used.
- **`"winter_max"`**: when the bin-mean MLD exceeds 100 m (a proxy for deep winter convection), the bin-maximum MLD is used instead of the mean. This better captures the deepest convective mixing that sets the baseline nutrient inventory before the productive season.

#### 5d. NCP integration depth
For each pair of consecutive time bins *t* and *t+1*, the NCP integration depth is:

```
NCP_integration_depth(t) = max(MLD(t), MLD(t−1), Zeu(t), Zeu(t−1))
```

where Zeu is the euphotic zone depth (set to a constant `zeu_default = 40 m` across the North Atlantic in the current configuration). Taking the maximum across the current and previous MLD and Zeu ensures that the integration layer spans the full depth range that was ever within the mixed layer or lit zone during the interval bounded by the two profiles. This avoids underestimating NCP when the mixed layer shallows rapidly.

A 3-point centered rolling mean (`zoo::rollapply`, width = 3) is applied to the integration depth time series to reduce noise before integration.

#### 5e. Nitrate integration
For each time bin, nitrate is linearly interpolated on a 1-m grid from 0 to the integration depth using `approxfun`, then integrated as:

```
int_N = mean(nitrate(0 → depth)) × depth   [mmol m⁻²]
```

This is equivalent to trapezoidal integration on the 1-m grid. A minimum of 2 valid nitrate points within the layer is required; otherwise NA is returned.

A 3-point centered rolling mean is applied to the integrated nitrate time series before differencing.

#### 5f. NCP formula
NCP is computed as the decrease in vertically integrated nitrate between consecutive bins, converted to carbon using the Redfield C:N ratio and normalized by the elapsed time:

```
NCP = (int_N(t−1) − int_N(t)) × 6.625 / dt   [mmol C m⁻² d⁻¹]
```

where `6.625` is the Redfield C:N molar ratio (106 C : 16 N) and `dt` is the number of days between bin mid-points. A positive NCP indicates net autotrophy (nitrate being consumed by phytoplankton); a negative value indicates net heterotrophy or remineralization.

A soft outlier correction is applied: if NCP < −60 mmol C m⁻² d⁻¹, the value is halved. This damps extreme apparent remineralization events that likely reflect sampling artefacts rather than true biological signal.

#### 5g. Sensitivity test
The computation is repeated for all time steps listed in `ncp_time_steps` (typically 10, 15, 20, 30 days). Results are overlaid in a sensitivity plot to assess robustness to the choice of temporal binning. Shorter time steps resolve intra-seasonal variability but are noisier due to sparse sampling; longer time steps smooth the signal but may average over bloom events.

### Step 6 — Uncertainty Estimation (`compute_ncp_uncertainty.R`)

Uncertainty is propagated via a Monte Carlo approach with `n_mc` iterations (typically 100). Two sources of uncertainty are jointly sampled in each iteration:

1. **Profile bootstrap**: within each time bin, profiles are resampled with replacement (block bootstrap preserving the temporal bin structure). This propagates spatial and inter-float variability in nitrate concentrations.

2. **CANYON-B neural network uncertainty**: Gaussian noise with standard deviation equal to `canyon_nitrate_ci` (the per-level CI returned by CANYON-B, defaulting to 1.2 mmol m⁻³ when missing) is added independently to each depth level of each resampled profile.

After resampling and noise injection, the full NCP computation (binning, integration, differencing) is repeated identically to Step 5. The distribution of NCP estimates across all iterations is summarized as:
- `NCP_mean`, `NCP_sd`: mean and standard deviation
- `NCP_q05`, `NCP_q95`: 5th and 95th percentiles (90% credible interval)

### Step 7 — Annual and Seasonal Summaries (`summary_plots.R`)

From the Monte Carlo output, the following annual metrics are derived:

- **ANCP** (Annual NCP): `ANCP = Σ (NCP_mean × dt) / 1000` over a calendar year (mol C m⁻² yr⁻¹). Negative NCP periods contribute negatively.
- **sANCP** (Seasonal/productive NCP): same as ANCP but with NCP clipped to zero before summing (`pmax(NCP_mean, 0)`). Captures only the productive (autotrophic) signal, ignoring apparent remineralization periods.
- **ANCP from climatology**: monthly mean NCP × 30.5 d, summed over 12 months. Provides an estimate of average annual export independent of year-to-year variability in sampling coverage.
- **ANCP partition**: annual positive and negative NCP components are plotted separately to show the balance between production and apparent respiration/remineralization periods.
- **Monthly climatology**: mean NCP per calendar month, averaged across all years.

Uncertainty on ANCP is propagated assuming independence between time steps: `σ_ANCP = sqrt(Σ (σ_NCP × dt)²) / 1000`.

---

## Single-Float Mode (`compute_ncp_float.R`)

For individual float tracking, the basin-mode spatial pooling is bypassed. Two sub-modes exist:

- **Time-windowed mode** (`float_time_window` set, e.g. "15 days"): profiles within each window are bin-averaged before NCP is computed, exactly as in basin mode but for a single float.
- **Per-profile mode** (`float_time_window = null`): NCP is computed from consecutive pairs of raw profiles, without any binning. The integration depth for each pair is `max(MLD_i, MLD_{i-1}, Zeu_i, Zeu_{i-1})`.

---

## Key Parameters (North Atlantic Configuration)

| Parameter | Value | Notes |
|---|---|---|
| `zeu_default` | 40 m | Constant euphotic zone depth; no satellite Zeu used |
| `ncp_time_steps` | 10, 15, 20, 30 days | Sensitivity range |
| `mld_method` | `"winter_max"` | Uses bin-max MLD when mean MLD > 100 m |
| `n_mc` | 100 | Monte Carlo iterations for uncertainty |
| `date_start` / `date_end` | 2020-01-01 / 2026-01-01 | Data window |
| MLD threshold | Δσ₀ = 0.03 kg m⁻³ | Referenced to 0–5 m |
| Redfield ratio | 6.625 | C:N = 106:16 |
| NCP outlier cap | −60 mmol C m⁻² d⁻¹ | Values below this are halved |
| MLD IQR fence | 1.5 × IQR | Per-basin outlier rejection |

---

## Software and Library Versions

### Pipeline Orchestration
- **Snakemake** — workflow management; version used at time of analysis (conda-based installation)

### Python (data download and NPP)
- **Python** 3.10.14 (conda-forge, MSC v.1938 64-bit Windows)
- **pandas** 2.3.3 — data handling for manifests and download bookkeeping
- **argopy** — BGC-Argo GDAC interface (`argo_gdac` utility)
- **netCDF4** — reading Sprof NetCDF files (Python scripts; R uses `ncdf4`)
- **numpy** — numerical operations in Python scripts

### R (core NCP computation)
R version used is from the Positron IDE installation. All packages from CRAN unless noted.

| Package | Role |
|---|---|
| `tidyverse` (dplyr, ggplot2, purrr, readr, tidyr, lubridate, stringr) | Data manipulation, visualization, I/O |
| `zoo` | `na.approx` (linear interpolation along time), `rollapply` (rolling mean) |
| `sf` | Spatial filtering — polygon intersection, `st_within`, WGS84 CRS |
| `data.table` / `fread` | Fast streaming read of the large (~15 GB) merged CSV with column selection |
| `gsw` | TEOS-10 thermodynamic equations: SA, CT, σ₀ |
| `castr` | MLD computation (`castr::mld()`) |
| `seacarb` | Required internally by the CANYON-B implementation for pCO₂ calculations (not used for NCP directly) |
| `mgcv` | Required internally by CANYON-B (`mgcv::in.out` for Arctic polar shift correction) |
| `ncdf4` | Reading Argo Sprof NetCDF files in R |
| `furrr` | Parallel execution of CANYON-B across floats (multisession `future_map`) |
| `here` | Path management |

### CANYON-B
- **CANYON-B** neural network (Bittig et al. 2018, *Frontiers in Marine Science*): in-house vectorized R implementation (`scripts/utils/fastr_canyon.R`) using pre-trained weight files distributed with the original MATLAB/Python code. The implementation reproduces the published network architecture (committee of multi-layer perceptrons, 1–2 hidden layers, tanh activation) and full uncertainty decomposition (committee uncertainty, noise variance, parameter uncertainty, input propagation uncertainty, measurement uncertainty).

---

## Key Methodological Choices for the Paper

1. **No direct nitrate measurements**: Nitrate is entirely CANYON-B-predicted from temperature, salinity, and oxygen. This is both a strength (uniform coverage, no sensor calibration drift) and a limitation (CANYON-B uncertainty ≈ 1–2 mmol m⁻³ per level; systematic biases possible in undersampled regions).

2. **Integration depth = max(MLD, Zeu) across consecutive profiles**: Ensures the budget is closed when the mixed layer deepens (entraining nutrients) or when Zeu exceeds the MLD. Using only the instantaneous MLD would underestimate the depth over which production integrates, especially in the North Atlantic where MLD varies from <20 m in summer to >500 m in winter.

3. **Winter_max MLD method**: When bin-mean MLD > 100 m, the bin-maximum is used. This captures the true depth of winter convective mixing that pre-conditions the nutrient inventory for the spring bloom, rather than an average over a period when the mixed layer may be evolving rapidly.

4. **Basin pooling of multiple floats**: All BGC-Argo floats within each named basin polygon are pooled at each time step. This increases spatial representativeness and reduces sampling noise, at the cost of spatial detail within the basin.

5. **Constant Zeu = 40 m**: A single constant is used for the North Atlantic (as opposed to satellite-derived Zeu varying in space and time). This simplification is conservative and appropriate when MLD regularly exceeds Zeu (which is common in this region), making the integration depth MLD-dominated throughout most of the year.

6. **NCP expressed in carbon units**: The nitrate budget (mmol N m⁻² d⁻¹) is converted to carbon via the Redfield ratio (C:N = 6.625), yielding mmol C m⁻² d⁻¹. Annual integrals are expressed in mol C m⁻² yr⁻¹.

7. **sANCP vs ANCP**: sANCP (NCP clipped at 0) is reported alongside ANCP to distinguish between gross production and apparent respiration. Negative NCP periods in winter could reflect true heterotrophy, but may also arise from noise in the nitrate budget at a time when seasonal signals are small relative to CANYON-B uncertainty.

---

## References to Cite

- **CANYON-B**: Bittig, H.C., Steinhoff, T., Claustre, H., Fiedler, B., Williams, N.L., Sauzède, R., Körtzinger, A., Gattuso, J.-P. (2018). An alternative to static climatologies: robust estimation of open ocean CO₂ variables and nutrient concentrations from T, S, and O₂ data using Bayesian neural networks. *Frontiers in Marine Science*, 5, 328. https://doi.org/10.3389/fmars.2018.00328
- **TEOS-10 / gsw**: IOC, SCOR and IAPSO (2010). The international thermodynamic equation of seawater – 2010: Calculation and use of thermodynamic properties. *Intergovernmental Oceanographic Commission, Manuals and Guides* No. 56, UNESCO.
- **MLD criterion**: de Boyer Montégut, C., Madec, G., Fischer, A.S., Lazar, A., Iudicone, D. (2004). Mixed layer depth over the global ocean: An examination of profile data and a profile-based climatology. *Journal of Geophysical Research: Oceans*, 109, C12003.
- **BGC-Argo data**: Biogeochemical-Argo Planning Group (2016). The scientific rationale, design, and implementation plan for a Biogeochemical-Argo float array. *Biogeochemical-Argo Planning Group*. https://doi.org/10.13155/46601
