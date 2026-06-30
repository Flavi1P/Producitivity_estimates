# Hands-off: add an air-sea O2 flux correction to the 3902681 oxygen budget

This is a **self-contained implementation spec**. You do not need prior
conversation context. It tells you exactly what to add to
`Scripts/o2_budget_3902681.R`, the equations and constants to use, the units and
sign conventions, and how to verify. All paths are relative to the repo root
`C:\Users\petit\Documents\Producitivity_estimates`. **Run everything from the
repo root.** R is **4.5.2 only**:

```powershell
& "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/o2_budget_3902681.R
```

Read `.claude/o2_budget_3902681_explained.md` first — it explains the existing
script in full. This document only adds the **air-sea diffusive flux** term on
top of it. Do **not** modify `Scripts/ncp_o2budget_float_3902681.R` (a separate
older script).

---

## 0. Goal in one paragraph

The existing script computes the rate of change of the float's 0–40 m O2
inventory (`rate_mmol_o2_m2_d`, signed dO2/dt) for three pairing variants
(`prev`, `night`, `net`). That rate currently mixes biology and gas exchange.
Add the **diffusive air-sea O2 flux** `F_as` (Wanninkhof 2014 piston velocity ×
surface O2 disequilibrium, driven by ERA5 hourly wind + sea-level pressure) and
report a **flux-corrected** rate `rate_corr = rate − F_as`. Keep the raw `rate`
column unchanged. Diffusive flux only — **no bubble injection** in this pass.

---

## 1. New inputs

### 1a. ERA5 NetCDF — `Data/Raw/ERA5/era5_wind_slp_3902681.nc`
Download with a **new** script `Scripts/download_era5_wind.py` (template in §5).
Hourly ERA5 single-levels, variables `10m_u_component_of_wind` (`u10`),
`10m_v_component_of_wind` (`v10`), `mean_sea_level_pressure` (`msl`, Pa). Area
`[65, -40, 58, -10]` (N, W, S, E), dates **2024-10-01 → 2026-02-11** (covers the
float record 2024-10-02 → 2026-02-10). `Data/Raw/` is **gitignored** — never
commit the NetCDF.

- **User prerequisite:** a `~/.cdsapirc` with Copernicus CDS credentials. The
  download is a manual step (like the PAR manifest in `download_par.py`); do not
  automate credential setup.
- Python env: conda `cmts_learn_olci` (the repo's download env, see
  `clearsky_daily_par.py`): `C:/Users/petit/anaconda3/envs/cmts_learn_olci/python.exe`.
- ERA5 `time` units are **"hours since 1900-01-01"** → POSIXct:
  `as.POSIXct(time_hours*3600, origin="1900-01-01", tz="UTC")`.

### 1b. Everything else is unchanged
Float NetCDF `Data/Raw/Floats/3902681_Sprof.nc`, reference CSVs, all as in the
existing script. No new R dependency — `ncdf4`, `tidyverse`, `lubridate`, `gsw`,
`zoo` are already loaded. **`gsw` 1.2.0 on this machine has `gsw_O2sol`**
(Garcia & Gordon 1992 solubility) — use it; do not hand-roll a solubility
polynomial.

---

## 2. Equations, constants, units, signs (read carefully)

All concentrations in **mmol m⁻³** (same as the existing inventory). Flux and
rates in **mmol O2 m⁻² d⁻¹**. Sign convention: **`F_as` positive = O2 INTO the
ocean.**

**Surface O2 saturation (per profile).** Using surface `SA`, `CT` (TEOS-10) and
surface density `rho_surf`:
```
o2_eq_umolkg = gsw_O2sol(SA, CT, p = 0, longitude = lon, latitude = lat)  # umol/kg
o2_eq        = o2_eq_umolkg * rho_surf / 1000                              # mmol m-3, at 1 atm
```
SLP correction to actual pressure (first-order; water-vapour term omitted — note
in a comment):
```
o2_eq_p = o2_eq * (msl / 101325)      # msl in Pa, 101325 Pa = 1 atm
```

**Schmidt number for O2 (Wanninkhof 2014, Table 1), T in °C:**
```
Sc = 1920.4 - 135.6*T + 5.2122*T^2 - 0.10939*T^3 + 0.00093777*T^4
```

**Piston velocity (Wanninkhof 2014), cm h⁻¹ → m d⁻¹** (factor 0.24 = 24/100):
```
k = 0.251 * <U2> * (Sc/660)^(-0.5) * 0.24      # m d-1
```
where `<U2>` is the **mean of (u10² + v10²)** over the budget interval (mean of
squared wind, not square of mean — this preserves the gustiness that drives
quadratic gas exchange).

**Diffusive flux and corrected rate:**
```
dO2  = o2_surf_mean - o2_eq_p_mean       # mmol m-3 ; >0 = supersaturated
F_as = -k * dO2                          # mmol O2 m-2 d-1 ; >0 = into ocean
rate_corr = rate - F_as                  # biological + physical residual
```
So supersaturation (`dO2>0`) gives `F_as<0` (outgassing); undersaturation gives
`F_as>0` (ingassing). Both `o2_surf_mean` and `o2_eq_p_mean` are the **average
of the pair's two profiles** (start and end), to represent the interval.

---

## 3. Code changes to `Scripts/o2_budget_3902681.R`

### 3a. Read ERA5 once, after the float `nc_close` (near line 46)
```r
era5_path <- "Data/Raw/ERA5/era5_wind_slp_3902681.nc"
e_nc  <- nc_open(era5_path)
e_lon <- ncvar_get(e_nc, "longitude")
e_lat <- ncvar_get(e_nc, "latitude")
e_t   <- as.POSIXct(ncvar_get(e_nc, "time") * 3600, origin = "1900-01-01", tz = "UTC")
e_u10 <- ncvar_get(e_nc, "u10")      # [lon x lat x time]
e_v10 <- ncvar_get(e_nc, "v10")
e_msl <- ncvar_get(e_nc, "msl")      # Pa
nc_close(e_nc)
# check dimension order with dim(e_u10) vs length(e_lon)/length(e_lat)/length(e_t)
```
> Confirm the array dim order (ncdf4 returns `[lon, lat, time]` for ERA5 here,
> but verify with `dim()` and `e_nc$var$u10$dim` and index accordingly.)

`era5_at()` helper — nearest grid cell to the midpoint position, interval means:
```r
era5_at <- function(t0, t1, lon0, lat0) {
  i  <- which.min(abs(e_lon - lon0))
  j  <- which.min(abs(e_lat - lat0))
  kt <- which(e_t >= t0 & e_t <= t1)
  if (length(kt) == 0) kt <- which.min(abs(as.numeric(e_t - (t0 + (t1 - t0)/2))))
  u <- e_u10[i, j, kt]; v <- e_v10[i, j, kt]; m <- e_msl[i, j, kt]
  list(u2 = mean(u^2 + v^2, na.rm = TRUE), msl = mean(m, na.rm = TRUE))
}
```

### 3b. Capture surface state inside the existing inventory loop
The loop (existing lines ~65–98) already computes `SA`, `CT`, `rho`, `o2_vol`,
`depth`, and `o2_grid` on `grid <- 0:40`. Add per-profile vectors before the
loop:
```r
o2_surf_vec  <- rep(NA_real_, n_prof)   # mean O2 0-10 m, mmol m-3
o2_eq_vec    <- rep(NA_real_, n_prof)   # equilibrium O2 0-10 m, mmol m-3 (1 atm)
```
Inside the loop, right after `o2_grid` is built and validated (after existing
line ~95):
```r
top <- depth <= 10
SA_s <- mean(SA[top]); CT_s <- mean(CT[top]); rho_s <- mean(rho[top])
o2_surf_vec[p] <- mean(o2_grid[grid <= 10])                       # mmol m-3
o2_eq_umolkg   <- gsw_O2sol(SA_s, CT_s, 0, lon[p], lat[p])        # umol/kg
o2_eq_vec[p]   <- o2_eq_umolkg * rho_s / 1000                     # mmol m-3 (1 atm)
```
> `SA`, `CT`, `rho` are currently computed on the **filtered+reordered** levels
> (existing lines 76–88 filter by `keep` then `order(depth)`). Make sure `top`
> uses the same reordered `depth`. Simplest: compute `SA_s`/`CT_s`/`rho_s` from
> the reordered vectors (apply the same `ord` used for `depth`/`o2_vol`). Double-
> check the reorder so surface T/S aren't taken from a deep level.

Add the two new columns to the `prof` tibble (existing lines ~100–106) so they
survive `arrange(time)` and the head/tail `slice` trim:
```r
o2_surf = o2_surf_vec,
o2_eq   = o2_eq_vec
```

### 3c. Compute flux in `make_pair()` (existing lines ~129–143)
Extend `make_pair(d, i, j, type)` to also return flux columns. Keep the existing
`rate_mmol_o2_m2_d` exactly as is.
```r
make_pair <- function(d, i, j, type) {
  dt_days <- as.numeric(difftime(d$time[j], d$time[i], units = "days"))
  rate    <- (d$inventory[j] - d$inventory[i]) / dt_days

  T_s     <- mean(c(d$? ... ))   # need surface temperature; see note below
  Sc      <- 1920.4 - 135.6*T_s + 5.2122*T_s^2 - 0.10939*T_s^3 + 0.00093777*T_s^4
  lon_mid <- mean(c(d$lon[i], d$lon[j]))
  lat_mid <- mean(c(d$lat[i], d$lat[j]))
  w       <- era5_at(d$time[i], d$time[j], lon_mid, lat_mid)
  k       <- 0.251 * w$u2 * (Sc/660)^(-0.5) * 0.24             # m d-1

  o2_eq_p_i <- d$o2_eq[i] * (w$msl/101325)
  o2_eq_p_j <- d$o2_eq[j] * (w$msl/101325)
  dO2  <- mean(c(d$o2_surf[i], d$o2_surf[j])) - mean(c(o2_eq_p_i, o2_eq_p_j))
  Fas  <- -k * dO2                                             # mmol O2 m-2 d-1, +into ocean

  tibble(
    type = type,
    mtime = d$time[i] + (d$time[j] - d$time[i]) / 2,
    t_start = d$time[i], t_end = d$time[j], dt_days = dt_days,
    phase_start = d$phase[i], phase_end = d$phase[j],
    inv_start_mmol_m2 = d$inventory[i], inv_end_mmol_m2 = d$inventory[j],
    rate_mmol_o2_m2_d = rate,
    k_m_d = k, u2_m2_s2 = w$u2, o2_eq_mmol_m3 = mean(c(o2_eq_p_i, o2_eq_p_j)),
    delta_o2_mmol_m3 = dO2, Fas_mmol_o2_m2_d = Fas,
    rate_corr_mmol_o2_m2_d = rate - Fas
  )
}
```
> **Surface temperature for `Sc`:** the loop currently does not store surface
> temperature into `prof`. Add a third per-profile vector `sst_vec[p] <-
> mean(tt[top])` (using in-situ `tt` on the reordered levels) and a `sst` column
> in `prof`, then `T_s <- mean(c(d$sst[i], d$sst[j]))`. Salinity for Schmidt
> isn't needed (W2014 Sc is for 35 psu seawater); using float SST is sufficient.

`make_pair` is called from three places (existing lines ~148, ~159) — they all
flow through this one function, so all three budget tables get the flux columns
automatically.

### 3d. CSV (existing Step 3, lines ~168–173)
Add the new columns to the `select()` so they're written:
```r
budget_out <- budget |>
  select(type, mtime, t_start, t_end, dt_days, phase_start, phase_end,
         inv_start_mmol_m2, inv_end_mmol_m2, rate_mmol_o2_m2_d,
         k_m_d, u2_m2_s2, o2_eq_mmol_m3, delta_o2_mmol_m3,
         Fas_mmol_o2_m2_d, rate_corr_mmol_o2_m2_d)
```

### 3e. Plot (existing Step 4)
Add flux-corrected `net` and `night` series so raw vs corrected vs reference are
visible together. Mirror the existing `mine_net`/`mine_night` blocks (lines
~257–269) with `value = rate_corr_mmol_o2_m2_d` and a new `source` level
`"computed (flux-corr)"`; add it to `src_levels` (line ~291) and to the
`scale_linetype_manual` / `scale_shape_manual` maps. For night, remember the
existing **loss convention flip** (`-rate`) — apply the same negation to
`rate_corr` for the night series. Keep `coord_cartesian(ylim = c(-200, 400))`
and the 8-pt `roll_smooth`. Optionally also overlay a smoothed `Fas` series to
show seasonal gas-exchange magnitude.

---

## 4. Critical gotchas

- **Raw `rate` column is sacred.** The verification block (`report_fit`,
  existing lines ~214–244) compares raw `net`/`night` to the reference series and
  the slopes must stay ~1.0. `rate_corr` is *expected* to diverge from the
  references (they carry no flux correction) — that is NOT a regression. Do not
  "correct" the references.
- **MLD > 40 m caveat (document, don't fix).** `F_as` is a surface boundary flux
  applied to the whole 0–40 m inventory. This is exact only when MLD ≤ 40 m. In
  Iceland Basin / Irminger winter, MLD is far deeper than 40 m, so gas exchange
  ventilates water below 40 m and the simple correction over-attributes flux to
  the layer in winter. Add a comment. An MLD-aware version (using `mlotst` from
  `Data/Raw/cmems_obs-mob_glo_phy_my_0.125deg_P1M-m_*.nc`) is a later extension.
- **ERA5 array dimension order** — verify with `dim()` before indexing
  `[i, j, kt]`; ERA5-from-CDS is usually `[lon, lat, time]` but confirm.
- **Surface bin reorder** — `SA/CT/rho/tt` must be indexed on the same
  reordered (`order(depth)`) levels as `depth`, or surface T/S will be wrong.
- **Sign**: `F_as` positive into ocean; `rate_corr = rate − F_as`. Double-check
  against `delta_o2`: `delta_o2 > 0` ⇒ `F_as < 0`.
- **Keep units in O2** — no O2→C `/2` anywhere (consistent with the existing
  script).

---

## 5. New file: `Scripts/download_era5_wind.py`

```python
# Downloads ERA5 hourly 10 m wind + MSL pressure over the float-3902681 box.
# Requires a ~/.cdsapirc with Copernicus CDS credentials (manual setup).
# Run with the project download env:
#   C:/Users/petit/anaconda3/envs/cmts_learn_olci/python.exe Scripts/download_era5_wind.py
import cdsapi, os

os.makedirs("Data/Raw/ERA5", exist_ok=True)
c = cdsapi.Client()
c.retrieve(
    "reanalysis-era5-single-levels",
    {
        "product_type": "reanalysis",
        "variable": [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "mean_sea_level_pressure",
        ],
        "date": "2024-10-01/2026-02-11",
        "time": [f"{h:02d}:00" for h in range(24)],
        "area": [65, -40, 58, -10],   # N, W, S, E
        "format": "netcdf",
    },
    "Data/Raw/ERA5/era5_wind_slp_3902681.nc",
)
print("wrote Data/Raw/ERA5/era5_wind_slp_3902681.nc")
```
> Note: the current CDS backend may return `.nc` named variables `u10`, `v10`,
> `msl` with dims `valid_time`/`latitude`/`longitude`. If the variable/time names
> differ from §3a, adjust the `ncvar_get` names and the time-origin accordingly
> (newer CDS files sometimes use `valid_time` = "seconds since 1970-01-01").

---

## 6. Verification checklist

1. Run `download_era5_wind.py`; confirm `Data/Raw/ERA5/era5_wind_slp_3902681.nc`
   exists and ERA5 time range brackets the float JULD range (print both).
2. Run the budget script (command at top). It must complete and write
   `Data/Processed/o2_budget_3902681.csv` (now with the 6 new columns) and
   `Output/o2_budget_3902681_comparison.png`.
3. Physical sanity (print summaries):
   - `k_m_d` ≈ 1–10 m d⁻¹ for 5–15 m s⁻¹ winds.
   - `Fas_mmol_o2_m2_d` order ±1…±100; negative (outgassing) during summer
     supersaturation, positive (ingassing) in winter undersaturation; sign
     opposite to `delta_o2_mmol_m3`.
   - `o2_eq_mmol_m3` ≈ 250–330 mmol m⁻³ for cold subpolar surface water.
4. Regression guard: existing `report_fit()` slopes for raw `net`/`night` still
   ~1.0 (raw `rate` untouched).
5. Visual: corrected curve sits above raw in summer (removing outgassing) and
   below raw in winter (removing ingassing).
6. After implementing, update `.claude/o2_budget_3902681_explained.md` — change
   the §0 "no air-sea flux term" statement and document the new term, columns,
   sign convention, ERA5 dependency, and the MLD>40 m caveat.
```
```

---

## 7. Out of scope (do not do unless asked)

- Bubble injection (Liang 2013 / Bittig) — diffusive only here.
- Reuer (2007) gas-residence-time weighted piston velocity — interval-mean k is
  sufficient for this pass.
- MLD-aware partitioning of the flux below 40 m.
- Modifying `Scripts/ncp_o2budget_float_3902681.R` or the reference CSVs.
