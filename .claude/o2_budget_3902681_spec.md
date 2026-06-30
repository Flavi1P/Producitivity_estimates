# Implementation spec — profile-by-profile O2 budget for float 3902681

This is a **hands-off, self-contained** specification. You (the implementing
agent) do not need prior conversation context. Build exactly what is described
here. All paths are relative to the repo root
`C:\Users\petit\Documents\Producitivity_estimates`. Run everything **from the
repo root**.

---

## 0. Goal

From BGC-Argo float **3902681**, compute a profile-by-profile **dissolved-oxygen
budget**: the rate of change of the 0–40 m O2 inventory between profiles,
expressed in **mmol O2 m⁻² d⁻¹**. Produce three budget variants, write them to a
CSV, and make one overlay plot comparing them against two pre-existing reference
series.

Deliverables:
1. New script `Scripts/o2_budget_3902681.R`.
2. Output CSV `Data/Processed/o2_budget_3902681.csv`.
3. Output plot `Output/o2_budget_3902681_comparison.png`.

**Do not** modify the existing `Scripts/ncp_o2budget_float_3902681.R`.

---

## 1. Confirmed methodology decisions (do not deviate)

| Decision | Value |
|---|---|
| Integration depth | **0–40 m**, fixed |
| O2 unit conversion | µmol/kg → mmol/m³ via **gsw in-situ density** |
| Budget timestamp | **midpoint** of the two profiles in each pair |
| Air-sea flux | **none** (omit) |
| Output units | **mmol O2 m⁻² d⁻¹**, signed (negative = O2 loss). **No O2→C /2 conversion.** |
| Data source | read `Data/Raw/Floats/3902681_Sprof.nc` directly |

---

## 2. Environment / how to run

R **4.5.2** only (not 4.1.3):

```powershell
& "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/o2_budget_3902681.R
```

Required packages (all already used elsewhere in this repo): `ncdf4`,
`tidyverse`, `lubridate`, `gsw`, `zoo`. Do not add new dependencies.

---

## 3. Input data

### 3.1 Sprof NetCDF — `Data/Raw/Floats/3902681_Sprof.nc`

2D variables are `[N_LEVELS × N_PROF]` (level × profile). Read with
`ncdf4::ncvar_get`. Pattern to copy: `Scripts/sprof_to_csv.R` (flatten 2D arrays,
per-profile interpolation with `zoo::na.approx`).

Variables to read:

| Variable | Meaning | Notes |
|---|---|---|
| `JULD` | profile time | days since **1950-01-01 00:00 UTC**. `time = as.POSIXct(JULD*86400, origin="1950-01-01", tz="UTC")` |
| `LATITUDE`, `LONGITUDE` | per-profile position | length `N_PROF` |
| `PRES` | pressure | decibar; 2D |
| `TEMP` | in-situ temperature | °C; 2D |
| `PSAL` | practical salinity | 2D |
| `DOXY_ADJUSTED` | dissolved O2 | **µmol/kg**; 2D |
| `DOXY_ADJUSTED_QC` | QC flags | char array, 2D; keep flags `1,2,5,8` |

Confirmed facts about this file: ~493 profiles, 2024-10-02 → 2026-02-10,
DOXY_ADJUSTED surface values ~205–358 µmol/kg, PRES in decibar.

`DOXY_ADJUSTED_QC` is read as a character matrix (one char per level). Parse each
char to integer; treat a level as good if QC ∈ {1,2,5,8}. (If QC parsing is
awkward, an acceptable fallback is to keep all finite `DOXY_ADJUSTED` values, but
prefer QC filtering.)

### 3.2 Reference CSVs (comma-separated, `read_csv`)

- `Data/Processed/O2_float_net_change.csv` — columns:
  `mtime, Date, Int_dO2dt_40m_mmol_m2_d`  (195 rows incl. header → 194 data)
- `Data/Processed/O2_float_night_loss.csv` — columns:
  `mtime, Date, Int_O2_loss_40m_mmol_m2_d`  (197 rows)

`mtime` is a **MATLAB datenum**. Convert:
`time = as.POSIXct((mtime - 719529) * 86400, origin = "1970-01-01", tz = "UTC")`.
The value columns are already in mmol O2 m⁻² d⁻¹. **Use them as-is (do NOT divide
by 2).** These are the "net change" and "night loss" series respectively. Their
exact derivation is unknown — they are for comparison only.

---

## 4. Step-by-step algorithm

### Step 1 — per-profile O2 inventory (0–40 m)

For each profile column `p` (1..N_PROF):

1. Extract that profile's vectors: `PRES[,p]`, `TEMP[,p]`, `PSAL[,p]`,
   `DOXY_ADJUSTED[,p]`, QC[,p].
2. Drop levels where DOXY_ADJUSTED is NA or QC not in {1,2,5,8}, or T/S NA.
3. Compute **in-situ density** with gsw (pressure in dbar, lon/lat scalars):
   ```r
   SA  <- gsw_SA_from_SP(PSAL, PRES, lon[p], lat[p])
   CT  <- gsw_CT_from_t(SA, TEMP, PRES)
   rho <- gsw_rho(SA, CT, PRES)          # kg m-3
   o2_vol <- DOXY_ADJUSTED * rho / 1000  # mmol O2 m-3
   ```
4. Convert pressure to depth: `depth <- gsw_z_from_p(PRES, lat[p]) * -1` (gives
   positive metres down), or simply treat `PRES` ≈ depth (difference is <0.5 % in
   top 40 m — using PRES directly is acceptable; state which you used).
5. Regrid `o2_vol` onto a regular **0–40 m at 1 m** grid:
   ```r
   grid <- 0:40
   o2_grid <- approx(x = depth, y = o2_vol, xout = grid, rule = 2)$y
   ```
   `rule = 2` extends the shallowest/deepest available value to the grid ends
   (handles the small near-surface gap). **Require** that the profile actually
   reaches ≥ 40 m (`max(depth) >= 40`) and has ≥ ~5 valid levels in 0–60 m;
   otherwise set inventory = NA and drop the profile.
6. **Inventory** = trapezoidal integral over 0–40 m:
   ```r
   inventory <- sum(diff(grid) * (head(o2_grid,-1) + tail(o2_grid,-1)) / 2)
   ```
   Units: mmol O2 m⁻².

Collect a tidy table `prof` with one row per profile:
`prof_index, time, lat, lon, lst, phase, inventory`.

- Local solar time: `lst <- (hour(time) + minute(time)/60 + lon/15) %% 24`.
- Phase: `phase <- ifelse(lst < 12, "dawn", "dusk")` (same rule as
  `Scripts/float_3902681_trajectory.py`).

Sort `prof` by `time` ascending and drop NA-inventory rows.

> Sampling pattern for sanity (do not hard-code, just expect it): the float runs
> a ~48 h cycle — **dusk** profile (LST ~18.4 h) → ~11 h overnight → **pre-dawn**
> profile (LST ~5.6 h) → ~37 h → next dusk. The first ~5 profiles (early Oct
> 2024) are near midday and irregular; the regular regime starts ~2024-10-12.

### Step 2 — build the three budgets

Define a helper that, given an ordered profile pair (row i, row j = later),
returns:
```
dt_days = as.numeric(difftime(time_j, time_i, units = "days"))
rate    = (inventory_j - inventory_i) / dt_days     # mmol O2 m-2 d-1, signed
mtime   = time_i + (time_j - time_i)/2              # midpoint
```

Compute over the time-sorted `prof` table:

1. **`prev`** — every consecutive pair `(i, i+1)`. Tag the transition
   `phase_i -> phase_j`. This is the broadest budget ("with the previous
   profile").
2. **`night`** (night loss) — the subset of consecutive pairs where
   `phase_i == "dusk"` AND `phase_j == "dawn"` AND `dt_days < 0.7` (the ~11 h
   dark interval). Expected: O2 loss (mostly negative) = respiration.
3. **`net`** (net change) — consecutive **same-phase** pairs: for the subseries
   of `dusk` profiles take each consecutive `dusk[k] -> dusk[k+1]`, and likewise
   for `dawn`. (~48 h apart; the prompt's "two dusk or two dawn"). Net community
   O2 change over a full diel cycle.

Combine into one long table with a `type` column (`prev` / `night` / `net`).

### Step 3 — write CSV

`Data/Processed/o2_budget_3902681.csv` (comma-separated, `write_csv`), columns:
```
type, mtime, t_start, t_end, dt_days, phase_start, phase_end,
inv_start_mmol_m2, inv_end_mmol_m2, rate_mmol_o2_m2_d
```

### Step 4 — comparison plot

`Output/o2_budget_3902681_comparison.png` (ggplot2, `theme_bw`, width 11 ×
height 6, dpi 200).

- X axis: date. Y axis: `mmol O2 m⁻² d⁻¹`. Add `geom_hline(yintercept = 0)`.
- Plot, distinguishing **mine vs reference** by linetype/shape and **quantity**
  by colour:
  - my **net change** (`type == "net"`) vs reference net change
    (`Int_dO2dt_40m_mmol_m2_d`),
  - my **night loss** (`type == "night"`) vs reference night loss
    (`Int_O2_loss_40m_mmol_m2_d`).
  - Optionally include my `prev` series faintly.
- Raw diel rates are noisy: overlay points **plus** a light smoother (loess
  `span ≈ 0.2–0.3`, or a centered rolling mean ~6–10 points via
  `zoo::rollapply`) so trends are legible. Keep both mine and reference on the
  same smoothing so the comparison is fair.
- Title: `Float 3902681: O2 budget (0–40 m) — computed vs reference`.
- `ggsave(...)` and `cat()` the saved path.

---

## 5. Verification (run after building)

1. Script runs clean and writes the CSV + PNG.
2. Counts: expect on the order of ~190–195 `net` and ~190–195 `night` estimates
   spanning 2024-10 → 2026-02 (reference files have 194 / 196 data rows). A wildly
   different count means the phase/pairing logic is wrong.
3. Quantitative agreement: join my `net` (and `night`) estimates to the reference
   on nearest `mtime` (tolerance ~1 day), report Pearson correlation and median
   difference. Strong correlation + near-zero bias validates the method. A
   constant multiplicative offset would implicate the density/units step; a
   constant additive offset or scatter would implicate pairing/timestamps.
4. Eyeball: night loss predominantly negative; net change oscillating around 0
   with plausible seasonal structure (larger magnitude in productive season).

Report the correlation/bias numbers and the two output paths back to the user.

---

## 6. Pitfalls

- The interp CSV `Data/Intermediate/Floats/argo_3902681_interp.csv` has **no
  time-of-day** (date only) — it cannot separate dusk/dawn. Do **not** use it;
  read the Sprof.
- Keep units in O2 throughout; the *existing* plotting script divides O2 by 2 to
  get carbon — we do **not**.
- `gsw` functions are vectorised; pass whole profile vectors. lon/lat are scalars
  per profile (recycle).
- QC array is character; parse carefully. If in doubt, log how many levels pass
  QC for the first few profiles.
- Guard against profiles that don't reach 40 m (drop them; expect a few gaps).
