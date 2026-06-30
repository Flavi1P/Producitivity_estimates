# Handoff — vertically-resolved (10 m) nighttime respiration field for float 3902681

> This is a **self-contained implementation handoff**. You do not need prior
> conversation context. It tells you exactly what to add to
> `Scripts/o2_budget_3902681.R`, why, the sign conventions, and how to verify.
> Read `.claude/o2_budget_3902681_explained.md` (esp. §2, §7, §8) for the full
> background on the existing script — this builds directly on it. All paths are
> relative to the repo root `C:\Users\petit\Documents\Producitivity_estimates`.
> Run everything **from the repo root**, R **4.5.2** only.
>
> (When implementing, you may copy this file to
> `.claude/o2_budget_respiration_profile_handoff.md` to sit alongside the other
> budget docs.)

---

## 0. Goal

The existing script already produces a **mixed-layer-integrated** nighttime
community respiration index `R` (Step 5g, §8 of the explainer): one number per
night pair, areal + volumetric. The task is to **vertically resolve** it into a
**depth (10 m bins, 0–200 m) × time** field of respiration rate, in **volumetric
units (mmol O2 m⁻³ d⁻¹)**, with physical fluxes removed as far as is defensible
from the float alone.

These four design choices are **already decided** (do not re-ask):
1. **Window:** dark night only — the existing `dusk → dawn`, `dt < 0.7 d` pairs.
2. **Physical removal:** air-sea exchange **plus an explicit vertical-diffusion
   term** (`Kz · ∂²O2/∂z²`). Entrainment and advection are neglected (justified
   below).
3. **Depth range:** 0–200 m, 10 m bins (20 bins).
4. **Units:** volumetric (mmol O2 m⁻³ d⁻¹) only.

Deliverable: a new section **Step 5h** appended to `Scripts/o2_budget_3902681.R`.
**Do not** modify `Scripts/ncp_o2budget_float_3902681.R` (a different, older
script). **Do not** alter any existing step — Step 5h only appends, and the
existing verification block (NET slope ≈ 1.015, NIGHT slope ≈ 1.020) must stay
identical (your regression guard).

---

## 1. Physical basis & the equation to implement

In a **fixed 10 m bin** during the dark window photosynthesis is zero, so the
bin's O2 tendency is respiration plus physical exchange:
```
dO2/dt(z) = −R(z) + S_airsea(z) + S_diff(z)
```
Solve for respiration (positive = O2 consumed):
```
R(z) = S_airsea(z) + S_diff(z) − dO2/dt(z)        [mmol O2 m⁻³ d⁻¹]
```

Terms:
- **`dO2/dt(z)`** — `(o2_end − o2_start) / dt_days` of the bin-mean O2
  concentration (mmol m⁻³), from the two casts of the pair.
- **`S_airsea(z)`** — air-sea gas exchange enters only the mixed layer and is
  spread over it by mixing. Use `F_as / MLD̄` for bins whose midpoint is **inside**
  the mixed layer (`z_mid < MLD̄`), and **0** below. This is the same `F_as/MLD`
  convention `make_pair_mld` already uses. `F_as` (mmol O2 m⁻² d⁻¹) comes from the
  existing `airsea_flux(d, i, j)` (Wanninkhof 2014; `F_as > 0` = into ocean).
- **`S_diff(z) = Kz · ∂²O2/∂z²`** — flux-convergence of vertical turbulent
  diffusion. `∂²O2/∂z²` is the **centered second difference** of the 10 m-binned
  O2 profile, `(o2[k−1] − 2·o2[k] + o2[k+1]) / Δz²` (Δz = 10 m), computed per cast
  then **averaged over the pair's two casts**. Constant background
  `Kz = 1e-5 m² s⁻¹` (≈ 0.864 m² d⁻¹). In a well-mixed layer the profile is flat
  ⇒ curvature ≈ 0 ⇒ `S_diff` self-vanishes; it activates at the pycnocline where
  the gradient is sharp. Bin ends (top/bottom, no neighbour) get `S_diff = NA → 0`
  or one-sided — set to 0 and note it.

**Why entrainment and advection are dropped (don't add them):**
- Entrainment is **negligible over an 11 h window** — §8 of the explainer: median
  overnight entrainment correction ≈ 0.00 mmol m⁻² d⁻¹ (the MLD barely moves
  overnight).
- Advection: the displacement-vs-rate diagnostic (Step 5c) found Spearman
  `|rate|` vs inter-profile displacement ≈ 0 — no detectable lateral-gradient
  signal.

---

## 2. What already exists in the script to reuse (do NOT rewrite)

All of these are defined **before** Step 5g, so a section appended after Step 5g
sees them:

| Symbol | What it is |
|---|---|
| `o2_prof_list[[prof_index]]` | list with `$depth`, `$o2_vol` (mmol O2 m⁻³), depth-ordered, per cast. **Indexed by original `prof_index`, NOT row number** (§7 gotcha). |
| `layer_mean_o2(pi, z_top, z_bot)` | trapezoidal **bin-mean** O2 over `[z_top, z_bot]` for stored profile `pi`. Returns **NA** if the cast doesn't reach `z_bot` (no extrapolation). **This is your per-bin O2 extractor — no new integration code needed.** |
| `airsea_flux(d, i, j)` | returns `list(k, u2, o2_eq, dO2, Fas)`; `Fas` in mmol O2 m⁻² d⁻¹, +into ocean. |
| `prof` (tibble) | one row per kept profile, columns `prof_index, time, lat, lon, phase ("dawn"/"dusk"), inventory, mld, o2_surf, o2_eq, sst`. Already time-sorted and head/tail trimmed. |
| `SIG_O2` | `1.0` mmol m⁻³ — adjusted-optode precision, for the per-pair uncertainty. |
| `roll_smooth_w(x, w, k)` | inverse-variance-weighted centered rolling mean (use for optional per-bin temporal smoothing; `w = 1/σ²`). |

Sampling reminder (§2 of explainer): night pairs are consecutive `prof` rows with
`phase[i] == "dusk"`, `phase[i+1] == "dawn"`, `dt < 0.7 d` (~11 h). There are
~234 such pairs. This is exactly how the existing `night_tbl` is filtered from
`prev_tbl` — mirror that filter so your night set is self-consistent.

---

## 3. Implementation (append as "Step 5h")

Place it **after Step 5g** (the existing respiration index) and before/after
Step 5c — anywhere after `o2_prof_list`, `prof`, `airsea_flux`, `layer_mean_o2`,
`roll_smooth_w` are defined. Suggested structure:

### 3.1 Parameters
```r
Z_EDGES   <- seq(0, 200, by = 10)              # 21 edges -> 20 bins
KZ_M2_S   <- 1e-5                               # background vertical diffusivity
KZ_M2_D   <- KZ_M2_S * 86400                    # -> m^2 d^-1 (~0.864)
out_csv_respz <- "Data/Processed/o2_budget_3902681_respiration_profile.csv"
out_png_respz <- "Output/o2_budget_3902681_respiration_profile.png"
```

### 3.2 `bin_respiration(d, i, j)` — one night pair → long tibble (20 rows)
For ordered rows `i` (earlier, dusk) and `j` (later, dawn) of `prof`:
1. `dt_days <- as.numeric(difftime(d$time[j], d$time[i], units = "days"))`.
2. Bin-mean O2 for each bin from both casts:
   `o2_i[k] <- layer_mean_o2(d$prof_index[i], Z_EDGES[k], Z_EDGES[k+1])`, same for
   `o2_j[k]`. (Vectorize over `k = 1..20`.)
3. `dO2dt[k] <- (o2_j[k] - o2_i[k]) / dt_days`. (NA where either bin-mean is NA —
   cast too shallow; keep as NA, do not extrapolate.)
4. Curvature per cast then averaged:
   `curv_i[k] <- (o2_i[k-1] - 2*o2_i[k] + o2_i[k+1]) / 100`  (Δz² = 100 m²), same
   for `curv_j`; `curv[k] <- (curv_i[k] + curv_j[k]) / 2`. Ends (`k = 1`, `k = 20`)
   → `0`. `S_diff[k] <- KZ_M2_D * curv[k]`.
5. `mld_bar <- mean(c(max(d$mld[i], 40), max(d$mld[j], 40)))`.
   `fx <- airsea_flux(d, i, j)`.
   `z_mid[k] <- (Z_EDGES[k] + Z_EDGES[k+1]) / 2`.
   `S_airsea[k] <- ifelse(z_mid[k] < mld_bar, fx$Fas / mld_bar, 0)`.
6. `R_vol[k] <- S_airsea[k] + S_diff[k] - dO2dt[k]`.
7. `sigma_R_vol[k] <- SIG_O2 * sqrt(2) / dt_days`  (depth-independent, Fix-3
   convention; same for all bins).
8. Return a tibble (long): `mtime` (pair midpoint = `d$time[i] + (d$time[j]-d$time[i])/2`),
   `t_start, t_end, dt_days, z_top, z_bot, z_mid, mld_bar, o2_start, o2_end,
   dO2dt, S_airsea, S_diff, R_vol, sigma_R_vol`.

### 3.3 Build the night set and map over it
```r
n <- nrow(prof)
night_idx <- which(prof$phase[-n] == "dusk" & prof$phase[-1] == "dawn" &
                   as.numeric(difftime(prof$time[-1], prof$time[-n], units="days")) < 0.7)
respz_tbl <- purrr::map_dfr(night_idx, ~ bin_respiration(prof, .x, .x + 1))
```
(Match the existing code's idioms — it uses `map_dfr` and base indexing already.)

### 3.4 CSV
`write_csv(respz_tbl, out_csv_respz)` — long form, one row per (night pair × bin).
`cat(sprintf("Wrote %s (%d rows)\n", out_csv_respz, nrow(respz_tbl)))`.

### 3.5 Plot — `R(z, t)` heatmap
- `geom_raster`/`geom_tile`: `x = mtime`, `y = z_mid`, `fill = R_vol`.
- `scale_y_reverse()` (depth down).
- Diverging fill centered at 0: `scale_fill_gradient2(low=…, mid="white", high=…,
  midpoint=0, limits=c(-6, 6), oob=scales::squish)` — the ±6 range matches the
  volumetric R scale in §8 so winter noise doesn't wash out the summer signal.
- Overlay MLD: a line of `mld_bar` (or `prof$mld`) vs `mtime` so the in-ML vs
  interior split is visible; dashed `geom_hline(yintercept = 40)` for the floor.
- `theme_bw()`, ~11×6, dpi 200. Title makes clear: dark-window respiration,
  flux + diffusion corrected, volumetric.
- Optional second panel or alpha-mask: grey out bins where `abs(R_vol) <
  sigma_R_vol` (winter is at the detection limit, §8) so the trustworthy
  (summer/interior) signal stands out.

### 3.6 Printed sanity (stdout, mirror Step 5g style)
- `n` night pairs, bins per pair.
- Per-season (`winter = month ∈ {11,12,1,2,3}`, else summer) median `R_vol`,
  fraction `R_vol > 0`.
- Median `|S_airsea|` (should be ~0 below the ML, non-zero in-ML) and median
  `|S_diff|` (should be ≪ `|R_vol|` in the stratified interior — confirms the
  diffusion correction is small there).
- **Kz sensitivity:** report median `|S_diff|` recomputed with `Kz = 1e-4 m² s⁻¹`
  (×10) so the reader sees how much the diffusion choice moves R.

---

## 4. Sign & units conventions (get these right)

| Quantity | Sign / units |
|---|---|
| `R_vol` | **positive = O2 consumed** (respiration). mmol O2 m⁻³ d⁻¹. |
| `S_airsea = Fas/MLD̄` | `Fas > 0` = O2 into ocean (a source). Enters as `+` in `R = S_airsea + S_diff − dO2dt`. |
| `S_diff = Kz·∂²O2/∂z²` | positive curvature ⇒ flux convergence ⇒ O2 source ⇒ `+`. |
| `dO2dt` | signed concentration tendency; O2 dropping ⇒ negative ⇒ raises `R`. |

- **Volumetric only.** No `z_int`/MLD multiplier anywhere (that is the
  depth-amplification artifact the volumetric form exists to avoid — §7a). If a
  reviewer wants areal per bin it is recoverable as `R_vol × 10` (bin thickness),
  but the deliverable is volumetric.
- **No O2 → carbon `/2`.** This script keeps everything in O2 (explainer §2.
  Gotchas). Do not divide by 2.
- **Never extrapolate below the deepest level.** `layer_mean_o2` already returns
  NA for bins the cast doesn't reach — propagate the NA, don't fill it.

---

## 5. Verification checklist

1. Run: `& "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/o2_budget_3902681.R`
   (from repo root). It should finish without error and print the existing blocks
   **plus** your new Step 5h summary.
2. **Regression guard:** the existing verification block must be unchanged —
   `NET … slope(mine~0+ref)=1.015`, `NIGHT … slope=1.020` (§4 of explainer). If
   these move, you accidentally touched shared state — revert and isolate Step 5h.
3. New CSV `Data/Processed/o2_budget_3902681_respiration_profile.csv` exists,
   ~234 pairs × 20 bins (minus NA-dropped deep bins for short casts), long form.
4. **Hand-check one row:** pick a pair+bin and confirm
   `R_vol == S_airsea + S_diff − dO2dt` and `dO2dt == (o2_end − o2_start)/dt_days`.
5. Open `Output/o2_budget_3902681_respiration_profile.png`: expect a coherent
   **subsurface respiration maximum in the productive season** (summer interior
   `R_vol` ~ +2…+6) beneath a near-zero/ noisy mixed layer; winter columns
   noise-dominated (matches §8 — trust summer/interior, not individual winter
   bins). MLD overlay should look sensible (deep in winter, shallow in summer).
6. Sanity stdout: in the stratified interior `|S_diff| ≪ |R_vol|`; `S_airsea` is
   zero below the ML and only acts in-ML.

---

## 6. Gotchas (carried from the explainer)

- `o2_prof_list` is indexed by **`prof_index`**, not the `prof` row number — always
  look up `d$prof_index[i]`.
- `castr` masks `stats::integrate`/`stats::smooth` on attach (harmless here).
- The `max(): no non-missing arguments` warning in the Step-1 loop is pre-existing
  and unrelated.
- Read the float from the **Sprof** (already loaded into `o2_prof_list`); do **not**
  use `Data/Intermediate/Floats/argo_3902681_interp.csv` (no time-of-day, can't
  separate dusk/dawn).
- `Kz = 1e-5 m² s⁻¹` is a **chosen constant background diffusivity** — document it
  as a parameter at the top of Step 5h and in the plot subtitle so it's obviously
  tunable. The §3.6 sensitivity line shows its leverage.
