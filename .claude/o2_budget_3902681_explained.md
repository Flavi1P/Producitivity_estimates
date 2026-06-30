# How the float-3902681 O2 budget works — handoff for the next agent

This is a **self-contained explainer**. You do not need prior conversation
context. It describes what `Scripts/o2_budget_3902681.R` does, why each step is
there, the sign conventions, and the validated findings. All paths are relative
to the repo root `C:\Users\petit\Documents\Producitivity_estimates`. Run
everything **from the repo root**.

The original build spec is `.claude/o2_budget_3902681_spec.md`. This document
reflects the **current, working** state of the script (which has since been
refined past that spec). Where they differ, **this document and the script
win**.

---

## 0. One-paragraph summary

From BGC-Argo float **3902681** the script computes a profile-by-profile
**dissolved-O2 budget**: the rate of change of the **0–40 m O2 inventory**
between profiles, in **mmol O2 m⁻² d⁻¹**. It produces three budget variants
(`prev`, `night`, `net`) and writes all of them (with the flux columns below) to
a CSV. It now **also** adds a **diffusive air-sea O2 flux** term (Wanninkhof
2014, ERA5-forced) and reports a flux-corrected rate (`rate_corr = rate − F_as`).
The comparison PNG overlays the **flux-corrected** `net`/`night` series against
the two pre-existing reference series (the raw and `/2` series stay in the CSV).
There is still **no O2→carbon `/2` conversion** — everything stays in O2.
**Bubble injection is now included** (Liang et al. 2013): the diffusive `F_as` is
left untouched and two bubble terms (`Fc` complete-dissolution, `Fp`
partial-dissolution) are **added** to form `F_as_total = F_as + Fc + Fp`, with a
separate `rate_corr_total = rate − F_as_total`. The diffusive-only `F_as` /
`rate_corr` columns and the validated §4 slopes are unchanged — bubbles are a
new, parallel set of columns. See §2c.

It **also** now builds a second, parallel budget — the **MLD-based budget** —
that integrates O2 not over the fixed 0–40 m slab but from the surface to the
**deeper of the two profiles' mixed-layer depths** for each pair (floored at
40 m). This deepens the inventory into the full mixed layer (so winter sub-40 m
ventilation is captured and the air-sea flux acts on the physically correct
depth). See §7.

Outputs:
- `Data/Processed/o2_budget_3902681.csv` (0–40 m budget, with 6 air-sea flux columns)
- `Data/Processed/o2_budget_3902681_mld.csv` (MLD-based budget, with `z_int_m`/`mld_*` columns)
- `Output/o2_budget_3902681_comparison.png`
- `Output/o2_budget_3902681_mld_comparison.png` (0–40 m vs MLD-deepest net change + integration depth)
- `Output/o2_budget_3902681_airsea_diagnostics.png` (the flux term, step by step)

Run:
```powershell
& "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/o2_budget_3902681.R
```
R **4.5.2** only. Packages (all already used in this repo): `ncdf4`,
`tidyverse`, `lubridate`, `gsw`, `zoo`, `castr` (for `mld()`), `patchwork` (for
the two-panel MLD figure). No genuinely new dependencies — `castr`/`patchwork`
are already used by `Scripts/format_for_ncp.R`.

---

## 1. Inputs

### Float NetCDF — `Data/Raw/Floats/3902681_Sprof.nc`
2D variables are `[N_LEVELS × N_PROF]` (level × profile). ~493 profiles,
2024-10-02 → 2026-02-10. Variables read: `JULD` (days since 1950-01-01),
`LATITUDE`, `LONGITUDE` (per-profile), `PRES` (dbar), `TEMP`, `PSAL`,
`DOXY_ADJUSTED` (**µmol/kg**), `DOXY_ADJUSTED_QC`.

`DOXY_ADJUSTED_QC` comes back as a **character vector of length N_PROF**, each
element an `N_LEVELS`-long string with one flag char per level. `parse_qc_col()`
splits it per character to integers; a level is kept only if its flag ∈
**{1, 2, 5, 8}** (and T/S/DOXY/PRES are all finite).

### Reference CSVs (comma-separated, read with `read_csv`)
- `Data/Processed/O2_float_net_change.csv` — `mtime, Date, Int_dO2dt_40m_mmol_m2_d`
- `Data/Processed/O2_float_night_loss.csv` — `mtime, Date, Int_O2_loss_40m_mmol_m2_d`

`mtime` is a **MATLAB datenum**:
`time = as.POSIXct((mtime - 719529) * 86400, origin = "1970-01-01", tz = "UTC")`.
The value columns are already in mmol O2 m⁻² d⁻¹. **Used as-is — never divided by
2.** Their exact derivation is unknown; they are for comparison only. (See §5 on
their apparent sampling.)

---

## 2. The algorithm, step by step

### Step 1 — per-profile 0–40 m O2 inventory
For each profile column `p`:
1. Keep levels passing QC (above) with finite T/S/DOXY/PRES; skip the profile if
   fewer than 5 levels survive.
2. **Convert µmol/kg → mmol m⁻³ via gsw in-situ density** (this is the only place
   units could go wrong):
   ```r
   SA  <- gsw_SA_from_SP(PSAL, PRES, lon[p], lat[p])
   CT  <- gsw_CT_from_t(SA, TEMP, PRES)
   rho <- gsw_rho(SA, CT, PRES)        # kg m-3
   o2_vol <- DOXY_ADJUSTED * rho / 1000  # mmol O2 m-3
   ```
3. Depth from pressure: `depth <- gsw_z_from_p(PRES, lat[p]) * -1` (positive m
   down).
4. **Require** the profile reaches ≥ 40 m and has ≥ 5 valid levels in 0–60 m;
   otherwise inventory = NA and the profile is dropped.
5. Regrid `o2_vol` onto `grid <- 0:40` (1 m) with `approx(..., rule = 2)` (rule 2
   extends the shallowest/deepest value to the grid ends, covering the small
   near-surface gap).
6. **Inventory** = trapezoidal integral over 0–40 m (units mmol O2 m⁻²):
   ```r
   inventory <- sum(diff(grid) * (head(o2_grid,-1) + tail(o2_grid,-1)) / 2)
   ```

In the same loop, two more per-profile quantities are computed and stored for
the MLD-based budget (§7):
- **MLD** via `castr::mld(sigma0, depth, ref.depths = 0:2, criteria = 0.03,
  default.depth = 300)` — the **same call as `Scripts/format_for_ncp.R`**. `sigma0`
  is `gsw_sigma0(SA, CT)` on the depth-ordered levels; `default.depth = 300` is
  returned when the 0.03 kg m⁻³ criterion is never met (near-homogeneous column or
  mixing deeper than the cast).
- the **full depth-ordered `(depth, o2_vol)` profile** is kept in a list
  (`o2_prof_list`, indexed by the original `prof_index` so it survives the
  prof-table filtering), so each cast can be re-integrated to an arbitrary bottom
  depth later.

This yields the `prof` table: one row per profile with
`prof_index, time, lat, lon, lst, phase, inventory, mld`.

- **Local solar time:** `lst <- (hour(time) + minute(time)/60 + lon/15) %% 24`.
- **Phase:** `phase <- ifelse(lst < 12, "dawn", "dusk")` (same rule as
  `Scripts/float_3902681_trajectory.py`).

`prof` is then sorted by `time`, NA-inventory rows dropped, and **the first and
last 5 profiles are trimmed** (`slice`): the early ones are irregular near-midday
casts that produce huge spurious diel rates; the tail is trimmed symmetrically.
After trimming, ~477 of ~494 profiles remain.

### Sampling pattern (so the pairing logic makes sense)
The float runs a ~48 h cycle: **dusk** profile (LST ~18 h) → ~11 h overnight →
**pre-dawn** profile (LST ~6 h) → ~37 h → next dusk. So per cycle there is one
dusk and one dawn profile, ~2 profiles/day.

### Step 2 — the three budgets
`make_pair(d, i, j, type)` takes an ordered pair (row `j` later than `i`) and
returns the signed rate plus a **midpoint timestamp**:
```
dt_days = difftime(time_j, time_i, "days")
rate    = (inventory_j - inventory_i) / dt_days   # signed dO2/dt, mmol O2 m-2 d-1
mtime   = time_i + (time_j - time_i)/2            # midpoint
```

The three variants:
1. **`prev`** — every consecutive pair `(i, i+1)`. The broadest budget. ~476
   pairs.
2. **`night`** — the subset of consecutive pairs with `phase_start == "dusk"`,
   `phase_end == "dawn"`, and `dt_days < 0.7` (the ~11 h dark interval). This is
   the overnight respiration signal. ~234 pairs.
3. **`net`** — consecutive **same-phase** pairs: dusk→dusk and dawn→dawn,
   computed on each phase's subseries then combined. ~48 h apart; net community
   O2 change over a full diel cycle. ~475 pairs (note: this is **both** same-phase
   series, so ~2× the per-cycle count — see §5).

All three are bound into one long `budget` table with a `type` column.

### Step 2b — air-sea diffusive O2 flux (Wanninkhof 2014)
Computed inside `make_pair()`, so all three budget types carry it. The float
loop additionally stores three per-profile surface quantities (means over the
top **0–10 m**, taken from the depth-reordered levels): `o2_surf` (mmol m⁻³),
`o2_eq` (Garcia & Gordon 1992 solubility via `gsw_O2sol`, µmol/kg × ρ/1000 →
mmol m⁻³ at 1 atm) and `sst` (in-situ °C). For each pair:

```
T_s = mean(sst_i, sst_j)
Sc  = 1920.4 − 135.6 T_s + 5.2122 T_s² − 0.10939 T_s³ + 0.00093777 T_s⁴   # W2014 Table 1
<U2> = interval mean of (u10² + v10²)            # ERA5 hourly, nearest grid cell to midpoint
k   = 0.251 · <U2> · (Sc/660)^(−0.5) · 0.24      # cm h⁻¹ → m d⁻¹  (0.24 = 24/100)
o2_eq_p = o2_eq · (msl / 101325)                 # SLP correction; water-vapour term omitted
dO2 = mean(o2_surf) − mean(o2_eq_p)              # >0 = supersaturated
F_as = −k · dO2                                  # mmol O2 m⁻² d⁻¹, +into ocean
rate_corr = rate − F_as
```

**Sign:** `F_as` positive = O2 **into** the ocean. Supersaturation (`dO2 > 0`)
⇒ `F_as < 0` (outgassing); undersaturation ⇒ `F_as > 0` (ingassing) — `F_as` is
the mirror of `dO2`. ERA5 forcing comes from
`Data/Raw/ERA5/era5_wind_slp_3902681.nc` (hourly `u10`/`v10`/`msl`,
`hours since 1900-01-01`, dims `[lon × lat × time]`), produced by
`Scripts/download_era5_wind.py`. That script's concatenation step **drops the
`expver`/`number` coords and writes NETCDF4_CLASSIC** — the raw CDS vlen-string
`expver` coordinate **segfaults R's `ncdf4`**, so never feed the raw chunks to R.

**MLD > 40 m caveat:** for the **0–40 m budget**, `F_as` is a surface boundary
flux applied to the whole 0–40 m inventory — exact only when MLD ≤ 40 m. In
Iceland Basin / Irminger winter the MLD is far deeper than 40 m, so gas exchange
ventilates water below 40 m and the simple correction over-attributes flux to
the 0–40 m layer in winter. The **MLD-based budget (§7)** addresses this by
integrating to the deeper of the two MLDs, so `F_as` acts on the whole mixed
layer — the depth a surface gas-exchange term physically belongs to.

The flux helper itself was refactored into a shared `airsea_flux(d, i, j)`
function (the surface term is depth-independent), called by both `make_pair`
(0–40 m) and `make_pair_mld` (§7). This is a pure refactor — the 0–40 m
verification slopes are unchanged (§4).

### Step 2c — bubble-mediated O2 flux (Liang et al. 2013)
Computed inside `airsea_flux()`, **added on top of** the diffusive `F_as` (which
is left exactly as in §2b — its value, the `Fas_mmol_o2_m2_d` column, and the §4
slopes are untouched). Follows the `fas_L13.m` reference implementation
(Nicholson `gas_toolbox`). Two terms, both **+ into the ocean**, evaluated at
**each hourly ERA5 wind step over the interval and then averaged** (the flux is
strongly nonlinear in wind, `~ u*_w^3.86`, so averaging the wind first would
under-estimate it):

```
Cd     = 0.0012 (u10≤11) | 4.9e-4 + 6.5e-5·u10 (11<u10<20) | 0.0018 (u10≥20)
u*_a   = u10·√Cd                                  # air-side friction velocity
u*_w   = u*_a / √(ρ_w/ρ_air)   (ρ_air = 1.225)    # water-side friction velocity
slpc   = (msl/101325 − pH2O) / (1 − pH2O)         # Weiss&Price 1980 vapour, RH=1
Kb     = 5.5·u*_w^2.76·(Sc/660)^(−2/3) · 86400    # bubble transfer velocity, m d⁻¹
ΔP     = 1.5244·u*_w^1.06                          # fractional bubble overpressure
Fc     = X_g·5.56·u*_w^3.86 ·1000·86400            # complete-dissolution, ALWAYS ≥0
Fp     = Kb·(O2_eq·(1+ΔP)·slpc − O2_surf)          # partial, saturation-dependent
F_as_total = F_as + Fc + Fp
rate_corr_total = rate − F_as_total
```

`X_g = 0.20946` (O2 mole fraction in dry air). `Sc` is the same W2014 O2 Schmidt
number used by the diffusive term; `O2_eq` is the stored 1-atm equilibrium
(`o2_eq_vec`); `O2_surf`, `ρ_w`, `S`, `T` are the pair-mean surface values (new
per-profile stores `sss_vec`, `rho_surf_vec`). `Fc` injects O2 regardless of
saturation (large bubbles fully collapse) so it is always ≥ 0; `Fp` exchanges
against an overpressure-raised target, so in undersaturated high-wind winter it
drives a large **ingassing**. For this float bubbles are a **material ~40 % of
the median total flux** (high-wind subpolar North Atlantic) — the reason they
could not be neglected. `airsea_flux()` returns `Fc`, `Fp`, `Fas_bub` (=Fc+Fp)
and `Fas_tot` alongside the diffusive `Fas`; both `make_pair` and `make_pair_mld`
expose them. A sign/magnitude summary (`Fc` all ≥ 0, bubble share of total) is
printed in the verification block.

**Scope caveat:** bubbles are wired into the **air-sea correction columns** of
both the 0–40 m and MLD CSVs (`rate_corr_total*`), but the **moving mixed-layer
residual** (§7a Fix 2, `rate_resid_*`) and the **nighttime respiration index R**
(§8) still subtract the **diffusive-only** `F_as` (they use `fx$Fas`, not
`fx$Fas_tot`). Folding bubbles into those derived biological residuals is a
deliberate next step, not yet done — if you do it, swap `fx$Fas` → `fx$Fas_tot`
in the `make_pair_mld` Fix-2 block and re-validate the §8 R seasonal line.

### Step 3 — CSV
`Data/Processed/o2_budget_3902681.csv`, columns:
```
type, mtime, t_start, t_end, dt_days, phase_start, phase_end,
inv_start_mmol_m2, inv_end_mmol_m2, rate_mmol_o2_m2_d,
k_m_d, u2_m2_s2, o2_eq_mmol_m3, delta_o2_mmol_m3,
Fas_mmol_o2_m2_d, rate_corr_mmol_o2_m2_d,
Fc_mmol_o2_m2_d, Fp_mmol_o2_m2_d, Fas_bubble_mmol_o2_m2_d,
Fas_total_mmol_o2_m2_d, rate_corr_total_mmol_o2_m2_d
```
The last five columns are the §2c bubble terms: `Fas_total = Fas + Fc + Fp` and
`rate_corr_total = rate − Fas_total` (diffusive + bubble). The diffusive-only
`Fas_mmol_o2_m2_d` / `rate_corr_mmol_o2_m2_d` are kept unchanged. The MLD CSV
(`o2_budget_3902681_mld.csv`) gets the same areal bubble columns plus volumetric
`Fas_total_vol_mmol_o2_m3_d` / `rate_corr_total_vol_mmol_o2_m3_d`.
**`rate_mmol_o2_m2_d` in the CSV is always the signed dO2/dt** (negative = O2
loss), for all three types, and is **left untouched** by the flux term — the
correction lives in the separate `rate_corr_mmol_o2_m2_d` column. The
loss-convention flip described below is applied only for the night
comparison/plot, **not** to the stored CSV.

### Step 4 — verification (printed to stdout)
`report_fit()` nearest-joins my estimates to a reference series on `mtime`
(tolerance 1 day) and prints **n**, **Pearson r**, **median(mine − ref)**, and
the **slope through the origin** (`lm(value ~ 0 + ref_val)`). The slope is the
key diagnostic: ~1 ⇒ 1:1 agreement; ~2 ⇒ mine is twice the reference (which would
imply a 0.5 coefficient on the reference side); ~0.5 ⇒ mine is half.

### Step 5 — comparison plot
`Output/o2_budget_3902681_comparison.png` (`theme_bw`, 11×6, dpi 200).
- **Smoothed lines only — no per-pair scatter.** Each series is the centered
  **8-point rolling mean** (`roll_smooth`) of the per-cycle estimates; the raw
  diel points and the faint grey `prev` cloud were intentionally dropped for
  legibility.
- **Only two source series are shown:** `computed (flux-corr)` (the
  `rate_corr = rate − F_as` estimate, **solid**) and `reference` (**dashed**).
  The raw `computed` and the `computed / 2` coefficient-check series are no
  longer plotted (they still live in the CSV / verification block).
- Colour = **quantity** (`net change` **blue** `#3690c0`, `night loss` **red**
  `#cb181d`); linetype = **source** (solid = flux-corr, dashed = reference).
- The flux-corrected `net change` sits **above** the raw inventory rate in
  summer (removing outgassing) and below it in winter (removing ingassing);
  `night loss` uses the same `−rate` loss flip applied to `rate_corr`.
- `coord_cartesian(ylim = c(-200, 400))` clips residual diel noise so the
  seasonal signal stays legible.
- `plot_df` still assembles all sources (raw `computed`, `computed / 2`,
  flux-corr, reference); the `src_levels` filter just before the `ggplot` call
  is what restricts the figure to flux-corr vs reference — widen it to bring the
  other series back.

### Step 6 — air-sea flux diagnostics plot
`Output/o2_budget_3902681_airsea_diagnostics.png` (`theme_bw`, 10×11, dpi 200).
Five stacked facets built from the dense `prev` series, showing the flux chain
end to end: wind speed `sqrt<U2>` → piston velocity `k` → O2 saturation (%) →
surface disequilibrium `dO2` → `F_as`. Blue line = 8-pt rolling mean; dashed
references at 100 % / 0. The script also prints a min/median/mean/max sanity
summary of `k`, `o2_eq`, `dO2`, `F_as` and a `dO2>0 ⇒ F_as<0` sign check.

---

## 3. Sign conventions — read this carefully

There are two conventions in play and the script deliberately keeps them
straight:

| Quantity | Reference column | Reference sign | My CSV `rate` | In plot/verification |
|---|---|---|---|---|
| **net change** | `Int_dO2dt_40m_mmol_m2_d` | signed dO2/dt | signed dO2/dt | used directly (same convention) |
| **night loss** | `Int_O2_loss_40m_mmol_m2_d` | **positive = O2 consumed** (a *loss* magnitude) | signed dO2/dt (negative when O2 drops) | my night rate is **negated** (`-rate`) to the loss convention so it lines up with the reference |

So: **net change** needs no flip — both sides are dO2/dt. **Night** does: the
reference stores loss as a positive number, my `rate` is a signed derivative, so
the script plots/compares `-rate` for night. The CSV itself stays in the signed
dO2/dt convention for every type — do not "fix" the CSV to match the plot.

If a future request says "show respiration" (negative) instead of "loss"
(positive), flip the night sign back; if it says "use the comparison's
convention" (the current state), keep night as a positive loss.

---

## 4. Validated findings (current outputs)

Nearest-join (≤ 1 day), against the two reference series:

| Comparison | n | Pearson r | slope (mine ~ 0 + ref) |
|---|---|---|---|
| NET (net change, 0–40 m) | 387 | 0.80 | **1.015** |
| NIGHT (loss convention, 0–40 m) | 194 | 0.90 | **1.020** |
| NIGHT ÷ 2 | 194 | 0.90 | 0.510 |
| MLD-NET (net change, MLD-deepest) | 382 | 0.29 | 2.46 |

Interpretation:
- Both 0–40 m `net` and `night` agree with the reference **≈1:1** (slope ~1.0)
  with near-zero bias → the density/units conversion (§2 step 1.2) and the
  pairing/timestamps are correct.
- **There is no hidden factor-of-2 (no 0.5 "respiration coefficient").** Halving
  my night estimate drives the slope to 0.51 and *worsens* RMSE — i.e. the
  agreement is already 1:1 without any coefficient.
- **MLD-NET deliberately does NOT match the 40 m reference** (slope 2.46, r 0.29)
  — that is the point: ~77 % of pairs integrate deeper than 40 m (median 143 m),
  so the deeper layer legitimately diverges in winter. The ~23 % of summer pairs
  that floor at 40 m coincide with the 0–40 m series. **Do not "fix" this slope
  toward 1** — it is comparing a different (deeper) layer, not the reference's.

---

## 5. Gotchas and known quirks

- **`net` count (~475) is ~2× the reference's (~194).** This is by design: `net`
  computes *both* the dusk→dusk and the dawn→dawn same-phase subseries, so there
  are ~2 estimates per ~48 h cycle. The reference apparently reports ~1 per
  cycle. The nearest-join still matches well (r = 0.80). Don't "fix" this unless
  asked.
- **The reference CSVs are not smooth.** Despite a verbal claim that the
  reference is an "18-day moving average", the stored values are noisy
  point-to-point (night-loss range −478…+1104, lag-1-difference SD > series SD).
  My raw per-cycle night estimate matches them 1:1; *smoothing my series to 18
  days drops the correlation to ~0.45*, which it only would if the reference is
  per-cycle. If the reference truly is an 18-day average elsewhere, the CSVs we
  read appear to hold the **pre-smoothing per-cycle values**. Worth confirming
  which file the averaging step writes to before trusting the "18-day" framing.
- **Do not use `Data/Intermediate/Floats/argo_3902681_interp.csv`** for this — it
  has date only, no time-of-day, so it cannot separate dusk/dawn. Read the Sprof.
- **PRES ≈ depth** in the top 40 m (difference < 0.5 %). The script uses
  `gsw_z_from_p` properly, but if you ever simplify, treating PRES as depth is
  acceptable here.
- **`gsw` functions are vectorised** — pass whole profile vectors; lon/lat are
  scalars per profile (recycled).
- **QC array is character**, parsed per char to integer. If parsing ever looks
  off, log how many levels pass QC for the first few profiles (the script keeps
  flags {1,2,5,8}).
- **Do not modify** `Scripts/ncp_o2budget_float_3902681.R` (a separate, older
  script). This budget lives entirely in `Scripts/o2_budget_3902681.R`.
- **No O2→C `/2`.** The *other* (NCP/respiration→carbon) scripts in this repo do
  divide O2 by 2; this one deliberately does not. Keep units in O2 throughout.
- **The MLD budget shares the trim/QC/flux of the 0–40 m budget** but drops a few
  extra pairs (~5/475 `net`) when a cast is shallower than the common MLD — the
  rate is `NA` there by design (no fabricated O2 below the cast). See §7.
- **`integrate_o2()` requires the cast to reach `z_bot`** before integrating. If
  you ever loosen this, you reintroduce `approx(rule = 2)` extrapolation that
  invents a flat O2 column below the deepest level — exactly the failure mode the
  guard prevents.

---

## 6. If you need to change something

- **Different integration depth (0–40 m budget):** change `grid <- 0:40` and the
  two coverage guards (`max(depth) >= 40`, `sum(depth <= 60) >= 5`) consistently.
- **Different night window:** the `dt_days < 0.7` filter defines "overnight".
- **Different trim:** `n_drop <- 5` controls the head/tail trim.
- **Plot y-window:** `coord_cartesian(ylim = ...)`. Widen if a request needs the
  full productive-season excursions; tighten for legibility.
- **MLD-budget floor / criterion:** the 40 m floor is `max(d$mld[i], d$mld[j], 40)`
  in `make_pair_mld`; the MLD criterion/`default.depth` live in the `castr::mld`
  call in the Step 1 loop (§7).
- After any change, re-run and re-read the printed verification block: the
  slopes should stay ~1.0 for the 0–40 m `net` and `night` (loss convention). A
  slope that jumps to ~2 or ~0.5 means a units/sign regression — check §2 step
  1.2 and §3 first. (The **MLD-NET** slope ~2.4 is expected — §4 — don't treat it
  as a regression.)

---

## 7. The MLD-based budget (integrate to the deeper of the two MLDs)

A second, parallel budget that integrates O2 to the **mixed-layer depth instead
of the fixed 0–40 m slab**. Motivation: in Iceland Basin / Irminger winter the
MLD is hundreds of metres, so the 0–40 m slab misses the bulk of the ventilated
column and the surface air-sea flux is applied to the wrong (too thin) layer.

### What it does
- **Integration depth per pair:** `z_int = max(MLD_start, MLD_end, 40)` — the
  **deeper of the two profiles' MLDs**, never shallower than 40 m. So if both
  MLDs are < 40 m (summer), the budget falls back to the 0–40 m depth and
  coincides with the original budget; otherwise it deepens.
- **Both profiles are integrated to the same `z_int`** (via `integrate_o2()`), so
  the inventory *difference* is a like-for-like comparison over a common layer.
  `integrate_o2()` trapezoidal-integrates the stored `(depth, o2_vol)` on a 1 m
  grid ending exactly at `z_int`, and returns `NA` if the cast does not reach
  `z_int` (so the pair is dropped rather than extrapolated).
- **Same three pairings** as the 0–40 m budget (`prev`/`night`/`net`), via
  `make_pair_mld()`, and the **same air-sea flux** (`airsea_flux()`) — but here
  `F_as` correctly acts on the full mixed layer (resolves the §2b caveat).

### Sign / units
Identical to the 0–40 m budget: `rate_mmol_o2_m2_d` is signed dO2/dt (negative =
loss), `rate_corr = rate − F_as`, units mmol O2 m⁻² d⁻¹, no O2→C `/2`. The night
loss-convention flip (§3) is **not** applied to the stored CSV.

### Outputs
- **CSV `Data/Processed/o2_budget_3902681_mld.csv`** — the common-`z_int` columns
  (`z_int_m`, `mld_start_m`, `mld_end_m`, `dmld_m`, `disp_km`, `rate_mmol_o2_m2_d`,
  `rate_corr_mmol_o2_m2_d`, …) **plus** the three fix groups below (Fix 1
  volumetric, Fix 2 moving-ML + entrainment, Fix 3 σ). All three `type`s.
- **PNG `Output/o2_budget_3902681_mld_comparison.png`** (`patchwork`, 11×8) —
  two stacked panels (see §7a for the post-fix version).
- **PNG `Output/o2_budget_3902681_mld_decomposition.png`** — the moving-ML areal
  budget split into `total ML tendency = F_as + entrainment + residual` (Fix 2).
- **PNG `Output/o2_budget_3902681_amplification.png`** — areal-vs-volumetric
  per-pair box spread by season; the winter blow-up collapses in volumetric (Fix 1).
- **PNG `Output/o2_budget_3902681_mld_volumetric.png`** — the net-change budget as
  a **concentration-unit timeseries** (mmol O2 m⁻³ d⁻¹): moving-ML residual
  (green, error-weighted, with raw per-pair grey points), MLD common-z and
  0-40 m flux-corr, and `reference / 40`. The amplification-free headline view —
  winter and summer on one comparable scale (Fix 1).

### Validated behaviour (current run)
- Integration depth: median **143 m**, max **735 m**; **23 %** of `net` pairs sit
  at the 40 m floor (summer), ~5/475 dropped for insufficient cast depth.
- MLD-NET vs the 0–40 m reference: r 0.29, slope 2.46 — **expected divergence**,
  the deeper layer is a different quantity (§4). Summer (floored at 40 m) tracks
  the 0–40 m series.

---

### 7a. The MLD-depth amplification fixes (Fix 1/2/3)

The common-`z_int` differencing above has a **depth-amplification artifact**:
within a pair both casts share `z_int`, so
`rate_areal = z_int · Δ⟨O2⟩ / dt` — in deep winter (`z_int` 100–700 m) a
few-mmol-m⁻³ concentration difference (signal *or* sensor noise) is mechanically
blown up to hundreds–thousands of mmol m⁻² d⁻¹. A displacement-vs-ΔMLD
diagnostic confirmed this is **not** lateral advection (Spearman `|rate|` vs
`disp_km` ≈ 0). Three fixes (`.claude/o2_budget_3902681_mld_fixes.md`):

- **Fix 1 — volumetric columns.** `rate_vol = rate / z_int`,
  `Fas_vol = Fas / z_int`, `rate_corr_vol = (rate − Fas) / z_int`
  (mmol O2 m⁻³ d⁻¹). The `z_int` factor is divided out, so winter/summer become
  comparable. Quantified: the winter/summer per-pair **SD ratio drops from
  ≈3.9 (areal) to ≈0.8 (volumetric)** — the amplification was entirely the units.
  `Fas / MLD` is the Cornec & Fassbender (2025) Eq. 5 mixed-layer convention.
- **Fix 2 — moving mixed-layer budget + explicit entrainment.** Instead of a
  fixed control volume that changes depth between pairs, track the
  **mixed-layer-mean O2 over each cast's OWN MLD** (`o2ml_start/end`, MLD floored
  at 40 m) and model the deepening explicitly, mirroring `compute_ncp()`:
  ```
  we              = max(0, ΔMLD) / dt                       (m d⁻¹, 0 if shoaling)
  entrain_vol     = (we / MLD̄) · (O2_below − ⟨O2⟩_ML)       (mmol m⁻³ d⁻¹, Cornec Eq. 6)
  rate_resid_vol  = rate_ml_vol − Fas_ml_vol − entrain_vol  (biological residual)
  ```
  `O2_below` = mean O2 in the 20 m just below the earlier MLD. Areal versions
  (`rate_ml_`, `entrain_`, `rate_resid_mmol_o2_m2_d`) = vol × `mld_bar_m`. Sign
  check passes: deepening into lower-O2 water ⇒ `entrain < 0` (ML O2 drops).
- **Fix 3 — uncertainty + weighted smoothing.** `sigma_rate_vol = SIG_O2·√2/dt`
  (`SIG_O2 = 1 mmol m⁻³`, the adjusted-optode precision; depth cancels in
  volumetric), `sigma_rate_resid = sigma_rate_vol · mld_bar` (areal). The
  headline moving-ML residual line is smoothed with `roll_smooth_w()`
  (inverse-variance weighted, `w = 1/σ²`), so noisy deep-winter pairs don't
  dominate the seasonal mean. **No clamps / winsorizing / median filters** — the
  residual winter spread is genuine concentration variability, represented not
  hidden.

**Headline = areal** (user/PI choice 2026-06-26): the plot and analysis quantity
is `rate_resid_mmol_o2_m2_d` (moving-ML, flux + entrainment corrected, in
mmol O2 m⁻² d⁻¹); the volumetric columns are kept in the CSV and are the
amplification-free diagnostic. Areal is recoverable from any vol column as
`vol · mld_bar_m`.

### 7a outputs
- **`o2_budget_3902681_mld_comparison.png`** top panel now shows three computed
  areal series vs the reference: `0-40 m (flux-corr)` (blue, unchanged),
  `MLD common-z (flux-corr)` (orange dotted — the pre-fix amplified series) and
  `MLD moving (flux+entrain-corr)` (green, the headline residual, error-weighted
  smoothed). Bottom panel = common `z_int` (unchanged).
- The validated 0–40 m budget is **untouched**: NET slope **1.015**, NIGHT slope
  **1.020** after the fixes (§4).

### Gotchas specific to the MLD budget
- **`o2_prof_list` is indexed by the original `prof_index`**, not the row number
  in the filtered/trimmed `prof` table. `make_pair_mld` looks up
  `d$prof_index[i]` — if you ever rebuild `prof`, keep the `prof_index` column.
- **`castr` masks `stats::integrate` and `stats::smooth`** on attach — harmless
  here, but note it if a later edit relies on those base functions.
- The **`max(): no non-missing arguments` warning** during the Step 1 loop is
  pre-existing (the 0–40 m coverage guard on an empty profile) and unrelated to
  the MLD code.

---

## 8. Nighttime respiration index R (Step 5g)

A community **respiration rate** `R` derived from the **clean dark window only**.
User/PI choice (2026-06-26): window = nighttime `dusk → dawn`; index = the
absolute corrected respiration rate (areal + volumetric).

### Why the dark window (and only it)
This float samples a clean **~11 h dark interval** (`dusk → dawn`, `dt ≈ 0.47 d`,
234 pairs) but has **no clean daytime window** — `dawn → dusk` is **~37 h** and
spans day+night+day (verified: `dawn→dusk` median `dt = 1.53 d`, zero pairs
< 0.7 d). In the dark photosynthesis is zero, so the moving mixed-layer O2 budget
reduces to `d⟨O2⟩/dt = −R + F_as/MLD + entrainment`, hence the biological
residual already computed by `make_pair_mld` **is** `−R`:
```
R = −rate_resid = −(rate − F_as − entrainment)     (positive = O2 consumed)
```
Reported areal (`R_mmol_o2_m2_d`) and volumetric (`R_vol_mmol_o2_m3_d`, the
amplification-free version). σ carries over from Fix 3; the seasonal line is the
inverse-variance-weighted 8-pt mean.

### Outputs
- **CSV `Data/Processed/o2_budget_3902681_respiration.csv`** (234 rows): `mtime`,
  interval, `mld_bar_m`, `we_m_d`, ML-mean O2, `Fas_`, `entrain_`, `R_mmol_o2_m2_d`,
  `R_vol_mmol_o2_m3_d`, and the σ columns.
- **PNG `Output/o2_budget_3902681_respiration.png`** — two panels: areal R (top)
  and volumetric R (bottom), grey raw per-pair points + red weighted seasonal line.

### Findings and limitations (read before using R)
- **Entrainment is negligible over 11 h** (median correction ≈ 0.00 mmol m⁻² d⁻¹):
  the MLD barely moves overnight, so "corrected for entrainment" is *included but
  immaterial at the night timescale* — the night correction is effectively
  air-sea only. (Entrainment matters for the multi-day net/seasonal budget, §7,
  not here.)
- **Seasonal structure is sensible: R is positive in the productive season**
  (summer ~ +2…+6 mmol O2 m⁻³ d⁻¹ volumetric — more biomass, more respiration)
  and **drops to ~0 / slightly negative in winter**.
- **Winter R is at the method's detection limit and not trustworthy.** Only ~45 %
  of night pairs give `R > 0`; the volumetric median is ≈ −0.5. Over an 11 h dark
  window the genuine respiration signal (~0.5 mmol m⁻³ d⁻¹) is **below the per-pair
  noise** (`σ_R_vol = SIG_O2·√2/dt ≈ 3 mmol m⁻³ d⁻¹`), and in winter the air-sea
  ingassing term dominates the budget, so a small over/under-estimate of `F_as`
  (large in high winds, deep ML) pushes the corrected residual negative. The
  *areal* winter values are additionally amplification-blown (range −12238 … +2876),
  so **use the volumetric R, and trust the summer/seasonal weighted line, not
  individual winter pairs.** This mirrors the whole analysis: winter is
  noise/flux-dominated, summer carries the biological signal.
