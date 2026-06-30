# Handoff: fixing the MLD-depth "multiplication" amplification in the O2 budget

This is a **self-contained implementation spec** for the next agent. You do not
need prior conversation context. It explains *why* the MLD-based O2 budget
produces huge, noisy winter rates, and *how* to fix it. All paths are relative to
the repo root `C:\Users\petit\Documents\Producitivity_estimates`; run everything
**from the repo root** with **R 4.5.2**:

```powershell
& "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/o2_budget_3902681.R
```

Read alongside:
- `.claude/o2_budget_3902681_explained.md` — how the whole script works (the
  MLD-based budget is its **§7**; the 0–40 m budget is the rest). **Do not break
  the validated 0–40 m budget** — its verification slopes must stay ~1.0 (see that
  doc's §4).
- `Scripts/o2_budget_3902681.R` — the only file to edit. The MLD budget lives in
  `make_pair_mld()`, the helper `integrate_o2()`, and the assembly block
  "Step 2c". The MLD outputs are `Data/Processed/o2_budget_3902681_mld.csv` and
  `Output/o2_budget_3902681_mld_comparison.png`.
- `Scripts/ncp_function.R` — `compute_ncp()` already implements an entrainment
  term for **nitrate** the same way Fix 2 needs it for O2. Mirror it.

---

## 0. TL;DR — what to implement, in order

1. **Fix 1 (do first, small, robust): volumetric units.** Add concentration-rate
   columns (`… / z_int`) to the MLD budget and make the comparison plot use them.
   This *is* the direct cure for the amplification. §3.
2. **Fix 2 (physical completion, larger): explicit entrainment term**, so the
   mixed-layer deepening is modelled instead of left in the residual. Compose with
   Fix 1. §4.
3. **Fix 3 (honesty): propagate per-pair uncertainty (∝ `z_int`) and
   error-weight the smoothing**, and/or compute winter budgets over multi-day
   windows. §5.
4. **Do NOT** apply cosmetic clamps / winsorizing / hard y-limits / median
   filters to hide the spikes. §6.

There is one **decision for the user/PI** before finalizing (areal vs volumetric
as the headline series) — see §7. If unconfirmed, default to: **keep areal in the
CSV, make volumetric the headline of the plot and the analysis quantity.**

---

## 1. What the MLD budget currently does (the starting point)

`make_pair_mld(d, i, j, type)` (in `Scripts/o2_budget_3902681.R`) takes an ordered
pair of profiles and:
- sets a common integration depth `z_int = max(d$mld[i], d$mld[j], 40)` (the
  deeper of the two mixed-layer depths, floored at 40 m);
- integrates **both** casts to that common `z_int` via `integrate_o2()` →
  `inv_i`, `inv_j` (mmol O2 m⁻²);
- `rate = (inv_j - inv_i) / dt_days` (signed dO2/dt, mmol O2 m⁻² d⁻¹);
- subtracts the Wanninkhof-2014 air-sea flux `fx$Fas` (from the shared
  `airsea_flux()` helper) → `rate_corr_mmol_o2_m2_d = rate - fx$Fas`.

It already stores `z_int_m`, `mld_start_m`, `mld_end_m`, `dmld_m`
(= `mld[j]-mld[i]`), and `disp_km` (great-circle inter-profile displacement). The
full depth-ordered O2 profile per cast is available in the global list
`o2_prof_list[[prof_index]]` as `list(depth = …, o2_vol = …)` (mmol O2 m⁻³) — you
will need this for Fix 2.

Median `z_int` ≈ 143 m, max ≈ 735 m; ~23 % of pairs sit at the 40 m floor
(summer), ~5/475 net pairs are dropped because a cast is shallower than `z_int`.

---

## 2. The problem — depth amplification (this is a units artifact, not a bug)

Within a pair both profiles use the **same** `z_int`, so the inventory difference
factorizes exactly:

```
rate_areal = z_int · (⟨O2⟩_j − ⟨O2⟩_i) / dt          [mmol m⁻² d⁻¹]
```

where `⟨O2⟩_k` is the **layer-mean** O2 concentration over 0…`z_int` for profile
`k`. All the biological/physical information is in `Δ⟨O2⟩` — a few mmol m⁻³. In
winter `z_int` is 100–700 m, so a 1 mmol m⁻³ mean difference (real signal, sensor
noise, or two casts ~11–37 h apart sampling slightly different water) becomes
~`z_int` mmol m⁻² → ~hundreds–thousands of mmol m⁻² d⁻¹. Winter just *has* the
big `z_int`, so winter gets the big, noisy numbers. **Nothing is computed wrong —
the multiplication is mechanical.**

A pairwise diagnostic already confirmed this is not lateral advection: Spearman
correlation of `|rate_corr|` vs `disp_km` ≈ 0 (even slightly negative in winter),
vs `|dmld_m|` only +0.12 within winter (+0.30 all-season, mostly seasonal
covariation). See `Output/o2_budget_3902681_disp_vs_dmld.png` and the
"Advection vs entrainment diagnostic" printout in the script. So: **don't chase a
horizontal-advection correction; fix the depth scaling and model entrainment.**

---

## 3. Fix 1 — do the budget in concentration (volumetric) units

**Idea:** stop carrying the `z_int` factor into the quantity you analyze and plot.
Report the **volumetric tendency**:

```
rate_vol      = rate_areal / z_int                   [mmol m⁻³ d⁻¹]
Fas_vol       = Fas / z_int                           (surface flux diluted over the layer)
rate_corr_vol = (rate_areal − Fas) / z_int = rate_corr_areal / z_int
```

This removes the `z_int` scaling entirely, so winter and summer become directly
comparable and the spurious winter excursions collapse to their genuine
concentration-level size. It also makes the gas term **consistent**: a surface
flux diluted over the mixed layer is `Fas / MLD` — exactly Cornec & Fassbender
(2025) Eq. 5 (`d[DIC]_Gas/dt = (k·…)/MLD`), the standard mixed-layer-budget
convention (everything in mmol m⁻³ d⁻¹; multiply by MLD only at the very end if an
areal NCP number is wanted).

### Implementation
In `make_pair_mld()`, after `fx <- airsea_flux(d, i, j)`, add to the returned
tibble (keep the existing areal columns — they're still useful and feed the 40 m
verification framing):

```r
    rate_vol_mmol_o2_m3_d      = rate / z_int,
    Fas_vol_mmol_o2_m3_d       = fx$Fas / z_int,
    rate_corr_vol_mmol_o2_m3_d = (rate - fx$Fas) / z_int,
```

Add those three names to the `select(...)` for `budget_mld_out` (the
`o2_budget_3902681_mld.csv` writer).

### Plot (`o2_budget_3902681_mld_comparison.png`)
Switch the top (rate) panel to **concentration units**. Two series:
`MLD-deepest (vol, flux-corr)` = `rate_corr_vol`, smoothed; and the reference.
**The reference (`Int_dO2dt_40m_mmol_m2_d`) is a 0–40 m areal rate — its
volumetric equivalent is `reference / 40`.** So plot `ref/40` to compare in
mmol m⁻³ d⁻¹. Update the y-axis label to `mmol O₂ m⁻³ d⁻¹` and widen/retune
`coord_cartesian(ylim = …)` (the volumetric values are ~1–2 orders of magnitude
smaller than the areal ones). Keep the bottom integration-depth panel as is.

### Expected outcome
Winter variance in the plotted series drops dramatically; the residual winter
noise is the genuine concentration variability (handled by Fix 3). The 0–40 m
budget and its verification block are untouched.

---

## 4. Fix 2 — explicit entrainment term (physical completion)

Fix 1 removes the mechanical scaling but the mixed-layer deepening is a **real**
process; right now it is bundled into the residual rate. Decompose the inventory
change with the product rule:

```
d(z·⟨O2⟩)/dt = z̄ · d⟨O2⟩/dt   +   ⟨O2⟩̄ · dz/dt
                 └ tendency you want   └ entrainment / deepening of low-O2 water
```

and add an explicit **entrainment tendency** (concentration units), mirroring
`compute_ncp()` in `Scripts/ncp_function.R` (which already does this for nitrate):

```
we              = max(0, ΔMLD) / dt                 (deepening rate; 0 if shoaling)
entrainment_vol = (we / MLD) · (O2_below − ⟨O2⟩_ML) [mmol m⁻³ d⁻¹]
```

where `⟨O2⟩_ML` is the mixed-layer mean O2 and `O2_below` is the mean O2 just
below the (earlier) MLD — both obtainable from `o2_prof_list[[prof_index]]`
(integrate/average over the relevant depth ranges, same pattern as
`integrate_o2()`). This is Cornec & Fassbender Eq. 6 in concentration form. The
biological/uncorrected residual then becomes, in concentration units:

```
rate_resid_vol = rate_vol − Fas_vol − entrainment_vol   [− diffusion_vol, optional]
```

Optionally also add diapycnal diffusion across the MLD base (Cornec Eq. 8,
`Kz·∂O2/∂z / MLD`, `Kz ≈ 1e-5 m² s⁻¹`) — lower priority; document if added.

### Conceptual choice to be aware of (flag in code comments)
The current code integrates both casts to a **common `z_int` and differences** —
that is a *fixed control-volume* budget, in which entrainment is internal
redistribution and partly cancels, but the **control volume changes between
consecutive pairs** (each pair has its own `z_int`), so the rates do not string
into one consistent inventory budget. The standard, cleaner alternative is a
**moving mixed-layer budget**: track `⟨O2⟩_ML(t)` over each cast's *own* MLD and
apply the entrainment term above. **Recommended target: the moving mixed-layer,
concentration-unit budget** (matches Cornec and `ncp_function.R`). If you keep the
common-`z_int` formulation, state in a comment that entrainment is then a
bottom-boundary flux term, not the Eq. 6 form, and don't double-count.

---

## 5. Fix 3 — uncertainty propagation + weighted smoothing (honesty, not hiding)

The amplification also amplifies noise, so represent it rather than suppress it:

- **Per-pair uncertainty:** `σ_rate_areal ≈ z_int · σ_⟨O2⟩ · √2 / dt`, i.e.
  `σ_rate_vol ≈ σ_⟨O2⟩ · √2 / dt` (the `z_int` cancels in volumetric — another
  reason to prefer Fix 1). Estimate `σ_⟨O2⟩` from the DOXY sensor precision and/or
  the within-layer scatter. Add a `sigma_rate_*` column.
- **Error-weight the smoothing:** replace the plain centered rolling mean
  (`roll_smooth()` in the script) with an inverse-variance-weighted rolling mean
  so noisy deep-winter pairs don't dominate the seasonal line.
- **Multi-day windows in winter (optional):** difference `⟨O2⟩` over ~10–15 day
  windows (the repo's standard NCP window is 15 d; see CLAUDE.md) instead of
  pair-to-pair, so you difference a *trend*, not high-frequency pair noise.

---

## 6. What NOT to do

- **No cosmetic clamps.** Do not winsorize, hard-clip, or median-filter the areal
  rate to make the plot look calm. That biases the seasonal mean and hides the
  physics. (The `NCP/2` clamp heuristic in `ncp_function.R` is exactly the kind of
  band-aid to avoid here — do not copy it.)
- **Do not "fix" the MLD-NET-vs-40 m-reference slope toward 1.** They are
  different layers; divergence is expected (explainer §4).
- **Do not touch** `make_pair()` / the 0–40 m budget logic, the QC parsing, the
  trim, or `Scripts/ncp_o2budget_float_3902681.R` (a separate older script).
- **No O2→C `/2` conversion** anywhere — this budget stays in O2 (explainer §5).

---

## 7. Decision point for the user/PI (ask before finalizing)

**Should the headline MLD series be areal (mmol m⁻² d⁻¹) or volumetric
(mmol m⁻³ d⁻¹)?**
- *Volumetric* → un-amplified concentration tendency, comparable across seasons,
  matches the literature; loses the areal-NCP magnitude.
- *Areal* → keeps the column-integrated magnitude (what you'd report as NCP) but
  carries the amplification.

Recommended default if unconfirmed: **CSV keeps both** (areal + volumetric);
**plot and analysis lead with volumetric**; areal recoverable as
`rate_vol · MLD` at the end.

---

## 8. Verification after implementing

Re-run the script and check the printed block + outputs:
- **0–40 m budget unchanged:** `report_fit` slopes still ~1.015 (NET) / ~1.020
  (NIGHT). If these move, you touched the wrong code path — revert and isolate.
- **Amplification gone:** the plotted volumetric winter series should be ~1–2
  orders of magnitude smaller and visibly less spiky than the areal one. Print
  `sd(rate_corr)` vs `sd(rate_corr_vol)` for winter (Nov–Mar) vs summer to
  quantify; the winter/summer SD ratio should drop sharply in volumetric.
- **Entrainment sign sanity:** deepening (`dmld_m > 0`) with lower O2 below should
  give a **negative** O2 entrainment tendency (ML O2 decreases). Print a sign
  check like the existing air-sea `delta>0 ⇒ Fas<0` check.
- **CSV columns present:** `rate_vol_mmol_o2_m3_d`, `Fas_vol_mmol_o2_m3_d`,
  `rate_corr_vol_mmol_o2_m3_d` (+ entrainment / sigma columns if added) in
  `o2_budget_3902681_mld.csv`.
- Update `.claude/o2_budget_3902681_explained.md` §7 (and §4 if numbers change) to
  document whatever you implement.

---

## 9. File / symbol map

| Thing | Where |
|---|---|
| Script (only file to edit) | `Scripts/o2_budget_3902681.R` |
| MLD pair budget | `make_pair_mld()` |
| Profile integrator (0→z) | `integrate_o2(pi, z_bot)` |
| Stored full O2 profiles | global `o2_prof_list[[prof_index]]` → `$depth`, `$o2_vol` (mmol m⁻³) |
| Air-sea flux helper | `airsea_flux(d, i, j)` → `list(k, u2, o2_eq, dO2, Fas)` |
| MLD per profile | `prof$mld` (castr `mld`, sigma0 0.03, default 300 m) |
| Smoother | `roll_smooth(x, k = 8)` |
| MLD CSV writer | `budget_mld_out` `select(...)` → `out_csv_mld` |
| MLD plot | "Step 5b" block → `out_png_mld` |
| Reference series | `ref_net` (`Int_dO2dt_40m_mmol_m2_d`, **0–40 m areal**; ÷40 for volumetric) |
| Entrainment template (nitrate) | `compute_ncp()` in `Scripts/ncp_function.R` |
