# Figure Review — Data-Processing & Plotting Guidelines

**Purpose:** Actionable guidance for the coding agent that owns the data-processing and
plotting code for the five main-text figures. The scientific message is **unchanged**:
GOP ≥ NPP ≥ NCP, seasonal decoupling of NCP and NPP, and an early e-ratio peak. These
notes make that message land more cleanly and defensibly — they do not alter the findings.

**Last updated:** 2026-07-01

---

## Cross-cutting issues (fix these first — they touch several figures)

1. **Harmonize smoothing windows.** GOP is an 18-day running mean, NPP 10-day, NCP 30-day
   (Fig 2), but 20-day in Fig 3, and other windows were tested per the handoff. When three
   curves are overlaid to argue an *ordering* and a *phase lag*, unequal smoothing is a
   confound — a lag can be partly an artefact of different window widths. Pick one window
   (or a small justified set) applied consistently, exposed as a single config parameter.
   If windows must differ (different native sampling), state why in the caption and verify
   the phase lag survives equal smoothing.

2. **Reconcile temporal coverage.** NCP spans 2015–2025 (11 yr, Fig 3); NPP spans 2015–2023
   (9 yr, Fig 4). Fig 5's e-ratio therefore divides two climatologies built from *different
   years*. Recompute the e-ratio on the **common year set** and label coverage consistently
   in every figure and caption.

3. **Adopt one winter-masking convention.** Winter is unreliable for two independent reasons
   (GOP diel method fails under deep mixing / low light; CbPM has no retrievals). They are
   currently handled differently (GOP left as raw noise; NPP splined to zero). Define a
   single `valid_season` mask per product, propagate it through every figure, and render
   masked periods identically (e.g. greyed band + no line, or dashed).

4. **Consistent uncertainty semantics.** Bands currently mix meanings (IQR, ±1 SD spatial,
   ±1 SD Monte Carlo, p10–p90 inter-annual). Keep whichever is scientifically correct per
   panel, but make the label wording explicit and the visual style consistent so readers
   can compare envelopes across figures.

5. **No silent axis clipping.** Assert that plotted data ⊆ axis limits (see Figs 2 and 3).

---

## Figure 1 — trajectory & sampling geometry
`fig1_float_3902681_trajectory.png`

Strong and clear. Small processing/render fixes:

- **Deploy marker colour mismatch:** the deployment is a *green* star, but on the date
  colourbar Oct 2024 maps to *dark purple*. Colour the start/end markers to match the
  colormap endpoints, or switch them to a neutral symbol so no date is inferred from the
  star's fill.
- **Overlapping annotations:** the "deploy 2024-10-02" / "last 2026-02-10" labels sit on top
  of the dense trajectory tangle near the Iceland shelf. Offset them with leader lines.
- **Schematic vs. text ordering:** the bottom-right schematic draws dawn→dusk→dawn, but the
  method brackets a *single night* (dusk, then pre-dawn ~11 h later). The 0–60 h schematic
  implies ~30 h dawn-to-dusk spacing, reading as inconsistent with "night-bracketing."
  Redraw so the dusk→pre-dawn pair that isolates nighttime O₂ loss is the obvious unit.

---

## Figure 2 — GOP/NPP/NCP synthesis (needs the most work)
`fig2_synthesis_GOP_NPP_NCP_areal_volumetric.png`

Message (three-way ordering + P1/P2 decoupling) is good, but panel (a) undercuts it:

- **Winter GOP is pure noise and visually dominates.** From Nov–Mar the red GOP line spikes
  between roughly +300 and −150 and the IQR band saturates the whole panel background. In
  that window GOP frequently drops *below* NPP/NCP, which visually **contradicts the caption's
  claim that GOP ≥ NPP ≥ NCP "at all times of year."** Apply the winter mask (cross-cutting
  item 3): suppress or grey out GOP where the diel method is unreliable, and restrict the
  "ceiling" claim to the productive season where it actually holds.
- **Data exceeds the axis.** The y-axis clips at 300 while the GOP IQR runs off the top. Set
  limits that contain the plotted data, or explicitly note clipping.
- **The IQR band is too wide to be informative.** An 18-day interquartile range that fills the
  panel tells the reader little. Consider a robust central estimate with a tighter interval
  (bootstrap CI of the mean, or median ± MAD), or a longer smoothing window — and align the
  band semantics with Figs 3–5.
- **Legend sits in the noisiest corner.** Move it to a clear region or add a solid background box.

---

## Figure 3 — basin NCP (decadal + synthetic)
`fig3_ncp_timeseries_publication.png`

Cleanest figure; message is crisp. Minor:

- Panel (b) grey individual-year lines are **clipped at the top** (spring peaks ~55 and a Dec
  excursion run past the axis). Raise `ylim` to contain them, or note the clip.
- Panel (a) markers are very dense — thinning them (or line-only) reduces ink without losing
  the "repeatable spring excursion" point.
- Confirm the bin width shown here (20-day) is the one standardized on per cross-cutting item 1.

---

## Figure 4 — basin NPP (decadal + synthetic)
`fig4_npp_timeseries_publication.png`

Reproducible and clear. Watch:

- **Winter spline forced to exactly 0.** The dashed interpolation drives NPP to a hard zero each
  winter. This almost certainly biases the annual NPP integral low and is physically
  questionable. Decide explicitly: a small non-zero winter floor, or clearly label the winter
  segment as "not estimated" rather than "= 0." Any annual-total NPP (~225 g C m⁻² yr⁻¹) must
  document how winter was treated.
- Coverage label (9 yr / 2015–2023) must be reconciled with NCP per cross-cutting item 2.

---

## Figure 5 — the e-ratio money figure
`fig5_eratio_synthetic_year.png`

Message lands, but the dual axis needs tightening:

- **Two competing zero lines.** There is a grey dotted line at production = 0 *and* another at
  e-ratio = 0 (which falls at production ≈ 35 on the left scale). A reader can misread which
  zero applies to which curve. Draw e-ratio gridlines/ticks only on the right axis, drop the
  second dotted line, or colour-code each reference line to its axis.
- **Verify the left↔right axis mapping in code.** The linear transform placing e-ratio
  1.0/0.0/−1.0 against the production scale should be a single documented function, not
  hand-tuned, so the horizontal reference lines are provably correct.
- **Recompute on common years** (cross-cutting item 2) — this is where the NCP/NPP year mismatch
  propagates into a quoted number.
- The e-ratio blows up wherever NPP is small; restricting to Feb–Oct is correct — make that
  mask the same `valid_season` object used elsewhere.

---

## Checklist for the coding agent

- [ ] One shared smoothing-window parameter; document any deviation.
- [ ] One `valid_season` mask per product, propagated to every figure and rendered identically.
- [ ] Common-year climatology for NCP and NPP before any ratio.
- [ ] No silent axis clipping — assert plotted data ⊆ axis limits.
- [ ] Explicit, documented winter treatment for NPP (no hard zero without justification).
- [ ] Uncertainty-band semantics consistent across figures (same statistic, same label wording).
- [ ] Fig 1 start/end markers coloured to match (or neutral to) the date colormap.
- [ ] Fig 2 GOP winter masked; IQR band tightened; y-limits contain data.
- [ ] Fig 5 dual-axis mapping factored into one documented function; single zero reference per axis.
