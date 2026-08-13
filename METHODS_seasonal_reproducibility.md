# Seasonal reproducibility of NPP and NCP — method note

Written to support the Materials & Methods section. Describes what
`figure_config.per_year_climatology_regression()` computes and why, as applied in
`Scripts/ncp_timeseries_figure.py` (Fig. 3c) and `Scripts/npp_timeseries_figure.py`
(Fig. 4c).

Implementation: `Scripts/figure_config.py` → `per_year_climatology_regression`.
Per-year outputs: `Output/npp_seasonal_reproducibility.csv`,
`Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_seasonal_reproducibility.csv`.

---

## 1. The question

Given a multi-year record of a seasonally-varying quantity, how closely does any
individual year follow the mean seasonal cycle? Two distinct ways a year can
depart from the norm are separated:

- a change in **shape or phase** (bloom earlier/later, different rise and decay
  form) → measured by R²;
- a change in **amplitude** (same shape, scaled up or down) → measured by the
  regression slope.

An amplitude-only anomaly leaves R² at 1 and moves the slope; a timing anomaly
drops R² while leaving the slope near 1. Reporting both separates them.

## 2. Climatology

The seasonal climatology is the inter-annual mean by calendar month, placed at
mid-month day-of-year and interpolated to a weekly grid with a cubic spline.

- **NCP**: 12 monthly means from the 20-day-bin record, interpolated with a
  **periodic** cubic spline so December wraps continuously into January
  (appropriate for a signed quantity that is defined year-round).
- **NPP**: 9 monthly means (February–October; the CbPM product has no output in
  November–January, when light is insufficient), interpolated with a natural
  cubic spline. The unobserved winter tail is not used in the analysis below.

Inter-annual spread is summarised by the 10th–90th percentile across years,
computed the same way.

## 3. Leave-one-out predictor

**The climatology used to score a given year is rebuilt from all *other* years.**

With only 9 (NPP) to 11 (NCP) years, each year contributes ~10% of the mean, so
scoring a year against a climatology that includes it makes each year partly
predict itself and biases every score upward — most strongly for the years that
are most anomalous, which is exactly where the bias does the most damage. The
climatology *plotted* in Fig. 3b / 4b is the conventional all-years one; only
the regression predictor is leave-one-out. This should be stated explicitly, as
the two differ slightly.

## 4. Regression model

For year *Y*, using that year's own observations at their native sampling times
*t* (20-day bins for NCP, monthly for NPP — no re-gridding or interpolation of
the observations):

```
obs_Y(t) = a · clim_(−Y)(doy(t)) + b + ε
```

fitted by least squares, where `clim_(−Y)` is the leave-one-out climatology
evaluated at the same day-of-year. Reported quantities:

| symbol | name | meaning |
|---|---|---|
| `R²` | coefficient of determination of the fit | fraction of that year's variance explained by the climatological **shape and phase**; invariant to amplitude and offset |
| `a` | slope | **amplitude anomaly**; a > 1 = seasonal cycle stronger than the norm |
| `b` | intercept | constant offset (units of the variable) |
| `NSE` | Nash–Sutcliffe efficiency | skill of the **raw** climatology as a direct predictor, i.e. against the 1:1 line with no fitted a, b. Always ≤ R²; negative means the climatology predicts that year worse than that year's own mean |
| `RMSE`, `bias` | | deviation from the raw climatology, in data units |

R² and NSE together are informative: R² high with NSE much lower means the year
had the right *shape* but the wrong *amplitude or offset*.

## 5. Weighting — and why NCP and NPP are treated differently

### NCP: inverse-variance weighted (w = 1/σ²)

The nitrate-budget NCP carries a per-bin Monte Carlo uncertainty σ (CANYON-B
prediction CI + profile bootstrap, 200 iterations). **σ is strongly seasonal**
(median across 2015–2025, mmol C m⁻² d⁻¹):

| Jun–Sep | Oct | Nov–May |
|---|---|---|
| 1.7 – 2.7 | 5.0 | 10.2 – 15.7 (individual bins to 31) |

In winter σ is as large as the signal itself. Unweighted, the score is then
dominated by bins that are not significantly different from zero — several are
physically impossible, e.g. 2023-01-20 = **+32.4 ± 21.8** mmol C m⁻² d⁻¹, a
positive January NCP under negligible light.

Weighting the fit and all statistics by 1/σ² changes the result materially:

| | 2022 | 2023 | mean, 11 yr |
|---|---|---|---|
| unweighted | 0.49 | **0.21** | 0.727 |
| 1/σ² weighted | 0.85 | **0.88** | **0.894** |

Unweighted, 2022 and 2023 appear to be strongly anomalous years. Weighted, both
are unremarkable (2023 slope moves from 0.49 to 1.12). **The apparent anomaly is
a winter-precision artefact, not an anomalous seasonal cycle** — this is the
substantive reason the weighted version is used as the primary result. The
unweighted values are retained in the CSV and plotted as black dashes in
Fig. 3c as a transparency/sensitivity check rather than being dropped.

Implementation detail: σ is floored at 5% of each year's median σ before
inverting, so that a single near-zero uncertainty cannot dominate the fit.
Weights are normalised to sum to 1 (irrelevant to the fit, but it makes the
weighted RMSE and bias directly comparable to their unweighted counterparts).

### NPP: deliberately unweighted

`npp_int_std` in the CMEMS/CbPM product is the **spatial** standard deviation of
NPP across the domain's grid cells within a month, *not* an uncertainty on the
monthly mean. Weighting by it would down-weight spatially heterogeneous months
rather than poorly-constrained ones — a different and unjustified operation.
The NPP record also has no analogous precision problem to correct: its winter
months are simply **absent** (`n_cells = 0`), not noisy.

**Consequence for the paper:** the NPP and NCP R² are not strictly like-for-like
and a direct numerical comparison should be qualified.

## 6. Results

**NPP** (2015–2023, 9 yr, Feb–Oct months, unweighted):
mean R² = **0.878**, range 0.834–0.965; slopes 0.81–1.06.
Every year is within ±20% of the climatological amplitude, and no year departs
materially in shape or phase.

**NCP** (2015–2025, 11 yr, 20-day bins, 1/σ² weighted):
mean R² = **0.894**, range 0.809–0.945; slopes 0.81–1.13.

So on the comparable measure both records are reproducible to a similar degree,
with no outlier year in either. Note that the *unweighted* NCP mean (0.727, with
2023 at 0.21) would support the opposite claim — that NCP is far less
reproducible than NPP. That contrast is a statement about winter measurement
precision in the nitrate budget, not about interannual variability in
productivity, and should not be presented as the latter.

## 7. Sensitivity: excluding anomalous years

Removing 2022 and 2023 from the NCP climatology was considered and **rejected**:

1. It is circular — defining the norm by excluding the years that disagree with
   it, then reporting that the remainder agrees.
2. The effect is negligible: dropping both years shifts the climatology by at
   most 4.25 mmol C m⁻² d⁻¹ against a seasonal amplitude of 53.6 (< 8%,
   concentrated in November–December).
3. It is unnecessary once σ is used: the two years are not anomalous under the
   weighted metric (§5).

## 8. Suggested Methods text

> Seasonal reproducibility was quantified by regressing each year's
> observations on the seasonal climatology evaluated at the same day-of-year,
> obs_Y(t) = a·clim(doy(t)) + b. To keep the reference independent of the year
> being scored, the climatology was recomputed for each year with that year
> excluded (leave-one-out); with 9–11 years, retaining the target year
> otherwise inflates the agreement. The coefficient of determination R² measures
> the fraction of the year's variance explained by the climatological seasonal
> shape and phase (invariant to amplitude and offset), while the slope a
> measures the amplitude anomaly relative to the climatological cycle.
>
> For NCP the regression and all associated statistics were weighted by the
> inverse variance of the Monte Carlo uncertainty (w = 1/σ²). This is necessary
> because σ is strongly seasonal — median 1.7–2.7 mmol C m⁻² d⁻¹ in June–
> September but 10–16 in November–May, comparable to the signal — so that
> unweighted scores are dominated by winter bins that are not significantly
> different from zero. NPP was left unweighted because the available spread
> (`npp_int_std`) is a spatial standard deviation across grid cells rather than
> an uncertainty on the monthly mean; the NPP and NCP R² are therefore not
> strictly equivalent measures.
>
> Over 2015–2023 the NPP seasonal cycle was highly reproducible (mean R² = 0.88,
> range 0.83–0.97; slopes 0.81–1.06). Over 2015–2025 the NCP seasonal cycle was
> reproducible to a comparable degree (weighted mean R² = 0.89, range 0.81–0.95;
> slopes 0.81–1.13). Unweighted, the NCP mean falls to 0.73 with two apparent
> outlier years (2022, 2023); this difference is attributable to winter
> measurement precision rather than to interannual variability in production.
