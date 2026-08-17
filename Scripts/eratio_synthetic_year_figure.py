"""
Synthetic-year NPP, NCP and e-ratio (NCP:NPP) — Iceland Basin & Irminger Sea.

The NPP and NCP curves on the primary axis are the exact synthetic-year
builders from npp_timeseries_figure.py and ncp_timeseries_figure.py, so the
climatologies are identical to the ones in the standalone NPP and NCP figures:

  NPP : CbPM synthetic year (Feb-Oct observed, Nov-Jan spline-interpolated),
        converted mg C -> mmol C m-2 d-1 (/12) to share NCP's units.
  NCP : nitrate-drawdown synthetic year (full 12-month periodic climatology).

The e-ratio is NOT the ratio of those two climatologies. It is computed
year by year on the full timeseries and only then averaged:

    e_y(doy) = NCP_y(doy) / NPP_y(doy)      for each year y = 2015..2023
    e(doy)   = mean_y e_y(doy)   +/-  SD_y e_y(doy)

  * NCP_y comes from a cubic spline through the *whole* 20-day-bin record
    (not per-calendar-year splines), so winter values are interpolated across
    the year boundary from real neighbouring bins.
  * NPP_y is that year's own Fig-4 spline through its 9 observed monthly means
    (Feb-Oct), zero-anchored at DOY 1/365 — i.e. Feb-Oct is observation-driven
    and Nov-Jan is that year's winter interpolation.

Averaging the ratios rather than ratioing the averages is what keeps the
inter-annual uncertainty: the +/-1 SD ribbon on the e-ratio is the year-to-year
spread of the ratio itself, and the Monte-Carlo NCP uncertainty is propagated
alongside it (reported in the CSV and in the seasonal summary).

Only years present in BOTH records are used (2015-2023; the CbPM/CMEMS NPP
product stops in 2023 — see figure_config.COMMON_YEAR_*). Days where a year's
NPP falls below NPP_FLOOR_FRAC of that year's seasonal peak are dropped for
that year, and a DOY is only plotted where at least MIN_YEARS years survive.
That floor is what sets the drawn span (Feb-Oct at the default 0.15): the
point-wise ratio of per-year curves genuinely diverges towards the winter,
where NPP -> 0 while NCP stays negative, so the deep-winter e-ratio is not
estimable from these data. The unfloored per-year ratios are still written to
eratio_per_year.csv.

The headline number is the season-integrated ratio, sum(NCP)/sum(NPP) over the
observed Feb-Oct window, computed per year and then averaged +/- 1 SD: being an
integral it is insensitive to the near-zero-NPP shoulder days.

Input:  Output/cmems_npp_timeseries_domain_mean.csv
        Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_uncertainty.xlsx
Output: Output/IcelandIrminger_2015_2025/eratio_synthetic_year.png
        Output/IcelandIrminger_2015_2025/eratio_doy_stats.csv
        Output/IcelandIrminger_2015_2025/eratio_per_year.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parent))
import npp_timeseries_figure as npp_mod   # noqa: E402
import ncp_timeseries_figure as ncp_mod   # noqa: E402
from figure_config import (               # noqa: E402
    COMMON_YEAR_START, COMMON_YEAR_END,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR   = REPO_ROOT / "Output" / "IcelandIrminger_2015_2025"
OUT_PNG   = OUT_DIR / "eratio_synthetic_year.png"
OUT_DOY   = OUT_DIR / "eratio_doy_stats.csv"
OUT_YEAR  = OUT_DIR / "eratio_per_year.csv"

MGC_TO_MMOLC = 1.0 / 12.0   # matches convention used elsewhere in this repo

GREEN        = npp_mod.GREEN
GREEN_LIGHT  = npp_mod.GREEN_LIGHT
BLUE         = ncp_mod.BLUE
BLUE_LIGHT   = ncp_mod.BLUE_LIGHT
ERATIO_COLOR = "#7a2048"
ERATIO_LIGHT = "#d9b6c4"
MONTH_LABELS = npp_mod.MONTH_LABELS
OBS_MONTHS   = set(range(2, 11))   # Feb-Oct: real CbPM NPP; else winter spline
# Where NPP -> 0 the point-wise ratio NCP:NPP diverges (each year's spline is
# anchored to 0 at DOY 1/365), and — unlike the old ratio-of-climatologies —
# averaging per-year ratios lets a single small-NPP year blow the ensemble up.
# A year-DOY therefore only enters the ensemble where that year's NPP exceeds
# a floor set as a fraction of its own seasonal peak (self-scaling), with a
# small absolute backstop. Both tunable; lowering NPP_FLOOR_FRAC extends the
# curve towards the winter at the cost of a rapidly widening SD ribbon
# (0.15 -> Feb-Oct, ribbon within +/-1.5; 0.05 -> near-full year, ribbon > 4).
NPP_FLOOR_FRAC = 0.15              # fraction of that year's peak NPP
NPP_FLOOR_ABS  = 5.0               # mmol C m-2 d-1
# A DOY is only drawn where this many years survive the floor + coverage cuts.
MIN_YEARS    = 5
STEP_DAYS    = 7


def _plot_segments(ax, x, y, mask, connect, **plot_kw):
    """Plot y over contiguous runs of `mask`, extending each run by one point
    into a neighbour that is `connect`-eligible so adjacent solid/dashed
    segments meet with no visible gap. Only the first run carries the label."""
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return
    label  = plot_kw.pop("label", None)
    breaks = np.where(np.diff(idx) > 1)[0] + 1
    for k, seg in enumerate(np.split(idx, breaks)):
        i0, i1 = seg[0], seg[-1]
        if i0 - 1 >= 0 and connect[i0 - 1]:
            i0 -= 1
        if i1 + 1 < len(x) and connect[i1 + 1]:
            i1 += 1
        ax.plot(x[i0:i1 + 1], y[i0:i1 + 1],
                label=label if k == 0 else None, **plot_kw)


def _fill_segments(ax, x, lo, hi, mask, **fill_kw):
    """fill_between over contiguous runs of `mask` only."""
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return
    label  = fill_kw.pop("label", None)
    breaks = np.where(np.diff(idx) > 1)[0] + 1
    for k, seg in enumerate(np.split(idx, breaks)):
        sl = slice(seg[0], seg[-1] + 1)
        ax.fill_between(x[sl], lo[sl], hi[sl],
                        label=label if k == 0 else None, **fill_kw)


# ---------------------------------------------------------------------------
# Per-year e-ratio on the full timeseries
# ---------------------------------------------------------------------------

def ncp_full_record_splines(ncp_df: pd.DataFrame):
    """Cubic splines of NCP mean and MC sigma through the *entire* 20-day-bin
    record (no per-year fitting), so a given year's winter is interpolated from
    the real bins on either side of the year boundary.

    Returns (t0, t_max, cs_mean, cs_sd) with t in days since t0.
    """
    data = ncp_df.dropna(subset=["ncp_mean"]).sort_values("date")
    t0 = data["date"].min()
    t  = (data["date"] - t0).dt.days.to_numpy(dtype=float)
    return (t0, float(t.max()),
            CubicSpline(t, data["ncp_mean"].to_numpy(dtype=float)),
            CubicSpline(t, data["ncp_sd"].to_numpy(dtype=float)))


def ncp_year_on_doy(ncp_splines, year: int, doy: np.ndarray):
    """Evaluate the full-record NCP spline (mean, sigma) at `doy` of `year`.
    NaN outside the record's time span."""
    t0, t_max, cs_m, cs_s = ncp_splines
    dates = pd.Timestamp(year, 1, 1) + pd.to_timedelta(doy - 1, unit="D")
    t = (dates - t0).days.to_numpy(dtype=float)
    inside = (t >= 0) & (t <= t_max)
    m = np.where(inside, cs_m(t), np.nan)
    s = np.where(inside, cs_s(t), np.nan)
    return m, s


def npp_year_on_doy(npp_df: pd.DataFrame, year: int, doy: np.ndarray):
    """That year's own Fig-4 NPP spline (zero-anchored at DOY 1/365) evaluated
    at `doy`, converted mg C -> mmol C m-2 d-1. NaN if the year has < 3 months.
    """
    yr_df = npp_df[npp_df["date"].dt.year == year]
    d, y = npp_mod.individual_year_spline(yr_df, step_days=STEP_DAYS)
    if d is None:
        return np.full_like(doy, np.nan, dtype=float)
    assert np.allclose(d, doy), "per-year NPP spline grid must match doy grid"
    return y * MGC_TO_MMOLC


def build_eratio_ensemble(npp_df: pd.DataFrame, ncp_df: pd.DataFrame,
                          doy: np.ndarray):
    """Per-year e-ratio on the common DOY grid, then the across-year statistics.

    Returns (stats, per_year) where

      stats : dict of ndarray on `doy`
          er_mean, er_sd   ensemble mean and 1 SD of the *ratio* across years
          er_sem           standard error of the mean (er_sd / sqrt(n))
          er_mc            MC-propagated sigma on the ensemble mean, i.e.
                           sqrt(mean_y (sd_NCP_y / NPP_y)^2 / n)
          er_tot           sqrt(er_sem^2 + er_mc^2)
          n                number of contributing years
      per_year : DataFrame  long-form year x doy table (the full timeseries the
          statistics are built from, including the unfloored ratios), written
          to CSV for traceability.
    """
    ncp_spl = ncp_full_record_splines(ncp_df)
    years = [y for y in range(COMMON_YEAR_START, COMMON_YEAR_END + 1)
             if (npp_df["date"].dt.year == y).any()]

    rows = []
    er_cols, mc_cols = [], []
    for yr in years:
        npp_y = npp_year_on_doy(npp_df, yr, doy)
        ncp_y, ncp_sd_y = ncp_year_on_doy(ncp_spl, yr, doy)

        floor = max(NPP_FLOOR_ABS, NPP_FLOOR_FRAC * np.nanmax(npp_y))
        finite = np.isfinite(npp_y) & np.isfinite(ncp_y)
        ok = finite & (npp_y >= floor)
        # unfloored ratio (kept in the CSV so nothing is silently discarded)
        er_all = ncp_y / np.where(finite & (npp_y > 0), npp_y, np.nan)
        denom = np.where(ok, npp_y, np.nan)     # avoid 0-division outside `ok`
        er = ncp_y / denom
        mc = np.abs(ncp_sd_y) / denom

        er_cols.append(er)
        mc_cols.append(mc)
        rows.append(pd.DataFrame({
            "year": yr, "doy": doy,
            "npp_mmolC": npp_y, "ncp_mmolC": ncp_y, "ncp_sd_mmolC": ncp_sd_y,
            "npp_floor": floor, "above_floor": ok,
            "eratio_unfloored": er_all,
            "eratio": er, "eratio_sigma_mc": mc,
        }))

    E  = np.vstack(er_cols)          # (n_years, n_doy)
    MC = np.vstack(mc_cols)

    # Column-wise nan-aware statistics done by hand (rather than nanmean/nanstd)
    # so all-NaN DOYs stay silent instead of raising "mean of empty slice".
    finite = np.isfinite(E)
    n = finite.sum(axis=0)
    enough = n >= MIN_YEARS
    n_safe = np.maximum(n, 1)

    er_mean = np.where(enough,
                       np.where(finite, E, 0.0).sum(axis=0) / n_safe, np.nan)
    # ddof=1: sample SD of the year-to-year spread of the ratio
    var = (np.where(finite, (E - er_mean) ** 2, 0.0).sum(axis=0)
           / np.maximum(n - 1, 1))
    er_sd = np.where(enough, np.sqrt(var), np.nan)
    er_mc = np.where(enough,
                     np.sqrt(np.where(finite, MC ** 2, 0.0).sum(axis=0)
                             / n_safe / n_safe),
                     np.nan)
    er_sem = er_sd / np.sqrt(n_safe)
    er_tot = np.sqrt(er_sem ** 2 + er_mc ** 2)

    stats = {"er_mean": er_mean, "er_sd": er_sd, "er_sem": er_sem,
             "er_mc": er_mc, "er_tot": er_tot, "n": n, "years": years}
    return stats, pd.concat(rows, ignore_index=True)


def seasonal_summary(per_year: pd.DataFrame):
    """Per-year season-integrated e-ratio (sum NCP / sum NPP over the observed
    Feb-Oct window) and its across-year mean +/- SD.

    This is the headline number: it is an integral ratio, so it is not
    sensitive to the near-zero-NPP shoulder days the point-wise ratio is.
    """
    d = per_year.copy()
    d["month"] = [(pd.Timestamp(2021, 1, 1)
                   + pd.Timedelta(days=int(x) - 1)).month for x in d["doy"]]
    # The NPP floor is a guard on the *point-wise* ratio only; the integral is
    # stable at low NPP, so it runs over the whole observed Feb-Oct window.
    obs = d[d["month"].isin(OBS_MONTHS)
            & d["npp_mmolC"].notna() & d["ncp_mmolC"].notna()]
    g = (obs.groupby("year")
            .agg(n_doy=("doy", "size"),
                 npp_sum=("npp_mmolC", "sum"),
                 ncp_sum=("ncp_mmolC", "sum"),
                 eratio_pointwise_mean=("eratio", "mean"))
            .reset_index())
    g["eratio_int"] = g["ncp_sum"] / g["npp_sum"]
    return g[["year", "n_doy", "npp_sum", "ncp_sum",
              "eratio_int", "eratio_pointwise_mean"]]


def main() -> int:
    # ---- NPP synthetic year (mg C -> mmol C) ----
    npp_df = pd.read_csv(npp_mod.CSV_FILE, parse_dates=["date"])
    doy, npp_mean, npp_p10, npp_p90, npp_clim = npp_mod.build_synthetic_year(
        npp_df, step_days=STEP_DAYS)
    npp_mean_c = npp_mean * MGC_TO_MMOLC
    npp_p10_c  = npp_p10  * MGC_TO_MMOLC
    npp_p90_c  = npp_p90  * MGC_TO_MMOLC

    # ---- NCP synthetic year (already mmol C) ----
    # The NCP record runs to 2025 but NPP stops in 2023: the plotted NCP
    # climatology is restricted to the common years so the two curves and the
    # title describe the same window. The per-year e-ratio below still splines
    # the *full* record (2015-2025) — that only supplies interpolation support
    # either side of a common year, never an extra year to average over.
    ncp_df_full = ncp_mod.load_data()
    ncp_df = ncp_df_full[
        ncp_df_full["date"].dt.year.between(COMMON_YEAR_START, COMMON_YEAR_END)
    ].reset_index(drop=True)
    doy_ncp, ncp_mean, ncp_p10, ncp_p90, ncp_clim = ncp_mod.build_synthetic_year(
        ncp_df, step_days=STEP_DAYS)
    assert np.allclose(doy, doy_ncp), "NPP/NCP synthetic-year grids must match"

    # ---- e-ratio: computed per year on the full timeseries, then averaged ----
    st, per_year = build_eratio_ensemble(npp_df, ncp_df_full, doy)
    eratio, er_sd, er_n = st["er_mean"], st["er_sd"], st["n"]
    years = st["years"]

    doy_month = np.array(
        [(pd.Timestamp(2021, 1, 1) + pd.Timedelta(days=int(d) - 1)).month for d in doy]
    )
    is_obs = np.isin(doy_month, list(OBS_MONTHS))   # Feb-Oct: real CbPM NPP
    is_win = ~is_obs                                # Nov-Jan: NPP from spline
    have   = np.isfinite(eratio)

    # ---- tables ----
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "doy": doy, "month": doy_month,
        "eratio_mean": eratio, "eratio_sd": er_sd,
        "eratio_sem": st["er_sem"], "eratio_sigma_mc": st["er_mc"],
        "eratio_sigma_total": st["er_tot"], "n_years": er_n,
        "npp_clim_mmolC": npp_mean_c, "ncp_clim_mmolC": ncp_mean,
    }).to_csv(OUT_DOY, index=False)
    per_year.to_csv(OUT_YEAR, index=False)
    print(f"Wrote {OUT_DOY}")
    print(f"Wrote {OUT_YEAR}")

    seas = seasonal_summary(per_year)
    er_int_m, er_int_s = seas["eratio_int"].mean(), seas["eratio_int"].std(ddof=1)
    print(f"\nYears used: {years[0]}-{years[-1]} (n = {len(years)})")
    print(seas.round(3).to_string(index=False))
    print(f"\nFeb-Oct integrated e-ratio = {er_int_m:.3f} +/- {er_int_s:.3f} "
          f"(mean +/- 1 SD across {len(seas)} years)")
    j_peak = int(np.nanargmax(np.where(have, eratio, -np.inf)))
    print(f"Peak point-wise e-ratio = {eratio[j_peak]:.3f} "
          f"+/- {er_sd[j_peak]:.3f} (1 SD across years) at DOY {doy[j_peak]:.0f}"
          f", n = {er_n[j_peak]} yr")
    print(f"e-ratio drawn over DOY {doy[have].min():.0f}-{doy[have].max():.0f} "
          f"(NPP floor = max({NPP_FLOOR_ABS:.0f}, "
          f"{NPP_FLOOR_FRAC:.2f} x that year's peak NPP))")

    # ---- figure ----
    mpl.rcParams.update({"axes.spines.right": True})   # need the right spine here

    fig, ax1 = plt.subplots(figsize=(9.5, 5.8), dpi=300)
    ax2 = ax1.twinx()

    ax1.axhline(0, color="0.65", lw=0.8, ls=":", zorder=1)

    ax1.fill_between(doy, npp_p10_c, npp_p90_c, color=GREEN_LIGHT, alpha=0.55,
                     linewidth=0, zorder=2, label="NPP p10–p90 (inter-annual)")
    ax1.plot(doy, npp_mean_c, color=GREEN, lw=2.0, zorder=4,
             label="NPP (CbPM, 0-200 m)")

    ax1.fill_between(doy, ncp_p10, ncp_p90, color=BLUE_LIGHT, alpha=0.55,
                     linewidth=0, zorder=2, label="NCP p10–p90 (inter-annual)")
    ax1.plot(doy, ncp_mean, color=BLUE, lw=2.0, zorder=4,
             label="NCP (nitrate budget)")

    ax1.set_ylabel("Production (mmol C m$^{-2}$ d$^{-1}$)", fontsize=11)
    # Zero must sit at the same height on both axes. The e-ratio axis is fixed
    # to (-1, 1), i.e. zero at mid-height, so the production axis is made
    # symmetric about zero as well.
    prod_lim = 1.05 * np.nanmax(
        np.abs(np.concatenate([npp_p10_c, npp_p90_c, ncp_p10, ncp_p90]))
    )
    ax1.set_ylim(-prod_lim, prod_lim)
    ax1.set_xlim(1, 365)
    mid_doys = [pd.Timestamp(2021, m, 15).timetuple().tm_yday for m in range(1, 13)]
    ax1.set_xticks(mid_doys)
    ax1.set_xticklabels(MONTH_LABELS)
    ax1.set_xlabel("Month", fontsize=11)
    ax1.legend(loc="upper left", fontsize=8.5, ncol=1)

    # e-ratio, overlaid on secondary axis
    ax2.axhline(0, color=ERATIO_COLOR, lw=0.7, ls=":", alpha=0.4, zorder=3)
    ax2.axhline(1, color=ERATIO_COLOR, lw=0.8, ls="--", alpha=0.5, zorder=3)

    # Inter-annual +/-1 SD of the per-year ratios
    er_lo_band = eratio - er_sd
    er_hi_band = eratio + er_sd
    _fill_segments(ax2, doy, er_lo_band, er_hi_band, have,
                   color=ERATIO_LIGHT, alpha=0.45, linewidth=0, zorder=5,
                   label="e-ratio ±1 SD (inter-annual)")

    # Observed-NPP window: solid line + markers
    _plot_segments(ax2, doy, eratio, is_obs & have, connect=have,
                   color=ERATIO_COLOR, lw=2.4, ls="-",
                   marker="o", ms=4.5, markevery=3, zorder=7,
                   label="e-ratio (NCP:NPP), observed NPP")
    # Winter: interpolated-NPP window, dashed (matches Fig 4 winter styling)
    _plot_segments(ax2, doy, eratio, is_win & have, connect=have,
                   color=ERATIO_COLOR, lw=1.8, ls="--", zorder=7,
                   label="e-ratio, winter (interpolated NPP)")

    ax2.set_ylabel("e-ratio  (NCP : NPP)", fontsize=11, color=ERATIO_COLOR)
    # Fixed to (-1, 1) so zero sits at mid-height, matching the symmetric
    # production axis. The +/-1 SD ribbon is clipped where it runs outside
    # that range; the unclipped values are in eratio_doy_stats.csv.
    ax2.set_ylim(-1.0, 1.0)
    ax2.tick_params(axis="y", colors=ERATIO_COLOR)
    ax2.spines["right"].set_color(ERATIO_COLOR)

    h2, l2 = ax2.get_legend_handles_labels()
    ax2.legend(h2, l2, loc="upper right", fontsize=8.5)

    fig.suptitle(
        "Synthetic-year NPP, NCP and e-ratio\n"
        f"Iceland Basin & Irminger Sea (pooled, {years[0]}–{years[-1]})",
        fontsize=11.5, y=0.99,
    )
    # Caption below the axes — it used to sit inside ax1 and overplot the
    # autumn NCP curve.
    caption = [
        f"e-ratio computed per year on the full {years[0]}–{years[-1]} "
        f"timeseries (n = {len(years)} yr) and then averaged, not as the ratio "
        "of the two climatologies; band = ±1 SD across years.",
    ]
    if bool((is_win & have).any()):
        caption.append("Solid = observed CbPM NPP (Feb–Oct); dashed = winter "
                       "spline-interpolated NPP (Nov–Jan).")
    caption.append(
        f"A year is dropped where its NPP < {NPP_FLOOR_FRAC:.2f} × its seasonal "
        f"peak (the point-wise ratio diverges as NPP → 0) and a DOY is drawn "
        f"where ≥ {MIN_YEARS} yr remain, hence the "
        f"DOY {doy[have].min():.0f}–{doy[have].max():.0f} span."
    )
    caption.append(
        f"Season-integrated e-ratio (∑NCP/∑NPP over Feb–Oct, per year) = "
        f"{er_int_m:.2f} ± {er_int_s:.2f} (mean ± 1 SD across years; "
        f"range {seas['eratio_int'].min():.2f}–{seas['eratio_int'].max():.2f})."
    )

    fig.tight_layout()
    fig.text(0.5, -0.005, "\n".join(caption), ha="center", va="top",
             fontsize=7.5, color="0.45")
    fig.savefig(OUT_PNG, bbox_inches="tight", dpi=300)
    print(f"Wrote {OUT_PNG}")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
