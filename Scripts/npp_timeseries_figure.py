"""
Publication-quality NPP figure — three panels.

Panel A  Full 2015-2023 monthly timeseries with cubic-spline interpolation
         across winter (Nov-Jan) data gaps.  Interpolated segments are shown
         as dashed lines; spatial ±1 SD ribbon covers the observed season.

Panel B  Weekly synthetic-year climatology derived by averaging each calendar
         month across all 9 years (inter-annual p10/p90 ribbon + mean), with
         individual annual cycles overlaid as thin lines coloured by how
         closely each year tracks the climatology (R^2).

Panel C  Seasonal reproducibility. Each year's observed months (Feb-Oct) are
         regressed on a *leave-one-out* climatology (the year under test
         excluded from the multi-year mean) evaluated at the same day-of-year:
             NPP_year(t) = slope * NPP_clim(doy(t)) + intercept
         Bars are the R^2 of that fit — the fraction of the year's variance
         explained by the climatological seasonal shape and phase, i.e. how
         "normal" the year's seasonal cycle was. Black diamonds (right axis)
         are the slope: the amplitude anomaly, >1 = stronger seasonal cycle
         than the norm.  Colour scale is shared with the NCP figure, so the
         two are directly comparable.

Input:  Output/cmems_npp_timeseries_domain_mean.csv
Output: Output/npp_timeseries_publication.png        (all three panels)
        Output/npp_synthetic_year.png                 (panel B only)
        Output/npp_seasonal_reproducibility.png       (panel C only)
        Output/npp_seasonal_reproducibility.csv       (per-year metrics)
"""

from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import (            # noqa: E402
    GREEN, GREEN_LIGHT,
    R2_CMAP, R2_VMIN, R2_VMAX,
    assert_within_axis,
    per_year_climatology_regression,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
CSV_FILE  = REPO_ROOT / "Output" / "cmems_npp_timeseries_domain_mean.csv"
OUT_BOTH  = REPO_ROOT / "Output" / "npp_timeseries_publication.png"
OUT_SYNTH = REPO_ROOT / "Output" / "npp_synthetic_year.png"
OUT_REPRO = REPO_ROOT / "Output" / "npp_seasonal_reproducibility.png"
OUT_STATS = REPO_ROOT / "Output" / "npp_seasonal_reproducibility.csv"

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
# GREEN / GREEN_LIGHT now come from figure_config (single source of truth);
# only the figure-local greys stay here.
GREY_LINE   = "#aaaaaa"
WINTER_FACE = "#f0f0f0"   # subtle shading for interpolated winters

MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
WINTER_MONTHS = {11, 12, 1}     # months with no CbPM output

# Fraction of the reproducibility panel's height reserved above the data for
# the legend (applied to both the R^2 and the slope axis so they stay aligned).
HEADROOM = 0.13

mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size":       10,
    "axes.linewidth":  0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "legend.fontsize":  9,
})


# ---------------------------------------------------------------------------
# Spline helpers
# ---------------------------------------------------------------------------

def _augment_winter_anchors(df: pd.DataFrame) -> pd.DataFrame:
    """Add Dec-15 = 0 anchors (and start/end = 0) so the spline dips to
    near-zero in winter rather than over-shooting between seasons."""
    yr_min = int(df["date"].dt.year.min())
    yr_max = int(df["date"].dt.year.max())
    extra = []
    extra.append({"date": pd.Timestamp(yr_min, 1, 1),
                  "npp_int_mean": 0.0, "npp_int_std": 0.0})
    for yr in range(yr_min, yr_max + 1):
        extra.append({"date": pd.Timestamp(yr, 12, 15),
                      "npp_int_mean": 0.0, "npp_int_std": 0.0})
    extra.append({"date": pd.Timestamp(yr_max + 1, 1, 1),
                  "npp_int_mean": 0.0, "npp_int_std": 0.0})
    out = (pd.concat([df, pd.DataFrame(extra)], ignore_index=True)
             .sort_values("date")
             .reset_index(drop=True))
    return out


def build_full_spline(df: pd.DataFrame, step_days: int = 7):
    """Cubic splines (mean, lo = mean−std, hi = mean+std) on a weekly grid.
    Returns (dates_dense, y_mean, y_lo, y_hi, is_winter_mask)."""
    df_aug = _augment_winter_anchors(df)
    data = df_aug.dropna(subset=["npp_int_mean"])

    t0     = data["date"].min()
    t_days = (data["date"] - t0).dt.days.values.astype(float)
    y_m    = data["npp_int_mean"].values
    y_lo   = (data["npp_int_mean"] - data["npp_int_std"]).values.clip(0)
    y_hi   = (data["npp_int_mean"] + data["npp_int_std"]).values

    cs_m  = CubicSpline(t_days, y_m)
    cs_lo = CubicSpline(t_days, y_lo)
    cs_hi = CubicSpline(t_days, y_hi)

    t_max   = (df_aug["date"].max() - t0).days
    t_dense = np.arange(0, t_max + step_days, step_days, dtype=float)
    dates   = pd.to_datetime(t0) + pd.to_timedelta(t_dense, unit="D")

    y_mean = cs_m(t_dense).clip(0)
    y_lo_d = cs_lo(t_dense).clip(0)
    y_hi_d = cs_hi(t_dense).clip(0)

    is_winter = np.array([d.month in WINTER_MONTHS for d in dates])
    return dates, y_mean, y_lo_d, y_hi_d, is_winter


def build_climatology(df: pd.DataFrame):
    """Per-calendar-month inter-annual statistics → cubic splines.

    Factored out of ``build_synthetic_year`` so the same climatology can be
    re-evaluated at arbitrary days-of-year by the seasonal-reproducibility
    regression (including on year-subsets, for the leave-one-out norm).

    Returns (cs_mean, cs_p10, cs_p90, clim_df).
    """
    df = df.copy()
    df["month"] = df["date"].dt.month
    valid = df.dropna(subset=["npp_int_mean"])

    clim = (valid.groupby("month")["npp_int_mean"]
                 .agg(
                     clim_mean="mean",
                     clim_std="std",
                     clim_p10=lambda x: np.percentile(x, 10),
                     clim_p90=lambda x: np.percentile(x, 90),
                 )
                 .reset_index()
                 .sort_values("month"))
    clim["doy"] = clim["month"].apply(
        lambda m: pd.Timestamp(2020, m, 15).timetuple().tm_yday
    )

    # NB the [1] / [365] zero anchors only shape the (unobserved) winter tails;
    # the reproducibility regression is evaluated on the observed Feb-Oct
    # months only, so it is unaffected by them.
    doy_pts  = [1]   + clim["doy"].tolist()           + [365]
    mean_pts = [0.0] + clim["clim_mean"].tolist()     + [0.0]
    p10_pts  = [0.0] + clim["clim_p10"].tolist()      + [0.0]
    p90_pts  = [0.0] + clim["clim_p90"].tolist()      + [0.0]

    return (CubicSpline(doy_pts, mean_pts),
            CubicSpline(doy_pts, p10_pts),
            CubicSpline(doy_pts, p90_pts),
            clim)


def build_synthetic_year(df: pd.DataFrame, step_days: int = 7):
    """Inter-annual climatology per calendar month → weekly spline.
    Returns (doy_dense, y_mean, y_p10, y_p90, clim_df)."""
    cs_m, cs_p10, cs_p90, clim = build_climatology(df)
    doy_dense = np.linspace(1, 365, int(365 / step_days) + 1)
    return (doy_dense,
            cs_m(doy_dense).clip(0),
            cs_p10(doy_dense).clip(0),
            cs_p90(doy_dense).clip(0),
            clim)


def compute_reproducibility(df: pd.DataFrame, leave_one_out: bool = True):
    """Per-year regression of the observed monthly NPP on the climatology.

    With ``leave_one_out`` (default) the climatology used as the predictor for
    a given year is rebuilt from the *other* years only, so each year is scored
    against a norm it did not help define. See
    ``figure_config.per_year_climatology_regression`` for the returned columns.

    Deliberately **unweighted**, unlike the NCP figure. ``npp_int_std`` is the
    *spatial* standard deviation of NPP across the domain's grid cells within a
    month, not an uncertainty on the monthly mean, so 1/sigma^2 weighting would
    down-weight spatially heterogeneous months rather than poorly-constrained
    ones. The NPP product also has no winter-precision problem to correct for:
    the winter months are simply absent (n_cells = 0), not noisy. Because of
    this the NPP and NCP R^2 are not strictly like-for-like -- state that in
    any direct comparison.
    """
    def clim_fn(year, doy):
        ref = df[df["date"].dt.year != year] if leave_one_out else df
        cs_m, *_ = build_climatology(ref)
        return cs_m(doy)

    return per_year_climatology_regression(df["date"], df["npp_int_mean"],
                                           clim_fn, min_points=4)


def individual_year_spline(df_year: pd.DataFrame, step_days: int = 7):
    """Spline for a single year's monthly data + zero anchors."""
    data = df_year.dropna(subset=["npp_int_mean"])
    if len(data) < 3:
        return None, None
    doy_pts = ([1]
               + [d.timetuple().tm_yday for d in data["date"]]
               + [365])
    y_pts   = [0.0] + data["npp_int_mean"].tolist() + [0.0]
    cs = CubicSpline(doy_pts, y_pts)
    doy_dense = np.linspace(1, 365, int(365 / step_days) + 1)
    return doy_dense, cs(doy_dense).clip(0)


# ---------------------------------------------------------------------------
# Panel-drawing functions
# ---------------------------------------------------------------------------

def _winter_spans(ax, yr_min: int, yr_max: int):
    """Shade Nov 15 → Feb 15 winters as very light grey bands."""
    for yr in range(yr_min, yr_max + 2):
        start = pd.Timestamp(yr - 1, 11, 15)
        end   = pd.Timestamp(yr, 2, 15)
        ax.axvspan(start, end, color=WINTER_FACE, linewidth=0, zorder=0)


def draw_panel_a(ax, df: pd.DataFrame):
    dates, y_m, y_lo, y_hi, is_win = build_full_spline(df)

    yr_min = int(df["date"].dt.year.min())
    yr_max = int(df["date"].dt.year.max())
    _winter_spans(ax, yr_min, yr_max)

    # ±1 SD ribbon — observed season only
    obs = ~is_win
    # fill gaps between observed blocks for ribbon continuity
    ax.fill_between(dates, y_lo, y_hi,
                    where=obs,
                    color=GREEN_LIGHT, linewidth=0, alpha=0.65,
                    label="±1 SD (spatial)")

    # Spline line — split observed / interpolated
    # Build continuous segments for each category
    def _plot_segments(mask, lw, ls, color, label, **kw):
        """Plot contiguous masked segments; only the first gets a legend label."""
        idx  = np.where(mask)[0]
        if len(idx) == 0:
            return
        breaks = np.where(np.diff(idx) > 1)[0] + 1
        segs   = np.split(idx, breaks)
        for k, seg in enumerate(segs):
            i0  = max(seg[0] - 1, 0)
            i1  = min(seg[-1] + 2, len(dates))
            lbl = label if k == 0 else "_nolegend_"
            ax.plot(dates[i0:i1], y_m[i0:i1], lw=lw, ls=ls,
                    color=color, label=lbl, **kw)

    _plot_segments(~is_win, lw=2.0, ls="-",  color=GREEN,
                   label="NPP (CbPM, 0-200 m)", zorder=4)
    _plot_segments(is_win,  lw=1.2, ls="--", color="#6aaa85",
                   label="Spline interpolation (winter)", zorder=3)

    # Observed monthly data points
    valid = df.dropna(subset=["npp_int_mean"])
    ax.scatter(valid["date"], valid["npp_int_mean"],
               s=18, color=GREEN, zorder=5, linewidths=0.5,
               edgecolors="white")

    ax.set_xlim(pd.Timestamp(yr_min, 1, 1), pd.Timestamp(yr_max + 1, 1, 1))
    ax.set_ylim(-50, 1850)
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(100))

    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4, 7, 10]))

    ax.set_ylabel("NPP (mg C m$^{-2}$ d$^{-1}$)", fontsize=11)
    ax.set_xlabel("")

    # Legend — custom order
    handles, labels = ax.get_legend_handles_labels()
    # Add a grey patch for winter shading
    win_patch = mpatches.Patch(facecolor=WINTER_FACE, edgecolor="0.7",
                               linewidth=0.5, label="Winter (no data)")
    ax.legend(handles + [win_patch],
              labels  + ["Winter (no data)"],
              loc="upper left", ncol=2, fontsize=9,
              handlelength=2.0, columnspacing=1.5)
    ax.text(-0.06, 1.03, "(a)", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top")


def draw_panel_b(ax, df: pd.DataFrame, stats: pd.DataFrame | None = None,
                 panel_label: str = "(b)"):
    doy, y_m, y_p10, y_p90, clim = build_synthetic_year(df)

    cmap = plt.get_cmap(R2_CMAP)
    norm = mpl.colors.Normalize(vmin=R2_VMIN, vmax=R2_VMAX)
    r2_by_year = ({} if stats is None
                  else dict(zip(stats["year"], stats["r2"])))

    # Individual year splines — coloured by their R^2 against the climatology
    years = sorted(df["date"].dt.year.unique())
    ymaxs = []
    for yr in years:
        yr_df = df[df["date"].dt.year == yr].copy()
        d, y  = individual_year_spline(yr_df)
        if d is None:
            continue
        r2 = r2_by_year.get(int(yr))
        if r2 is None:
            col, alpha = GREY_LINE, 0.55
        else:
            col, alpha = cmap(norm(r2)), 0.9
        ax.plot(d, y, lw=1.0, color=col, alpha=alpha, zorder=2)
        ymaxs.append(np.nanmax(y))

    # p10-p90 ribbon
    ax.fill_between(doy, y_p10, y_p90,
                    color=GREEN_LIGHT, linewidth=0, alpha=0.6,
                    label="p10–p90 (inter-annual)", zorder=3)

    # Mean line — split into observed / interpolated DOY
    doy_month = np.array(
        [(pd.Timestamp(2020, 1, 1) + pd.Timedelta(days=int(d) - 1)).month
         for d in doy]
    )
    is_obs  = np.isin(doy_month, list(range(2, 11)))
    is_win  = ~is_obs

    def _plot_seg_doy(mask, lw, ls, color, label, **kw):
        idx    = np.where(mask)[0]
        if len(idx) == 0:
            return
        breaks = np.where(np.diff(idx) > 1)[0] + 1
        segs   = np.split(idx, breaks)
        for k, seg in enumerate(segs):
            i0  = max(seg[0] - 1, 0)
            i1  = min(seg[-1] + 2, len(doy))
            lbl = label if k == 0 else "_nolegend_"
            ax.plot(doy[i0:i1], y_m[i0:i1], lw=lw, ls=ls,
                    color=color, label=lbl, **kw)

    # White halo keeps the climatology readable over the coloured year lines.
    _plot_seg_doy(is_obs, lw=2.4, ls="-",  color=GREEN,
                  label="Mean NPP (CbPM)", zorder=5,
                  path_effects=[pe.Stroke(linewidth=4.2, foreground="white"),
                                pe.Normal()])
    _plot_seg_doy(is_win, lw=1.2, ls="--", color="#6aaa85",
                  label="Spline interpolation (winter)", zorder=4)

    # Observed monthly climatology dots
    ax.scatter(clim["doy"], clim["clim_mean"],
               s=32, color=GREEN, zorder=6, linewidths=0.6,
               edgecolors="white")

    # x-axis: month labels at mid-month DOY
    mid_doys   = [pd.Timestamp(2020, m, 15).timetuple().tm_yday
                  for m in range(1, 13)]
    ax.set_xticks(mid_doys)
    ax.set_xticklabels(MONTH_LABELS, fontsize=9)
    ax.set_xlim(1, 365)
    # y-limits sized to contain every individual-year spline rather than a
    # hard-coded ceiling the coloured year lines could run through.
    hi = max(ymaxs + [float(np.nanmax(y_p90))])
    ax.set_ylim(-50, hi * 1.12)
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(100))

    ax.set_ylabel("NPP (mg C m$^{-2}$ d$^{-1}$)", fontsize=11)
    ax.set_xlabel("Month", fontsize=11)

    # Year labels in top-right corner
    ax.text(0.98, 0.97, f"n = {len(years)} yr  ({years[0]}–{years[-1]})",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, color="0.45")

    handles, labels = ax.get_legend_handles_labels()
    if stats is not None:
        handles.append(Line2D([], [], color=GREY_LINE, lw=1.0))
        labels.append("Individual years (colour = R$^2$ vs climatology)")
    ax.legend(handles, labels, loc="upper left", ncol=1, fontsize=9)
    if panel_label:
        # -0.075 (not -0.06) clears the 4-digit y tick labels
        ax.text(-0.075, 1.03, panel_label, transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top")

    assert_within_axis(ax, y_m, y_p10, y_p90)


def draw_panel_c(ax, stats: pd.DataFrame, panel_label: str = "(c)"):
    """Seasonal reproducibility: per-year R^2 against the leave-one-out
    climatology (bars) plus the regression slope / amplitude anomaly
    (diamonds, right axis)."""
    cmap = plt.get_cmap(R2_CMAP)
    norm = mpl.colors.Normalize(vmin=R2_VMIN, vmax=R2_VMAX)

    x = np.arange(len(stats))
    colors = [cmap(norm(v)) for v in stats["r2"]]

    ax.bar(x, stats["r2"], width=0.68, color=colors,
           edgecolor="white", linewidth=0.6, zorder=3,
           label="R$^2$ vs climatology (unweighted)")

    mean_r2 = float(stats["r2"].mean())
    ax.axhline(mean_r2, color="0.35", lw=1.0, ls="--", zorder=4,
               label=f"mean R$^2$ = {mean_r2:.2f}")

    ax.set_xticks(x)
    ax.set_xticklabels(stats["year"].astype(str), fontsize=9)
    ax.set_xlim(-0.7, len(stats) - 0.3)
    # HEADROOM of the top 13% is kept free on BOTH axes so the legend never
    # sits on top of a bar or a slope marker.
    ax.set_ylim(0, 1.0 / (1.0 - HEADROOM))
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
    ax.set_ylabel("R$^2$ vs climatology", fontsize=11)
    ax.set_xlabel("Year", fontsize=11)

    # Slope = amplitude anomaly, on a secondary axis. Everything belonging to
    # this axis (markers, reference line, ticks, label) is drawn in GREEN so
    # the axis a symbol belongs to is readable without consulting the legend.
    ax2 = ax.twinx()
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_color(GREEN)
    ax2.axhline(1.0, color=GREEN, lw=0.8, ls=":", alpha=0.6, zorder=2)
    ax2.plot(x, stats["slope"], ls="none", marker="D", ms=5.5,
             color=GREEN, mec="white", mew=0.7, zorder=6,
             label="Slope (amplitude vs norm)")
    smin, smax = float(stats["slope"].min()), float(stats["slope"].max())
    spread = max(smax - smin, 0.2)
    lo2 = min(smin, 1.0) - 0.25 * spread
    hi2 = max(smax, 1.0) + 0.25 * spread
    ax2.set_ylim(lo2, lo2 + (hi2 - lo2) / (1.0 - HEADROOM))
    ax2.set_ylabel("Regression slope", fontsize=11, color=GREEN)
    ax2.tick_params(axis="y", labelsize=9, colors=GREEN)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper center", ncol=3, fontsize=9,
              framealpha=1.0, frameon=True, edgecolor="0.85",
              borderpad=0.4, columnspacing=1.4)

    if panel_label:
        ax.text(-0.075, 1.04, panel_label, transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top")
    ax.set_axisbelow(True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    df = pd.read_csv(CSV_FILE, parse_dates=["date"])

    stats = compute_reproducibility(df)
    OUT_STATS.parent.mkdir(parents=True, exist_ok=True)
    stats.to_csv(OUT_STATS, index=False)
    print(f"Wrote {OUT_STATS}")
    print(stats.round(3).to_string(index=False))

    yr_min = int(df["date"].dt.year.min())
    yr_max = int(df["date"].dt.year.max())
    n_yr   = yr_max - yr_min + 1

    # ---- 3-panel combined figure ----
    fig, axes = plt.subplots(3, 1, figsize=(10, 11.5), dpi=300,
                             gridspec_kw={"hspace": 0.42,
                                          "height_ratios": [1, 1, 0.72]})

    draw_panel_a(axes[0], df)
    draw_panel_b(axes[1], df, stats)
    draw_panel_c(axes[2], stats)

    fig.suptitle(
        "Column-integrated NPP 0–200 m — Iceland Basin & Irminger Sea\n"
        f"CbPM applied to CMEMS 3D BGC product ({yr_min}–{yr_max})",
        fontsize=11, y=0.975,
    )
    fig.text(
        0.5, 0.055,
        "(c) Each year's observed months (Feb–Oct) regressed on a "
        "leave-one-out monthly climatology at the same day-of-year:\n"
        "R$^2$ = variance explained by the climatological seasonal shape and "
        "phase; slope = amplitude relative to the norm. Unweighted — "
        "npp_int_std is a spatial\nSD across grid cells, not an uncertainty on "
        "the monthly mean (the NCP panel is 1/$\\sigma^2$ weighted; the two "
        "R$^2$ are therefore not strictly like-for-like).",
        ha="center", va="top", fontsize=8.5, color="0.35",
    )

    OUT_BOTH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BOTH, bbox_inches="tight", dpi=300)
    print(f"Wrote {OUT_BOTH}")
    plt.close(fig)

    # ---- Synthetic year — standalone ----
    fig2, ax2 = plt.subplots(figsize=(7.5, 4.5), dpi=300)
    draw_panel_b(ax2, df, stats, panel_label="")
    ax2.set_title(
        "Synthetic-year NPP — Iceland Basin & Irminger Sea\n"
        f"CbPM / CMEMS BGC 3D ({n_yr}-yr weekly climatology, "
        f"{yr_min}–{yr_max})",
        fontsize=10,
    )
    fig2.tight_layout()
    fig2.savefig(OUT_SYNTH, bbox_inches="tight", dpi=300)
    print(f"Wrote {OUT_SYNTH}")
    plt.close(fig2)

    # ---- Seasonal reproducibility — standalone ----
    fig3, ax3 = plt.subplots(figsize=(7.5, 3.8), dpi=300)
    draw_panel_c(ax3, stats, panel_label="")
    ax3.set_title(
        "Seasonal reproducibility of basin NPP\n"
        "each year vs leave-one-out climatology "
        f"({yr_min}–{yr_max})",
        fontsize=10,
    )
    fig3.tight_layout()
    fig3.savefig(OUT_REPRO, bbox_inches="tight", dpi=300)
    print(f"Wrote {OUT_REPRO}")
    plt.close(fig3)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
