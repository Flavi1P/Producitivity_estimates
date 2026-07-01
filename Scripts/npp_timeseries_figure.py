"""
Publication-quality NPP figure — two panels.

Panel A  Full 2015-2023 monthly timeseries with cubic-spline interpolation
         across winter (Nov-Jan) data gaps.  Interpolated segments are shown
         as dashed lines; spatial ±1 SD ribbon covers the observed season.

Panel B  Weekly synthetic-year climatology derived by averaging each calendar
         month across all 9 years (inter-annual p10/p90 ribbon + mean), with
         individual annual cycles overlaid as thin grey lines.

Input:  Output/cmems_npp_timeseries_domain_mean.csv
Output: Output/npp_timeseries_publication.png   (both panels)
        Output/npp_synthetic_year.png            (panel B only)
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
from scipy.interpolate import CubicSpline

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
CSV_FILE  = REPO_ROOT / "Output" / "cmems_npp_timeseries_domain_mean.csv"
OUT_BOTH  = REPO_ROOT / "Output" / "npp_timeseries_publication.png"
OUT_SYNTH = REPO_ROOT / "Output" / "npp_synthetic_year.png"

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
GREEN       = "#1a6e3c"
GREEN_LIGHT = "#c0deca"
GREY_LINE   = "#aaaaaa"
WINTER_FACE = "#f0f0f0"   # subtle shading for interpolated winters

MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
WINTER_MONTHS = {11, 12, 1}     # months with no CbPM output

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


def build_synthetic_year(df: pd.DataFrame, step_days: int = 7):
    """Inter-annual climatology per calendar month → weekly spline.
    Returns (doy_dense, y_mean, y_p10, y_p90, clim_df)."""
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
                 .reset_index())
    clim["doy"] = clim["month"].apply(
        lambda m: pd.Timestamp(2020, m, 15).timetuple().tm_yday
    )

    doy_pts  = [1]   + clim["doy"].tolist()           + [365]
    mean_pts = [0.0] + clim["clim_mean"].tolist()     + [0.0]
    p10_pts  = [0.0] + clim["clim_p10"].tolist()      + [0.0]
    p90_pts  = [0.0] + clim["clim_p90"].tolist()      + [0.0]

    cs_m   = CubicSpline(doy_pts, mean_pts)
    cs_p10 = CubicSpline(doy_pts, p10_pts)
    cs_p90 = CubicSpline(doy_pts, p90_pts)

    doy_dense = np.linspace(1, 365, int(365 / step_days) + 1)
    return (doy_dense,
            cs_m(doy_dense).clip(0),
            cs_p10(doy_dense).clip(0),
            cs_p90(doy_dense).clip(0),
            clim)


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


def draw_panel_b(ax, df: pd.DataFrame):
    doy, y_m, y_p10, y_p90, clim = build_synthetic_year(df)

    # Individual year splines — thin grey
    years = sorted(df["date"].dt.year.unique())
    for yr in years:
        yr_df = df[df["date"].dt.year == yr].copy()
        d, y  = individual_year_spline(yr_df)
        if d is not None:
            ax.plot(d, y, lw=0.8, color=GREY_LINE, alpha=0.55, zorder=2)

    # p10-p90 ribbon
    ax.fill_between(doy, y_p10, y_p90,
                    color=GREEN_LIGHT, linewidth=0, alpha=0.65,
                    label="p10–p90 (inter-annual)", zorder=3)

    # Mean line — split into observed / interpolated DOY
    valid_months = set(range(2, 11))   # Feb-Oct
    is_obs = np.array([
        any(pd.Timestamp(2020, 1, 1) + pd.Timedelta(days=int(d) - 1) is not None
            and (pd.Timestamp(2020, 1, 1) + pd.Timedelta(days=int(d) - 1)).month
            in valid_months
            for _ in [0])
        for d in doy
    ])
    # Simpler: compute from DOY directly
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

    _plot_seg_doy(is_obs, lw=2.2, ls="-",  color=GREEN,
                  label="Mean NPP (CbPM)", zorder=5)
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
    ax.set_ylim(-50, 1850)
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(100))

    ax.set_ylabel("NPP (mg C m$^{-2}$ d$^{-1}$)", fontsize=11)
    ax.set_xlabel("Month", fontsize=11)

    # Year labels in top-right corner
    ax.text(0.98, 0.97, f"n = {len(years)} yr  ({years[0]}–{years[-1]})",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, color="0.45")

    ax.legend(loc="upper left", ncol=1, fontsize=9)
    ax.text(-0.06, 1.03, "(b)", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    df = pd.read_csv(CSV_FILE, parse_dates=["date"])

    # ---- 2-panel combined figure ----
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), dpi=300,
                             gridspec_kw={"hspace": 0.38})

    draw_panel_a(axes[0], df)
    draw_panel_b(axes[1], df)

    fig.suptitle(
        "Column-integrated NPP 0–200 m — Iceland Basin & Irminger Sea\n"
        "CbPM applied to CMEMS 3D BGC product (2015–2023)",
        fontsize=11, y=0.98,
    )

    OUT_BOTH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BOTH, bbox_inches="tight", dpi=300)
    print(f"Wrote {OUT_BOTH}")
    plt.close(fig)

    # ---- Synthetic year — standalone ----
    fig2, ax2 = plt.subplots(figsize=(7.5, 4.5), dpi=300)
    draw_panel_b(ax2, df)
    ax2.set_title(
        "Synthetic-year NPP — Iceland Basin & Irminger Sea\n"
        "CbPM / CMEMS BGC 3D (9-yr weekly climatology, 2015–2023)",
        fontsize=10,
    )
    fig2.tight_layout()
    fig2.savefig(OUT_SYNTH, bbox_inches="tight", dpi=300)
    print(f"Wrote {OUT_SYNTH}")
    plt.close(fig2)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
