"""
Publication-quality NCP figure — two panels. Mirrors the style of
npp_timeseries_figure.py but for the nitrate-drawdown NCP budget, with the
shaded ribbon propagating the Monte Carlo simulation uncertainty (CANYON-B
prediction CI + profile bootstrap, 200 iterations) rather than a spatial SD.

Panel A  Full 2015-2025, 20-day-bin timeseries (BGC-Argo nitrate budget,
         Iceland Basin + Irminger Sea pooled). Cubic-spline smoothing through
         the actual bins (the underlying record is continuous, unlike the
         CbPM/NPP product there is no winter data gap to interpolate across).
         Ribbon = Monte Carlo mean ± 1 SD.

Panel B  Weekly synthetic-year climatology built from a periodic cubic
         spline through the per-calendar-month mean (inter-annual p10/p90
         ribbon + mean), individual annual cycles overlaid as thin grey
         lines.

Input:  Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_uncertainty.xlsx
Output: Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_timeseries_publication.png
        Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_synthetic_year.png
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.interpolate import CubicSpline

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT  = Path(__file__).resolve().parents[1]
BASIN_DIR  = REPO_ROOT / "Output" / "IcelandIrminger_2015_2025" / "ncp" / "IcelandIrminger"
XLSX_FILE  = BASIN_DIR / "ncp_uncertainty.xlsx"
OUT_BOTH   = BASIN_DIR / "ncp_timeseries_publication.png"
OUT_SYNTH  = BASIN_DIR / "ncp_synthetic_year.png"

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
BLUE       = "#1f4e79"
BLUE_LIGHT = "#b8cbe0"
GREY_LINE  = "#aaaaaa"
ZERO_LINE  = "#888888"

MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

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
# Data loading
# ---------------------------------------------------------------------------

def load_data() -> pd.DataFrame:
    df = pd.read_excel(XLSX_FILE)
    df = df.rename(columns={
        "date_grid": "date",
        "NCP_mean":  "ncp_mean",
        "NCP_sd":    "ncp_sd",
        "NCP_q05":   "ncp_q05",
        "NCP_q95":   "ncp_q95",
    })
    df["date"] = pd.to_datetime(df["date"])
    return df.sort_values("date").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Spline helpers
# ---------------------------------------------------------------------------

def build_full_spline(df: pd.DataFrame, step_days: int = 7):
    """Cubic splines (mean, lo = mean-sd, hi = mean+sd) on a weekly grid.
    No zero-anchoring: NCP is signed (negative = net heterotrophy) and the
    record is already continuous year-round, so there is no gap to bridge.
    Returns (dates_dense, y_mean, y_lo, y_hi)."""
    data = df.dropna(subset=["ncp_mean"])

    t0     = data["date"].min()
    t_days = (data["date"] - t0).dt.days.values.astype(float)
    y_m    = data["ncp_mean"].values
    y_lo   = (data["ncp_mean"] - data["ncp_sd"]).values
    y_hi   = (data["ncp_mean"] + data["ncp_sd"]).values

    cs_m  = CubicSpline(t_days, y_m)
    cs_lo = CubicSpline(t_days, y_lo)
    cs_hi = CubicSpline(t_days, y_hi)

    t_max   = t_days.max()
    t_dense = np.arange(0, t_max + step_days, step_days, dtype=float)
    dates   = pd.to_datetime(t0) + pd.to_timedelta(t_dense, unit="D")

    return dates, cs_m(t_dense), cs_lo(t_dense), cs_hi(t_dense)


def build_synthetic_year(df: pd.DataFrame, step_days: int = 7):
    """Inter-annual climatology per calendar month -> periodic weekly spline
    (Dec wraps smoothly into Jan, appropriate for a signed, full-year series).
    Returns (doy_dense, y_mean, y_p10, y_p90, clim_df)."""
    d = df.copy()
    d["month"] = d["date"].dt.month
    valid = d.dropna(subset=["ncp_mean"])

    clim = (valid.groupby("month")["ncp_mean"]
                 .agg(
                     clim_mean="mean",
                     clim_std="std",
                     clim_p10=lambda x: np.percentile(x, 10),
                     clim_p90=lambda x: np.percentile(x, 90),
                 )
                 .reset_index()
                 .sort_values("month"))
    clim["doy"] = clim["month"].apply(
        lambda m: pd.Timestamp(2021, m, 15).timetuple().tm_yday
    )

    # Close the annual cycle: replicate the first (January) anchor one full
    # period (365 d) later so the periodic spline wraps Dec -> Jan smoothly.
    doy_pts  = clim["doy"].tolist() + [clim["doy"].iloc[0] + 365]
    mean_pts = clim["clim_mean"].tolist() + [clim["clim_mean"].iloc[0]]
    p10_pts  = clim["clim_p10"].tolist()  + [clim["clim_p10"].iloc[0]]
    p90_pts  = clim["clim_p90"].tolist()  + [clim["clim_p90"].iloc[0]]

    cs_m   = CubicSpline(doy_pts, mean_pts, bc_type="periodic")
    cs_p10 = CubicSpline(doy_pts, p10_pts,  bc_type="periodic")
    cs_p90 = CubicSpline(doy_pts, p90_pts,  bc_type="periodic")

    doy_dense = np.linspace(1, 365, int(365 / step_days) + 1)
    return doy_dense, cs_m(doy_dense), cs_p10(doy_dense), cs_p90(doy_dense), clim


def individual_year_spline(df_year: pd.DataFrame, step_days: int = 7):
    """Natural spline for a single year's bins, evaluated only over that
    year's own DOY range (no zero anchors, no extrapolation)."""
    data = df_year.dropna(subset=["ncp_mean"]).sort_values("date")
    if len(data) < 4:
        return None, None
    doy_pts = [d.timetuple().tm_yday for d in data["date"]]
    y_pts   = data["ncp_mean"].tolist()
    cs = CubicSpline(doy_pts, y_pts)
    doy_dense = np.linspace(min(doy_pts), max(doy_pts),
                            int((max(doy_pts) - min(doy_pts)) / step_days) + 1)
    return doy_dense, cs(doy_dense)


# ---------------------------------------------------------------------------
# Panel-drawing functions
# ---------------------------------------------------------------------------

def draw_panel_a(ax, df: pd.DataFrame):
    dates, y_m, y_lo, y_hi = build_full_spline(df)

    yr_min = int(df["date"].dt.year.min())
    yr_max = int(df["date"].dt.year.max())

    ax.axhline(0, color=ZERO_LINE, lw=0.8, ls=":", zorder=1)

    ax.fill_between(dates, y_lo, y_hi,
                    color=BLUE_LIGHT, linewidth=0, alpha=0.65,
                    label="±1 SD (Monte Carlo, n=200)")

    ax.plot(dates, y_m, lw=2.0, ls="-", color=BLUE,
            label="NCP (nitrate budget, 20-day bins)", zorder=4)

    valid = df.dropna(subset=["ncp_mean"])
    ax.scatter(valid["date"], valid["ncp_mean"],
               s=16, color=BLUE, zorder=5, linewidths=0.5,
               edgecolors="white")

    ax.set_xlim(pd.Timestamp(yr_min, 1, 1), pd.Timestamp(yr_max + 1, 1, 1))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))

    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4, 7, 10]))

    ax.set_ylabel("NCP (mmol C m$^{-2}$ d$^{-1}$)", fontsize=11)
    ax.set_xlabel("")

    ax.legend(loc="upper left", ncol=1, fontsize=9)
    ax.text(-0.06, 1.03, "(a)", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top")


def draw_panel_b(ax, df: pd.DataFrame):
    doy, y_m, y_p10, y_p90, clim = build_synthetic_year(df)

    years = sorted(df["date"].dt.year.unique())
    for yr in years:
        yr_df = df[df["date"].dt.year == yr].copy()
        d, y  = individual_year_spline(yr_df)
        if d is not None:
            ax.plot(d, y, lw=0.8, color=GREY_LINE, alpha=0.55, zorder=2)

    ax.axhline(0, color=ZERO_LINE, lw=0.8, ls=":", zorder=1)

    ax.fill_between(doy, y_p10, y_p90,
                    color=BLUE_LIGHT, linewidth=0, alpha=0.65,
                    label="p10–p90 (inter-annual)", zorder=3)

    ax.plot(doy, y_m, lw=2.2, ls="-", color=BLUE,
            label="Mean NCP (nitrate budget)", zorder=5)

    ax.scatter(clim["doy"], clim["clim_mean"],
               s=32, color=BLUE, zorder=6, linewidths=0.6,
               edgecolors="white")

    mid_doys = [pd.Timestamp(2021, m, 15).timetuple().tm_yday
                for m in range(1, 13)]
    ax.set_xticks(mid_doys)
    ax.set_xticklabels(MONTH_LABELS, fontsize=9)
    ax.set_xlim(1, 365)
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))

    ax.set_ylabel("NCP (mmol C m$^{-2}$ d$^{-1}$)", fontsize=11)
    ax.set_xlabel("Month", fontsize=11)

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
    df = load_data()

    # ---- 2-panel combined figure ----
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), dpi=300,
                             gridspec_kw={"hspace": 0.38})

    draw_panel_a(axes[0], df)
    draw_panel_b(axes[1], df)

    fig.suptitle(
        "Basin-scale NCP — Iceland Basin & Irminger Sea (pooled)\n"
        "BGC-Argo nitrate-drawdown budget, CANYON-B (2015–2025)",
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
        "Synthetic-year NCP — Iceland Basin & Irminger Sea (pooled)\n"
        "BGC-Argo nitrate budget (weekly climatology, 2015–2025)",
        fontsize=10,
    )
    fig2.tight_layout()
    fig2.savefig(OUT_SYNTH, bbox_inches="tight", dpi=300)
    print(f"Wrote {OUT_SYNTH}")
    plt.close(fig2)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
