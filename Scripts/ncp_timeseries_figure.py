"""
Publication-quality NCP figure — three panels. Mirrors the style of
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
         ribbon + mean), individual annual cycles overlaid as thin lines
         coloured by how closely each year tracks the climatology (R^2).

Panel C  Seasonal reproducibility. Each year's 20-day bins are regressed on a
         *leave-one-out* climatology (the year under test excluded from the
         multi-year mean) evaluated at the same day-of-year:
             NCP_year(t) = slope * NCP_clim(doy(t)) + intercept
         The fit is inverse-variance weighted by the Monte Carlo sigma, which
         is 2-5 mmol C m-2 d-1 in the productive season but 15-32 in winter;
         unweighted, a few unconstrained winter bins dominate the score.
         Bars are the R^2 of that fit — the fraction of the year's variance
         explained by the climatological seasonal shape and phase, i.e. how
         "normal" the year's seasonal cycle was. Black diamonds (right axis)
         are the slope: the amplitude anomaly, >1 = stronger seasonal cycle
         than the norm. Open circles are the unweighted R^2 (sensitivity).

Input:  Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_uncertainty.xlsx
Output: Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_timeseries_publication.png
        Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_synthetic_year.png
        Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_seasonal_reproducibility.png
        Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_seasonal_reproducibility.csv
"""

from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import (            # noqa: E402
    BLUE, BLUE_LIGHT, ZERO_LINE,
    R2_CMAP, R2_VMIN, R2_VMAX,
    assert_within_axis,
    per_year_climatology_regression,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT  = Path(__file__).resolve().parents[1]
BASIN_DIR  = REPO_ROOT / "Output" / "IcelandIrminger_2015_2025" / "ncp" / "IcelandIrminger"
XLSX_FILE  = BASIN_DIR / "ncp_uncertainty.xlsx"
OUT_BOTH   = BASIN_DIR / "ncp_timeseries_publication.png"
OUT_SYNTH  = BASIN_DIR / "ncp_synthetic_year.png"
OUT_REPRO  = BASIN_DIR / "ncp_seasonal_reproducibility.png"
OUT_STATS  = BASIN_DIR / "ncp_seasonal_reproducibility.csv"

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
# BLUE / BLUE_LIGHT / ZERO_LINE now come from figure_config (single source of
# truth); only the figure-local greys stay here.
GREY_LINE  = "#aaaaaa"

MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# Fraction of the reproducibility panel's height reserved above the data for
# the legend (applied to both the R^2 and the slope axis so they stay aligned).
# 0.22 because the NCP legend runs to two rows (weighted + unweighted entries).
HEADROOM = 0.22

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


def build_climatology(df: pd.DataFrame):
    """Per-calendar-month inter-annual statistics -> periodic cubic splines.

    Factored out of ``build_synthetic_year`` so the same climatology can be
    re-evaluated at arbitrary days-of-year by the seasonal-reproducibility
    regression (including on year-subsets, for the leave-one-out norm).

    Returns (cs_mean, cs_p10, cs_p90, clim_df).
    """
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

    # bc_type="periodic" also makes evaluation outside the knot range wrap
    # periodically, so early-January DOYs (< Jan 15) are interpolated, not
    # extrapolated.
    cs_m   = CubicSpline(doy_pts, mean_pts, bc_type="periodic")
    cs_p10 = CubicSpline(doy_pts, p10_pts,  bc_type="periodic")
    cs_p90 = CubicSpline(doy_pts, p90_pts,  bc_type="periodic")

    return cs_m, cs_p10, cs_p90, clim


def build_synthetic_year(df: pd.DataFrame, step_days: int = 7):
    """Inter-annual climatology per calendar month -> periodic weekly spline
    (Dec wraps smoothly into Jan, appropriate for a signed, full-year series).
    Returns (doy_dense, y_mean, y_p10, y_p90, clim_df)."""
    cs_m, cs_p10, cs_p90, clim = build_climatology(df)
    doy_dense = np.linspace(1, 365, int(365 / step_days) + 1)
    return doy_dense, cs_m(doy_dense), cs_p10(doy_dense), cs_p90(doy_dense), clim


def compute_reproducibility(df: pd.DataFrame, leave_one_out: bool = True,
                            weighted: bool = True):
    """Per-year regression of the 20-day NCP bins on the seasonal climatology.

    With ``leave_one_out`` (default) the climatology used as the predictor for
    a given year is rebuilt from the *other* years only, so each year is scored
    against a norm it did not help define.

    With ``weighted`` (default) the fit is inverse-variance weighted by the
    Monte Carlo ``ncp_sd``. This is essential for this record: sigma is
    2-5 mmol C m-2 d-1 in the productive season but 15-32 in winter, so
    unweighted the score is dominated by a handful of winter bins that are not
    significantly different from zero (e.g. 2023-01-20 = +32.4 +/- 21.8, a
    positive January NCP, which is not physically possible). Unweighted those
    bins drove 2022/2023 to R^2 = 0.49 / 0.21; weighted they are 0.85 / 0.88,
    i.e. the apparent anomaly was a precision artefact, not an anomalous year.
    ``weighted=False`` reproduces the unweighted sensitivity test.

    See ``figure_config.per_year_climatology_regression`` for the columns.
    """
    def clim_fn(year, doy):
        ref = df[df["date"].dt.year != year] if leave_one_out else df
        cs_m, *_ = build_climatology(ref)
        return cs_m(doy)

    return per_year_climatology_regression(
        df["date"], df["ncp_mean"], clim_fn,
        sigma=df["ncp_sd"] if weighted else None,
    )


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


def draw_panel_b(ax, df: pd.DataFrame, stats: pd.DataFrame | None = None,
                 panel_label: str = "(b)"):
    doy, y_m, y_p10, y_p90, clim = build_synthetic_year(df)

    cmap = plt.get_cmap(R2_CMAP)
    norm = mpl.colors.Normalize(vmin=R2_VMIN, vmax=R2_VMAX)
    r2_by_year = ({} if stats is None
                  else dict(zip(stats["year"], stats["r2"])))

    years = sorted(df["date"].dt.year.unique())
    ymins, ymaxs = [], []
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
        ymins.append(np.nanmin(y))
        ymaxs.append(np.nanmax(y))

    ax.axhline(0, color=ZERO_LINE, lw=0.8, ls=":", zorder=1)

    ax.fill_between(doy, y_p10, y_p90,
                    color=BLUE_LIGHT, linewidth=0, alpha=0.55,
                    label="p10–p90 (inter-annual)", zorder=3)

    # White halo keeps the climatology readable over the coloured year lines
    # (the dark end of the R^2 colour ramp is close to the NCP navy).
    ax.plot(doy, y_m, lw=2.4, ls="-", color=BLUE,
            label="Mean NCP (nitrate budget)", zorder=5,
            path_effects=[pe.Stroke(linewidth=4.2, foreground="white"),
                          pe.Normal()])

    ax.scatter(clim["doy"], clim["clim_mean"],
               s=32, color=BLUE, zorder=6, linewidths=0.6,
               edgecolors="white")

    mid_doys = [pd.Timestamp(2021, m, 15).timetuple().tm_yday
                for m in range(1, 13)]
    ax.set_xticks(mid_doys)
    ax.set_xticklabels(MONTH_LABELS, fontsize=9)
    ax.set_xlim(1, 365)
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))

    # y-limits sized to contain every individual-year spline (Task 4: the
    # grey year lines used to run off the top of the panel).
    lo = min(ymins + [np.nanmin(y_p10)])
    hi = max(ymaxs + [np.nanmax(y_p90)])
    pad = 0.08 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)

    ax.set_ylabel("NCP (mmol C m$^{-2}$ d$^{-1}$)", fontsize=11)
    ax.set_xlabel("Month", fontsize=11)

    ax.text(0.98, 0.97, f"n = {len(years)} yr  ({years[0]}–{years[-1]})",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, color="0.45")

    handles, labels = ax.get_legend_handles_labels()
    if stats is not None:
        handles.append(Line2D([], [], color=GREY_LINE, lw=1.0))
        labels.append("Individual years (colour = R$^2$ vs climatology)")
    ax.legend(handles, labels, loc="upper left", ncol=1, fontsize=9)
    ax.text(-0.06, 1.03, panel_label, transform=ax.transAxes,
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
           label="R$^2$ vs climatology (1/$\\sigma^2$ weighted)")

    # Unweighted sensitivity, if carried. Drawn as a dash spanning the bar
    # width (not a free-floating marker) so it unambiguously reads as "this
    # bar's other value" and cannot be confused with the slope markers, which
    # live on the other axis.
    if "r2_unweighted" in stats.columns:
        ax.hlines(stats["r2_unweighted"], x - 0.34, x + 0.34,
                  color="0.15", lw=1.8, zorder=5,
                  label="R$^2$ unweighted (sensitivity)")

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
    # this axis (markers, reference line, ticks, label) is drawn in BLUE so the
    # axis a symbol belongs to is readable without consulting the legend.
    ax2 = ax.twinx()
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_color(BLUE)
    ax2.axhline(1.0, color=BLUE, lw=0.8, ls=":", alpha=0.6, zorder=2)
    ax2.plot(x, stats["slope"], ls="none", marker="D", ms=5.5,
             color=BLUE, mec="white", mew=0.7, zorder=6,
             label="Slope (amplitude vs norm)")
    smin, smax = float(stats["slope"].min()), float(stats["slope"].max())
    spread = max(smax - smin, 0.2)
    lo2 = min(smin, 1.0) - 0.25 * spread
    hi2 = max(smax, 1.0) + 0.25 * spread
    ax2.set_ylim(lo2, lo2 + (hi2 - lo2) / (1.0 - HEADROOM))
    ax2.set_ylabel("Regression slope", fontsize=11, color=BLUE)
    ax2.tick_params(axis="y", labelsize=9, colors=BLUE)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper center", ncol=2, fontsize=8.5,
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
    df = load_data()

    stats = compute_reproducibility(df, weighted=True)

    # Unweighted sensitivity test, carried in the CSV so the effect of the
    # noisy winter bins stays visible rather than hidden (see docstring of
    # compute_reproducibility).
    unw = compute_reproducibility(df, weighted=False)
    stats = stats.merge(
        unw[["year", "slope", "r2"]].rename(columns={"slope": "slope_unweighted",
                                                     "r2": "r2_unweighted"}),
        on="year", how="left",
    )

    OUT_STATS.parent.mkdir(parents=True, exist_ok=True)
    stats.to_csv(OUT_STATS, index=False)
    print(f"Wrote {OUT_STATS}")
    print(stats.round(3).to_string(index=False))
    print(f"  mean R2 (1/sigma^2 weighted) = {stats['r2'].mean():.3f}")
    print(f"  mean R2 (unweighted)         = {stats['r2_unweighted'].mean():.3f}")

    yr_min = int(df["date"].dt.year.min())
    yr_max = int(df["date"].dt.year.max())

    # ---- 3-panel combined figure ----
    fig, axes = plt.subplots(3, 1, figsize=(10, 11.5), dpi=300,
                             gridspec_kw={"hspace": 0.42,
                                          "height_ratios": [1, 1, 0.72]})

    draw_panel_a(axes[0], df)
    draw_panel_b(axes[1], df, stats)
    draw_panel_c(axes[2], stats)

    fig.suptitle(
        "Basin-scale NCP — Iceland Basin & Irminger Sea (pooled)\n"
        f"BGC-Argo nitrate-drawdown budget, CANYON-B ({yr_min}–{yr_max})",
        fontsize=11, y=0.975,
    )
    fig.text(
        0.5, 0.055,
        "(c) Each year's 20-day bins regressed on a leave-one-out monthly "
        "climatology at the same day-of-year, weighted by 1/$\\sigma^2$ "
        "(Monte Carlo).\n"
        "R$^2$ = variance explained by the climatological seasonal shape and "
        "phase; slope = amplitude relative to the norm. Open circles show the "
        "unweighted R$^2$: the\n"
        "difference is set by winter bins whose $\\sigma$ (15–32) is as large "
        "as the signal, not by anomalous seasonal cycles.",
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
        "Synthetic-year NCP — Iceland Basin & Irminger Sea (pooled)\n"
        f"BGC-Argo nitrate budget (weekly climatology, {yr_min}–{yr_max})",
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
        "Seasonal reproducibility of basin NCP\n"
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
