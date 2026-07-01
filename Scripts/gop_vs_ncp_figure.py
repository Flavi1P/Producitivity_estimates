"""
Publication-quality figure: Gross Oxygen Production (GOP) from the
high-frequency profiling float 3902681, alongside NCP from that same float's
own nitrate-drawdown budget (not the pooled basin product). Mirrors the
house style of npp_timeseries_figure.py and ncp_timeseries_figure.py (Arial
rcParams, spline-smoothed lines, consistent colour scheme).

GOP is now integrated over 0-200 m rather than the fixed 0-40 m slab used in
the original version of this figure. Rationale (see
Scripts/o2_budget_3902681_multidepth.R, which recomputes the diel O2 budget
integrated to several fixed candidate depths -- 40/100/150/200 m -- using the
same profile ingest and diffusive-only Wanninkhof 2014 air-sea correction as
the parent 0-40 m budget, and Scripts/gop_depth_sensitivity_figure.py, which
produces the comparison below):
  - 0-40 m under-samples the mixed layer most of the year: over the plot
    window the median MLD is ~88 m and MLD <= 40 m only ~35% of the time, so
    a 40 m slab misses a large fraction of the productive layer's O2
    inventory outside summer.
  - The decisive test is PHASE alignment with NCP's spring bloom rise, not
    whole-record correlation or noise level (both roughly flat across depth
    and dominated by the physiologically uninteresting deep-winter months).
    Cross-correlating smoothed GOP against NCP restricted to the spring
    bloom window (robust to the exact window: Feb-Jun, Jan-Jul, Mar-May all
    agree) shows GOP visibly LAGGING NCP's April rise at shallow depths --
    by +15 d at 40 m, +13 d at 100 m, +9 d at 150 m -- collapsing to 0-day
    lag at 200 m, where zero-lag r also peaks (0.76 vs. 0.34-0.62 at
    40/100/150 m). This is physically expected: in March-April the mixed
    layer here is frequently deeper than 40-150 m, so a shallow O2 budget
    only "sees" the production signal once the mixed layer shoals later in
    the season -- it is an integration-depth artifact, not a real
    biological lag. 200 m is deep enough to capture the O2 inventory change
    in near real time with the nitrate-budget NCP during the bloom onset,
    which is the comparison this figure exists to make.
  - The trade-off is a noisier smoothed line outside the bloom (larger IQR
    envelope, occasional large winter excursions) -- an acceptable cost for
    getting the bloom-onset timing right. See
  Output/float_3902681_GOP_depth_sensitivity.png for the 4-panel comparison
  (correlation, SNR, and bloom-window lag are annotated per depth).

GOP (mat_and_meth.md S2.2.1 describes the ORIGINAL 0-40 m convention; this
figure now uses the 0-200 m variant above): the flux-corrected net diel O2
change and the flux-corrected nighttime O2 loss (Wanninkhof 2014 diffusive
air-sea correction) are each smoothed with an 18-day centered running window,
then summed:
    GOP = smooth18(net diel change, flux-corr) + smooth18(nighttime loss, flux-corr)
The 18-day window statistic is a median (with the interquartile range shown
as the envelope) rather than a mean/SD: the raw per-cycle diel rate is noisy
and has occasional extreme outliers that a mean would not be robust to.
GOP stays in O2 units (mmol O2 m-2 d-1) -- no O2->C conversion (per
o2_budget_3902681_explained.md, this script's convention is deliberately
O2-only for the float budget).

NCP: float 3902681's own nitrate-drawdown NCP (Data/Processed/
ncp_float_3902681_30d.xlsx), a single-float 30-day-window product -- distinct
from the pooled Iceland Basin + Irminger Sea basin-scale product (S2.2.4)
used in an earlier version of this figure. This product carries no
per-estimate uncertainty, so only the smoothed line (no ribbon) is shown.

Because GOP (O2) and NCP (C) are different currencies, NCP is shown on its
own twin y-axis (both symmetric about zero so the zero-lines coincide)
rather than assuming an unstated photosynthetic quotient.

Plot window: 2024-11-01 to 2025-11-01.

Input:  Data/Processed/o2_budget_3902681_multidepth.xlsx (depth_m == 200)
        Data/Processed/ncp_float_3902681_30d.xlsx
Output: Output/float_3902681_GOP_vs_basin_NCP.png
        Output/float_3902681_GOP_vs_basin_NCP.pdf
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
REPO_ROOT   = Path(__file__).resolve().parents[1]
O2_XLSX     = REPO_ROOT / "Data" / "Processed" / "o2_budget_3902681_multidepth.xlsx"
FLOAT_NCP_XLSX = REPO_ROOT / "Data" / "Processed" / "ncp_float_3902681_30d.xlsx"
OUT_PNG     = REPO_ROOT / "Output" / "float_3902681_GOP_vs_basin_NCP.png"
OUT_PDF     = REPO_ROOT / "Output" / "float_3902681_GOP_vs_basin_NCP.pdf"

# Chosen from the 40/100/150/200 m depth-sensitivity test described in the
# module docstring (Scripts/o2_budget_3902681_multidepth.R produces all four):
# 200 m is the only depth at which GOP's spring rise is in phase (0-day lag)
# with the nitrate-budget NCP's spring rise.
INTEGRATION_DEPTH_M = 200

WINDOW_DAYS = 18     # per mat_and_meth.md S2.2.1
MIN_N       = 3      # minimum raw points required in a smoothing window

PLOT_START  = pd.Timestamp("2024-11-01")
PLOT_END    = pd.Timestamp("2025-11-01")

# ---------------------------------------------------------------------------
# Style (matches npp_timeseries_figure.py / ncp_timeseries_figure.py)
# ---------------------------------------------------------------------------
GOP_COLOR   = "#c1440e"
GOP_LIGHT   = "#f0c4a8"
NCP_COLOR   = "#1f4e79"
NCP_LIGHT   = "#b8cbe0"
ZERO_LINE   = "#888888"

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
    "legend.frameon": False,
    "legend.fontsize":  9,
})


# ---------------------------------------------------------------------------
# GOP: centered time-window smoothing of the flux-corrected O2 budget
# ---------------------------------------------------------------------------

def centered_window_stats(times: np.ndarray, values: np.ndarray,
                           grid: np.ndarray, window_days: float):
    """Centered time-window median/IQR/count of `values` at each `grid` time,
    over all `times` within +/- window_days/2. Median/IQR (rather than
    mean/SD) because the raw per-cycle diel rate is noisy and occasionally
    has extreme outliers (e.g. a single night pair at -1360 mmol O2 m-2 d-1
    in May 2025) that would otherwise dominate the smoothed estimate."""
    half = np.timedelta64(int(window_days / 2 * 86400), "s")
    n = len(grid)
    med = np.full(n, np.nan)
    p25 = np.full(n, np.nan)
    p75 = np.full(n, np.nan)
    count = np.zeros(n, dtype=int)
    for i, t in enumerate(grid):
        mask = (times >= t - half) & (times <= t + half)
        v = values[mask]
        count[i] = v.size
        if v.size >= MIN_N:
            med[i] = np.median(v)
            p25[i], p75[i] = np.percentile(v, [25, 75])
    return med, p25, p75, count


def load_gop() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_excel(O2_XLSX)
    df["mtime"] = pd.to_datetime(df["mtime"], utc=True).dt.tz_localize(None)
    df = df[df["depth_m"] == INTEGRATION_DEPTH_M]

    net = df.loc[df["type"] == "net", ["mtime", "rate_corr_mmol_o2_m2_d"]].dropna()
    net = net.sort_values("mtime")

    night = df.loc[df["type"] == "night", ["mtime", "rate_corr_mmol_o2_m2_d"]].dropna()
    night = night.sort_values("mtime")
    night_resp = -night["rate_corr_mmol_o2_m2_d"].to_numpy()   # +ve = O2 consumed (respiration)

    t0 = min(net["mtime"].min(), night["mtime"].min())
    t1 = max(net["mtime"].max(), night["mtime"].max())
    grid = pd.date_range(t0, t1, freq="1D").to_numpy()

    net_m, net_p25, net_p75, net_n = centered_window_stats(
        net["mtime"].to_numpy(), net["rate_corr_mmol_o2_m2_d"].to_numpy(),
        grid, WINDOW_DAYS)
    resp_m, resp_p25, resp_p75, resp_n = centered_window_stats(
        night["mtime"].to_numpy(), night_resp, grid, WINDOW_DAYS)

    ok = (net_n >= MIN_N) & (resp_n >= MIN_N)
    gop = pd.DataFrame({
        "date":      pd.to_datetime(grid),
        "gop_mean":  net_m + resp_m,          # sum of the two window medians
        "gop_lo":    net_p25 + resp_p25,       # approx. envelope (IQR sum)
        "gop_hi":    net_p75 + resp_p75,
    })[ok].reset_index(drop=True)

    print(f"[GOP] 0-{INTEGRATION_DEPTH_M} m, {WINDOW_DAYS}-day centered smoothing of "
          f"flux-corrected net diel change ({len(net)} pairs) + nighttime loss "
          f"({len(night)} pairs) -> {len(gop)} daily GOP estimates "
          f"({gop['date'].min().date()} to {gop['date'].max().date()})")

    return gop, net


# ---------------------------------------------------------------------------
# Float 3902681's own nitrate-drawdown NCP (single-float product, 30-day
# window; distinct from the pooled Iceland Basin + Irminger Sea basin-scale
# product used in the earlier version of this figure)
# ---------------------------------------------------------------------------

def load_float_ncp(t0: pd.Timestamp, t1: pd.Timestamp) -> pd.DataFrame:
    df = pd.read_excel(FLOAT_NCP_XLSX)
    df["date"] = pd.to_datetime(df["date"])
    # collapse duplicate-date profile pairs (same 30-day NCP estimate) to one row/date
    df = df.groupby("date", as_index=False)["NCP"].mean().sort_values("date")

    pad = pd.Timedelta(days=15)
    sub = df[(df["date"] >= t0 - pad) & (df["date"] <= t1 + pad)].copy()
    print(f"[NCP] float 3902681 nitrate-drawdown budget (30-day window): "
          f"{len(sub)} estimates overlapping the plot window "
          f"({sub['date'].min().date() if len(sub) else 'n/a'} to "
          f"{sub['date'].max().date() if len(sub) else 'n/a'})")
    return sub


def spline_ncp(sub: pd.DataFrame, step_days: int = 3):
    t0 = sub["date"].min()
    t_days = (sub["date"] - t0).dt.days.to_numpy(dtype=float)
    y_m = sub["NCP"].to_numpy()

    cs_m = CubicSpline(t_days, y_m)

    t_dense = np.arange(0, t_days.max() + step_days, step_days)
    dates = t0 + pd.to_timedelta(t_dense, unit="D")
    return dates, cs_m(t_dense)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    gop, net_raw = load_gop()

    # ---- restrict to the requested plot window ----
    gop = gop[(gop["date"] >= PLOT_START) & (gop["date"] <= PLOT_END)].reset_index(drop=True)
    net_raw = net_raw[(net_raw["mtime"] >= PLOT_START) & (net_raw["mtime"] <= PLOT_END)]

    ncp_sub = load_float_ncp(PLOT_START, PLOT_END)
    ncp_dates, ncp_m = spline_ncp(ncp_sub)
    ncp_mask = (ncp_dates >= PLOT_START) & (ncp_dates <= PLOT_END)
    ncp_dates, ncp_m = (a[ncp_mask] for a in (ncp_dates, ncp_m))
    ncp_sub_plot = ncp_sub[(ncp_sub["date"] >= PLOT_START) & (ncp_sub["date"] <= PLOT_END)]

    fig, ax1 = plt.subplots(figsize=(11, 6), dpi=300)
    ax2 = ax1.twinx()

    # ---- symmetric y-limits so both zero-lines coincide at plot centre ----
    # O2 axis (GOP) is set from the actual IQR envelope range (not a hardcoded
    # clip): at the chosen 0-200 m integration depth the envelope is far wider
    # than the 0-40 m case (see module docstring), so a fixed +/-200/400 clip
    # would truncate real winter excursions.
    m1 = np.nanmax(np.abs(np.concatenate([gop["gop_lo"], gop["gop_hi"]]))) * 1.15
    m2 = np.nanmax(np.abs(ncp_m)) * 1.15
    ax1.set_ylim(-m1, m1)
    ax2.set_ylim(-m2, m2)
    ax1.axhline(0, color=ZERO_LINE, lw=0.8, ls=":", zorder=1)

    # ---- faint raw net diel change (context / noise level), clipped to the
    # display range so isolated outliers don't stretch the axis further ----
    ax1.plot(net_raw["mtime"],
              net_raw["rate_corr_mmol_o2_m2_d"].clip(-m1, m1),
              lw=0.5, color=GOP_COLOR, alpha=0.18, zorder=2)

    # ---- GOP (float O2 budget) ----
    ax1.fill_between(gop["date"], gop["gop_lo"], gop["gop_hi"],
                      color=GOP_LIGHT, linewidth=0, alpha=0.7, zorder=3,
                      label="GOP IQR (18-day window)")
    ax1.plot(gop["date"], gop["gop_mean"], lw=2.0, color=GOP_COLOR, zorder=5,
              label="GOP (flux-corrected O2 budget, float 3902681)")

    # ---- NCP (float 3902681's own nitrate-drawdown budget) ----
    ax2.plot(ncp_dates, ncp_m, lw=2.0, color=NCP_COLOR, zorder=4,
              label="NCP (float 3902681 nitrate budget, 30-day window)")
    ax2.scatter(ncp_sub_plot["date"], ncp_sub_plot["NCP"], s=16,
                color=NCP_COLOR, zorder=6, linewidths=0.5, edgecolors="white")

    # ---- axes cosmetics ----
    ax1.set_xlim(PLOT_START, PLOT_END)
    ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax1.xaxis.set_minor_locator(mdates.MonthLocator())

    ax1.set_ylabel(r"GOP, 0-200 m (mmol O$_2$ m$^{-2}$ d$^{-1}$)", fontsize=11,
                    color=GOP_COLOR)
    ax2.set_ylabel(r"NCP (mmol C m$^{-2}$ d$^{-1}$)", fontsize=11,
                    color=NCP_COLOR)
    ax1.tick_params(axis="y", labelcolor=GOP_COLOR)
    ax2.tick_params(axis="y", labelcolor=NCP_COLOR)
    ax2.spines["top"].set_visible(False)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8.5, ncol=1)

    fig.suptitle(
        "Float 3902681 (Iceland Basin, Nov 2024–Nov 2025)\n"
        "Gross Oxygen Production (0-200 m) vs. the float's own nitrate-budget NCP",
        fontsize=12, y=0.98,
    )
    ax1.set_title(
        "GOP = flux-corrected net diel O$_2$ change + flux-corrected nighttime "
        "loss (Wanninkhof 2014), 0-200 m integration, 18-day running median",
        fontsize=9, color="0.35", pad=10,
    )

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, bbox_inches="tight", dpi=300)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
