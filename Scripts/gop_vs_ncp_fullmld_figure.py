"""
Variant of Scripts/gop_vs_ncp_movingmld_figure.py using an integration depth that
tracks the FULL mixed-layer depth (MLD), floored at 50 m with NO upper cap,
instead of the [50, 200] m clamp of the sibling moving-MLD figure.

  z_bot(pair) = max(50, max(MLD_i, MLD_j))

Rationale: the 200 m cap in gop_vs_ncp_movingmld_figure.py keeps the smoothed
summer line clean but truncates the deep-convective winter/spring mixed layer, so
the O2 inventory change of the ventilating layer is under-sampled exactly when the
ML is deepest. Removing the cap integrates over the actual mixed layer at every
pair (matching the parent MLD budget in Scripts/o2_budget_3902681.R, which floors
at 40 m with no cap). The 50 m floor is retained to avoid integrating over an
unrealistically thin summer ML. See Scripts/o2_budget_3902681_fullmld.R for the
O2-budget recomputation; everything else (18-day centered median smoothing, IQR
envelope, float 3902681's own nitrate-drawdown NCP on a twin axis) is unchanged
from gop_vs_ncp_movingmld_figure.py.

Input:  Data/Processed/o2_budget_3902681_fullmld.csv
        Data/Processed/ncp_float_3902681_30d.xlsx
Output: Output/float_3902681_GOP_fullMLD_vs_basin_NCP.png / .pdf
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
REPO_ROOT      = Path(__file__).resolve().parents[1]
O2_CSV         = REPO_ROOT / "Data" / "Processed" / "o2_budget_3902681_fullmld.csv"
FLOAT_NCP_XLSX = REPO_ROOT / "Data" / "Processed" / "ncp_float_3902681_30d.xlsx"
NPP_CSV        = REPO_ROOT / "Data" / "Processed" / "icb_npp_2025.csv"
OUT_PNG        = REPO_ROOT / "Output" / "float_3902681_GOP_fullMLD_vs_basin_NCP.png"
OUT_PDF        = REPO_ROOT / "Output" / "float_3902681_GOP_fullMLD_vs_basin_NCP.pdf"

MG_C_PER_MMOL = 12.0   # mg C -> mmol C (repo convention, ncp_o2budget_float_3902681.R)

MLD_FLOOR = 50    # m  (no upper cap: integration bottom follows the full MLD)

WINDOW_DAYS = 18     # per mat_and_meth.md S2.2.1
MIN_N       = 3      # minimum raw points required in a smoothing window

PLOT_START  = pd.Timestamp("2024-11-01")
PLOT_END    = pd.Timestamp("2025-11-01")

# ---------------------------------------------------------------------------
# Style (matches the parent figure)
# ---------------------------------------------------------------------------
GOP_COLOR   = "#c1440e"
GOP_LIGHT   = "#f0c4a8"
NCP_COLOR   = "#1f4e79"
NCP_LIGHT   = "#b8cbe0"
NPP_COLOR   = "#2e7d32"
ZERO_LINE   = "#888888"
MLD_COLOR   = "#5a5a5a"

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
    over all `times` within +/- window_days/2 (median/IQR for robustness to
    occasional extreme per-cycle diel outliers)."""
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


def load_gop() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(O2_CSV)
    df["mtime"] = pd.to_datetime(df["mtime"], utc=True).dt.tz_localize(None)

    net = df.loc[df["type"] == "net", ["mtime", "rate_corr_mmol_o2_m2_d"]].dropna()
    net = net.sort_values("mtime")

    night = df.loc[df["type"] == "night", ["mtime", "rate_corr_mmol_o2_m2_d"]].dropna()
    night = night.sort_values("mtime")
    night_resp = -night["rate_corr_mmol_o2_m2_d"].to_numpy()   # +ve = O2 consumed

    # moving integration depth per pair (for the context sub-panel)
    zbot = df.loc[df["type"] == "net", ["mtime", "z_bot_m"]].dropna().sort_values("mtime")

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
        "gop_mean":  net_m + resp_m,
        "gop_lo":    net_p25 + resp_p25,
        "gop_hi":    net_p75 + resp_p75,
    })[ok].reset_index(drop=True)

    print(f"[GOP] full MLD bottom (floor {MLD_FLOOR} m, no cap), {WINDOW_DAYS}-day "
          f"centered smoothing of flux-corrected net diel change ({len(net)} pairs) "
          f"+ nighttime loss ({len(night)} pairs) -> {len(gop)} daily GOP estimates "
          f"({gop['date'].min().date()} to {gop['date'].max().date()})")

    return gop, net, zbot


# ---------------------------------------------------------------------------
# Float 3902681's own nitrate-drawdown NCP (single-float, 30-day window)
# ---------------------------------------------------------------------------

def load_float_ncp(t0: pd.Timestamp, t1: pd.Timestamp) -> pd.DataFrame:
    df = pd.read_excel(FLOAT_NCP_XLSX)
    df["date"] = pd.to_datetime(df["date"])
    df = df.groupby("date", as_index=False)["NCP"].mean().sort_values("date")

    pad = pd.Timedelta(days=15)
    sub = df[(df["date"] >= t0 - pad) & (df["date"] <= t1 + pad)].copy()
    print(f"[NCP] float 3902681 nitrate-drawdown budget (30-day window): "
          f"{len(sub)} estimates overlapping the plot window")
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
# NPP (CbPM, Iceland Basin, 10-day binned, already depth-integrated / areal)
# ---------------------------------------------------------------------------

def load_npp(t0: pd.Timestamp, t1: pd.Timestamp) -> pd.DataFrame:
    df = pd.read_csv(NPP_CSV)
    df["date"] = pd.to_datetime(df["date_10day"])
    df = df.dropna(subset=["npp"]).sort_values("date")
    df["npp_mmol"] = df["npp"] / MG_C_PER_MMOL     # mg C -> mmol C m-2 d-1

    pad = pd.Timedelta(days=15)
    sub = df[(df["date"] >= t0 - pad) & (df["date"] <= t1 + pad)].copy()
    print(f"[NPP] CbPM Iceland Basin (10-day, already areal): {len(sub)} bins "
          f"overlapping the plot window; peak {sub['npp_mmol'].max():.0f} mmol C m-2 d-1")
    return sub


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    gop, net_raw, zbot = load_gop()

    gop = gop[(gop["date"] >= PLOT_START) & (gop["date"] <= PLOT_END)].reset_index(drop=True)
    net_raw = net_raw[(net_raw["mtime"] >= PLOT_START) & (net_raw["mtime"] <= PLOT_END)]
    zbot = zbot[(zbot["mtime"] >= PLOT_START) & (zbot["mtime"] <= PLOT_END)]

    ncp_sub = load_float_ncp(PLOT_START, PLOT_END)
    ncp_dates, ncp_m = spline_ncp(ncp_sub)
    ncp_mask = (ncp_dates >= PLOT_START) & (ncp_dates <= PLOT_END)
    ncp_dates, ncp_m = (a[ncp_mask] for a in (ncp_dates, ncp_m))
    ncp_sub_plot = ncp_sub[(ncp_sub["date"] >= PLOT_START) & (ncp_sub["date"] <= PLOT_END)]

    npp_sub = load_npp(PLOT_START, PLOT_END)
    npp_plot = npp_sub[(npp_sub["date"] >= PLOT_START) & (npp_sub["date"] <= PLOT_END)]

    # two stacked panels: main GOP/NCP comparison + a thin integration-depth strip
    fig, (ax1, axd) = plt.subplots(
        2, 1, figsize=(11, 6.8), dpi=300,
        gridspec_kw={"height_ratios": [4, 1], "hspace": 0.08}, sharex=True)
    ax2 = ax1.twinx()

    # ---- symmetric y-limits so both zero-lines coincide ----
    m1 = np.nanmax(np.abs(np.concatenate([gop["gop_lo"], gop["gop_hi"]]))) * 1.15
    m2 = np.nanmax(np.abs(ncp_m)) * 1.15
    ax1.set_ylim(-250, 500)
    ax2.set_ylim(-250, 500)
    ax1.axhline(0, color=ZERO_LINE, lw=0.8, ls=":", zorder=1)

    ax1.plot(net_raw["mtime"],
              net_raw["rate_corr_mmol_o2_m2_d"].clip(-m1, m1),
              lw=0.5, color=GOP_COLOR, alpha=0.18, zorder=2)

    ax1.fill_between(gop["date"], gop["gop_lo"], gop["gop_hi"],
                      color=GOP_LIGHT, linewidth=0, alpha=0.7, zorder=3,
                      label="GOP IQR (18-day window)")
    ax1.plot(gop["date"], gop["gop_mean"], lw=2.0, color=GOP_COLOR, zorder=5,
              label="GOP (flux-corrected O2 budget, full MLD, floor 50 m)")

    ax2.plot(ncp_dates, ncp_m, lw=2.0, color=NCP_COLOR, zorder=4,
              label="NCP (float 3902681 nitrate budget, 30-day window)")
    ax2.scatter(ncp_sub_plot["date"], ncp_sub_plot["NCP"], s=16,
                color=NCP_COLOR, zorder=6, linewidths=0.5, edgecolors="white")

    ax2.plot(npp_plot["date"], npp_plot["npp_mmol"], lw=1.8, color=NPP_COLOR,
              marker="o", ms=3.5, mfc="white", mec=NPP_COLOR, mew=0.8, zorder=4,
              label="NPP (CbPM, Iceland Basin, 10-day)")

    # ---- integration-depth strip (per-pair moving z_bot, auto-scaled, no cap) ----
    axd.plot(zbot["mtime"], zbot["z_bot_m"], lw=0.8, color=MLD_COLOR, alpha=0.9)
    axd.axhline(MLD_FLOOR, color=ZERO_LINE, lw=0.6, ls="--")
    zbot_max = float(zbot["z_bot_m"].max()) if len(zbot) else MLD_FLOOR
    axd.set_ylim(zbot_max + 25, MLD_FLOOR - 15)   # depth increases downward, auto to deepest
    axd.set_ylabel("z$_{int}$ (m)", fontsize=9)
    axd.tick_params(axis="y", labelsize=8)

    # ---- axes cosmetics ----
    ax1.set_xlim(PLOT_START, PLOT_END)
    axd.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    axd.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    axd.xaxis.set_minor_locator(mdates.MonthLocator())

    ax1.set_ylabel(r"GOP, full MLD (mmol O$_2$ m$^{-2}$ d$^{-1}$)", fontsize=11,
                    color=GOP_COLOR)
    ax2.set_ylabel(r"NCP / NPP (mmol C m$^{-2}$ d$^{-1}$)", fontsize=11, color=NCP_COLOR)
    ax1.tick_params(axis="y", labelcolor=GOP_COLOR)
    ax2.tick_params(axis="y", labelcolor=NCP_COLOR)
    ax2.spines["top"].set_visible(False)
    axd.spines["top"].set_visible(False)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8.5, ncol=1)

    fig.suptitle(
        "Float 3902681 (Iceland Basin, Nov 2024-Nov 2025)\n"
        "Gross Oxygen Production (full MLD, floor 50 m) vs. the float's own nitrate-budget NCP",
        fontsize=12, y=0.99,
    )
    ax1.set_title(
        "GOP = flux-corrected net diel O$_2$ change + flux-corrected nighttime "
        "loss (Wanninkhof 2014); integration bottom = MLD floored at 50 m (no cap), "
        "18-day running median",
        fontsize=9, color="0.35", pad=10,
    )

    fig.tight_layout(rect=(0, 0, 1, 0.93))
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, bbox_inches="tight", dpi=300)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
