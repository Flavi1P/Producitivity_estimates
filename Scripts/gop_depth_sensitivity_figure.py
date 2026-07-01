"""
Depth-sensitivity diagnostic (not the main figure): compares GOP from float
3902681, integrated to each of four fixed candidate depths (40/100/150/200 m,
Data/Processed/o2_budget_3902681_multidepth.xlsx), against the float's own
nitrate-drawdown NCP. Produced to justify the 0-200 m integration depth used
in Scripts/gop_vs_ncp_figure.py.

Two diagnostics are reported per depth:
  - whole-record Pearson r and smoothed-GOP signal-to-noise ratio (18-day
    median std / median IQR envelope width) -- both are roughly flat across
    depth and dominated by the physiologically uninteresting deep-winter
    months, so neither discriminates well between candidates.
  - the DECISIVE one: cross-correlation lag between smoothed GOP and NCP
    restricted to the spring bloom window (Feb-Jun), reporting the best-fit
    lag in days (positive = GOP lags NCP) and the zero-lag correlation. This
    lag shrinks monotonically with integration depth (+15, +13, +9, +0 days
    at 40/100/150/200 m) because in March-April the mixed layer here is
    frequently deeper than 40-150 m, so shallow integration only picks up
    the production signal once the mixed layer shoals later in the season.
    200 m is the only candidate with 0-day lag (and the highest zero-lag r,
    0.76) -- i.e. the only depth where GOP's spring rise is in phase with
    NCP's spring rise, which is the comparison this figure exists to make.

Input:  Data/Processed/o2_budget_3902681_multidepth.xlsx
        Data/Processed/ncp_float_3902681_30d.xlsx
Output: Output/float_3902681_GOP_depth_sensitivity.png
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.interpolate import CubicSpline

REPO_ROOT = Path(__file__).resolve().parents[1]
MD_XLSX   = REPO_ROOT / "Data" / "Processed" / "o2_budget_3902681_multidepth.xlsx"
NCP_XLSX  = REPO_ROOT / "Data" / "Processed" / "ncp_float_3902681_30d.xlsx"
OUT_PNG   = REPO_ROOT / "Output" / "float_3902681_GOP_depth_sensitivity.png"

WINDOW_DAYS = 18
MIN_N = 3
PLOT_START = pd.Timestamp("2024-11-01")
PLOT_END = pd.Timestamp("2025-11-01")

# Spring bloom onset window used for the lag test (the decisive diagnostic --
# see module docstring). Result is not sensitive to the exact bounds: Jan-Jul
# and Mar-May give the same qualitative ranking.
BLOOM_START = pd.Timestamp("2025-02-15")
BLOOM_END   = pd.Timestamp("2025-06-15")


def bloom_lag(gop: pd.DataFrame, ncp_dates: np.ndarray, ncp_vals: np.ndarray):
    """Best-fit lag (in whole days) and zero-lag Pearson r between smoothed
    GOP and NCP, restricted to the bloom window. Positive lag = GOP lags
    NCP. Both series are resampled onto a common 1-day grid first so that
    one index step in the cross-correlation is exactly one day (the NCP
    spline elsewhere in this script is sampled every 3 days, which would
    otherwise silently rescale the lag axis)."""
    day_grid = mdates.date2num(pd.date_range(BLOOM_START, BLOOM_END, freq="1D"))

    ncp_num = mdates.date2num(pd.DatetimeIndex(ncp_dates))
    ncp_1d = np.interp(day_grid, ncp_num, ncp_vals)

    gmask = (gop["date"] >= BLOOM_START - pd.Timedelta(days=30)) & \
            (gop["date"] <= BLOOM_END + pd.Timedelta(days=30))
    gsub = gop[gmask]
    gop_1d = np.interp(day_grid, mdates.date2num(gsub["date"]), gsub["gop_mean"])

    a = (gop_1d - gop_1d.mean()) / gop_1d.std()
    b = (ncp_1d - ncp_1d.mean()) / ncp_1d.std()
    xcorr = np.correlate(a, b, mode="full") / len(a)
    lags = np.arange(-len(a) + 1, len(a))          # now in whole days
    best_lag = int(lags[np.argmax(xcorr)])
    r0 = np.corrcoef(a, b)[0, 1]
    return best_lag, r0


def centered_window_stats(times, values, grid, window_days):
    half = np.timedelta64(int(window_days / 2 * 86400), "s")
    n = len(grid)
    med = np.full(n, np.nan); p25 = np.full(n, np.nan); p75 = np.full(n, np.nan)
    count = np.zeros(n, dtype=int)
    for i, t in enumerate(grid):
        mask = (times >= t - half) & (times <= t + half)
        v = values[mask]
        count[i] = v.size
        if v.size >= MIN_N:
            med[i] = np.median(v)
            p25[i], p75[i] = np.percentile(v, [25, 75])
    return med, p25, p75, count


def main() -> int:
    df = pd.read_excel(MD_XLSX)
    df["mtime"] = pd.to_datetime(df["mtime"], utc=True).dt.tz_localize(None)

    ncp = pd.read_excel(NCP_XLSX)
    ncp["date"] = pd.to_datetime(ncp["date"])
    ncp = ncp.groupby("date", as_index=False)["NCP"].mean().sort_values("date")
    pad = pd.Timedelta(days=15)
    ncp_sub = ncp[(ncp["date"] >= PLOT_START - pad) & (ncp["date"] <= PLOT_END + pad)].copy()
    t0 = ncp_sub["date"].min()
    t_days = (ncp_sub["date"] - t0).dt.days.to_numpy(dtype=float)
    cs = CubicSpline(t_days, ncp_sub["NCP"].to_numpy())
    t_dense = np.arange(0, t_days.max() + 3, 3)
    ncp_dates = t0 + pd.to_timedelta(t_dense, unit="D")
    ncp_vals = cs(t_dense)
    mask = (ncp_dates >= PLOT_START) & (ncp_dates <= PLOT_END)
    ncp_dates, ncp_vals = ncp_dates[mask], ncp_vals[mask]

    depths = sorted(df["depth_m"].unique())
    fig, axes = plt.subplots(len(depths), 1, figsize=(11, 3 * len(depths)),
                              sharex=True, dpi=150)

    print(f"{'depth_m':>8} {'corr(GOP,NCP)':>14} {'median_IQR_width':>17} "
          f"{'signal_std':>11} {'SNR':>6} {'bloom_lag_d':>12} {'bloom_r0':>9}")
    for ax, z in zip(axes, depths):
        sub = df[df["depth_m"] == z]
        net = sub[sub["type"] == "net"].sort_values("mtime")
        night = sub[sub["type"] == "night"].sort_values("mtime")
        night_resp = -night["rate_corr_mmol_o2_m2_d"].to_numpy()

        tt0 = min(net["mtime"].min(), night["mtime"].min())
        tt1 = max(net["mtime"].max(), night["mtime"].max())
        grid = pd.date_range(tt0, tt1, freq="1D").to_numpy()

        net_m, net_p25, net_p75, net_n = centered_window_stats(
            net["mtime"].to_numpy(), net["rate_corr_mmol_o2_m2_d"].to_numpy(),
            grid, WINDOW_DAYS)
        resp_m, resp_p25, resp_p75, resp_n = centered_window_stats(
            night["mtime"].to_numpy(), night_resp, grid, WINDOW_DAYS)

        ok = (net_n >= MIN_N) & (resp_n >= MIN_N)
        gop = pd.DataFrame({
            "date": pd.to_datetime(grid),
            "gop_mean": net_m + resp_m,
            "gop_lo": net_p25 + resp_p25,
            "gop_hi": net_p75 + resp_p75,
        })[ok].reset_index(drop=True)
        gop = gop[(gop["date"] >= PLOT_START) & (gop["date"] <= PLOT_END)]

        ax2 = ax.twinx()
        m1 = np.nanmax(np.abs(np.concatenate([gop["gop_lo"], gop["gop_hi"]]))) * 1.15
        m2 = np.nanmax(np.abs(ncp_vals)) * 1.15
        ax.set_ylim(-m1, m1); ax2.set_ylim(-m2, m2)
        ax.axvspan(BLOOM_START, BLOOM_END, color="#238b45", alpha=0.08, lw=0, zorder=0)
        ax.axhline(0, color="grey", lw=0.8, ls=":")
        ax.fill_between(gop["date"], gop["gop_lo"], gop["gop_hi"],
                         color="#f0c4a8", alpha=0.7, lw=0)
        ax.plot(gop["date"], gop["gop_mean"], color="#c1440e", lw=2)
        ax2.plot(ncp_dates, ncp_vals, color="#1f4e79", lw=2)

        gop_interp = np.interp(mdates.date2num(ncp_dates),
                                mdates.date2num(gop["date"]), gop["gop_mean"])
        r = np.corrcoef(gop_interp, ncp_vals)[0, 1]
        width = gop["gop_hi"] - gop["gop_lo"]
        snr = gop["gop_mean"].std() / width.median()
        best_lag, r0 = bloom_lag(gop, ncp_dates, ncp_vals)

        print(f"{z:8d} {r:14.3f} {width.median():17.1f} "
              f"{gop['gop_mean'].std():11.1f} {snr:6.3f} {best_lag:12d} {r0:9.3f}")

        ax.set_title(
            f"Integration depth = {z} m   |   whole-record r = {r:.2f}, SNR = {snr:.2f}   |   "
            f"bloom window (shaded): lag = {best_lag:+d} d, r$_0$ = {r0:.2f}",
            fontsize=10)
        ax.set_ylabel("GOP (mmol O$_2$ m$^{-2}$ d$^{-1}$)", color="#c1440e", fontsize=8)
        ax2.set_ylabel("NCP (mmol C m$^{-2}$ d$^{-1}$)", color="#1f4e79", fontsize=8)

    axes[-1].xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    fig.suptitle(
        "Float 3902681: GOP depth-sensitivity test (orange, left axis) vs. "
        "nitrate-budget NCP (blue, right axis)", fontsize=12, y=1.0)
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    print(f"Wrote {OUT_PNG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
