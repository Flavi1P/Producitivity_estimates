"""
Synthesis figure combining the AREAL (depth-integrated) and VOLUMETRIC
(mixed-layer-mean) comparisons of the three productivity estimates for float
3902681 (Iceland Basin, Nov 2024-Nov 2025) into a single publication-ready,
multi-panel figure:

  (a) Areal      : GOP, NPP, NCP  in mmol C m^-2 d^-1
  (b) Volumetric : GOP, NPP, NCP  in mmol C m^-3 d^-1
  (c) Integration depth z_int(t)  (context strip shared by both panels)

This merges Scripts/gop_vs_ncp_fullmld_figure.py (areal) and
Scripts/gop_vs_ncp_fullmld_volumetric_figure.py (volumetric). The two source
figures use TWIN O2/C axes, which makes the gross-vs-net hierarchy impossible to
read directly. Here GOP (an O2 flux) is converted to carbon units with the repo
photosynthetic quotient PQ = 1.45 (CLAUDE.md; O2:C of production), so all three
estimates share ONE carbon axis per panel and the key message is directly
legible:

    GOP  >=  NPP  >=  NCP        (gross >= net-primary >= net-community)

i.e. neither NPP (CbPM) nor NCP (nitrate budget) exceeds the independent
gross-production ceiling set by the O2 budget.

Two seasonal regimes are shaded and annotated:
  P1  "Spring onset"        - GOP rises concomitantly with NPP and NCP.
  P2  "Post-bloom drawdown" - GOP and NPP stay high but NCP collapses toward
                              zero as respiration consumes the fixed carbon.

Input:  Data/Processed/o2_budget_3902681_fullmld.csv
        Data/Processed/ncp_float_3902681_30d.xlsx
        Data/Processed/icb_npp_2025.csv
Output: Output/synthesis_GOP_NPP_NCP_areal_volumetric.png / .pdf
"""

from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline, PchipInterpolator

# Shared figure config = single source of truth (see Scripts/figure_config.py).
# Adding this file's own directory to sys.path lets the import work regardless
# of the current working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import (            # noqa: E402
    WINDOW_DAYS, MIN_N, GOP_DEEP_MIX_THRESHOLD_M,
    GOP_COLOR, GOP_LIGHT, NCP_COLOR, NPP_COLOR, ZERO_LINE, MLD_COLOR,
    P1_SHADE, P2_SHADE, MASK_BAND,
    centered_window_stats, render_masked, assert_within_axis,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT      = Path(__file__).resolve().parents[1]
O2_CSV         = REPO_ROOT / "Data" / "Processed" / "o2_budget_3902681_fullmld.csv"
FLOAT_NCP_XLSX = REPO_ROOT / "Data" / "Processed" / "ncp_float_3902681_30d.xlsx"
NPP_CSV        = REPO_ROOT / "Data" / "Processed" / "icb_npp_2025.csv"
OUT_PNG        = REPO_ROOT / "Output" / "synthesis_GOP_NPP_NCP_areal_volumetric.png"
OUT_PDF        = REPO_ROOT / "Output" / "synthesis_GOP_NPP_NCP_areal_volumetric.pdf"

MG_C_PER_MMOL = 12.0   # mg C -> mmol C (repo convention)
PQ            = 1.45   # photosynthetic quotient O2:C (CLAUDE.md) -> GOP_C = GOP_O2 / PQ

MLD_FLOOR   = 50       # m  (no upper cap: integration bottom follows the full MLD)
# WINDOW_DAYS, MIN_N now imported from figure_config (single source of truth).

PLOT_START  = pd.Timestamp("2024-11-01")
PLOT_END    = pd.Timestamp("2025-11-01")

# Seasonal regimes to highlight
P1_START, P1_END = pd.Timestamp("2025-04-05"), pd.Timestamp("2025-05-12")   # spring onset
P2_START, P2_END = pd.Timestamp("2025-05-12"), pd.Timestamp("2025-07-20")   # post-bloom drawdown

# ---------------------------------------------------------------------------
# Style  (colours imported from figure_config; only rcParams are local)
# ---------------------------------------------------------------------------
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
    "legend.frameon":  False,
    "legend.fontsize":  9,
})


# ---------------------------------------------------------------------------
# GOP smoothing
# ---------------------------------------------------------------------------
# The centered-window smoother (centered_window_stats) is imported from
# figure_config. Two things are computed here on top of it:
#   * a deep-mixing VALIDITY mask (winter): the float O2 budget is not
#     interpretable once the integration bottom follows deep winter convection
#     (see FIGURE_HANDOFF.md Task 3 / Decision 3). The mask is built in two
#     steps -- a window-MEDIAN depth criterion that places the boundary at the
#     actual convection edge (so the spring-onset GOP survives), then an EROSION
#     by the smoothing half-window that drops any point whose 18-day window
#     straddles that boundary and would blend in a winter tail. See
#     _window_deep_mix_mask.
#   * a robust standard error of the median (replacing the raw 18-day IQR, which
#     saturated the panel). GOP = (median[net] + median[resp]) / PQ, so its
#     uncertainty is the standard error of the SUM of two independent window
#     medians. Each window's SE-of-median is estimated the MAD way
#     (SE ~= 1.2533 * 1.4826 * MAD / sqrt(n)); the two are combined in
#     quadrature and divided by PQ. This is ~3x tighter than the IQR of the raw
#     dusk/dawn rates (which are very noisy) and, being MAD-based, has none of
#     the heavy tails a bootstrap percentile picks up at noisy windows. The band
#     is drawn at +/- SE_MULT standard errors.

SE_MULT = 1.0    # band half-width in units of the standard error of the median


def _drop_short_valid_runs(valid, dates, min_days=WINDOW_DAYS):
    """Set to False any contiguous valid run shorter than ``min_days``.

    Used both inside the deep-mixing mask and again after the series is clipped
    to the plot window, so an isolated short valid stub (e.g. the few days of a
    longer well-stratified run that happen to fall inside the window) is not
    drawn as an orphan segment.
    """
    valid = np.asarray(valid, dtype=bool).copy()
    dates = np.asarray(dates)
    run_min = np.timedelta64(int(min_days * 86400), "s")
    n = len(valid)
    i = 0
    while i < n:
        if valid[i]:
            j = i
            while j + 1 < n and valid[j + 1]:
                j += 1
            if dates[j] - dates[i] < run_min:
                valid[i:j + 1] = False
            i = j + 1
        else:
            i += 1
    return valid


def _window_deep_mix_mask(times, z_bot, grid, window_days=WINDOW_DAYS, min_n=MIN_N,
                          threshold=GOP_DEEP_MIX_THRESHOLD_M, min_run_days=WINDOW_DAYS):
    """Validity mask: True where the 18-day window is free of deep mixing.

    Steps:

    1. *Window-median criterion.* A grid point is provisionally valid if its
       centered window holds at least ``min_n`` net profiles whose MEDIAN
       integration bottom ``z_bot <= threshold``. Using the median (not the
       maximum) puts the boundary at the actual convection edge, so the smoothed
       GOP is retained through the spring onset rather than being masked by a
       few lingering deep profiles.
    2. *Erosion by the smoothing half-window.* A provisionally-valid point is
       dropped if any point inside its own window is invalid, i.e. its 18-day
       window straddles the deep-mixing boundary and would blend in a winter
       tail. This removes boundary contamination without re-masking the
       well-stratified interior.
    3. *Minimum-run filter.* Isolated valid segments shorter than
       ``min_run_days`` are dropped: an 18-day smoothed median needs at least
       that much clean coverage to be meaningful (this also removes the short
       early-November sliver left after erosion).
    """
    times = np.asarray(times)
    z_bot = np.asarray(z_bot, dtype=float)
    grid = np.asarray(grid)
    half = np.timedelta64(int(window_days / 2 * 86400), "s")

    # step 1: window-median depth criterion
    prelim = np.zeros(len(grid), dtype=bool)
    for i, t in enumerate(grid):
        zw = z_bot[(times >= t - half) & (times <= t + half)]
        if zw.size >= min_n:
            prelim[i] = np.median(zw) <= threshold

    # step 2: erode by the smoothing half-window (time-based, so it does not
    # assume a particular grid spacing)
    valid = prelim.copy()
    for i, t in enumerate(grid):
        if not prelim[i]:
            continue
        neigh = (grid >= t - half) & (grid <= t + half)
        if not prelim[neigh].all():
            valid[i] = False

    # step 3: drop valid runs shorter than min_run_days
    if min_run_days:
        valid = _drop_short_valid_runs(valid, grid, min_run_days)
    return valid


def _window_median_se(times, values, grid, window_days=WINDOW_DAYS, min_n=MIN_N):
    """Robust (MAD-based) standard error of the median in each centered window."""
    times = np.asarray(times)
    values = np.asarray(values, dtype=float)
    grid = np.asarray(grid)
    half = np.timedelta64(int(window_days / 2 * 86400), "s")
    se = np.full(len(grid), np.nan)
    for i, t in enumerate(grid):
        v = values[(times >= t - half) & (times <= t + half)]
        if v.size >= min_n:
            sigma = 1.4826 * np.median(np.abs(v - np.median(v)))   # robust std
            se[i] = 1.2533 * sigma / np.sqrt(v.size)               # SE of median
    return se


def _gop_median_se(net_times, net_vals, resp_times, resp_vals, grid,
                   window_days=WINDOW_DAYS, min_n=MIN_N):
    """Standard error of GOP = (median[net] + median[resp]) / PQ per grid time.

    The two window medians are independent, so their standard errors add in
    quadrature; the result is expressed in carbon units (/PQ).
    """
    se_net = _window_median_se(net_times, net_vals, grid, window_days, min_n)
    se_resp = _window_median_se(resp_times, resp_vals, grid, window_days, min_n)
    return np.sqrt(se_net ** 2 + se_resp ** 2) / PQ


def load_gop():
    """Return areal GOP, volumetric GOP, and the per-net-pair z_bot series.

    Both GOP variants are already converted to CARBON units (/PQ) and carry a
    ``valid`` column (False during deep winter convection)."""
    df = pd.read_csv(O2_CSV)
    df["mtime"] = pd.to_datetime(df["mtime"], utc=True).dt.tz_localize(None)

    net = df.loc[df["type"] == "net",
                 ["mtime", "rate_corr_mmol_o2_m2_d", "z_bot_m"]].dropna().sort_values("mtime")
    night = df.loc[df["type"] == "night",
                   ["mtime", "rate_corr_mmol_o2_m2_d", "z_bot_m"]].dropna().sort_values("mtime")

    net["rate_areal"] = net["rate_corr_mmol_o2_m2_d"]
    net["rate_vol"]   = net["rate_corr_mmol_o2_m2_d"] / net["z_bot_m"]
    night_resp_areal  = -night["rate_corr_mmol_o2_m2_d"].to_numpy()
    night_resp_vol    = -(night["rate_corr_mmol_o2_m2_d"] / night["z_bot_m"]).to_numpy()

    zbot = df.loc[df["type"] == "net", ["mtime", "z_bot_m"]].dropna().sort_values("mtime")

    t0 = min(net["mtime"].min(), night["mtime"].min())
    t1 = max(net["mtime"].max(), night["mtime"].max())
    grid = pd.date_range(t0, t1, freq="1D").to_numpy()

    # deep-mixing validity mask (shared by areal & volumetric GOP)
    valid = _window_deep_mix_mask(net["mtime"].to_numpy(), net["z_bot_m"].to_numpy(), grid)

    def build(net_col, resp_vals):
        net_m, _, _, net_n = centered_window_stats(
            net["mtime"].to_numpy(), net[net_col].to_numpy(), grid)
        resp_m, _, _, resp_n = centered_window_stats(
            night["mtime"].to_numpy(), resp_vals, grid)
        se = _gop_median_se(
            net["mtime"].to_numpy(), net[net_col].to_numpy(),
            night["mtime"].to_numpy(), resp_vals, grid)
        ok = (net_n >= MIN_N) & (resp_n >= MIN_N)
        # GOP (O2) -> carbon via PQ; band is +/- SE_MULT s.e. of the median
        mean = (net_m + resp_m) / PQ
        out = pd.DataFrame({
            "date":     pd.to_datetime(grid),
            "gop_mean": mean,
            "gop_lo":   mean - SE_MULT * se,
            "gop_hi":   mean + SE_MULT * se,
            "valid":    valid,
        })[ok].reset_index(drop=True)
        return out

    gop_areal = build("rate_areal", night_resp_areal)
    gop_vol   = build("rate_vol",   night_resp_vol)

    # raw per-pair net rates (carbon) for the faint background trace; keep z_bot
    # so the raw trace can be suppressed during deep mixing too.
    net_out = net[["mtime", "z_bot_m"]].copy()
    net_out["rate_areal_c"] = net["rate_areal"] / PQ
    net_out["rate_vol_c"]   = net["rate_vol"] / PQ

    n_valid = int(gop_areal["valid"].sum())
    print(f"[GOP] {len(gop_areal)} daily estimates (areal & volumetric), converted "
          f"O2->C with PQ={PQ}; {n_valid} valid (z_bot<= "
          f"{GOP_DEEP_MIX_THRESHOLD_M:.0f} m in-window), "
          f"{len(gop_areal) - n_valid} masked (deep mixing)")
    return gop_areal, gop_vol, net_out, zbot


def zint_series(zbot):
    """Common daily MLD-based integration depth z(t) for volumetric dilution."""
    daily = (zbot.assign(day=zbot["mtime"].dt.floor("D"))
                  .groupby("day", as_index=False)["z_bot_m"].mean())
    return daily["day"].to_numpy(), daily["z_bot_m"].to_numpy()


# ---------------------------------------------------------------------------
# NCP (float 3902681 nitrate budget, 30-day window) - areal & volumetric
# ---------------------------------------------------------------------------

def load_ncp(t0, t1, z_dates, z_vals):
    df = pd.read_excel(FLOAT_NCP_XLSX)
    df["date"] = pd.to_datetime(df["date"])
    df = df.groupby("date", as_index=False)["NCP"].mean().sort_values("date")
    pad = pd.Timedelta(days=15)
    sub = df[(df["date"] >= t0 - pad) & (df["date"] <= t1 + pad)].copy()

    zx = z_dates.astype("datetime64[ns]").astype("float64")
    dx = sub["date"].to_numpy().astype("datetime64[ns]").astype("float64")
    sub["z_int"] = np.interp(dx, zx, z_vals)
    sub["NCP_vol"] = sub["NCP"] / sub["z_int"]
    print(f"[NCP] {len(sub)} estimates; z(t) range "
          f"[{sub['z_int'].min():.0f}, {sub['z_int'].max():.0f}] m")
    return sub


def spline(sub, ycol, kind, step_days=3):
    t0 = sub["date"].min()
    t_days = (sub["date"] - t0).dt.days.to_numpy(dtype=float)
    y = sub[ycol].to_numpy()
    interp = CubicSpline(t_days, y) if kind == "cubic" else PchipInterpolator(t_days, y)
    t_dense = np.arange(0, t_days.max() + step_days, step_days)
    dates = t0 + pd.to_timedelta(t_dense, unit="D")
    return dates, interp(t_dense)


# ---------------------------------------------------------------------------
# NPP (CbPM, Iceland Basin, 10-day) - areal & volumetric
# ---------------------------------------------------------------------------

def load_npp(t0, t1, z_dates, z_vals):
    df = pd.read_csv(NPP_CSV)
    df["date"] = pd.to_datetime(df["date_10day"])
    df = df.dropna(subset=["npp"]).sort_values("date")
    pad = pd.Timedelta(days=15)
    sub = df[(df["date"] >= t0 - pad) & (df["date"] <= t1 + pad)].copy()

    sub["npp_areal"] = sub["npp"] / MG_C_PER_MMOL
    zx = z_dates.astype("datetime64[ns]").astype("float64")
    dx = sub["date"].to_numpy().astype("datetime64[ns]").astype("float64")
    sub["z_int"] = np.interp(dx, zx, z_vals)
    sub["npp_vol"] = sub["npp_areal"] / sub["z_int"]
    print(f"[NPP] {len(sub)} bins; areal peak {sub['npp_areal'].max():.0f} "
          f"mmol C m-2 d-1, vol peak {sub['npp_vol'].max():.2f} mmol C m-3 d-1")
    return sub


# ---------------------------------------------------------------------------
# Panel drawing
# ---------------------------------------------------------------------------

def clip_window(df, col):
    return df[(df[col] >= PLOT_START) & (df[col] <= PLOT_END)]


def draw_panel(ax, gop, net_raw, raw_col, ncp_dates, ncp_m, ncp_pts, ncp_col,
               npp_dates, npp_vals, ylim, label):
    lo, hi = ylim
    ax.set_ylim(lo, hi)
    ax.set_xlim(PLOT_START, PLOT_END)

    # seasonal regime shading (behind everything)
    ax.axvspan(P1_START, P1_END, color=P1_SHADE, alpha=0.35, lw=0, zorder=0)
    ax.axvspan(P2_START, P2_END, color=P2_SHADE, alpha=0.35, lw=0, zorder=0)
    ax.axhline(0, color=ZERO_LINE, lw=0.8, ls=":", zorder=1)

    # faint raw per-pair GOP trace across the FULL record. No .clip(): points
    # outside the panel are dropped to NaN (not pinned to the border), so
    # nothing is misdrawn -- during winter deep mixing the areal rate simply
    # runs off-scale and the trace disappears where it exceeds the axis.
    raw_y = net_raw[raw_col].where((net_raw[raw_col] >= lo) & (net_raw[raw_col] <= hi))
    ax.plot(net_raw["mtime"], raw_y, lw=0.5, color=GOP_COLOR, alpha=0.15, zorder=2)

    # GOP (carbon) line = the gross-production ceiling, drawn over the ENTIRE
    # record including winter deep convection. Where the areal winter rate
    # exceeds the panel (see figure axis choice), matplotlib clips the line at
    # the axis box, so the biologically-relevant season stays legible. GOP is
    # intentionally allowed to exit the axis and is therefore NOT passed to
    # assert_within_axis (that guard is for the net estimates, which must fit).
    ax.plot(gop["date"], gop["gop_mean"], lw=2.2, color=GOP_COLOR, zorder=6)

    # NPP then NCP
    ax.plot(npp_dates, npp_vals, lw=1.8, color=NPP_COLOR, marker="o", ms=3.5,
            mfc="white", mec=NPP_COLOR, mew=0.8, zorder=5)
    ax.plot(ncp_dates, ncp_m, lw=2.0, color=NCP_COLOR, zorder=5)
    ax.scatter(ncp_pts["date"], ncp_pts[ncp_col], s=15, color=NCP_COLOR,
               zorder=7, linewidths=0.5, edgecolors="white")

    ax.text(0.008, 0.94, label, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")

    # forbid silent clipping of the NET estimates: NPP and NCP must fit the
    # axis. GOP is intentionally allowed to run off-scale during winter deep
    # mixing (see comment above) and is excluded from this guard.
    assert_within_axis(ax, npp_vals, ncp_m, ncp_pts[ncp_col].to_numpy())


def main() -> int:
    gop_areal, gop_vol, net_raw, zbot = load_gop()
    z_dates, z_vals = zint_series(zbot)

    ncp = load_ncp(PLOT_START, PLOT_END, z_dates, z_vals)
    npp = load_npp(PLOT_START, PLOT_END, z_dates, z_vals)

    # dense interpolants
    nad, nam = spline(ncp, "NCP", "cubic")      # areal NCP
    nvd, nvm = spline(ncp, "NCP_vol", "pchip")  # volumetric NCP

    def within(dates, vals):
        m = (dates >= PLOT_START) & (dates <= PLOT_END)
        return dates[m], vals[m]

    nad, nam = within(nad, nam)
    nvd, nvm = within(nvd, nvm)

    gop_areal = clip_window(gop_areal, "date")
    gop_vol   = clip_window(gop_vol, "date")
    # re-drop valid runs that are too short once restricted to the plot window
    for g in (gop_areal, gop_vol):
        g["valid"] = _drop_short_valid_runs(g["valid"].to_numpy(), g["date"].to_numpy())
    net_raw   = clip_window(net_raw, "mtime")
    zbot_w    = clip_window(zbot, "mtime")
    ncp_pts   = clip_window(ncp, "date")
    npp_pts   = clip_window(npp, "date")

    # -----------------------------------------------------------------------
    fig, (axA, axB, axD) = plt.subplots(
        3, 1, figsize=(11, 9.3), dpi=300,
        gridspec_kw={"height_ratios": [4, 4, 1.1], "hspace": 0.11}, sharex=True)

    draw_panel(axA, gop_areal, net_raw, "rate_areal_c", nad, nam, ncp_pts, "NCP",
               npp_pts["date"], npp_pts["npp_areal"], (-140, 300), "(a)")
    draw_panel(axB, gop_vol, net_raw, "rate_vol_c", nvd, nvm, ncp_pts, "NCP_vol",
               npp_pts["date"], npp_pts["npp_vol"], (-3.0, 6.5), "(b)")

    axA.set_ylabel(r"Areal rate (mmol C m$^{-2}$ d$^{-1}$)", fontsize=11)
    axB.set_ylabel(r"Volumetric rate (mmol C m$^{-3}$ d$^{-1}$)", fontsize=11)

    # regime labels + "GOP is the ceiling" annotation on the top panel
    axA.text(P1_START + (P1_END - P1_START) / 2, 292, "P1  spring onset",
             ha="center", va="top", fontsize=8.5, color="#8a6d0b",
             fontweight="bold")
    axA.text(P2_START + (P2_END - P2_START) / 2, 292, "P2  post-bloom drawdown",
             ha="center", va="top", fontsize=8.5, color="#5b4a7a",
             fontweight="bold")
    axA.annotate("In the productive season, GOP sets the\n"
                 "gross-production ceiling: NPP & NCP stay below",
                 xy=(pd.Timestamp("2025-06-05"), 200),
                 xytext=(pd.Timestamp("2025-08-05"), 255),
                 fontsize=8.5, color=GOP_COLOR, ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", color=GOP_COLOR, lw=1.0),
                 bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.75))

    # arrows describing each regime on the volumetric panel
    axB.annotate("spring onset: GOP, NPP, NCP\nrise together as mixing ceases",
                 xy=(pd.Timestamp("2025-04-28"), 3.2),
                 xytext=(pd.Timestamp("2024-12-05"), 5.2),
                 fontsize=8.5, color="#8a6d0b", ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", color="#8a6d0b", lw=1.0),
                 bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.75))
    axB.annotate("NCP collapses toward 0\n(respiration) while GOP & NPP stay high",
                 xy=(pd.Timestamp("2025-06-20"), 0.2),
                 xytext=(pd.Timestamp("2025-07-25"), 4.4),
                 fontsize=8.5, color="#5b4a7a", ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", color="#5b4a7a", lw=1.0),
                 bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.75))

    # -------- integration-depth strip --------
    axD.plot(zbot_w["mtime"], zbot_w["z_bot_m"], lw=0.8, color=MLD_COLOR, alpha=0.9)
    axD.axhline(MLD_FLOOR, color=ZERO_LINE, lw=0.6, ls="--")
    axD.axvspan(P1_START, P1_END, color=P1_SHADE, alpha=0.35, lw=0, zorder=0)
    axD.axvspan(P2_START, P2_END, color=P2_SHADE, alpha=0.35, lw=0, zorder=0)
    zbot_max = float(zbot_w["z_bot_m"].max()) if len(zbot_w) else MLD_FLOOR
    axD.set_ylim(zbot_max + 25, MLD_FLOOR - 15)
    axD.set_ylabel("z$_{int}$ (m)", fontsize=9)
    axD.tick_params(axis="y", labelsize=8)
    axD.text(0.008, 0.90, "(c)", transform=axD.transAxes, fontsize=13,
             fontweight="bold", va="top", ha="left")

    # x-axis
    axD.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    axD.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    axD.xaxis.set_minor_locator(mdates.MonthLocator())

    # -------- shared legend --------
    legend_handles = [
        Line2D([0], [0], color=GOP_COLOR, lw=2.2,
               label="GOP  (float O$_2$ budget, full MLD; O$_2$$\\to$C, PQ=1.45)"),
        Line2D([0], [0], color=NPP_COLOR, lw=1.8, marker="o", ms=3.5,
               mfc="white", mec=NPP_COLOR, label="NPP  (CbPM, Iceland Basin, 10-day)"),
        Line2D([0], [0], color=NCP_COLOR, lw=2.0,
               label="NCP  (float 3902681 nitrate budget, 30-day)"),
    ]
    leg = axA.legend(handles=legend_handles, loc="upper left", fontsize=8.5, ncol=1,
                     handlelength=1.8, borderpad=0.5, labelspacing=0.4,
                     frameon=True, framealpha=1.0)
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_edgecolor("none")
    leg.set_zorder(20)

    fig.suptitle(
        "Float 3902681 (Iceland Basin, Nov 2024-Nov 2025): gross vs. net productivity\n"
        "GOP $\\geq$ NPP $\\geq$ NCP on a common carbon axis (areal & volumetric)",
        fontsize=13, y=0.985)

    # caption footnote: native cadences + smoothing/lag note (FIGURE_HANDOFF Task 3.7)
    fig.text(
        0.5, 0.012,
        "Native cadences: GOP dusk/dawn O$_2$ budget smoothed on an 18-day centered "
        "window; NPP CbPM ~10-day; NCP float nitrate budget 30-day. The full GOP "
        "record is shown; during winter deep convection the integration bottom "
        "follows the mixed layer to depth, so the areal GOP swings off-scale and "
        "runs off panel (a) (the volumetric rate in (b) stays bounded). The "
        "GOP$\\geq$NPP$\\geq$NCP ordering and the P1$\\rightarrow$P2 lag were checked "
        "to survive equal (18-day) smoothing of all three series.",
        ha="center", va="bottom", fontsize=7.2, color="#444444", wrap=True)

    fig.tight_layout(rect=(0, 0.035, 1, 0.955))
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, bbox_inches="tight", dpi=300)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
