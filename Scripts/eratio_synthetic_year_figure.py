"""
Synthetic-year NPP, NCP and e-ratio (NCP:NPP) — Iceland Basin & Irminger Sea.

Reuses the exact synthetic-year builders from npp_timeseries_figure.py and
ncp_timeseries_figure.py so the climatologies are identical to the ones in
the standalone NPP and NCP figures:

  NPP : CbPM synthetic year (Feb-Oct observed, Nov-Jan spline-interpolated),
        converted mg C -> mmol C m-2 d-1 (/12) to share NCP's units.
  NCP : nitrate-drawdown synthetic year (full 12-month periodic climatology).

Both are plotted on the same primary axis (mmol C m-2 d-1). The e-ratio
(NCP:NPP) is overlaid on a secondary axis and shown for the whole year:
Feb-Oct uses the actual CbPM NPP estimate (solid line), while Nov-Jan uses
the winter cubic-spline interpolation of NPP (dashed line, exactly the curve
drawn in Fig 4). The only points dropped are the deep-winter days where the
interpolated NPP falls below NPP_FLOOR (mmol C m-2 d-1) and the ratio would
diverge as NPP -> 0.

Input:  Output/cmems_npp_timeseries_domain_mean.csv
        Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_uncertainty.xlsx
Output: Output/IcelandIrminger_2015_2025/eratio_synthetic_year.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import npp_timeseries_figure as npp_mod   # noqa: E402
import ncp_timeseries_figure as ncp_mod   # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_PNG   = REPO_ROOT / "Output" / "IcelandIrminger_2015_2025" / "eratio_synthetic_year.png"

MGC_TO_MMOLC = 1.0 / 12.0   # matches convention used elsewhere in this repo

GREEN        = npp_mod.GREEN
GREEN_LIGHT  = npp_mod.GREEN_LIGHT
BLUE         = ncp_mod.BLUE
BLUE_LIGHT   = ncp_mod.BLUE_LIGHT
ERATIO_COLOR = "#7a2048"
MONTH_LABELS = npp_mod.MONTH_LABELS
OBS_MONTHS   = set(range(2, 11))   # Feb-Oct: real CbPM NPP; else winter spline
# Below this interpolated-NPP value the e-ratio (NCP:NPP) diverges as NPP -> 0
# (the Fig-4 spline is anchored to 0 at DOY 1/365), so those deep-winter days
# are dropped rather than plotted as a spike. Tunable.
NPP_FLOOR    = 5.0                 # mmol C m-2 d-1


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


def main() -> int:
    # ---- NPP synthetic year (mg C -> mmol C) ----
    npp_df = pd.read_csv(npp_mod.CSV_FILE, parse_dates=["date"])
    doy, npp_mean, npp_p10, npp_p90, npp_clim = npp_mod.build_synthetic_year(npp_df)
    npp_mean_c = npp_mean * MGC_TO_MMOLC
    npp_p10_c  = npp_p10  * MGC_TO_MMOLC
    npp_p90_c  = npp_p90  * MGC_TO_MMOLC

    # ---- NCP synthetic year (already mmol C) ----
    ncp_df = ncp_mod.load_data()
    doy_ncp, ncp_mean, ncp_p10, ncp_p90, ncp_clim = ncp_mod.build_synthetic_year(ncp_df)
    assert np.allclose(doy, doy_ncp), "NPP/NCP synthetic-year grids must match"

    # ---- e-ratio, full year (winter uses the Fig-4 spline NPP) ----
    doy_month = np.array(
        [(pd.Timestamp(2021, 1, 1) + pd.Timedelta(days=int(d) - 1)).month for d in doy]
    )
    is_obs = np.isin(doy_month, list(OBS_MONTHS))   # Feb-Oct: real CbPM NPP
    is_win = ~is_obs                                # Nov-Jan: NPP from spline
    # Only define the ratio where NPP is above the floor (avoid the NPP -> 0
    # singularity where the spline is anchored to 0).
    npp_ok = npp_mean_c >= NPP_FLOOR
    eratio = np.full_like(doy, np.nan, dtype=float)
    eratio[npp_ok] = ncp_mean[npp_ok] / npp_mean_c[npp_ok]

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
    ax1.set_xlim(1, 365)
    mid_doys = [pd.Timestamp(2021, m, 15).timetuple().tm_yday for m in range(1, 13)]
    ax1.set_xticks(mid_doys)
    ax1.set_xticklabels(MONTH_LABELS)
    ax1.set_xlabel("Month", fontsize=11)
    ax1.legend(loc="upper left", fontsize=8.5, ncol=1)

    # e-ratio, overlaid on secondary axis
    ax2.axhline(0, color=ERATIO_COLOR, lw=0.7, ls=":", alpha=0.4, zorder=3)
    ax2.axhline(1, color=ERATIO_COLOR, lw=0.8, ls="--", alpha=0.5, zorder=3)

    # Observed-NPP window: solid line + markers
    _plot_segments(ax2, doy, eratio, is_obs & npp_ok, connect=npp_ok,
                   color=ERATIO_COLOR, lw=2.4, ls="-",
                   marker="o", ms=4.5, markevery=3, zorder=6,
                   label="e-ratio (NCP:NPP), observed NPP")
    # Winter: interpolated-NPP window, dashed (matches Fig 4 winter styling)
    _plot_segments(ax2, doy, eratio, is_win & npp_ok, connect=npp_ok,
                   color=ERATIO_COLOR, lw=1.8, ls="--", zorder=6,
                   label="e-ratio, winter (interpolated NPP)")

    ax2.set_ylabel("e-ratio  (NCP : NPP)", fontsize=11, color=ERATIO_COLOR)
    # Contain the (bounded) winter dip, keeping the 1.0 reference visible.
    er_lo = np.nanmin(eratio)
    er_hi = np.nanmax(eratio)
    ax2.set_ylim(min(-1.0, np.floor((er_lo - 0.2) * 2) / 2),
                 max(1.3, np.ceil((er_hi + 0.2) * 2) / 2))
    ax2.tick_params(axis="y", colors=ERATIO_COLOR)
    ax2.spines["right"].set_color(ERATIO_COLOR)

    h2, l2 = ax2.get_legend_handles_labels()
    ax2.legend(h2, l2, loc="upper right", fontsize=8.5)

    fig.suptitle(
        "Synthetic-year NPP, NCP and e-ratio\n"
        "Iceland Basin & Irminger Sea (pooled, 2015–2025)",
        fontsize=11.5, y=0.99,
    )
    ax1.text(0.985, 0.02,
             "e-ratio solid = observed CbPM NPP (Feb–Oct)\n"
             "dashed = winter spline-interpolated NPP (Nov–Jan)\n"
             f"(undefined where interpolated NPP < {NPP_FLOOR:.0f} mmol C m$^{{-2}}$ d$^{{-1}}$)",
             transform=ax1.transAxes, ha="right", va="bottom",
             fontsize=7.0, color="0.45")

    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, bbox_inches="tight", dpi=300)
    print(f"Wrote {OUT_PNG}")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
