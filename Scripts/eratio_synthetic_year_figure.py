"""
Synthetic-year NPP, NCP and e-ratio (NCP:NPP) — Iceland Basin & Irminger Sea.

Reuses the exact synthetic-year builders from npp_timeseries_figure.py and
ncp_timeseries_figure.py so the climatologies are identical to the ones in
the standalone NPP and NCP figures:

  NPP : CbPM synthetic year (Feb-Oct observed, Nov-Jan spline-interpolated),
        converted mg C -> mmol C m-2 d-1 (/12) to share NCP's units.
  NCP : nitrate-drawdown synthetic year (full 12-month periodic climatology).

Both are plotted on the same primary axis (mmol C m-2 d-1). The e-ratio
(NCP:NPP) is overlaid on a secondary axis, restricted to the Feb-Oct window
where NPP is an actual satellite-derived estimate (not a winter spline
extrapolation), avoiding divide-by-near-zero artefacts.

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
OBS_MONTHS   = set(range(2, 11))   # Feb-Oct: real CbPM NPP, not winter spline


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

    # ---- e-ratio, restricted to the observed NPP window ----
    doy_month = np.array(
        [(pd.Timestamp(2021, 1, 1) + pd.Timedelta(days=int(d) - 1)).month for d in doy]
    )
    valid = np.isin(doy_month, list(OBS_MONTHS))
    eratio = np.full_like(doy, np.nan, dtype=float)
    eratio[valid] = ncp_mean[valid] / npp_mean_c[valid]

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
    ax2.plot(doy[valid], eratio[valid], color=ERATIO_COLOR, lw=2.4, ls="-",
             marker="o", ms=4.5, markevery=3, zorder=6,
             label="e-ratio (NCP:NPP)")

    ax2.set_ylabel("e-ratio  (NCP : NPP)", fontsize=11, color=ERATIO_COLOR)
    ax2.set_ylim(-1.0, 1.3)
    ax2.tick_params(axis="y", colors=ERATIO_COLOR)
    ax2.spines["right"].set_color(ERATIO_COLOR)

    h2, l2 = ax2.get_legend_handles_labels()
    ax2.legend(h2, l2, loc="upper right", fontsize=9)

    fig.suptitle(
        "Synthetic-year NPP, NCP and e-ratio\n"
        "Iceland Basin & Irminger Sea (pooled, 2015–2025)",
        fontsize=11.5, y=0.99,
    )
    ax1.text(0.985, 0.02,
             "e-ratio shown Feb–Oct only\n(no CbPM NPP estimate in winter)",
             transform=ax1.transAxes, ha="right", va="bottom",
             fontsize=7.5, color="0.45")

    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, bbox_inches="tight", dpi=300)
    print(f"Wrote {OUT_PNG}")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
