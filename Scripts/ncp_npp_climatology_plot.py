"""
Property-to-property seasonal plot: NCP vs NPP climatology (Irminger / ICB).

Mimics Data/Processed/npq_to_npp.png — one marker per calendar month,
coloured by month on a cyclic (twilight) colourbar, joined by grey arrows in
calendar order to trace the seasonal trajectory, each point labelled Jan..Dec.

X: column-integrated CbPM NPP climatology (mg C m-2 d-1)
       from Output/cmems_npp_climatology_domain_mean.csv
Y: nitrate-drawdown NCP climatology, 30-day window (mmol C m-2 d-1)
       from Output/ncp_climatology_icb_30d.csv

Months without an NPP value (winter light gap: Jan, Nov, Dec) are dropped,
since the plot needs both properties.

Run:
    & "C:/Users/petit/miniconda3/envs/UVP6/python.exe" Scripts/ncp_npp_climatology_plot.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import Normalize

REPO_ROOT = Path(__file__).resolve().parents[1]
NPP_CSV = REPO_ROOT / "Output" / "cmems_npp_climatology_domain_mean.csv"
NCP_CSV = REPO_ROOT / "Output" / "ncp_climatology_icb_30d.csv"
OUT_PNG = REPO_ROOT / "Output" / "ncp_to_npp_climatology.png"

MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def main() -> int:
    npp = pd.read_csv(NPP_CSV)[["month", "npp_int_mean"]]
    ncp = pd.read_csv(NCP_CSV)[["month", "mean_ncp"]]

    df = (npp.merge(ncp, on="month", how="outer")
            .sort_values("month")
            .dropna(subset=["npp_int_mean", "mean_ncp"]))
    months = df["month"].to_numpy()
    x = df["npp_int_mean"].to_numpy() / 12.0   # mg C -> mmol C
    y = df["mean_ncp"].to_numpy()
    print(f"Plotting {len(df)} months: "
          f"{', '.join(MONTH_NAMES[m - 1] for m in months)}")

    cmap = cm.twilight
    norm = Normalize(vmin=1, vmax=12)

    fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=130)

    # grey arrows joining consecutive (available) months
    for i in range(len(df) - 1):
        ax.annotate(
            "", xy=(x[i + 1], y[i + 1]), xytext=(x[i], y[i]),
            arrowprops=dict(arrowstyle="-|>", color="0.55", lw=1.6,
                            shrinkA=11, shrinkB=11),
            zorder=1,
        )

    ax.scatter(x, y, c=months, cmap=cmap, norm=norm,
               s=190, edgecolor="black", linewidth=1.1, zorder=2)

    # month labels, nudged up-left like the reference
    for xi, yi, m in zip(x, y, months):
        ax.annotate(MONTH_NAMES[m - 1], (xi, yi),
                    textcoords="offset points", xytext=(8, 8),
                    fontsize=10, fontweight="bold", zorder=3)

    ax.axhline(0, color="0.7", lw=0.8, ls="--", zorder=0)
    ax.set_xlabel("NPP (mmol C m$^{-2}$ d$^{-1}$)", fontsize=11)
    ax.set_ylabel("NCP (mmol C m$^{-2}$ d$^{-1}$)", fontsize=11)
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)
    ax.margins(0.12)

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, ticks=range(1, 13), pad=0.02)
    cbar.ax.set_yticklabels(MONTH_NAMES)
    cbar.set_label("Month", fontsize=11)

    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, bbox_inches="tight")
    print(f"Wrote {OUT_PNG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
