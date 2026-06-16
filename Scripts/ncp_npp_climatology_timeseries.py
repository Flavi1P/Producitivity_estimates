"""
Seasonal (climatological) time series of NPP and NCP, mean +/- SD ribbons.

Both productivity terms on a single mmol C m-2 d-1 axis vs calendar month:

  NPP : column-integrated CbPM NPP climatology, converted mg C -> mmol C (/12)
        from Output/cmems_npp_climatology_domain_mean.csv
        (mean = npp_int_mean, SD = npp_int_std = spatial spread of the field)
  NCP : nitrate-drawdown NCP climatology, 30-day window
        from Output/ncp_climatology_icb_30d.csv
        (mean = mean_ncp, SD = sd_ncp = inter-annual spread)

Run:
    & "C:/Users/petit/miniconda3/envs/UVP6/python.exe" Scripts/ncp_npp_climatology_timeseries.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
NPP_CSV = REPO_ROOT / "Output" / "cmems_npp_climatology_domain_mean.csv"
NCP_CSV = REPO_ROOT / "Output" / "ncp_climatology_icb_30d.csv"
OUT_PNG = REPO_ROOT / "Output" / "ncp_npp_climatology_timeseries.png"
OUT_PNG_NCP = REPO_ROOT / "Output" / "ncp_climatology_timeseries.png"

MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

NPP_COLOR = "#1f6e3f"
NCP_COLOR = "#1f4e79"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncp-only", action="store_true",
                    help="Plot only the NCP climatology (no NPP series).")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    ncp = pd.read_csv(NCP_CSV)

    fig, ax = plt.subplots(figsize=(9, 5.2), dpi=130)

    # NPP (skipped in --ncp-only mode)
    if not args.ncp_only:
        npp = pd.read_csv(NPP_CSV)
        npp["npp_mean"] = npp["npp_int_mean"] / 12.0    # mg C -> mmol C
        npp["npp_sd"] = npp["npp_int_std"] / 12.0
        n = npp.dropna(subset=["npp_mean"])
        ax.fill_between(n["month"], n["npp_mean"] - n["npp_sd"],
                        n["npp_mean"] + n["npp_sd"],
                        color=NPP_COLOR, alpha=0.20, linewidth=0)
        ax.plot(n["month"], n["npp_mean"], color=NPP_COLOR, lw=1.8,
                marker="o", ms=5, label="NPP (CbPM, satellite/CMEMS)")

    # NCP
    c = ncp.dropna(subset=["mean_ncp"])
    ax.fill_between(c["month"], c["mean_ncp"] - c["sd_ncp"],
                    c["mean_ncp"] + c["sd_ncp"],
                    color=NCP_COLOR, alpha=0.20, linewidth=0)
    ax.plot(c["month"], c["mean_ncp"], color=NCP_COLOR, lw=1.8,
            marker="s", ms=5, label="NCP (nitrate, 30-day window)")

    ax.axhline(0, color="0.6", lw=0.8, ls="--", zorder=0)
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(MONTH_NAMES)
    ax.set_xlim(0.7, 12.3)
    ax.set_xlabel("Month", fontsize=11)
    ax.set_ylabel("Production (mmol C m$^{-2}$ d$^{-1}$)", fontsize=11)
    if args.ncp_only:
        title = ("Iceland Basin / Irminger Sea — NCP climatology "
                 "(30-day window, mean $\\pm$ SD)")
        out_png = args.out or OUT_PNG_NCP
    else:
        title = ("Iceland Basin / Irminger Sea — NPP & NCP climatology "
                 "(mean $\\pm$ SD)")
        out_png = args.out or OUT_PNG
    ax.set_title(title, fontsize=11)
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)
    ax.legend(frameon=False, fontsize=10, loc="upper left")

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    print(f"Wrote {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
