"""
Domain-mean timeseries of column-integrated CbPM NPP (0-200 m).

Reads Output/cmems_cbpm_monthly_0p25deg.nc (produced by cmems_3d_cbpm.py),
clips to the study domain, computes the spatial mean of npp_int for every
calendar month in the record (2015-2023), and writes:

  Output/cmems_npp_timeseries_domain_mean.csv
      date, npp_int_mean, npp_int_std, npp_int_p10, npp_int_p90, n_cells
      (npp_int in mg C m-2 d-1)

  Output/cmems_npp_timeseries.png
      Time series plot of the domain mean +/- 1 SD ribbon.

Run:
    python Scripts/cmems_npp_timeseries.py
    python Scripts/cmems_npp_timeseries.py --lon-min -40 --lon-max -10 --lat-min 58 --lat-max 65
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import xarray as xr

REPO_ROOT = Path(__file__).resolve().parents[1]
NC_FILE  = REPO_ROOT / "Output" / "cmems_cbpm_monthly_0p25deg.nc"
OUT_CSV  = REPO_ROOT / "Output" / "cmems_npp_timeseries_domain_mean.csv"
OUT_PNG  = REPO_ROOT / "Output" / "cmems_npp_timeseries.png"

# Default domain: Iceland Basin + Irminger Sea (matches cmems_npp_climatology.py)
LON_MIN, LON_MAX = -40.0, -10.0
LAT_MIN, LAT_MAX =  58.0,  65.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in",      dest="nc_in",  type=Path, default=NC_FILE)
    ap.add_argument("--out-csv", dest="out_csv", type=Path, default=OUT_CSV)
    ap.add_argument("--out-png", dest="out_png", type=Path, default=OUT_PNG)
    ap.add_argument("--lon-min", type=float, default=LON_MIN)
    ap.add_argument("--lon-max", type=float, default=LON_MAX)
    ap.add_argument("--lat-min", type=float, default=LAT_MIN)
    ap.add_argument("--lat-max", type=float, default=LAT_MAX)
    args = ap.parse_args()

    ds = xr.open_dataset(args.nc_in)
    npp = ds["npp_int"].sel(
        longitude=slice(args.lon_min, args.lon_max),
        latitude=slice(args.lat_min,  args.lat_max),
    )
    print(
        f"Domain: lon {args.lon_min}..{args.lon_max}, "
        f"lat {args.lat_min}..{args.lat_max}  "
        f"-> {npp.sizes['latitude']} lat x {npp.sizes['longitude']} lon, "
        f"{npp.sizes['time']} months"
    )

    # --- spatial statistics per timestep ---
    rows = []
    for t in npp.time.values:
        field = npp.sel(time=t).values          # (lat, lon)
        vals  = field[np.isfinite(field)]
        if vals.size == 0:
            rows.append(dict(date=pd.Timestamp(t),
                             npp_int_mean=np.nan, npp_int_std=np.nan,
                             npp_int_p10=np.nan, npp_int_p90=np.nan,
                             n_cells=0))
            continue
        rows.append(dict(
            date=pd.Timestamp(t),
            npp_int_mean=round(float(np.mean(vals)),       3),
            npp_int_std= round(float(np.std(vals)),        3),
            npp_int_p10= round(float(np.percentile(vals, 10)), 3),
            npp_int_p90= round(float(np.percentile(vals, 90)), 3),
            n_cells=int(vals.size),
        ))

    df = pd.DataFrame(rows)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv}  ({len(df)} rows)")

    # --- plot ---
    fig, ax = plt.subplots(figsize=(11, 4.5), dpi=130)

    valid = df.dropna(subset=["npp_int_mean"])
    dates = valid["date"]
    mean  = valid["npp_int_mean"]
    sd    = valid["npp_int_std"]

    ax.fill_between(dates, mean - sd, mean + sd,
                    color="#1f6e3f", alpha=0.20, linewidth=0, label="±1 SD (spatial)")
    ax.plot(dates, mean, color="#1f6e3f", lw=1.8, marker="o", ms=3.5,
            label="NPP mean (CbPM, 0-200 m)")

    ax.axhline(0, color="0.6", lw=0.8, ls="--", zorder=0)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4, 7, 10]))
    ax.set_xlim(dates.min(), dates.max())
    ax.set_xlabel("Date", fontsize=11)
    ax.set_ylabel("NPP (mg C m$^{-2}$ d$^{-1}$)", fontsize=11)
    ax.set_title(
        "Iceland Basin / Irminger Sea — column-integrated NPP 0-200 m "
        "(CbPM × CMEMS 3D BGC, monthly mean ± 1 SD)",
        fontsize=10,
    )
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)
    ax.legend(frameon=False, fontsize=10)

    fig.tight_layout()
    fig.savefig(args.out_png, bbox_inches="tight")
    print(f"Wrote {args.out_png}")

    print(df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
