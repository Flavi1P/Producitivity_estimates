"""
Monthly climatology of column-integrated CbPM NPP over the Iceland Basin box.

Reads the gridded CbPM output produced by Scripts/cmems_3d_cbpm.py
(Output/cmems_cbpm_monthly_0p25deg.nc, monthly 0.25 deg, 2015-2023),
clips to lon -40..-10 / lat 58..64, collapses the 9-year monthly series to
a 12-month climatology (mean over years per calendar month) of npp_int, and
writes two CSVs:

  Output/cmems_npp_climatology_gridded.csv
      month, latitude, longitude, npp_int_mean, npp_int_std, n_years
      (full spatial field of the climatology, one row per cell per month)

  Output/cmems_npp_climatology_domain_mean.csv
      month, npp_int_mean, npp_int_std, npp_int_p10, npp_int_p90, n_cells
      (domain-averaged seasonal cycle; mean over all valid cells of the
       monthly-climatology field)

npp_int is column-integrated NPP 0-200 m in mg C m-2 d-1.

Run:
    & "C:/Users/petit/miniconda3/envs/UVP6/python.exe" Scripts/cmems_npp_climatology.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

REPO_ROOT = Path(__file__).resolve().parents[1]
NC_FILE = REPO_ROOT / "Output" / "cmems_cbpm_monthly_0p25deg.nc"
OUT_GRID = REPO_ROOT / "Output" / "cmems_npp_climatology_gridded.csv"
OUT_DOMAIN = REPO_ROOT / "Output" / "cmems_npp_climatology_domain_mean.csv"

# Requested box (inclusive).
LON_MIN, LON_MAX = -40.0, -10.0
LAT_MIN, LAT_MAX = 58.0, 64.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="nc_in", type=Path, default=NC_FILE)
    ap.add_argument("--out-grid", type=Path, default=OUT_GRID)
    ap.add_argument("--out-domain", type=Path, default=OUT_DOMAIN)
    args = ap.parse_args()

    ds = xr.open_dataset(args.nc_in)
    npp = ds["npp_int"].sel(
        longitude=slice(LON_MIN, LON_MAX),
        latitude=slice(LAT_MIN, LAT_MAX),
    )
    print(f"Box: lon {LON_MIN}..{LON_MAX}, lat {LAT_MIN}..{LAT_MAX}  "
          f"-> {npp.sizes['latitude']} lat x {npp.sizes['longitude']} lon, "
          f"{npp.sizes['time']} months "
          f"({str(ds.time.values.min())[:7]}..{str(ds.time.values.max())[:7]})")

    grp = npp.groupby("time.month")
    clim_mean = grp.mean("time", skipna=True)          # (month, lat, lon)
    clim_std = grp.std("time", skipna=True)
    n_years = grp.count("time")                         # valid years per cell

    # --- gridded climatology CSV ---
    grid = xr.Dataset(
        {
            "npp_int_mean": clim_mean,
            "npp_int_std": clim_std,
            "n_years": n_years,
        }
    ).to_dataframe().reset_index()
    grid = grid.dropna(subset=["npp_int_mean"])
    grid = grid[["month", "latitude", "longitude",
                 "npp_int_mean", "npp_int_std", "n_years"]]
    grid = grid.sort_values(["month", "latitude", "longitude"])
    grid["npp_int_mean"] = grid["npp_int_mean"].round(3)
    grid["npp_int_std"] = grid["npp_int_std"].round(3)
    args.out_grid.parent.mkdir(parents=True, exist_ok=True)
    grid.to_csv(args.out_grid, index=False)
    print(f"Wrote {args.out_grid}  ({len(grid)} rows)")

    # --- domain-mean seasonal cycle CSV ---
    # Spatial stats over the climatological field for each calendar month.
    rows = []
    for m in range(1, 13):
        field = clim_mean.sel(month=m).values
        vals = field[np.isfinite(field)]
        if vals.size == 0:
            rows.append(dict(month=m, npp_int_mean=np.nan, npp_int_std=np.nan,
                             npp_int_p10=np.nan, npp_int_p90=np.nan, n_cells=0))
            continue
        rows.append(dict(
            month=m,
            npp_int_mean=round(float(np.mean(vals)), 3),
            npp_int_std=round(float(np.std(vals)), 3),
            npp_int_p10=round(float(np.percentile(vals, 10)), 3),
            npp_int_p90=round(float(np.percentile(vals, 90)), 3),
            n_cells=int(vals.size),
        ))
    domain = pd.DataFrame(rows)
    domain.to_csv(args.out_domain, index=False)
    print(f"Wrote {args.out_domain}")
    print(domain.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
