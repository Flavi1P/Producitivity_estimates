"""
Materialise monthly chl(z) and Cphyto(z) on the CMEMS 0.25 deg grid
at native depths (22 levels in 0..200 m), matching the time axis of
Output/cmems_cbpm_monthly_0p25deg.nc.

Companion to Scripts/cmems_3d_cbpm.py: that script already computes
both fields in-memory per cell but only writes the CbPM outputs.  This
helper persists chl(z) and Cphyto(z) so downstream FDA analyses
(Scripts/npp_bioregions_depth.R) can read them without rerunning CbPM.

Run:
    & "C:/Users/petit/miniconda3/envs/UVP6/python.exe" Scripts/cmems_chl_cphyto_monthly.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr

REPO_ROOT = Path(__file__).resolve().parents[1]
BGC_FILE = (
    REPO_ROOT
    / "Data"
    / "Raw"
    / "cmems_obs-mob_glo_bgc-chl-poc_my_0.25deg_P7D-m_multi-vars_"
      "39.88W-10.12W_58.12N-64.88N_0.00-200.00m_2015-01-07-2023-12-27.nc"
)
OUT_FILE = REPO_ROOT / "Output" / "cmems_chl_cphyto_monthly_0p25deg.nc"


def bbp_to_cphyto(bbp: xr.DataArray) -> xr.DataArray:
    """Same coefficients as Scripts/cmems_3d_cbpm.py bbp_to_cphyto."""
    return (bbp / (470.0 / 440.0)) * 12128.0 + 0.59


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="nc_in", type=Path, default=BGC_FILE)
    ap.add_argument("--out", type=Path, default=OUT_FILE)
    args = ap.parse_args()

    print(f"Loading {args.nc_in.name}")
    ds = xr.open_dataset(args.nc_in)[["chl", "bbp"]]

    print("Aggregating weekly -> monthly mean")
    ds_m = ds.resample(time="1MS").mean(skipna=True)

    chl = ds_m["chl"].astype(np.float32)
    cphyto = bbp_to_cphyto(ds_m["bbp"]).astype(np.float32)
    cphyto.attrs.clear()

    out = xr.Dataset(
        data_vars=dict(
            chl=(("time", "depth", "latitude", "longitude"), chl.data,
                 {"units": "mg m-3",
                  "long_name": "Chlorophyll-a concentration"}),
            cphyto=(("time", "depth", "latitude", "longitude"), cphyto.data,
                    {"units": "mg C m-3",
                     "long_name": "Phytoplankton carbon (from bbp)"}),
        ),
        coords=dict(
            time=("time", ds_m["time"].values),
            depth=("depth", ds_m["depth"].values.astype(np.float32),
                   {"units": "m", "positive": "down"}),
            latitude=("latitude", ds_m["latitude"].values.astype(np.float32)),
            longitude=("longitude",
                       ds_m["longitude"].values.astype(np.float32)),
        ),
        attrs=dict(
            title=("Monthly 0.25-deg Chla and Cphyto at native depths "
                   "(CMEMS MULTIOBS BGC)"),
            source=str(args.nc_in.name),
            cphyto_formula=(
                "Cphyto = bbp / (470/440) * 12128 + 0.59 "
                "(matches Scripts/cmems_3d_cbpm.py)"),
        ),
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_netcdf(args.out)
    print(f"Wrote {args.out}  "
          f"({out['chl'].shape} chl, {out['cphyto'].shape} cphyto)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
