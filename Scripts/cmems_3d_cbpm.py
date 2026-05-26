"""
Gridded CbPM applied to the CMEMS 3D BGC product (chl, bbp, PAR).

Inputs
------
Data/Raw/cmems_obs-mob_glo_bgc-chl-poc_my_0.25deg_P7D-m_multi-vars_*.nc
    Weekly, 0.25 deg, 22 depth levels (0..200 m), 2015-01-07..2023-12-27.
    Variables used: chl, bbp, PAR.

What it does
------------
1. Resample weekly -> monthly mean (keeps the same 0.25 deg grid).
2. Convert bbp -> Cphyto using the same formula as Scripts/format_for_ncp.R:
       Cphyto = bbp / (470/440) * 12128 + 0.59       [mg C m-3]
3. Interpolate chl(z) and Cphyto(z) from the 22 native depths onto a
   200-point 1-m grid (z = 0..199 m).  Required by cbpm_argo().
4. Convert each cell's surface PAR (mid-month, monthly mean) to a daily
   light integral (mol photons m-2 d-1) using a pvlib clear-sky scaling
   at noon UTC.  The clear-sky curve depends only on (date, lat) so we
   cache one curve per (year, month, lat).
5. Call Scripts/cbpm_py/cbpm_argo for every (time, lat, lon) cell,
   collect pp_z, mu_z, mzeu, integrate the column.
6. Write Output/cmems_cbpm_monthly_0p25deg.nc.

Run
---
    # quick sanity pass (one month, decimated grid):
    python Scripts/cmems_3d_cbpm.py --test

    # full run:
    python Scripts/cmems_3d_cbpm.py
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pvlib
import xarray as xr

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "Scripts"))
from cbpm_py.cbpm_argo import cbpm_argo  # noqa: E402

BGC_FILE = (
    REPO_ROOT
    / "Data"
    / "Raw"
    / "cmems_obs-mob_glo_bgc-chl-poc_my_0.25deg_P7D-m_multi-vars_"
      "39.88W-10.12W_58.12N-64.88N_0.00-200.00m_2015-01-07-2023-12-27.nc"
)
OUT_FILE = REPO_ROOT / "Output" / "cmems_cbpm_monthly_0p25deg.nc"

DEPTH_GRID = np.arange(200, dtype=np.float32)        # 0..199 m, 1 m step
GHI_TO_PAR = 4.57                                    # umol J-1
MIN_PAR_CLEAR = 50.0                                 # umol m-2 s-1 (gate)


def bbp_to_cphyto(bbp: xr.DataArray) -> xr.DataArray:
    """Cphyto [mg C m-3] from bbp using the same coefficients as
    Scripts/format_for_ncp.R (Behrenfeld-style)."""
    return (bbp / (470.0 / 440.0)) * 12128.0 + 0.59


def interp_profile_to_1m(values: np.ndarray, src_depth: np.ndarray) -> np.ndarray:
    """Linear interp of a length-N profile onto the 0..199 m / 1 m grid.

    NaNs at the source are ignored.  Returns NaNs where extrapolation
    would be required (above the shallowest valid sample or below the
    deepest valid sample)."""
    mask = np.isfinite(values)
    if mask.sum() < 2:
        return np.full(DEPTH_GRID.size, np.nan, dtype=np.float32)
    return np.interp(
        DEPTH_GRID,
        src_depth[mask],
        values[mask],
        left=np.nan,
        right=np.nan,
    ).astype(np.float32)


def daily_par_from_noon(par_meas: float, lat: float, date: pd.Timestamp,
                        clearsky_cache: dict) -> float:
    """Convert a (monthly-mean) surface PAR snapshot into a daily light
    integral (mol photons m-2 d-1) by scaling a pvlib clear-sky curve so
    its noon-UTC value equals par_meas, then integrating.

    The clear-sky curve depends only on (date, lat); cache it.
    """
    if not np.isfinite(par_meas) or par_meas <= 0:
        return np.nan

    key = (date.year, date.month, round(float(lat), 4))
    cached = clearsky_cache.get(key)
    if cached is None:
        loc = pvlib.location.Location(float(lat), 0.0, tz="UTC")
        # use the 15th of the month as a representative day
        d = pd.Timestamp(year=date.year, month=date.month, day=15, tz="UTC")
        times = pd.date_range(
            f"{d:%Y-%m-%d} 00:00", f"{d:%Y-%m-%d} 23:59",
            freq="1min", tz="UTC",
        )
        ghi = loc.get_clearsky(times)["ghi"].to_numpy()
        par_clear = ghi * GHI_TO_PAR                # umol m-2 s-1
        par_clear_noon = float(par_clear[12 * 60])  # 12:00 sample
        cached = (par_clear, par_clear_noon)
        clearsky_cache[key] = cached

    par_clear, par_clear_noon = cached
    if not np.isfinite(par_clear_noon) or par_clear_noon < MIN_PAR_CLEAR:
        return np.nan
    scale = par_meas / par_clear_noon
    return float(par_clear.sum() * scale * 60.0 / 1e6)   # mol m-2 d-1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true",
                    help="Process only the first month and a 4x decimated grid.")
    ap.add_argument("--start", type=str, default=None,
                    help="Start date (YYYY-MM-DD), inclusive.")
    ap.add_argument("--end", type=str, default=None,
                    help="End date (YYYY-MM-DD), inclusive.")
    ap.add_argument("--out", type=Path, default=OUT_FILE)
    args = ap.parse_args()

    print(f"Loading {BGC_FILE.name}")
    ds = xr.open_dataset(BGC_FILE)[["chl", "bbp", "PAR"]]
    src_depth = ds["depth"].values.astype(np.float32)

    print("Aggregating weekly -> monthly mean")
    ds_m = ds.resample(time="1MS").mean(skipna=True)

    if args.test:
        ds_m = ds_m.isel(time=slice(5, 6),                   # one mid-year month
                         latitude=slice(None, None, 4),
                         longitude=slice(None, None, 4))
        print(f"--test: subset to time={ds_m.sizes['time']}, "
              f"lat={ds_m.sizes['latitude']}, lon={ds_m.sizes['longitude']}")
    elif args.start or args.end:
        ds_m = ds_m.sel(time=slice(args.start, args.end))
        print(f"Time-sliced to {args.start} .. {args.end}: "
              f"{ds_m.sizes['time']} months")

    chl = ds_m["chl"].astype(np.float32)
    bbp = ds_m["bbp"].astype(np.float32)
    cphyto = bbp_to_cphyto(bbp).astype(np.float32)
    par_surf = ds_m["PAR"].isel(depth=0).astype(np.float32)   # umol m-2 s-1

    times = ds_m["time"].values
    lats = ds_m["latitude"].values.astype(np.float32)
    lons = ds_m["longitude"].values.astype(np.float32)
    nt, nlat, nlon = len(times), len(lats), len(lons)

    npp = np.full((nt, DEPTH_GRID.size, nlat, nlon), np.nan, dtype=np.float32)
    mu = np.full_like(npp, np.nan)
    zeu = np.full((nt, nlat, nlon), np.nan, dtype=np.float32)
    dpar = np.full((nt, nlat, nlon), np.nan, dtype=np.float32)
    npp_int = np.full((nt, nlat, nlon), np.nan, dtype=np.float32)

    clearsky_cache: dict = {}
    n_cells = nt * nlat * nlon
    n_done = 0
    n_ok = 0
    t0 = time.time()

    print(f"Running CbPM on {n_cells} cells "
          f"({nt} months x {nlat} lat x {nlon} lon)")

    chl_arr = chl.values
    cphyto_arr = cphyto.values
    par_surf_arr = par_surf.values

    for it in range(nt):
        t = pd.Timestamp(times[it])
        for ilat in range(nlat):
            lat = float(lats[ilat])
            for ilon in range(nlon):
                n_done += 1

                par_meas = float(par_surf_arr[it, ilat, ilon])
                irr = daily_par_from_noon(par_meas, lat, t, clearsky_cache)
                if not np.isfinite(irr):
                    continue

                chl_z = interp_profile_to_1m(
                    chl_arr[it, :, ilat, ilon], src_depth)
                if not np.isfinite(chl_z[0]):
                    continue

                cphyto_z = interp_profile_to_1m(
                    cphyto_arr[it, :, ilat, ilon], src_depth)
                if not np.isfinite(cphyto_z[0]):
                    continue

                # cbpm_argo mutates chl_z in place (clips negatives) - pass a copy
                pp_z, mu_z, _par_z, _prcnt, _ntf, _igf, mzeu = cbpm_argo(
                    chl_z.copy(),
                    cphyto_z.copy(),
                    irr,
                    int(t.year),
                    int(t.month),
                    int(t.day),
                    lat,
                )

                npp[it, :, ilat, ilon] = pp_z
                mu[it, :, ilat, ilon] = mu_z
                zeu[it, ilat, ilon] = mzeu if np.isfinite(mzeu) else np.nan
                dpar[it, ilat, ilon] = irr
                # column integral (mg C m-3 d-1 * 1 m -> mg C m-2 d-1)
                npp_int[it, ilat, ilon] = np.nansum(pp_z)
                n_ok += 1

        elapsed = time.time() - t0
        rate = n_done / max(elapsed, 1e-6)
        eta = (n_cells - n_done) / max(rate, 1e-6)
        print(f"  t={pd.Timestamp(times[it]):%Y-%m}  "
              f"{n_done}/{n_cells}  ok={n_ok}  "
              f"{rate:.0f} cells/s  ETA {eta/60:.1f} min")

    print(f"Finished {n_ok}/{n_cells} cells with valid CbPM output "
          f"({100*n_ok/n_cells:.1f}%)")

    out = xr.Dataset(
        data_vars=dict(
            npp=(("time", "depth", "latitude", "longitude"), npp,
                 {"units": "mg C m-3 d-1",
                  "long_name": "Depth-resolved net primary production (CbPM)"}),
            mu=(("time", "depth", "latitude", "longitude"), mu,
                {"units": "d-1",
                 "long_name": "Phytoplankton division rate"}),
            zeu=(("time", "latitude", "longitude"), zeu,
                 {"units": "m",
                  "long_name": "1% surface PAR depth"}),
            daily_par=(("time", "latitude", "longitude"), dpar,
                       {"units": "mol photons m-2 d-1",
                        "long_name":
                            "Surface daily light integral (clear-sky scaled)"}),
            npp_int=(("time", "latitude", "longitude"), npp_int,
                     {"units": "mg C m-2 d-1",
                      "long_name":
                          "Column-integrated NPP 0-200 m (sum of pp_z * 1 m)"}),
        ),
        coords=dict(
            time=("time", times),
            depth=("depth", DEPTH_GRID,
                   {"units": "m", "positive": "down"}),
            latitude=("latitude", lats),
            longitude=("longitude", lons),
        ),
        attrs=dict(
            title="Monthly 0.25-deg CbPM NPP from CMEMS BGC 3D product",
            source=str(BGC_FILE.name),
            cbpm_reference=(
                "Westberry et al. 2008; Arteaga's CbPM-Argo adaptation "
                "in Scripts/cbpm_py/cbpm_argo.py"),
            cphyto_formula=(
                "Cphyto = bbp / (470/440) * 12128 + 0.59 "
                "(matches Scripts/format_for_ncp.R)"),
            par_handling=(
                "Surface PAR scaled to daily integral via pvlib clear-sky "
                "GHI*4.57, anchored at noon UTC on mid-month day"),
        ),
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_netcdf(args.out)
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
