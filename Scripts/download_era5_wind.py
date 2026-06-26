# Downloads ERA5 hourly 10 m wind + MSL pressure over the float-3902681 box.
#
# The Copernicus CDS rejects a single 17-month hourly request as "too large",
# so this downloads ONE FILE PER MONTH into Data/Raw/ERA5/chunks/ and then
# concatenates them along time into the single file the R budget script reads:
#   Data/Raw/ERA5/era5_wind_slp_3902681.nc
#
# Requires a ~/.cdsapirc with Copernicus CDS credentials (manual setup):
#   url: https://cds.climate.copernicus.eu/api
#   key: <YOUR-API-KEY>
# and you must accept the dataset licence once on the CDS website:
#   https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels
#
# Run with the project download env (from the repo root):
#   C:/Users/petit/anaconda3/envs/cmts_learn_olci/python.exe Scripts/download_era5_wind.py
import os
import cdsapi
import xarray as xr

OUT_DIR = "Data/Raw/ERA5"
CHUNK_DIR = os.path.join(OUT_DIR, "chunks")
OUT_FILE = os.path.join(OUT_DIR, "era5_wind_slp_3902681.nc")
os.makedirs(CHUNK_DIR, exist_ok=True)

# (year, month) covering the float record 2024-10-02 -> 2026-02-10,
# padded by one month at each end.
MONTHS = (
    [(2024, m) for m in range(10, 13)]
    + [(2025, m) for m in range(1, 13)]
    + [(2026, m) for m in (1, 2)]
)

VARIABLES = [
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "mean_sea_level_pressure",
]
AREA = [65, -40, 58, -10]   # N, W, S, E
TIMES = [f"{h:02d}:00" for h in range(24)]

c = cdsapi.Client()
chunk_files = []
for year, month in MONTHS:
    fn = os.path.join(CHUNK_DIR, f"era5_{year}_{month:02d}.nc")
    chunk_files.append(fn)
    if os.path.exists(fn):
        print(f"skip (exists): {fn}")
        continue
    ndays = 31  # CDS silently ignores days that don't exist in the month
    print(f"downloading {year}-{month:02d} ...")
    c.retrieve(
        "reanalysis-era5-single-levels",
        {
            "product_type": "reanalysis",
            "variable": VARIABLES,
            "year": str(year),
            "month": f"{month:02d}",
            "day": [f"{d:02d}" for d in range(1, ndays + 1)],
            "time": TIMES,
            "area": AREA,
            "format": "netcdf",
        },
        fn,
    )
    print(f"  wrote {fn}")

# ---- concatenate monthly chunks along the time dimension --------------------
# NOTE: the raw CDS chunks carry an `expver` variable-length string coordinate
# (and a scalar `number`). R's ncdf4 SEGFAULTS opening a file that contains a
# vlen-string variable, so we drop those auxiliary coords and write a clean
# NETCDF4_CLASSIC file with only u10/v10/msl + numeric lon/lat/time. The time
# axis is written as float64 "hours since 1900-01-01" (classic ERA5 encoding)
# so the R side can decode it with a fixed origin.
import numpy as np

print("concatenating monthly chunks ...")
ds = xr.open_mfdataset(chunk_files, combine="by_coords")
tname = "valid_time" if "valid_time" in ds.coords else "time"
ds = ds.sortby(tname)
nsteps = ds.sizes.get(tname, "?")

# strip the auxiliary coords/vars that crash ncdf4
for aux in ("expver", "number"):
    if aux in ds.coords:
        ds = ds.reset_coords(aux, drop=True)
    if aux in ds.variables:
        ds = ds.drop_vars(aux)

# rename the time dim to `time` and force a clean CF numeric encoding
ds = ds.rename({tname: "time"})
ds = ds[["u10", "v10", "msl"]]
enc = {
    "time": {"dtype": "float64", "units": "hours since 1900-01-01",
             "calendar": "gregorian"},
    "u10": {"dtype": "float32", "zlib": False},
    "v10": {"dtype": "float32", "zlib": False},
    "msl": {"dtype": "float32", "zlib": False},
}
ds.to_netcdf(OUT_FILE, format="NETCDF4_CLASSIC", encoding=enc)
ds.close()
print(f"wrote {OUT_FILE} (time = {nsteps} steps)")
