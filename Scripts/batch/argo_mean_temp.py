import xarray as xr
import numpy as np
from argopy import DataFetcher as ArgoDataFetcher
from argopy import set_options
import matplotlib.dates as mdates

# Enable caching (important on JASMIN)
set_options(parallel=True, src='gdac')

print("Fetching Argo data...")

argo = ArgoDataFetcher(progress = True)

# North Atlantic bounding box
ds = (
    argo.region([-40, -10, 58, 65, 0, 2000, "2019-01", "2019-02"])
        .to_xarray()
)

print(ds)

# Basic QC (keep only good TEMP)
if "TEMP_QC" in ds:
    ds["TEMP"] = ds.TEMP.where(ds.TEMP_QC == 1)

# Convert time
ds["TIME"] = xr.decode_cf(ds).TIME

# Define depth bins (example)
depth_bins = np.arange(0, 2001, 100)

print("Binning data...")

# Group by year + depth bins
ds_out = (
    ds.argo.groupby_pressure_bins(depth_bins)["TEMP"]
    .groupby(["TIME.year", "STD_PRES_BINS"]).mean()
)

print(ds_out)

print("Saving output...")

ds_out.to_netcdf("Output/argo_mean_temp.nc")

print("Done.")
