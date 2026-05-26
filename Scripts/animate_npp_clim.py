"""
Animate the climatological year of column-integrated NPP from the
gridded CbPM output.

Reads Output/cmems_cbpm_monthly_0p25deg.nc, collapses 2015-2023 monthly
fields to a 12-month climatology (mean per calendar month), and writes
a 12-frame GIF.

Run:
    & "C:/Users/petit/miniconda3/envs/UVP6/python.exe" Scripts/animate_npp_clim.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.animation import FuncAnimation, PillowWriter

REPO_ROOT = Path(__file__).resolve().parents[1]
NC_FILE = REPO_ROOT / "Output" / "cmems_cbpm_monthly_0p25deg.nc"
GIF_OUT = REPO_ROOT / "Output" / "cbpm_npp_climatology.gif"

MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="nc_in", type=Path, default=NC_FILE)
    ap.add_argument("--out", type=Path, default=GIF_OUT)
    ap.add_argument("--fps", type=int, default=2)
    ap.add_argument("--vmax", type=float, default=1600.0,
                    help="Color scale cap (mg C m-2 d-1). Default 1600 ~= p95.")
    args = ap.parse_args()

    ds = xr.open_dataset(args.nc_in)
    # Climatology: mean over years for each calendar month -> (12, lat, lon)
    clim = ds["npp_int"].groupby("time.month").mean("time", skipna=True)
    months = clim["month"].values
    lats = ds["latitude"].values
    lons = ds["longitude"].values

    # Domain-mean seasonal cycle for the bottom panel
    seasonal_mean = clim.mean(dim=("latitude", "longitude"),
                              skipna=True).values

    cmap = cmocean.cm.algae.copy()
    cmap.set_bad("#d0d0d0", alpha=1.0)

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(9, 5.5), dpi=110)
    gs = fig.add_gridspec(2, 1, height_ratios=[4.0, 1.0], hspace=0.30)
    ax = fig.add_subplot(gs[0], projection=proj)
    ax_ts = fig.add_subplot(gs[1])

    extent = [lons.min() - 0.125, lons.max() + 0.125,
              lats.min() - 0.125, lats.max() + 0.125]
    ax.set_extent(extent, crs=proj)
    ax.add_feature(cfeature.LAND.with_scale("50m"),
                   facecolor="#efe8d6", edgecolor="black", linewidth=0.4,
                   zorder=2)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"),
                   linewidth=0.4, zorder=3)
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="grey",
                      alpha=0.4, linestyle="--")
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 8}

    img = ax.pcolormesh(
        lons, lats, clim.isel(month=0).values,
        cmap=cmap, vmin=0.0, vmax=args.vmax,
        shading="nearest", transform=proj, zorder=1,
    )
    cbar = fig.colorbar(img, ax=ax, orientation="vertical",
                        pad=0.02, shrink=0.85, extend="max")
    cbar.set_label("Column-integrated NPP (mg C m$^{-2}$ d$^{-1}$)",
                   fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    title = ax.set_title("", fontsize=11)

    ax_ts.plot(months, seasonal_mean, color="#1f6e3f", linewidth=1.4,
               marker="o", markersize=4)
    ax_ts.set_ylabel("Domain mean\n(mg C m$^{-2}$ d$^{-1}$)", fontsize=8)
    ax_ts.set_xticks(months)
    ax_ts.set_xticklabels(MONTH_NAMES, fontsize=8)
    ax_ts.set_xlim(0.5, 12.5)
    ax_ts.set_ylim(0, np.nanmax(seasonal_mean) * 1.15)
    ax_ts.tick_params(labelsize=8)
    ax_ts.grid(alpha=0.25)
    cursor = ax_ts.axvline(months[0], color="#c0392b", linewidth=1.5)

    fig.suptitle("CMEMS 3D BGC -> CbPM   |   Iceland Basin   |   "
                 "climatological year (2015-2023)",
                 fontsize=11, y=0.98)

    def update(frame_idx: int):
        arr = clim.isel(month=frame_idx).values
        masked = np.ma.masked_invalid(arr)
        img.set_array(masked.ravel())
        title.set_text(MONTH_NAMES[frame_idx])
        cursor.set_xdata([months[frame_idx], months[frame_idx]])
        return img, title, cursor

    n_frames = len(months)
    print(f"Rendering {n_frames} frames -> {args.out}")
    anim = FuncAnimation(fig, update, frames=n_frames,
                         interval=1000 / args.fps, blit=False)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    anim.save(args.out, writer=PillowWriter(fps=args.fps))
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
