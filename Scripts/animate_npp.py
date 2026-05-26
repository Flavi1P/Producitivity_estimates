"""
Animate monthly column-integrated NPP from the gridded CbPM output.

Reads Output/cmems_cbpm_monthly_0p25deg.nc and writes a GIF of npp_int
month by month (2015-2023, 108 frames).

Run:
    & "C:/Users/petit/miniconda3/envs/UVP6/python.exe" Scripts/animate_npp.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.animation import FuncAnimation, PillowWriter

REPO_ROOT = Path(__file__).resolve().parents[1]
NC_FILE = REPO_ROOT / "Output" / "cmems_cbpm_monthly_0p25deg.nc"
GIF_OUT = REPO_ROOT / "Output" / "cbpm_npp_monthly.gif"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="nc_in", type=Path, default=NC_FILE)
    ap.add_argument("--out", type=Path, default=GIF_OUT)
    ap.add_argument("--fps", type=int, default=4)
    ap.add_argument("--vmax", type=float, default=1600.0,
                    help="Color scale cap (mg C m-2 d-1). Default 1600 ~= p95.")
    args = ap.parse_args()

    ds = xr.open_dataset(args.nc_in)
    npp = ds["npp_int"]                       # (time, lat, lon)
    times = pd.to_datetime(ds["time"].values)
    lats = ds["latitude"].values
    lons = ds["longitude"].values

    # Domain-mean for the inset timeline
    domain_mean = npp.mean(dim=("latitude", "longitude"), skipna=True).values

    # Colormap: cmocean.algae (sequential green for productivity); NaN -> grey
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

    # Use first month for initial image; will update with set_array each frame
    img = ax.pcolormesh(
        lons, lats, npp.isel(time=0).values,
        cmap=cmap, vmin=0.0, vmax=args.vmax,
        shading="nearest", transform=proj, zorder=1,
    )
    cbar = fig.colorbar(img, ax=ax, orientation="vertical",
                        pad=0.02, shrink=0.85, extend="max")
    cbar.set_label("Column-integrated NPP (mg C m$^{-2}$ d$^{-1}$)",
                   fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    title = ax.set_title("", fontsize=11)

    # Bottom timeline: domain-mean NPP with a moving cursor
    ax_ts.plot(times, domain_mean, color="#1f6e3f", linewidth=1.2)
    ax_ts.set_ylabel("Domain mean\n(mg C m$^{-2}$ d$^{-1}$)", fontsize=8)
    ax_ts.set_xlim(times[0], times[-1])
    ax_ts.set_ylim(0, np.nanmax(domain_mean) * 1.1)
    ax_ts.tick_params(labelsize=8)
    ax_ts.grid(alpha=0.25)
    cursor = ax_ts.axvline(times[0], color="#c0392b", linewidth=1.5)

    fig.suptitle("CMEMS 3D BGC → CbPM   |   Iceland Basin   |   0.25° monthly",
                 fontsize=11, y=0.98)

    def update(frame_idx: int):
        arr = npp.isel(time=frame_idx).values
        # pcolormesh with shading='nearest' uses a flat array, masked for NaN
        masked = np.ma.masked_invalid(arr)
        img.set_array(masked.ravel())
        title.set_text(f"{pd.Timestamp(times[frame_idx]):%B %Y}")
        cursor.set_xdata([times[frame_idx], times[frame_idx]])
        return img, title, cursor

    n_frames = len(times)
    print(f"Rendering {n_frames} frames -> {args.out}")
    anim = FuncAnimation(fig, update, frames=n_frames, interval=1000 / args.fps,
                         blit=False)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    anim.save(args.out, writer=PillowWriter(fps=args.fps))
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
