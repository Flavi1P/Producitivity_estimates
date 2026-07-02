"""
Trajectory + dawn/dusk sampling figure for BGC-Argo float 3902681.

Panels:
  A  trajectory map, coloured by date (Iceland Basin / Irminger Sea).
  C  schematic of the dusk -> pre-dawn night bracket: float depth vs time,
     with the day/night light cycle, isolating the ~11 h night O2-loss unit.

Run from repo root:
  C:/Users/flapet/AppData/Local/miniforge3/envs/productivity_py/python.exe Scripts/float_3902681_trajectory.py
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection

FLOAT = "3902681"
SPROF = f"Data/Raw/Floats/{FLOAT}_Sprof.nc"
OUT = f"Output/float_{FLOAT}_trajectory.png"

import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
ds = xr.open_dataset(SPROF)
juld = pd.to_datetime(ds["JULD"].values)
lat = ds["LATITUDE"].values.astype(float)
lon = ds["LONGITUDE"].values.astype(float)
cyc = ds["CYCLE_NUMBER"].values.astype(float)

df = pd.DataFrame({"time": juld, "lat": lat, "lon": lon, "cyc": cyc})
df = df.dropna(subset=["time", "lat", "lon"]).sort_values("time").reset_index(drop=True)

# local solar time (clock) -> dawn (<12) vs dusk (>=12)
lst = (df["time"].dt.hour + df["time"].dt.minute / 60 + df["lon"] / 15.0) % 24
df["lst"] = lst
df["phase"] = np.where(lst < 12, "dawn", "dusk")

t0, t1 = df["time"].min(), df["time"].max()
span_days = (t1 - t0).days
n = len(df)
ntrans = int((df["phase"].values[1:] != df["phase"].values[:-1]).sum())
print(f"{n} profiles | {t0.date()} -> {t1.date()} ({span_days} d)")
print(f"dawn={int((df.phase=='dawn').sum())} dusk={int((df.phase=='dusk').sum())} "
      f"| {ntrans}/{n-1} alternations")

tnum = mdates.date2num(df["time"])
cdawn, cdusk = "#f5a623", "#3b5b8c"

# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.edgecolor": "#444", "axes.linewidth": 0.8,
})
fig = plt.figure(figsize=(14, 5.8))
gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0],
                      wspace=0.16,
                      left=0.045, right=0.965, top=0.90, bottom=0.10)

proj = ccrs.LambertConformal(central_longitude=-22, central_latitude=61)
pc = ccrs.PlateCarree()

# === Panel A: trajectory map =================================================
axm = fig.add_subplot(gs[0], projection=proj)
pad_x, pad_y = 1.6, 0.6
axm.set_extent([df.lon.min() - pad_x, df.lon.max() + pad_x,
                df.lat.min() - pad_y, df.lat.max() + pad_y], crs=pc)

axm.add_feature(cfeature.NaturalEarthFeature(
    "physical", "land", "50m", facecolor="#e8e4dc", edgecolor="#9a9486",
    linewidth=0.5), zorder=1)
axm.add_feature(cfeature.NaturalEarthFeature(
    "physical", "ocean", "50m", facecolor="#eef4f7"), zorder=0)

gl = axm.gridlines(draw_labels=True, linewidth=0.4, color="#bbb",
                   alpha=0.7, linestyle=":")
gl.top_labels = gl.right_labels = False
gl.xlabel_style = gl.ylabel_style = {"size": 8, "color": "#555"}

pts = np.column_stack([df.lon, df.lat])
segs = np.concatenate([pts[:-1, None], pts[1:, None]], axis=1)
lc = LineCollection(segs, cmap="viridis", transform=pc, zorder=2,
                    linewidth=1.6, alpha=0.85)
lc.set_array(tnum[:-1])
axm.add_collection(lc)
axm.scatter(df.lon, df.lat, c=tnum, cmap="viridis", s=11, transform=pc,
            zorder=3, edgecolor="white", linewidth=0.25)
# Neutral start/end markers (white fill, black edge) so their colour implies no
# date on the viridis colourbar; shape alone distinguishes deploy (star) vs last
# (square).
axm.scatter(df.lon.iloc[0], df.lat.iloc[0], marker="*", s=340, c="white",
            transform=pc, zorder=5, edgecolor="k", linewidth=1.0)
axm.scatter(df.lon.iloc[-1], df.lat.iloc[-1], marker="s", s=95, c="white",
            transform=pc, zorder=5, edgecolor="k", linewidth=1.0)
# Labels offset into open water with leader lines so they don't sit on the dense
# Iceland-shelf track tangle.
_lbl_bbox = dict(boxstyle="round,pad=0.25", fc="white", ec="#888", lw=0.6)
_leader = dict(arrowstyle="-", color="#555", lw=0.8, shrinkA=2, shrinkB=7)
axm.annotate("deploy\n" + str(t0.date()), (df.lon.iloc[0], df.lat.iloc[0]),
             xytext=(-52, -46), textcoords="offset points", transform=pc,
             fontsize=8, color="#222", fontweight="bold",
             ha="center", va="center", bbox=_lbl_bbox, arrowprops=_leader)
axm.annotate("last\n" + str(t1.date()), (df.lon.iloc[-1], df.lat.iloc[-1]),
             xytext=(50, 40), textcoords="offset points", transform=pc,
             fontsize=8, color="#222", fontweight="bold",
             ha="center", va="center", bbox=_lbl_bbox, arrowprops=_leader)

cb = fig.colorbar(lc, ax=axm, orientation="horizontal", pad=0.05,
                  shrink=0.82, aspect=34)
cb.set_label("Profile date", fontsize=9)
cb.ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
cb.ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
cb.ax.tick_params(labelsize=7.5)

axm.set_title(f"BGC-Argo float {FLOAT}  ·  {n} profiles  ·  "
              f"Iceland Basin / Irminger Sea",
              fontsize=12, fontweight="bold", pad=8)

# === Panel C: the dusk -> pre-dawn night bracket =============================
# The respiration (night O2-loss) signal comes from a dusk profile paired with
# the immediately following pre-dawn profile ~11 h later (dusk->dawn, < 0.7 d;
# see night_tbl in o2_budget_3902681_fullmld.R). The O2 drawdown between those
# two casts, over the intervening night, IS night-time respiration. This panel
# makes that dusk->pre-dawn pair the obvious unit; daytime parking merely
# separates one night bracket from the next.
axs = fig.add_subplot(gs[1])
PARK = 1000.0          # parking depth (m)
H = 48.0               # window length (h) = 2 days of context

# profile surfacings (h): dusk (18 h) -> pre-dawn (29 h, ~11 h later) -> next
# dusk (42 h). The dusk1 -> dawn pair (18->29 h) is the highlighted night unit.
t_dusk1, t_dawn, t_dusk2 = 18.0, 29.0, 42.0

def vee(tc, half=1.1):
    """ascend-to-surface / descend-to-park V around a surfacing time tc."""
    return ([tc - half, tc, tc + half], [PARK, 0.0, PARK])

# build the depth path across the window
tx, dz = [0.0], [PARK]
for tc in (t_dusk1, t_dawn, t_dusk2):
    xs, ys = vee(tc)
    tx += xs; dz += ys
tx += [H]; dz += [PARK]

axs.set_xlim(0, H)
axs.set_ylim(PARK + 90, -230)        # headroom above 0 = "sky"

# day / night shading (sun up 06-18 local each day) across full height
for d in range(3):
    axs.axvspan(d * 24 + 6, d * 24 + 18, color="#ffe39a", alpha=0.55, zorder=0)
    axs.axvspan(d * 24 + 18, d * 24 + 30, color="#23355c", alpha=0.10, zorder=0)
axs.axvspan(0, 6, color="#23355c", alpha=0.10, zorder=0)

# emphasise the dusk -> pre-dawn NIGHT bracket (the respiration unit)
axs.axvspan(t_dusk1, t_dawn, color="#23355c", alpha=0.20, zorder=0.4)

# sun / moon glyphs in the sky band
for c in (12, 36):
    axs.text(c, -185, "☀", ha="center", va="center", fontsize=15,
             color="#e8a300", zorder=2)
axs.text(24, -185, "☾", ha="center", va="center", fontsize=12,
         color="#cdd6ea", zorder=2)

axs.plot(tx, dz, "-", color="#333", lw=2.0, zorder=3)
axs.axhline(0, color="#2e9e6f", lw=1.2, ls=":", zorder=1)   # sea surface

# surfacing markers + labels just above the surface
for tc, ph, lbl in [(t_dusk1, "dusk", "dusk\nprofile"),
                    (t_dawn, "dawn", "pre-dawn\nprofile"),
                    (t_dusk2, "dusk", "dusk\nprofile")]:
    col = cdawn if ph == "dawn" else cdusk
    axs.scatter([tc], [0], s=150, c=col, edgecolor="k", linewidth=0.5,
                zorder=5, clip_on=False)
    axs.annotate(lbl, (tc, 0), xytext=(0, 9),
                 textcoords="offset points", ha="center", va="bottom",
                 fontsize=8.5, fontweight="bold", color=col, zorder=6)

# what the bracket measures
axs.annotate(r"$\Delta$O$_2$ = night respiration",
             ((t_dusk1 + t_dawn) / 2, 175), ha="center", va="center",
             fontsize=8.5, fontweight="bold", color="#1b2a4a", zorder=6)

# park label in the clear day-time drift segment between the two nights
axs.text((t_dawn + t_dusk2) / 2, 905, "drift / park at ~1000 m", ha="center",
         fontsize=8, color="#555", style="italic", zorder=4)

# night-bracket span bar (~11 h)
axs.annotate("", xy=(t_dusk1, 1055), xytext=(t_dawn, 1055),
             arrowprops=dict(arrowstyle="<->", color="#1b2a4a", lw=1.5))
axs.text((t_dusk1 + t_dawn) / 2, 1085,
         "night bracket  ≈  11 h", ha="center", va="top",
         fontsize=9, fontweight="bold", color="#1b2a4a")

axs.set_ylabel("Depth  (m)")
axs.set_xlabel("Time  (hours)")
axs.set_xticks(np.arange(0, H + 1, 12))
axs.set_yticks([0, 250, 500, 750, 1000])
axs.set_title("Night O₂-loss bracket  (dusk → pre-dawn)",
              fontsize=11, fontweight="bold")
axs.tick_params(labelsize=8)

fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
print("saved ->", OUT)
