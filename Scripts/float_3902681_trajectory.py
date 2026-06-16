"""
Trajectory + dawn/dusk sampling figure for BGC-Argo float 3902681.

Panels:
  A  trajectory map, coloured by date (Iceland Basin / Irminger Sea).
  B  local-solar-time of profiles, zoomed on a typical ~2-week stretch,
     showing the cycle-to-cycle dawn <-> dusk flip in the real data.
  C  schematic of one dawn-dusk cycle: float depth vs time over 2.5 days,
     with the day/night light cycle, surfacing at dawn then dusk.

Run from repo root:
  C:/Users/petit/anaconda3/envs/cmts_learn_olci/python.exe Scripts/float_3902681_trajectory.py
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
fig = plt.figure(figsize=(14, 8.2))
gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1.0], height_ratios=[1.0, 1.0],
                      wspace=0.16, hspace=0.42,
                      left=0.045, right=0.965, top=0.90, bottom=0.085)

proj = ccrs.LambertConformal(central_longitude=-22, central_latitude=61)
pc = ccrs.PlateCarree()

# === Panel A: trajectory map =================================================
axm = fig.add_subplot(gs[:, 0], projection=proj)
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
axm.scatter(df.lon.iloc[0], df.lat.iloc[0], marker="*", s=320, c="#1a9850",
            transform=pc, zorder=5, edgecolor="k", linewidth=0.6)
axm.scatter(df.lon.iloc[-1], df.lat.iloc[-1], marker="s", s=90, c="#d73027",
            transform=pc, zorder=5, edgecolor="k", linewidth=0.6)
axm.annotate("deploy\n" + str(t0.date()), (df.lon.iloc[0], df.lat.iloc[0]),
             xytext=(8, -16), textcoords="offset points", transform=pc,
             fontsize=8, color="#1a6b35", fontweight="bold")
axm.annotate("last\n" + str(t1.date()), (df.lon.iloc[-1], df.lat.iloc[-1]),
             xytext=(8, 6), textcoords="offset points", transform=pc,
             fontsize=8, color="#a01e16", fontweight="bold")

cb = fig.colorbar(lc, ax=axm, orientation="horizontal", pad=0.05,
                  shrink=0.82, aspect=34)
cb.set_label("Profile date", fontsize=9)
cb.ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
cb.ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
cb.ax.tick_params(labelsize=7.5)

axm.set_title(f"BGC-Argo float {FLOAT}  ·  {n} profiles  ·  "
              f"Iceland Basin / Irminger Sea",
              fontsize=12, fontweight="bold", pad=8)

# === Panel B: zoom on a typical cycle (real data) ============================
axt = fig.add_subplot(gs[0, 1])
w0, w1 = pd.Timestamp("2025-04-21"), pd.Timestamp("2025-05-03")
win = df[(df.time >= w0) & (df.time <= w1)]
axt.axhspan(0, 12, color=cdawn, alpha=0.07)
axt.axhspan(12, 24, color=cdusk, alpha=0.07)
axt.axhline(12, color="#999", lw=0.7, ls="--")
axt.plot(win.time, win.lst, "-", color="#bbb", lw=1.0, zorder=2)
md = win.phase == "dawn"
axt.scatter(win.time[md], win.lst[md], s=55, c=cdawn, zorder=4,
            edgecolor="k", linewidth=0.4, label="dawn profile (~06 h)")
axt.scatter(win.time[~md], win.lst[~md], s=55, c=cdusk, zorder=4,
            edgecolor="k", linewidth=0.4, label="dusk profile (~18 h)")
axt.set_ylim(0, 24)
axt.set_yticks([0, 6, 12, 18, 24])
axt.set_ylabel("Local solar time  (h)")
axt.set_title("Profiles flip dawn ↔ dusk each cycle  (zoom)",
              fontsize=11, fontweight="bold")
axt.legend(loc="center left", fontsize=8, framealpha=0.92, handletextpad=0.3)
axt.xaxis.set_major_locator(mdates.DayLocator(interval=3))
axt.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
axt.tick_params(labelsize=8)
axt.grid(True, axis="y", alpha=0.25, lw=0.5)

# === Panel C: schematic of one dawn-dusk cycle ===============================
# Idealised: 2.5-day (60 h) window. Float parks at depth, surfaces to profile
# at dawn (light on) then at dusk (light off), matching the ~37 h / ~11 h gaps.
axs = fig.add_subplot(gs[1, 1])
PARK = 1000.0          # parking depth (m)
H = 60.0               # window length (h) = 2.5 days

# profile surfacings (h): dawn -> +36 h dusk -> +12 h dawn
t_dawn1, t_dusk, t_dawn2 = 6.0, 42.0, 54.0

def vee(tc, half=1.2):
    """ascend-to-surface / descend-to-park V around a surfacing time tc."""
    return ([tc - half, tc, tc + half], [PARK, 0.0, PARK])

# build the depth path across the window
tx, dz = [0.0], [PARK]
for tc in (t_dawn1, t_dusk, t_dawn2):
    xs, ys = vee(tc)
    tx += xs; dz += ys
tx += [H]; dz += [PARK]

axs.set_xlim(0, H)
axs.set_ylim(PARK + 90, -230)        # headroom above 0 = "sky"

# day / night shading (sun up 06-18 local each day) across full height
for d in range(3):
    axs.axvspan(d * 24 + 6, d * 24 + 18, color="#ffe39a", alpha=0.55, zorder=0)
    axs.axvspan(d * 24 + 18, d * 24 + 30, color="#23355c", alpha=0.12, zorder=0)
axs.axvspan(0, 6, color="#23355c", alpha=0.12, zorder=0)

# sun / moon glyphs in the sky band
for d in range(3):
    axs.text(d * 24 + 12, -185, "☀", ha="center", va="center", fontsize=15,
             color="#e8a300", zorder=2)
for c in (0, 24, 48):
    axs.text(c, -185, "☾", ha="center", va="center", fontsize=11,
             color="#5b6b8c", zorder=2)

axs.plot(tx, dz, "-", color="#333", lw=2.0, zorder=3)
axs.axhline(0, color="#2e9e6f", lw=1.2, ls=":", zorder=1)   # sea surface

# surfacing markers + labels just above the surface
for tc, ph in [(t_dawn1, "dawn"), (t_dusk, "dusk"), (t_dawn2, "dawn")]:
    col = cdawn if ph == "dawn" else cdusk
    axs.scatter([tc], [0], s=150, c=col, edgecolor="k", linewidth=0.5,
                zorder=5, clip_on=False)
    axs.annotate(f"{ph}\nprofile", (tc, 0), xytext=(0, 9),
                 textcoords="offset points", ha="center", va="bottom",
                 fontsize=8.5, fontweight="bold", color=col, zorder=6)

axs.text(30, 905, "drift / park at ~1000 m", ha="center", fontsize=8,
         color="#555", style="italic", zorder=4)

# 2.5-day span bar
axs.annotate("", xy=(0, 1055), xytext=(H, 1055),
             arrowprops=dict(arrowstyle="<->", color="#333", lw=1.3))
axs.text(H / 2, 1085, "one dawn–dusk cycle  ≈  2.5 days",
         ha="center", va="top", fontsize=9, fontweight="bold")

axs.set_ylabel("Depth  (m)")
axs.set_xlabel("Time  (hours)")
axs.set_xticks(np.arange(0, H + 1, 12))
axs.set_yticks([0, 250, 500, 750, 1000])
axs.set_title("Dawn–dusk cycle schematic", fontsize=11, fontweight="bold")
axs.tick_params(labelsize=8)

fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
print("saved ->", OUT)
