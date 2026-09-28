"""
CbPM NPP from the CMEMS 3D BGC product (2025), matched to float 3902681 profiles.

Inputs
------
Data/Raw/2025_NA/2025_NA/cmems_obs-mob_glo_bgc-chl-poc_my_0.25deg_P7D-m_*.nc
    Weekly, 0.25 deg, 36 depth levels (0..1000 m), 2025-02-05..2025-10-29.
    Same product (MULTIOBS_GLO_BGC_3D_MY_015_010) as the 2015-2023 gridded run.
    Read with h5py (see load_model) - use the cmts_learn_olci env.
Data/Processed/productivity_float_with_canyon.csv
    Float 3902681 profiles (semicolon / decimal comma); one position+date per
    prof_number.
Data/Processed/icb_npp_2025.csv
    Float NPP (10-day bins), for the comparison figure only.

What it does
------------
1. Match every float profile to the nearest model week / 0.25 deg cell.  If the
   nearest cell is empty (cloud gap), fall back to the nearest valid cell in
   the surrounding 3x3 box of the same week (`fallback` = True).
2. Per matched profile, same processing as Scripts/cmems_3d_cbpm.py:
   bbp -> Cphyto, chl & Cphyto -> 0..199 m grid, surface PAR -> daily light
   integral (pvlib clear-sky anchored at noon UTC, on the profile's own date),
   then cbpm_argo().
3. Write Data/Processed/cmems_npp_float3902681_2025.csv and
   Output/cmems_vs_float_npp_2025.png.

Caveats
-------
- Weekly-mean surface PAR -> daily dose via the clear-sky noon scaling: same
  assumption as the 2015-2023 run (comparable, not independent of it).
- Cphyto coefficients are the Argo bbp ones, kept identical to the gridded run.
- Model covers 2025-02-05..2025-10-29 only; profiles outside get NaN.

Run (from repo root)
---
    python Scripts/cmems_3d_cbpm_float_matchup.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "Scripts"))
from cbpm_py.cbpm_argo import cbpm_argo  # noqa: E402
from cmems_3d_cbpm import (  # noqa: E402
    bbp_to_cphyto, daily_par_from_noon, interp_profile_to_1m)
from figure_config import NPP_COLOR, GREEN_LIGHT, NCP_COLOR  # noqa: E402

MODEL_GLOB = "cmems_obs-mob_glo_bgc-chl-poc_my_0.25deg_P7D-m_*.nc"
MODEL_DIR = REPO_ROOT / "Data" / "Raw" / "2025_NA" / "2025_NA"
FLOAT_CSV = REPO_ROOT / "Data" / "Processed" / "productivity_float_with_canyon.csv"
FLOAT_NPP_CSV = REPO_ROOT / "Data" / "Processed" / "icb_npp_2025.csv"
OUT_CSV = REPO_ROOT / "Data" / "Processed" / "cmems_npp_float3902681_2025.csv"
OUT_PNG = REPO_ROOT / "Output" / "cmems_vs_float_npp_2025.png"

MAX_DT_DAYS = 4.0          # half a week + margin: beyond this, no match


def load_model(files):
    """Read the weekly NetCDF-4 files with h5py (not netCDF4/xarray): in the
    cmts_learn_olci env netCDF4 and h5py (pulled in by pvlib) ship clashing
    HDF5 DLLs and cannot be loaded in the same process.  Variables are plain
    float32 with _FillValue 99999, no scale/offset."""
    out = {k: [] for k in ("time", "chl", "bbp", "PAR")}
    for f in files:
        with h5py.File(f, "r") as h:
            depth = h["depth"][:].astype(np.float32)
            lat = h["latitude"][:]
            lon = h["longitude"][:]
            out["time"].append(h["time"][:].astype(float))   # h since 1950
            for k in ("chl", "bbp", "PAR"):
                a = h[k][:].astype(np.float32)
                a[a >= 9.9e4] = np.nan
                out[k].append(a)
    t = pd.Timestamp("1950-01-01") + pd.to_timedelta(
        np.concatenate(out["time"]), unit="h")
    return (pd.DatetimeIndex(t), depth, lat, lon,
            *(np.concatenate(out[k]) for k in ("chl", "bbp", "PAR")))


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = (np.sin((lat2 - lat1) / 2) ** 2
         + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2)
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def load_profiles() -> pd.DataFrame:
    d = pd.read_csv(FLOAT_CSV, sep=";", decimal=",",
                    usecols=["prof_number", "lon", "lat", "date"])
    p = d.drop_duplicates("prof_number").reset_index(drop=True)
    p["date"] = pd.to_datetime(p["date"])
    return p


def main() -> int:
    prof = load_profiles()
    files = sorted(MODEL_DIR.glob(MODEL_GLOB))
    print(f"{len(prof)} float profiles, {len(files)} model files")
    mtimes, src_depth, mlats, mlons, chl, bbp, par = load_model(files)

    keep = src_depth <= 200          # only need 0..199 m for CbPM
    src_depth = src_depth[keep]
    chl = chl[:, keep]                               # (t, z, lat, lon)
    cphyto = bbp_to_cphyto(bbp[:, keep])
    par0 = par[:, 0]                                 # (t, lat, lon)

    clearsky_cache: dict = {}
    rows = []
    for _, r in prof.iterrows():
        out = dict(prof_number=r.prof_number, date=r.date.date(),
                   lat=r.lat, lon=r.lon)
        it = int(np.abs(mtimes - r.date).argmin())
        dt = (mtimes[it] - r.date) / pd.Timedelta(days=1)
        ilat = int(np.abs(mlats - r.lat).argmin())
        ilon = int(np.abs(mlons - r.lon).argmin())

        # Nearest cell, else nearest valid cell of the 3x3 box (same week)
        cand = [(ilat + a, ilon + b) for a in (-1, 0, 1) for b in (-1, 0, 1)
                if 0 <= ilat + a < mlats.size and 0 <= ilon + b < mlons.size]
        cand = [c for c in cand
                if np.isfinite(chl[it, 0, c[0], c[1]])
                and np.isfinite(par0[it, c[0], c[1]])]
        if abs(dt) > MAX_DT_DAYS or not cand:
            rows.append(out)
            continue
        dist = [haversine_km(r.lat, r.lon, mlats[a], mlons[b]) for a, b in cand]
        jlat, jlon = cand[int(np.argmin(dist))]

        out.update(model_time=mtimes[it].date(), dt_days=dt,
                   model_lat=float(mlats[jlat]), model_lon=float(mlons[jlon]),
                   match_dist_km=float(np.min(dist)),
                   fallback=(jlat, jlon) != (ilat, ilon))

        irr = daily_par_from_noon(float(par0[it, jlat, jlon]), float(r.lat),
                                  r.date, clearsky_cache, day=None)
        chl_z = interp_profile_to_1m(chl[it, :, jlat, jlon], src_depth)
        cph_z = interp_profile_to_1m(cphyto[it, :, jlat, jlon], src_depth)
        out.update(daily_par=irr, chl_surf=float(chl_z[0]),
                   cphyto_surf=float(cph_z[0]))
        if np.isfinite(irr) and np.isfinite(chl_z[0]) and np.isfinite(cph_z[0]):
            # cbpm_argo mutates chl_z in place (clips negatives) - pass a copy
            pp_z, _mu, _par, _prc, _ntf, _igf, mzeu = cbpm_argo(
                chl_z.copy(), cph_z.copy(), irr,
                int(r.date.year), int(r.date.month), int(r.date.day),
                float(r.lat))
            zeu = float(mzeu) if np.isfinite(mzeu) else np.nan
            out.update(
                zeu=zeu,
                npp_int=float(np.nansum(pp_z)),        # mg C m-2 d-1, 0-200 m
                npp_int_zeu=(float(np.nansum(pp_z[:int(np.ceil(zeu)) + 1]))
                             if np.isfinite(zeu) else np.nan))
        rows.append(out)

    res = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT_CSV, index=False)
    n_ok = int(res["npp_int"].notna().sum())
    print(f"Wrote {OUT_CSV}: {n_ok}/{len(res)} profiles with model NPP "
          f"({int(res['fallback'].fillna(False).astype(bool).sum())} via 3x3 fallback)")

    # --- quick-look figure --------------------------------------------------
    fnpp = pd.read_csv(FLOAT_NPP_CSV, parse_dates=["date_10day"])
    res["date"] = pd.to_datetime(res["date"])
    ok = res.dropna(subset=["npp_int"])
    # model binned onto the float's 10-day bins (nearest bin centre)
    centres = fnpp["date_10day"].values
    ok = ok.assign(date_10day=centres[np.abs(
        ok["date"].values[:, None] - centres[None, :]).argmin(axis=1)])
    mbin = ok.groupby("date_10day")["npp_int"].agg(["mean", "std", "size"])

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.scatter(ok["date"], ok["npp_int"], s=10, color=GREEN_LIGHT,
               label="CMEMS 3D CbPM, per profile")
    ax.errorbar(mbin.index, mbin["mean"], yerr=mbin["std"], color=NPP_COLOR,
                marker="o", ms=4, lw=1.5, capsize=2,
                label="CMEMS 3D CbPM, 10-day mean ± sd")
    ax.plot(fnpp["date_10day"], fnpp["npp"], color=NCP_COLOR, marker="s",
            ms=4, lw=1.5, label="Float 3902681 CbPM (icb_npp_2025.csv)")
    ax.set_ylabel("NPP (mg C m$^{-2}$ d$^{-1}$)")
    ax.set_title("Float 3902681 vs CMEMS 3D BGC (2025) — CbPM NPP at float profiles")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=200)
    print(f"Wrote {OUT_PNG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
