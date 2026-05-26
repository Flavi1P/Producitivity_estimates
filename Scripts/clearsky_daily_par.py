"""
Convert instantaneous surface PAR measurements to daily PAR (mol photons m-2 d-1)
by scaling a pvlib clear-sky GHI curve to match the measured PAR at the
measurement time, then integrating over the full UTC day at 1-min resolution.

Algorithm (per row):
    1. Build a 1-min clear-sky GHI series for the measurement's UTC date at (lat, lon).
    2. Convert GHI (W m-2) -> PAR (umol photons m-2 s-1) via factor 4.57.
    3. Look up the clear-sky PAR at the measurement time, scale the curve so the
       value at that time equals the measured PAR.
    4. Sum the scaled curve and convert to daily light integral (mol m-2 d-1).

If the clear-sky PAR at the measurement time is below MIN_PAR_CLEAR
(i.e. measurement is near sunrise/sunset/night), scaling is undefined and
the row is returned as NaN.

CLI:
    python clearsky_daily_par.py --in <input.csv> --out <output.csv>

Input CSV columns:
    id            arbitrary unique row identifier (string)
    datetime_utc  ISO 8601 UTC timestamp of the PAR measurement
    lat           latitude (decimal degrees)
    lon           longitude (decimal degrees)
    par_meas      measured PAR (umol photons m-2 s-1)

Output CSV columns:
    id, daily_par_mol_m2_d
"""

import argparse
import sys

import numpy as np
import pandas as pd
import pvlib

GHI_TO_PAR = 4.57          # umol photons J-1 (Morel & Smith approximation)
MIN_PAR_CLEAR = 50.0       # umol m-2 s-1; below this, scaling is unreliable


def daily_par_from_instant(dt_utc: pd.Timestamp, lat: float, lon: float,
                           par_meas: float) -> float:
    if not np.isfinite(par_meas) or par_meas <= 0:
        return np.nan
    if not np.isfinite(lat) or not np.isfinite(lon):
        return np.nan

    location = pvlib.location.Location(lat, lon, tz="UTC")
    date_str = dt_utc.strftime("%Y-%m-%d")
    times = pd.date_range(f"{date_str} 00:00", f"{date_str} 23:59",
                          freq="1min", tz="UTC")

    clearsky = location.get_clearsky(times)
    par_clear = clearsky["ghi"] * GHI_TO_PAR  # umol m-2 s-1

    t_meas = times[times.get_indexer([dt_utc], method="nearest")[0]]
    par_clear_meas = float(par_clear.loc[t_meas])
    if not np.isfinite(par_clear_meas) or par_clear_meas < MIN_PAR_CLEAR:
        return np.nan

    scale = par_meas / par_clear_meas
    par_scaled = par_clear * scale
    # umol m-2 s-1 summed over 1-min steps -> umol m-2 day-1 -> mol m-2 day-1
    return float(par_scaled.sum() * 60.0 / 1e6)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_csv", required=True)
    ap.add_argument("--out", dest="out_csv", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)
    required = {"id", "datetime_utc", "lat", "lon", "par_meas"}
    missing = required - set(df.columns)
    if missing:
        print(f"Input CSV missing columns: {sorted(missing)}", file=sys.stderr)
        return 2

    df["datetime_utc"] = pd.to_datetime(df["datetime_utc"], utc=True)

    out = []
    for _, r in df.iterrows():
        dli = daily_par_from_instant(
            r["datetime_utc"], float(r["lat"]), float(r["lon"]),
            float(r["par_meas"]),
        )
        out.append({"id": r["id"], "daily_par_mol_m2_d": dli})

    pd.DataFrame(out).to_csv(args.out_csv, index=False)
    n_ok = sum(np.isfinite(o["daily_par_mol_m2_d"]) for o in out)
    print(f"clearsky_daily_par: {n_ok}/{len(out)} rows produced finite daily PAR",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
