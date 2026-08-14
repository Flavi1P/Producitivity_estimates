"""
What sets the April e-ratio, and why some years sit at the top of the +1 SD band.

Companion diagnostic to eratio_synthetic_year_figure.py. That figure shows the
April e-ratio peak (0.42 +/- 0.27 at DOY 99) as a 9-year mean with an
inter-annual SD ribbon; this script asks which years make the upper bound and
what physical/observational conditions they share.

Method
------
April is taken as DOY 92-120. Per year (2015-2023) the April means of the
per-year e-ratio, NPP and NCP come straight from eratio_per_year.csv (i.e. the
same per-year curves the figure averages). Candidate conditions are:

  zdepth_apr    the NCP integration depth in April, max(MLD, zeu) as used by
                compute_ncp(), from ncp_results.xlsx
  mld_apr_raw   median MLD of the raw float profiles in April - an independent
                check on zdepth_apr that does not go through the pipeline's
                smoothing spline
  frac_drawdown_by_apr
                fraction of the year's total nitrate drawdown (winter inventory
                maximum -> summer minimum) already completed by mid-April, i.e.
                how far into the seasonal cycle April falls that year
  intN_apr      standing nitrate inventory over the integration depth in April
  n_prof / n_float / lon / lat / frac_irminger
                April float sampling, to test the obvious confound that a year
                looks different only because the floats were somewhere else

Result (see the printed correlation table)
------------------------------------------
The April e-ratio is an NCP story: r = 0.97 with April NCP (CV 0.44) against
r = -0.35 with April NPP (CV 0.10). What both terms respond to is how deep the
mixed layer still is in April: deep April mixing depresses CbPM NPP
(r = -0.83, p = 0.005; -0.73 with the raw float MLD) and raises the
column-integrated nitrate drawdown (r = 0.47), so the ratio is pushed up from
both sides. Equivalently, high-e-ratio years are the years in which April still
falls early in the seasonal nitrate drawdown (rho = -0.77 with
frac_drawdown_by_apr). April float sampling explains nothing (|r| < 0.5, all
p > 0.2).

Input:  Output/IcelandIrminger_2015_2025/eratio_per_year.csv   (written by
            eratio_synthetic_year_figure.py - run that first)
        Output/IcelandIrminger_2015_2025/ncp/IcelandIrminger/ncp_results.xlsx
        Data/Processed/doxy_db_*_with_canyon.csv   (optional; float sampling)
Output: Output/IcelandIrminger_2015_2025/eratio_april_drivers.csv
        Output/IcelandIrminger_2015_2025/eratio_april_drivers.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import (            # noqa: E402
    GREEN, BLUE, ERATIO_COLOR, ZERO_LINE, MLD_COLOR,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR   = REPO_ROOT / "Output" / "IcelandIrminger_2015_2025"
ERATIO_CSV = OUT_DIR / "eratio_per_year.csv"
NCP_XLSX   = OUT_DIR / "ncp" / "IcelandIrminger" / "ncp_results.xlsx"
FLOAT_DBS  = sorted((REPO_ROOT / "Data" / "Processed").glob("doxy_db_*_with_canyon.csv"))
OUT_CSV    = OUT_DIR / "eratio_april_drivers.csv"
OUT_PNG    = OUT_DIR / "eratio_april_drivers.png"

APR = (92, 120)          # DOY window: 1-30 April
IRMINGER_LON = -28.0     # crude Irminger / Iceland Basin split

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})


def april_eratio() -> pd.DataFrame:
    er = pd.read_csv(ERATIO_CSV)
    a = er[(er.doy >= APR[0]) & (er.doy <= APR[1])]
    return (a.groupby("year")
             .agg(eratio_apr=("eratio", "mean"),
                  npp_apr=("npp_mmolC", "mean"),
                  ncp_apr=("ncp_mmolC", "mean"),
                  ncp_sd_apr=("ncp_sd_mmolC", "mean"))
             .reset_index())


def ncp_budget_terms(years) -> pd.DataFrame:
    d = pd.read_excel(NCP_XLSX)
    d["date_grid"] = pd.to_datetime(d["date_grid"])
    d["year"] = d.date_grid.dt.year
    d["doy"] = d.date_grid.dt.dayofyear
    # the mld column is mixed-type in the xlsx (floats, strings, and a handful
    # of strings Excel coerced to dates); take only what is genuinely numeric.
    # Every non-numeric cell is a shallow summer value, well after April.
    d["mld_num"] = pd.to_numeric(d["mld"], errors="coerce")
    d["zdepth"] = np.maximum(d["mld_num"], d["zeu"])

    rows = []
    for yr in years:
        y = d[d.year == yr].sort_values("doy")
        a = y[(y.doy >= APR[0]) & (y.doy <= APR[1] + 10)]
        win_max = y.loc[y.doy <= 100, "int_N_smooth"].max()
        sum_min = y.loc[(y.doy > 120) & (y.doy < 250), "int_N_smooth"].min()
        n_apr = a.int_N_smooth.mean()
        rows.append({
            "year": yr,
            "zdepth_apr": a.zdepth.mean(),
            "intN_apr": n_apr,
            "intN_winter_max": win_max,
            "frac_drawdown_by_apr": ((win_max - n_apr) / (win_max - sum_min)
                                     if np.isfinite(sum_min) else np.nan),
        })
    return pd.DataFrame(rows)


def float_sampling(years) -> pd.DataFrame:
    """April float sampling per year, from the CANYON-B float DB (one row per
    profile). Returns an empty frame if the (gitignored) DB is not present."""
    if not FLOAT_DBS:
        return pd.DataFrame({"year": list(years)})
    keep = ["float_wmo", "prof_number", "lon", "lat", "date", "depth", "MLD"]
    parts = []
    for f in FLOAT_DBS:
        for ch in pd.read_csv(f, sep=";", decimal=",", usecols=keep,
                              chunksize=2_000_000):
            ch = ch[ch.depth == 0]
            if len(ch):
                parts.append(ch)
    p = pd.concat(parts, ignore_index=True)
    p["date"] = pd.to_datetime(p["date"])
    p = p.drop_duplicates(["float_wmo", "prof_number", "date"])
    p["year"], p["doy"] = p.date.dt.year, p.date.dt.dayofyear
    a = p[(p.doy >= APR[0]) & (p.doy <= APR[1]) & p.year.isin(years)]
    return (a.groupby("year")
             .agg(n_prof=("prof_number", "size"),
                  n_float=("float_wmo", "nunique"),
                  lon_apr=("lon", "mean"),
                  lat_apr=("lat", "mean"),
                  mld_apr_raw=("MLD", "median"),
                  frac_irminger=("lon", lambda x: float((x < IRMINGER_LON).mean())))
             .reset_index())


def corr_table(t: pd.DataFrame, target: str = "eratio_apr") -> pd.DataFrame:
    rows = []
    for c in t.columns:
        if c in ("year", target) or not np.issubdtype(t[c].dtype, np.number):
            continue
        x, y = t[c].astype(float), t[target].astype(float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 4 or np.ptp(x[ok]) == 0:
            continue
        r, pr = pearsonr(x[ok], y[ok])
        s, ps = spearmanr(x[ok], y[ok])
        # leave-one-year-out range of r, to expose relations carried by one year
        jack = [pearsonr(x[ok].drop(i), y[ok].drop(i))[0] for i in x[ok].index]
        rows.append({"driver": c, "n": int(ok.sum()), "r": r, "p": pr,
                     "rho": s, "p_rho": ps,
                     "r_jack_min": min(jack), "r_jack_max": max(jack)})
    return (pd.DataFrame(rows)
              .reindex(pd.DataFrame(rows)["r"].abs().sort_values(ascending=False).index)
              .reset_index(drop=True))


def _z(v):
    v = np.asarray(v, dtype=float)
    return (v - v.mean()) / v.std(ddof=1)


def _fit_line(ax, x, y, color):
    ok = np.isfinite(x) & np.isfinite(y)
    b, a = np.polyfit(x[ok], y[ok], 1)
    xx = np.linspace(np.min(x[ok]), np.max(x[ok]), 10)
    ax.plot(xx, a + b * xx, color=color, lw=1.2, ls="--", alpha=0.55, zorder=2)
    return pearsonr(x[ok], y[ok])


def _scatter(ax, x, y, years, color, xlabel, ylabel, title, marker="o",
             label=None, annotate=True, fit=True):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ax.scatter(x, y, s=46, color=color, marker=marker, zorder=3,
               edgecolors="white", linewidths=0.8, label=label)
    if annotate:
        for xi, yi, yr in zip(x, y, years):
            ax.annotate(str(yr), (xi, yi), textcoords="offset points",
                        xytext=(6, 4), fontsize=7.5, color="0.35")
    r = p = np.nan
    if fit:
        r, p = _fit_line(ax, x, y, color)
    ax.set_xlabel(xlabel, fontsize=9.5)
    ax.set_ylabel(ylabel, fontsize=9.5)
    if title:
        ax.set_title(title, fontsize=9.5)
    return r, p


def main() -> int:
    t = april_eratio()
    years = t["year"].tolist()
    t = t.merge(ncp_budget_terms(years), on="year")
    samp = float_sampling(years)
    if len(samp.columns) > 1:
        t = t.merge(samp, on="year", how="left")
    t = t.sort_values("eratio_apr", ascending=False).reset_index(drop=True)

    mean_e = t.eratio_apr.mean()
    sd_e = t.eratio_apr.std(ddof=1)
    t["eratio_z"] = (t.eratio_apr - mean_e) / sd_e
    t.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}\n")

    print(f"=== April (DOY {APR[0]}-{APR[1]}) per year — "
          f"mean e-ratio {mean_e:.2f}, 1 SD {sd_e:.2f} ===")
    print(t.round(2).to_string(index=False))

    upper = t[t.eratio_z >= 1.0]
    print(f"\nAt or above the +1 SD bound (e-ratio >= {mean_e + sd_e:.2f}): "
          f"{', '.join(str(y) for y in upper.year) or 'none'}")

    print("\n=== correlation with the April e-ratio ===")
    ct = corr_table(t)
    print(ct.round(3).to_string(index=False))

    print("\n=== the two terms respond to April mixing in opposite directions ===")
    for c in ("ncp_apr", "npp_apr"):
        for z in ("zdepth_apr", "mld_apr_raw"):
            if z not in t:
                continue
            ok = np.isfinite(t[z]) & np.isfinite(t[c])
            r, p = pearsonr(t[z][ok], t[c][ok])
            print(f"  {z:>12} vs {c:>8}: r = {r:6.2f}  p = {p:.3f}")

    # ---------------------------------------------------------------- figure
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.0), dpi=300)
    fig.subplots_adjust(hspace=0.42, wspace=0.30)

    # (a) ranked April e-ratio with the mean / +-1 SD band
    ax = axes[0, 0]
    order = t.sort_values("eratio_apr", ascending=False)
    x = np.arange(len(order))
    cols = [ERATIO_COLOR if z >= 1 else "#c69ab0" for z in order.eratio_z]
    ax.axhspan(mean_e - sd_e, mean_e + sd_e, color=ERATIO_COLOR, alpha=0.10, lw=0)
    ax.axhline(mean_e, color=ERATIO_COLOR, lw=1.0, ls="--")
    ax.axhline(mean_e + sd_e, color=ERATIO_COLOR, lw=0.8, ls=":")
    ax.axhline(0, color=ZERO_LINE, lw=0.8, ls=":")
    ax.bar(x, order.eratio_apr, width=0.7, color=cols, edgecolor="white", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(order.year.astype(str), rotation=45, fontsize=8.5)
    ax.set_ylabel("April e-ratio (NCP : NPP)", fontsize=9.5)
    ax.set_title(f"April e-ratio by year  (mean {mean_e:.2f} ± {sd_e:.2f})",
                 fontsize=9.5)
    ax.text(0.98, 0.95, "dark = at/above +1 SD", transform=ax.transAxes,
            ha="right", va="top", fontsize=7.5, color="0.45")

    # (b) which term carries the variance. Both terms are standardised so they
    # share one axis — a twin axis here made the two point clouds impossible to
    # attribute, and the interesting quantity is the *relative* spread anyway.
    ax = axes[0, 1]
    rc, pc = _scatter(ax, _z(t.ncp_apr), t.eratio_apr, t.year, BLUE,
                      "April anomaly (SD units)", "April e-ratio", "")
    rn, pn = _scatter(ax, _z(t.npp_apr), t.eratio_apr, t.year, GREEN,
                      "April anomaly (SD units)", "April e-ratio", "",
                      marker="s", annotate=False)
    ax.set_title("The e-ratio follows NCP, not NPP", fontsize=9.5)
    ax.legend(handles=[
        mpl.lines.Line2D([], [], ls="none", marker="o", color=BLUE,
                         label=f"NCP  r = {rc:.2f}, p = {pc:.3f}   "
                               f"(CV {t.ncp_apr.std(ddof=1)/t.ncp_apr.mean():.2f})"),
        mpl.lines.Line2D([], [], ls="none", marker="s", color=GREEN,
                         label=f"NPP  r = {rn:.2f}, p = {pn:.3f}   "
                               f"(CV {t.npp_apr.std(ddof=1)/t.npp_apr.mean():.2f})"),
    ], loc="lower right", fontsize=8)

    # (c) both terms against the April integration depth, same standardisation
    ax = axes[1, 0]
    zd = t.zdepth_apr.to_numpy(float)
    rn2, pn2 = _scatter(ax, zd, _z(t.npp_apr), t.year, GREEN,
                        "April NCP integration depth, max(MLD, z$_{eu}$)  (m)",
                        "April anomaly (SD units)", "", marker="s")
    rc2, pc2 = _scatter(ax, zd, _z(t.ncp_apr), t.year, BLUE,
                        "April NCP integration depth, max(MLD, z$_{eu}$)  (m)",
                        "April anomaly (SD units)", "", annotate=False)
    ax.axhline(0, color=ZERO_LINE, lw=0.8, ls=":", zorder=1)
    ax.set_title("Deep April mixing: NPP down, NCP up", fontsize=9.5)
    ax.legend(handles=[
        mpl.lines.Line2D([], [], ls="none", marker="s", color=GREEN,
                         label=f"NPP  r = {rn2:.2f}, p = {pn2:.3f}"),
        mpl.lines.Line2D([], [], ls="none", marker="o", color=BLUE,
                         label=f"NCP  r = {rc2:.2f}, p = {pc2:.3f}"),
    ], loc="upper center", fontsize=8)

    # (d) seasonal phase: how far into the drawdown April falls
    ax = axes[1, 1]
    _scatter(ax, (t.frac_drawdown_by_apr * 100).to_numpy(),
             t.eratio_apr.to_numpy(), t.year, MLD_COLOR,
             "Annual nitrate drawdown completed by mid-April (%)",
             "April e-ratio",
             "High e-ratio = April early in the drawdown")
    rho, prho = spearmanr(t.frac_drawdown_by_apr, t.eratio_apr)
    ax.text(0.98, 0.92, f"Spearman ρ = {rho:.2f}, p = {prho:.3f}",
            transform=ax.transAxes, ha="right", fontsize=8, color="0.35")

    fig.suptitle(
        "What drives the April e-ratio peak — Iceland Basin & Irminger Sea, "
        f"{years[0]}–{years[-1]}",
        fontsize=11.5, y=0.975,
    )
    fig.text(
        0.5, 0.02,
        "April = DOY 92–120. Per-year e-ratio, NPP and NCP are the same "
        "per-year curves averaged in Fig 5; integration depth and nitrate "
        "inventory come from the NCP budget (ncp_results.xlsx).\n"
        "April float sampling (n profiles, n floats, mean position, Irminger "
        "fraction) shows no relation with the e-ratio (all |r| < 0.5, "
        "p > 0.2) and is carried in eratio_april_drivers.csv.",
        ha="center", va="top", fontsize=7.5, color="0.45",
    )
    fig.savefig(OUT_PNG, bbox_inches="tight", dpi=300)
    print(f"\nWrote {OUT_PNG}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
