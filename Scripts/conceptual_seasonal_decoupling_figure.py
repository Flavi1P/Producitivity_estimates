"""
Conceptual figure - seasonal decoupling of NPP and NCP controls the apparent
export efficiency (NCP:NPP) of the subpolar North Atlantic.

This is a SCHEMATIC, not a results figure: no error bars, no uncertainty, no
quantitative claim beyond the ratios quoted in the manuscript. The curve
*shapes* are nonetheless taken straight from the paper's own Iceland Basin /
Irminger climatology so the schematic cannot contradict the results figures:

    Output/IcelandIrminger_2015_2025/eratio_doy_stats.csv
        npp_clim_mmolC   CbPM seasonal climatology, maximum at DOY 169 (mid-June)
        ncp_clim_mmolC   nitrate-drawdown climatology, maximum at DOY 127
                         (early May) and negative in winter and autumn
        eratio_mean      year-by-year NCP:NPP averaged over 2015-2023, maximum
                         0.42 at DOY 99 (~9 April)

NPP and NCP are normalised by the NPP maximum, so the shared production axis is
dimensionless and carries no tick labels - only the *relative timing and
amplitude* of the two curves is asserted. The e-ratio panel keeps real numbers
because those are the manuscript's headline values.

The ratio is only drawn where NPP exceeds ERATIO_FLOOR_FRAC of its seasonal
maximum. Below that floor NPP -> 0 while NCP stays negative and the point-wise
ratio is not estimable - the same caveat the eratio_synthetic_year results
figure carries.

Layout
    row 0   seasonal-regime header band (3 ecological regimes)
    row 1   NPP and NCP on a shared, unlabelled production axis
    row 2   NCP:NPP (apparent export efficiency)
    row 3   water-column schematic per regime
    row 4   the two process chains (export-favourable vs recycling-dominated)

Colours are the repo's semantic palette (Scripts/figure_config.py), so the
schematic reads as part of the same figure set: NPP green, NCP blue, e-ratio
plum, nutrients teal, grazing/recycling orange, sinking particles dark neutral.

Run from the repo root:
    & "C:/Users/petit/anaconda3/envs/cmts_learn_olci/python.exe" Scripts/conceptual_seasonal_decoupling_figure.py

Output: Output/conceptual_seasonal_decoupling.png
        Output/conceptual_seasonal_decoupling.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Ellipse
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import (                      # noqa: E402
    NPP_COLOR, NCP_COLOR, ERATIO_COLOR, ZERO_LINE, MLD_COLOR,
    P1_SHADE, P2_SHADE,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "Output"
CLIM_CSV = OUT_DIR / "IcelandIrminger_2015_2025" / "eratio_doy_stats.csv"
OUT_PNG = OUT_DIR / "conceptual_seasonal_decoupling.png"
OUT_PDF = OUT_DIR / "conceptual_seasonal_decoupling.pdf"

# ---------------------------------------------------------------------------
# Palette (semantic; extends figure_config with the schematic-only elements)
# ---------------------------------------------------------------------------
NUTRIENT = "#2e7d6e"    # muted blue-green - nitrate
PHYTO = NPP_COLOR       # green - phytoplankton
GRAZER = "#c2551f"      # orange - zooplankton, grazing, respiration
PARTICLE = "#3f3f3f"    # dark neutral - sinking particles
WATER = "#eef4f8"       # water column below the mixed layer
MIXED = "#d6e6f1"       # mixed layer
INK = "#222222"
INK_SOFT = "#5b5b5b"
INK_MUTED = "#8a8a8a"
HALO = dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="none",
            alpha=0.78)

# ---------------------------------------------------------------------------
# Seasonal x-axis:  x = (doy - 1) * 12 / 365,  i.e. x = 0 is 1 January and
# x = 3.0 is 1 April. The figure spans January to the end of September.
# ---------------------------------------------------------------------------
X0, X1 = 0.30, 8.95
X = np.linspace(X0, X1, 1400)

MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"]

# Ecological regimes: (start, end, name, key characterisation, physics, shade,
# accent colour for the characterisation line)
REGIMES = [
    (X0, 2.60, "WINTER", "Net heterotrophy",
     "Deep mixing  -  nutrient accumulation", "#dfe6ec", "#4a5a66"),
    (2.60, 4.60, "EARLY SPRING", "Efficient export window",
     "Restratification  -  grazing lag", P1_SHADE, ERATIO_COLOR),
    (4.60, X1, "LATE SPRING & SUMMER", "Recycling-dominated",
     "Stratified  -  grazing and respiration increase", P2_SHADE, "#7a4a1f"),
]

# The point-wise ratio is only drawn where NPP exceeds this fraction of its
# seasonal maximum; below it NPP -> 0 and NCP:NPP is not estimable.
ERATIO_FLOOR_FRAC = 0.20
ERATIO_LEADIN = 0.30      # length of the dotted lead-in below the floor, in x


def doy_to_x(doy):
    return (np.asarray(doy, dtype=float) - 1.0) * 12.0 / 365.0


# ---------------------------------------------------------------------------
# Curves, read from the paper's own climatology
# ---------------------------------------------------------------------------
def load_curves():
    """Return smooth NPP / NCP / e-ratio callables on the seasonal x-axis.

    NPP and NCP are normalised by the NPP maximum. The e-ratio is the paper's
    year-by-year mean ratio, not the ratio of the two climatologies.
    """
    if not CLIM_CSV.exists():
        raise SystemExit(
            "missing %s - run Scripts/eratio_synthetic_year_figure.py first"
            % CLIM_CSV)

    d = pd.read_csv(CLIM_CSV)
    x = doy_to_x(d["doy"].to_numpy())
    npp = d["npp_clim_mmolC"].to_numpy(dtype=float)
    ncp = d["ncp_clim_mmolC"].to_numpy(dtype=float)
    npp_max = float(np.nanmax(npp))

    f_npp = PchipInterpolator(x, npp / npp_max)
    f_ncp = PchipInterpolator(x, ncp / npp_max)

    e = d["eratio_mean"].to_numpy(dtype=float)
    ok = np.isfinite(e)
    f_er = PchipInterpolator(x[ok], e[ok])

    # x below which the ratio stops being estimable (NPP below the floor)
    x_est = float(x[np.argmax(npp >= ERATIO_FLOOR_FRAC * npp_max)])
    return f_npp, f_ncp, f_er, x_est, npp_max


_f_npp, _f_ncp, _f_er, ERATIO_ESTIMABLE_FROM, NPP_MAX_MMOL = load_curves()


def npp_curve(x):
    return np.asarray(_f_npp(x), dtype=float)


def ncp_curve(x):
    return np.asarray(_f_ncp(x), dtype=float)


def eratio_curve(x):
    return np.asarray(_f_er(x), dtype=float)


def _argmax_x(f, lo, hi):
    xx = np.linspace(lo, hi, 8000)
    return float(xx[np.argmax(f(xx))])


NPP_PEAK_X = _argmax_x(npp_curve, 1.0, 8.5)
NCP_PEAK_X = _argmax_x(ncp_curve, 1.0, 8.5)
ERATIO_PEAK_X = _argmax_x(eratio_curve, 2.3, 6.0)
OFFSET_DAYS = (NPP_PEAK_X - NCP_PEAK_X) * 365.0 / 12.0


# ---------------------------------------------------------------------------
# Mixed-layer depth (conceptual, dimensionless: 0 = surface, -1 = "deep")
# ---------------------------------------------------------------------------
def mld(x):
    x = np.asarray(x, dtype=float)
    return -(0.09
             + 0.80 / (1.0 + np.exp((x - 2.95) / 0.62))
             + 0.34 / (1.0 + np.exp(-(x - 8.30) / 0.50)))


# ---------------------------------------------------------------------------
# Small drawing helpers
# ---------------------------------------------------------------------------
def season_shading(ax, alpha=0.30):
    for x0, x1, *_rest, shade, _accent in REGIMES:
        ax.axvspan(x0, x1, color=shade, alpha=alpha, lw=0, zorder=0)
    for _x0, x1, *_rest in REGIMES[:-1]:
        ax.axvline(x1, color="white", lw=1.6, zorder=1)
        ax.axvline(x1, color=INK_MUTED, lw=0.6, ls=(0, (4, 3)), zorder=1.1)


def month_axis(ax, labels=True):
    """Month ticks. Note the axes share x, so the tick *formatter* is shared -
    label visibility must be toggled per axis, not by blanking the labels."""
    ax.set_xlim(X0, X1)
    ax.set_xticks(np.arange(0.5, 9.0, 1.0))
    ax.set_xticklabels(MONTH_LABELS, fontsize=9.5, color=INK_SOFT)
    ax.tick_params(axis="x", length=3, color=INK_MUTED, pad=3,
                   labelbottom=labels)


def strip_spines(ax, keep=("left",)):
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)
        if side in keep:
            ax.spines[side].set_color(INK_MUTED)
            ax.spines[side].set_linewidth(0.8)


def curved_arrow(ax, xy_from, xy_to, color, rad=0.35, lw=1.3, ls="-",
                 alpha=1.0, zorder=6, mutation=8):
    ax.add_patch(FancyArrowPatch(
        xy_from, xy_to, connectionstyle="arc3,rad=%s" % rad,
        arrowstyle="-|>", mutation_scale=mutation, lw=lw, ls=ls,
        color=color, alpha=alpha, zorder=zorder, shrinkA=1, shrinkB=1))


def assert_within(ax, values, name):
    """Guard against silent clipping (same idea as figure_config's helper)."""
    lo, hi = ax.get_ylim()
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size and (v.min() < lo or v.max() > hi):
        raise SystemExit(
            "%s spans [%.3f, %.3f] outside the axis [%.3f, %.3f]"
            % (name, v.min(), v.max(), lo, hi))


# ===========================================================================
# Figure
# ===========================================================================
def build():
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.linewidth": 0.8,
        "axes.edgecolor": INK_MUTED,
        "text.color": INK,
        "axes.labelcolor": INK,
        "ytick.color": INK_SOFT,
        "xtick.color": INK_SOFT,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
    })

    fig = plt.figure(figsize=(12.5, 10.8), dpi=300)
    gs = GridSpec(
        5, 1, figure=fig,
        height_ratios=[0.85, 2.45, 1.55, 3.05, 0.70],
        hspace=0.13, left=0.075, right=0.982, top=0.925, bottom=0.082,
    )

    ax_band = fig.add_subplot(gs[0])
    ax_prod = fig.add_subplot(gs[1])
    ax_er = fig.add_subplot(gs[2], sharex=ax_prod)
    ax_col = fig.add_subplot(gs[3], sharex=ax_prod)
    ax_txt = fig.add_subplot(gs[4])

    # -- title ------------------------------------------------------------
    fig.text(0.075, 0.978,
             "Seasonal decoupling of primary and net community production "
             "controls apparent export efficiency",
             fontsize=13.5, fontweight="bold", ha="left", va="center",
             color=INK)
    fig.text(0.075, 0.955,
             "Subpolar North Atlantic (Iceland Basin / Irminger Sea)  -  "
             "conceptual schematic",
             fontsize=10, ha="left", va="center", color=INK_SOFT)

    # ===================================================================
    # Row 0 - regime header band
    # ===================================================================
    ax_band.set_xlim(X0, X1)
    ax_band.set_ylim(0, 1)
    ax_band.axis("off")
    for x0, x1, name, key, phys, shade, accent in REGIMES:
        ax_band.add_patch(FancyBboxPatch(
            (x0 + 0.04, 0.06), (x1 - x0) - 0.08, 0.88,
            boxstyle="round,pad=0,rounding_size=0.10",
            facecolor=shade, edgecolor="none", alpha=0.85, zorder=1))
        xc = 0.5 * (x0 + x1)
        ax_band.text(xc, 0.76, name, ha="center", va="center", fontsize=10.2,
                     fontweight="bold", color="#2b2b2b", zorder=2)
        ax_band.text(xc, 0.48, key, ha="center", va="center", fontsize=8.6,
                     fontweight="bold", color=accent, zorder=2)
        ax_band.text(xc, 0.20, phys, ha="center", va="center", fontsize=7.6,
                     color="#55606a", zorder=2)

    # ===================================================================
    # Row 1 - NPP and NCP (the temporal offset is the headline)
    # ===================================================================
    season_shading(ax_prod)
    strip_spines(ax_prod)
    month_axis(ax_prod, labels=False)

    npp = npp_curve(X)
    ncp = ncp_curve(X)

    ax_prod.set_ylim(-0.30, 1.36)
    ax_prod.axhline(0, color=ZERO_LINE, lw=0.9, ls=(0, (2, 3)), zorder=2)
    ax_prod.fill_between(X, 0, npp, color=NPP_COLOR, alpha=0.07, lw=0,
                         zorder=2.5)
    ax_prod.fill_between(X, 0, ncp, color=NCP_COLOR, alpha=0.11, lw=0,
                         zorder=2.6)
    ax_prod.plot(X, npp, color=NPP_COLOR, lw=2.7, zorder=4,
                 solid_capstyle="round")
    ax_prod.plot(X, ncp, color=NCP_COLOR, lw=2.7, zorder=4.1,
                 solid_capstyle="round")
    assert_within(ax_prod, np.concatenate([npp, ncp]), "NPP/NCP")

    ax_prod.set_yticks([])
    ax_prod.set_ylabel("Production\n(conceptual)", fontsize=10.5, color=INK,
                       labelpad=10)
    ax_prod.text(X0 + 0.09, 0.02, "0", fontsize=8.5, color=INK_SOFT,
                 ha="left", va="bottom", zorder=5)

    # direct labels - identity is never colour-alone. Both are right-aligned
    # on the same edge, each in the empty band just above its own curve.
    ax_prod.text(8.90, 0.83, "Net primary\nproduction (NPP)", color=NPP_COLOR,
                 fontsize=10.5, fontweight="bold", ha="right", va="center",
                 linespacing=1.25, zorder=6)
    ax_prod.text(8.90, 0.30, "Net community\nproduction (NCP)",
                 color=NCP_COLOR, fontsize=10.5, fontweight="bold",
                 ha="right", va="center", linespacing=1.25, zorder=6)
    ax_prod.plot([8.45, 8.45], [float(ncp_curve(8.45)) + 0.015, 0.21],
                 color=NCP_COLOR, lw=0.9, alpha=0.75, zorder=5)
    ax_prod.plot([8.45, 8.45], [float(npp_curve(8.45)) + 0.015, 0.74],
                 color=NPP_COLOR, lw=0.9, alpha=0.75, zorder=5)

    npp_pk = float(npp_curve(NPP_PEAK_X))
    ncp_pk = float(ncp_curve(NCP_PEAK_X))
    for xp, yp, col in ((NCP_PEAK_X, ncp_pk, NCP_COLOR),
                        (NPP_PEAK_X, npp_pk, NPP_COLOR)):
        ax_prod.plot([xp, xp], [-0.30, yp], color=col, lw=0.9,
                     ls=(0, (3, 3)), alpha=0.75, zorder=3)
        ax_prod.plot([xp], [yp], marker="o", ms=7.5, mfc="white", mec=col,
                     mew=2.1, zorder=7)

    off_y = 1.20
    ax_prod.annotate("", xy=(NPP_PEAK_X, off_y), xytext=(NCP_PEAK_X, off_y),
                     arrowprops=dict(arrowstyle="<|-|>", lw=1.5,
                                     color="#3a3a3a", mutation_scale=11,
                                     shrinkA=0, shrinkB=0), zorder=7)
    ax_prod.text(0.5 * (NCP_PEAK_X + NPP_PEAK_X), off_y + 0.05,
                 "NCP peaks ~6 weeks before NPP", ha="center", va="bottom",
                 fontsize=11.5, fontweight="bold", color="#2b2b2b", zorder=7)

    ax_prod.annotate("NCP maximum", xy=(NCP_PEAK_X, ncp_pk),
                     xytext=(NCP_PEAK_X - 0.20, 0.52),
                     fontsize=9.0, color=NCP_COLOR, ha="right", va="center",
                     arrowprops=dict(arrowstyle="-", lw=0.9, color=NCP_COLOR,
                                     shrinkA=2, shrinkB=6), zorder=7)
    ax_prod.text(NPP_PEAK_X + 0.13, 1.05, "NPP maximum", fontsize=9.0,
                 color=NPP_COLOR, ha="left", va="center", zorder=7)

    ax_prod.text(1.45, -0.21, "Net heterotrophy", fontsize=8.8,
                 style="italic", color=INK_SOFT, ha="center", va="center",
                 zorder=6)
    ax_prod.text(7.55, -0.21, "NPP still high, NCP below zero", fontsize=8.8,
                 style="italic", color=INK_SOFT, ha="center", va="center",
                 zorder=6)

    # ===================================================================
    # Row 2 - apparent export efficiency
    # ===================================================================
    season_shading(ax_er)
    strip_spines(ax_er)
    month_axis(ax_er, labels=True)

    er = eratio_curve(X)
    est = X >= ERATIO_ESTIMABLE_FROM
    lead = (X >= ERATIO_ESTIMABLE_FROM - ERATIO_LEADIN) & ~est

    ax_er.set_ylim(-0.64, 0.62)
    ax_er.axhline(0, color=ZERO_LINE, lw=0.9, ls=(0, (2, 3)), zorder=2)
    ax_er.fill_between(X[est], 0, er[est], where=er[est] > 0,
                       color=ERATIO_COLOR, alpha=0.14, lw=0, zorder=2.5)
    ax_er.fill_between(X[est], 0, er[est], where=er[est] < 0,
                       color=ERATIO_COLOR, alpha=0.07, lw=0, zorder=2.5)
    ax_er.plot(X[est], er[est], color=ERATIO_COLOR, lw=3.0, zorder=4,
               solid_capstyle="round")
    ax_er.plot(X[lead], er[lead], color=ERATIO_COLOR, lw=1.6,
               ls=(0, (1.5, 2.0)), alpha=0.5, zorder=4)
    assert_within(ax_er, np.concatenate([er[est], er[lead]]), "e-ratio")

    ax_er.set_yticks([-0.4, -0.2, 0.0, 0.2, 0.4])
    ax_er.set_yticklabels(["-0.4", "-0.2", "0", "0.2", "0.4"], fontsize=9)
    ax_er.set_ylabel("Apparent export\nefficiency (NCP:NPP)", fontsize=10.5,
                     color=ERATIO_COLOR, labelpad=10)
    ax_er.text(1.42, -0.17, "not estimable in winter\n(NPP $\\rightarrow$ 0)",
               fontsize=8.4, style="italic", color=INK_MUTED, ha="center",
               va="center", linespacing=1.3, zorder=5)

    er_pk = float(eratio_curve(ERATIO_PEAK_X))
    ax_er.plot([ERATIO_PEAK_X], [er_pk], marker="o", ms=8.5, mfc="white",
               mec=ERATIO_COLOR, mew=2.3, zorder=7)
    ax_er.annotate("Maximum apparent\nexport efficiency  ~0.40",
                   xy=(ERATIO_PEAK_X, er_pk),
                   xytext=(ERATIO_PEAK_X - 0.92, 0.58),
                   fontsize=10.5, fontweight="bold", color=ERATIO_COLOR,
                   ha="right", va="top", linespacing=1.3,
                   arrowprops=dict(arrowstyle="-|>", lw=1.6,
                                   color=ERATIO_COLOR, mutation_scale=11,
                                   shrinkA=4, shrinkB=5), zorder=7)
    ax_er.text(ERATIO_PEAK_X - 0.92, 0.31, "late March - mid-April",
               fontsize=8.8, color=ERATIO_COLOR, ha="right", va="top",
               style="italic", zorder=7)

    for xm, lab in ((4.5, "0.33"), (5.5, "0.21"), (6.5, "0.15")):
        ym = float(eratio_curve(xm))
        ax_er.plot([xm], [ym], marker="o", ms=4.5, color=ERATIO_COLOR,
                   zorder=6)
        ax_er.text(xm, ym - 0.07, lab, fontsize=8.8, color=ERATIO_COLOR,
                   ha="center", va="top", zorder=6)
    # NB: arc3 `rad` is a fraction of the chord in *display* space, so on a
    # wide, short panel even 0.15 bows the arrow far off the curve.
    curved_arrow(ax_er, (4.75, 0.36), (7.35, 0.04), "#5f5f5f", rad=-0.06,
                 lw=1.3, zorder=6, mutation=10)
    ax_er.text(6.05, 0.38, "progressive decline", fontsize=9.5,
               color="#4f4f4f", ha="center", va="bottom", style="italic",
               zorder=6)
    ax_er.text(7.70, -0.40, "becomes negative\nfrom August", fontsize=8.8,
               color=ERATIO_COLOR, ha="center", va="center", style="italic",
               linespacing=1.3, zorder=6, bbox=HALO)

    ax_er.text(8.86, 0.545, u"Season-integrated ratio  0.11 \u00b1 0.06",
               fontsize=8.8, color="#3f3f3f", ha="right", va="center",
               zorder=8,
               bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                         edgecolor="#c9c9c9", linewidth=0.7, alpha=0.95))

    # ===================================================================
    # Row 3 - water-column schematic
    # ===================================================================
    ax_col.set_ylim(-1.02, 0.14)
    ax_col.axis("off")
    month_axis(ax_col, labels=False)

    rng = np.random.default_rng(11)
    xs_fine = np.linspace(X0, X1, 900)
    mld_fine = mld(xs_fine)

    ax_col.fill_between(xs_fine, -1.02, 0.0, color=WATER, lw=0, zorder=0)
    ax_col.fill_between(xs_fine, mld_fine, 0.0, color=MIXED, lw=0, zorder=1)
    ax_col.plot(xs_fine, mld_fine, color=MLD_COLOR, lw=2.0, zorder=5,
                solid_capstyle="round")
    ax_col.axhline(0.0, color="#6f8fa5", lw=1.5, zorder=6)
    for _x0, x1, *_rest in REGIMES[:-1]:
        ax_col.plot([x1, x1], [-1.02, 0.0], color="white", lw=1.6, zorder=4)
        ax_col.plot([x1, x1], [-1.02, 0.0], color=INK_MUTED, lw=0.6,
                    ls=(0, (4, 3)), zorder=4.1)

    ax_col.annotate("Mixed layer depth", xy=(8.05, float(mld(8.05))),
                    xytext=(8.05, -0.47), fontsize=8.8, color=MLD_COLOR,
                    ha="center", va="center", zorder=10, bbox=HALO,
                    arrowprops=dict(arrowstyle="-", lw=0.9, color=MLD_COLOR,
                                    shrinkA=3, shrinkB=2))
    ax_col.annotate("", xy=(X0 + 0.13, -0.90), xytext=(X0 + 0.13, -0.14),
                    arrowprops=dict(arrowstyle="-|>", lw=1.0, color=INK_MUTED,
                                    mutation_scale=9), zorder=7)
    ax_col.text(X0 + 0.23, -0.52, "depth", fontsize=8.0, color=INK_MUTED,
                rotation=90, ha="left", va="center", zorder=7)

    def scatter_nutrients(x0, x1, n, ytop_frac, ybot=-0.98, size=15,
                          alpha=0.55):
        """Nitrate crosses between ytop_frac of the local MLD and ybot."""
        xs = rng.uniform(x0, x1, n)
        top = mld(xs) * ytop_frac
        ys = top + (ybot - top) * rng.random(n)
        ax_col.scatter(xs, ys, marker="+", s=size, c=NUTRIENT, alpha=alpha,
                       linewidths=0.9, zorder=3)

    def scatter_phyto(x0, x1, n, size, alpha=0.85, jitter_top=0.015):
        xs = rng.uniform(x0, x1, n)
        m = mld(xs)
        ys = -jitter_top + m * (0.05 + 0.85 * rng.random(n))
        ax_col.scatter(xs, ys, s=size, c=PHYTO, alpha=alpha,
                       edgecolors="white", linewidths=0.5, zorder=8)

    def draw_grazers(pos, w=0.070, h=0.036):
        for (gx, gy) in pos:
            ax_col.add_patch(Ellipse((gx, gy), w, h, angle=-18,
                                     facecolor=GRAZER, edgecolor="white",
                                     lw=0.6, alpha=0.95, zorder=9))

    # ---- winter --------------------------------------------------------
    scatter_nutrients(0.42, 2.50, 130, ytop_frac=0.0, size=17, alpha=0.62)
    scatter_phyto(0.65, 2.40, 9, size=11, alpha=0.55)
    for gx in (1.00, 1.95):
        curved_arrow(ax_col, (gx, -0.80), (gx, -0.28), NUTRIENT, rad=0.0,
                     lw=1.3, alpha=0.9, zorder=7, mutation=9)
    ax_col.text(1.47, -0.15, "NO$_3^-$ replenished\nby deep mixing",
                fontsize=9.0, color=NUTRIENT, fontweight="bold",
                ha="center", va="center", linespacing=1.3, zorder=10,
                bbox=HALO)
    ax_col.text(1.47, -0.94, "Low biomass, low light", fontsize=8.4,
                color=INK_SOFT, ha="center", va="center", style="italic",
                zorder=10, bbox=HALO)

    # ---- early spring (the hero panel) ---------------------------------
    scatter_nutrients(2.70, 4.50, 46, ytop_frac=0.55, size=16, alpha=0.60)
    scatter_nutrients(2.70, 4.50, 26, ytop_frac=0.0, ybot=-0.36, size=13,
                      alpha=0.34)
    scatter_phyto(2.75, 4.45, 34, size=52, alpha=0.88)
    for gx in (3.02, 3.95):
        curved_arrow(ax_col, (gx, -0.64), (gx, -0.40), NUTRIENT, rad=0.0,
                     lw=1.2, alpha=0.85, zorder=7, mutation=8)
    draw_grazers([(4.32, -0.24)])
    ax_col.text(4.42, -0.24, "few\ngrazers", fontsize=8.0, color=GRAZER,
                ha="left", va="center", style="italic", linespacing=1.2,
                zorder=10)
    ax_col.text(3.50, -0.06, "Large, fast-growing cells", fontsize=9.0,
                color=PHYTO, fontweight="bold", ha="center", va="center",
                zorder=10, bbox=HALO)

    # sinking particles - a conceptual implication, not a measured flux
    for sx in (2.95, 3.45, 3.95):
        for k in range(3):
            ax_col.plot([sx + 0.02 * k], [-0.52 - 0.075 * k], marker="s",
                        ms=2.6, color=PARTICLE, alpha=0.75 - 0.15 * k,
                        zorder=8)
        curved_arrow(ax_col, (sx, -0.48), (sx + 0.05, -0.84), PARTICLE,
                     rad=0.0, lw=1.1, alpha=0.7, zorder=7, mutation=8)
    ax_col.text(3.50, -0.94,
                "carbon accumulates and sinks\n(conceptual implication)",
                fontsize=8.0, color=PARTICLE, ha="center", va="center",
                style="italic", linespacing=1.25, zorder=10, bbox=HALO)

    # ---- late spring / summer ------------------------------------------
    scatter_nutrients(4.72, 8.88, 120, ytop_frac=1.02, size=16, alpha=0.55)
    scatter_phyto(4.75, 8.55, 46, size=13, alpha=0.80, jitter_top=0.008)
    draw_grazers([(4.95, -0.055), (5.40, -0.105), (5.85, -0.05),
                  (6.30, -0.11), (6.75, -0.05), (7.20, -0.105),
                  (7.65, -0.06), (8.10, -0.115), (8.55, -0.06)])
    for cx in (5.15, 6.05, 6.95):
        curved_arrow(ax_col, (cx, -0.028), (cx + 0.42, -0.028), GRAZER,
                     rad=-0.85, lw=1.2, alpha=0.9, zorder=9, mutation=8)
        curved_arrow(ax_col, (cx + 0.42, -0.128), (cx, -0.128), GRAZER,
                     rad=-0.85, lw=1.2, alpha=0.9, zorder=9, mutation=8)
    ax_col.text(6.60, 0.085, "Grazing + microbial recycling + respiration",
                fontsize=9.2, color=GRAZER, fontweight="bold", ha="center",
                va="center", zorder=10)

    curved_arrow(ax_col, (5.80, -0.26), (5.87, -0.50), PARTICLE, rad=0.0,
                 lw=1.0, ls=(0, (2, 2)), alpha=0.6, zorder=7, mutation=7)
    ax_col.text(5.80, -0.62, "less carbon leaves\nthe surface layer",
                fontsize=8.0, color=PARTICLE, ha="center", va="center",
                style="italic", linespacing=1.25, zorder=10, bbox=HALO)
    ax_col.text(7.15, -0.94, "Shallow, strongly stratified mixed layer",
                fontsize=8.4, color=INK_SOFT, ha="center", va="center",
                style="italic", zorder=10, bbox=HALO)

    # ===================================================================
    # Row 4 - the two process chains
    # ===================================================================
    ax_txt.set_xlim(X0, X1)
    ax_txt.set_ylim(0, 1)
    ax_txt.axis("off")

    chains = [
        (2.62, 4.58, P1_SHADE, "#5c2038",
         u"Production  \u2192  carbon retained as NCP\n"
         u"\u2192  HIGH apparent export efficiency"),
        (4.62, X1 - 0.02, P2_SHADE, "#4a3566",
         u"Production  \u2192  grazing, recycling and respiration  "
         u"\u2192  lower NCP\n\u2192  LOWER apparent export efficiency"),
    ]
    for x0, x1, shade, txt_col, label in chains:
        ax_txt.add_patch(FancyBboxPatch(
            (x0, 0.12), x1 - x0, 0.74,
            boxstyle="round,pad=0,rounding_size=0.12",
            facecolor=shade, edgecolor="none", alpha=0.60, zorder=1))
        ax_txt.text(0.5 * (x0 + x1), 0.49, label, ha="center", va="center",
                    fontsize=8.6, color=txt_col, fontweight="bold",
                    linespacing=1.4, zorder=2)

    fig.text(0.075, 0.010,
             "Conceptual schematic. The three curves follow the 2015-2023 "
             "Iceland Basin / Irminger Sea climatology, but the timing of the "
             "ecological processes drawn below them is illustrative, not "
             "measured.\n"
             "NCP:NPP is an apparent upper-ocean export-efficiency proxy - not "
             "a measurement of sinking particle flux; DOC accumulation, "
             "lateral transport, entrainment / detrainment and time lags all "
             "separate the two.\n"
             "The declining ratio is consistent with increasing recycling and "
             "a grazing / respiration lag; this study does not measure grazing "
             "directly.",
             fontsize=7.8, color=INK_SOFT, ha="left", va="bottom",
             linespacing=1.5)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300, facecolor="white")
    fig.savefig(OUT_PDF, facecolor="white")
    plt.close(fig)

    print("NPP climatology maximum : %.1f mmol C m-2 d-1" % NPP_MAX_MMOL)
    print("e-ratio maximum : x=%.2f (DOY %.0f)  value=%.3f"
          % (ERATIO_PEAK_X, ERATIO_PEAK_X * 365 / 12 + 1,
             float(eratio_curve(ERATIO_PEAK_X))))
    print("NCP maximum     : x=%.2f (DOY %.0f)"
          % (NCP_PEAK_X, NCP_PEAK_X * 365 / 12 + 1))
    print("NPP maximum     : x=%.2f (DOY %.0f)"
          % (NPP_PEAK_X, NPP_PEAK_X * 365 / 12 + 1))
    print("offset NCP->NPP : %.0f days (%.1f weeks)"
          % (OFFSET_DAYS, OFFSET_DAYS / 7.0))
    print("ratio drawn solid from x=%.2f (DOY %.0f)"
          % (ERATIO_ESTIMABLE_FROM, ERATIO_ESTIMABLE_FROM * 365 / 12 + 1))
    print("wrote %s" % OUT_PNG)
    print("wrote %s" % OUT_PDF)


if __name__ == "__main__":
    build()
