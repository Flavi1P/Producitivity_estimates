"""
Shared configuration and helpers for the productivity figure scripts.

This is the single source of truth for values that were previously copy-pasted
across the figure scripts (`synthesis_gop_npp_ncp_figure.py`, `gop_vs_ncp_*`,
`ncp_timeseries_figure.py`, `npp_timeseries_figure.py`,
`eratio_synthetic_year_figure.py`, `gop_depth_sensitivity_figure.py`):

  * the 18-day centered-window smoothing parameters,
  * the common 2015-2023 comparison window,
  * the estimate colours / shading colours,
  * the GOP deep-mixing (winter) validity threshold,

plus three helpers that let every figure render identically:

  * ``centered_window_stats`` - the one shared smoothing implementation,
  * ``render_masked``         - grey band + line suppression over invalid spans,
  * ``assert_within_axis``    - guard against silent y-axis clipping.

Importing this module has NO side effects on matplotlib state (it does not touch
``rcParams``) so each figure keeps its own styling.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Smoothing (per mat_and_meth.md S2.2.1)
# ---------------------------------------------------------------------------
WINDOW_DAYS = 18        # centered smoothing window width, days
MIN_N       = 3         # minimum raw points required inside a smoothing window

# ---------------------------------------------------------------------------
# Common comparison window
# ---------------------------------------------------------------------------
# NPP (reprocessed CMEMS CbPM) stops 2023-12-27, so basin comparisons and the
# e-ratio are restricted to the common 2015-2023 years (see FIGURE_HANDOFF.md).
COMMON_YEAR_START = 2015
COMMON_YEAR_END   = 2023

# ---------------------------------------------------------------------------
# GOP winter (deep-mixing) validity threshold
# ---------------------------------------------------------------------------
# The float O2 budget is not interpretable once the integration bottom follows
# deep winter convection; mask GOP where z_bot_m exceeds this depth.
GOP_DEEP_MIX_THRESHOLD_M = 200.0

# ---------------------------------------------------------------------------
# Colours (formerly redefined per file)
# ---------------------------------------------------------------------------
# GOP - float O2 budget (warm rust)
GOP_COLOR = "#c1440e"
GOP_LIGHT = "#f0c4a8"

# NCP - nitrate budget (blue). BLUE / BLUE_LIGHT are the names used in
# ncp_timeseries_figure.py; keep both as aliases.
NCP_COLOR  = "#1f4e79"
BLUE       = NCP_COLOR
BLUE_LIGHT = "#b8cbe0"

# NPP - CbPM (green). GREEN / GREEN_LIGHT are the names used in
# npp_timeseries_figure.py; keep both as aliases.
NPP_COLOR   = "#1a6e3c"
GREEN       = NPP_COLOR
GREEN_LIGHT = "#c0deca"

# e-ratio (NCP/NPP) - plum
ERATIO_COLOR = "#7a2048"

# Reference / structural colours
ZERO_LINE = "#888888"    # zero / reference dotted lines
MLD_COLOR = "#5a5a5a"    # mixed-layer / integration-depth trace
MASK_BAND = "#d9d9d9"    # grey band drawn over invalid ("not estimated") spans

# Seasonal-regime shading (single-float synthesis figure)
P1_SHADE = "#f6d98a"     # warm gold  - spring onset
P2_SHADE = "#c9b8dd"     # muted violet - post-bloom drawdown


# ---------------------------------------------------------------------------
# Shared smoothing implementation
# ---------------------------------------------------------------------------

def centered_window_stats(times, values, grid, window_days=WINDOW_DAYS, min_n=MIN_N):
    """Centered time-window median / IQR / count of ``values`` at each grid time.

    Parameters
    ----------
    times : array-like of datetime64
        Observation times of ``values``.
    values : array-like of float
        Values to summarise.
    grid : array-like of datetime64
        Times at which to evaluate the smoothed statistics.
    window_days : float
        Full width of the centered window (half is taken on each side).
    min_n : int
        Minimum number of raw points inside the window for a non-NaN result.

    Returns
    -------
    med, p25, p75 : ndarray of float
        Median and 25th/75th percentiles per grid point (NaN where ``count`` <
        ``min_n``).
    count : ndarray of int
        Number of raw points inside the window at each grid point.
    """
    times = np.asarray(times)
    values = np.asarray(values, dtype=float)
    grid = np.asarray(grid)

    half = np.timedelta64(int(window_days / 2 * 86400), "s")
    n = len(grid)
    med = np.full(n, np.nan)
    p25 = np.full(n, np.nan)
    p75 = np.full(n, np.nan)
    count = np.zeros(n, dtype=int)
    for i, t in enumerate(grid):
        mask = (times >= t - half) & (times <= t + half)
        v = values[mask]
        count[i] = v.size
        if v.size >= min_n:
            med[i] = np.median(v)
            p25[i], p75[i] = np.percentile(v, [25, 75])
    return med, p25, p75, count


# ---------------------------------------------------------------------------
# Masked-span rendering
# ---------------------------------------------------------------------------

def _contiguous_runs(flags):
    """Yield (start_idx, end_idx) inclusive index runs where ``flags`` is True."""
    flags = np.asarray(flags, dtype=bool)
    if flags.size == 0:
        return
    # indices where a run starts / ends
    padded = np.concatenate(([False], flags, [False]))
    diff = np.diff(padded.astype(int))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0] - 1
    for s, e in zip(starts, ends):
        yield int(s), int(e)


def render_masked(ax, x, valid_mask, y=None, *, line_kwargs=None,
                  band_color=MASK_BAND, band_alpha=0.55, band_zorder=0.5,
                  band_label=None):
    """Draw a grey band over invalid spans and (optionally) the valid-only line.

    Every figure calls this so masked ("not estimated" / deep-mixing) periods
    look identical: a light-grey ``axvspan`` over each contiguous invalid run,
    with the data line suppressed there.

    Parameters
    ----------
    ax : matplotlib Axes
    x : array-like
        X coordinates (datetime64 or numeric) aligned with ``valid_mask``.
    valid_mask : array-like of bool
        True where data is valid; False spans get the grey band / no line.
    y : array-like, optional
        If given, the line is plotted with values NaN-masked wherever
        ``valid_mask`` is False, so no line is drawn across invalid spans.
    line_kwargs : dict, optional
        Passed to ``ax.plot`` for the valid-only line (only used when ``y`` is
        given).
    band_color, band_alpha, band_zorder : styling for the invalid-span band.
    band_label : str, optional
        Legend label applied to the first invalid band only.

    Returns
    -------
    y_masked : ndarray or None
        ``y`` with invalid entries set to NaN (or None if ``y`` was None).
    """
    x = np.asarray(x)
    valid = np.asarray(valid_mask, dtype=bool)
    invalid = ~valid

    labelled = False
    for s, e in _contiguous_runs(invalid):
        # extend the span to the neighbouring valid samples so the band abuts
        # the drawn line rather than leaving a sliver gap
        left = x[s - 1] if s > 0 else x[s]
        right = x[e + 1] if e < len(x) - 1 else x[e]
        ax.axvspan(left, right, color=band_color, alpha=band_alpha,
                   lw=0, zorder=band_zorder,
                   label=(band_label if (band_label and not labelled) else None))
        labelled = labelled or (band_label is not None)

    y_masked = None
    if y is not None:
        y_masked = np.asarray(y, dtype=float).copy()
        y_masked[invalid] = np.nan
        if line_kwargs is not None:
            ax.plot(x, y_masked, **line_kwargs)
    return y_masked


# ---------------------------------------------------------------------------
# Clipping guard
# ---------------------------------------------------------------------------

def assert_within_axis(ax, *arrays, tol=1e-9):
    """Raise if any finite value in ``arrays`` falls outside ``ax.get_ylim()``.

    Used to forbid silent clipping: figures should set y-limits that actually
    contain their data rather than ``.clip()``-ing values into range.
    """
    lo, hi = ax.get_ylim()
    if lo > hi:                     # inverted axis (e.g. depth strip)
        lo, hi = hi, lo
    span = hi - lo
    pad = tol * (span if span > 0 else 1.0)

    for arr in arrays:
        a = np.asarray(arr, dtype=float)
        finite = a[np.isfinite(a)]
        if finite.size == 0:
            continue
        amin = float(finite.min())
        amax = float(finite.max())
        if amin < lo - pad or amax > hi + pad:
            raise ValueError(
                f"assert_within_axis: plotted values [{amin:.4g}, {amax:.4g}] "
                f"exceed axis y-limits [{lo:.4g}, {hi:.4g}]"
            )
