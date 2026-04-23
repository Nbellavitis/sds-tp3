"""Shared matplotlib style helpers for report-ready plots."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


# ARULA palette (ordered for maximum visual separation in line plots).
ARULA_COLORS = [
    "#0B1F3A",  # navy
    "#00A6A6",  # teal
    "#F26419",  # orange
    "#7B2CBF",  # violet
    "#2E8B57",  # green
    "#D7263D",  # red
    "#1D4ED8",  # blue
    "#F4B400",  # amber
    "#6D28D9",  # indigo
    "#0EA5E9",  # cyan
]


def apply_plot_style():
    """Apply a larger, readable default style for all static figures."""
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.titlesize": 18,
            "axes.labelsize": 17,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 13,
            "figure.titlesize": 18,
        }
    )


def get_distinct_series_styles(n):
    """Return clearly distinct (color, marker, linestyle) triplets for n series."""
    if n <= 0:
        return []

    color_count = len(ARULA_COLORS)
    markers = ["o", "s", "D", "^", "v", "P", "X", "*", "<", ">", "h", "8"]
    linestyles = ["-", "--", "-.", ":"]

    styles = []
    for idx in range(n):
        styles.append(
            {
                "color": ARULA_COLORS[idx % color_count],
                "marker": markers[idx % len(markers)],
                "linestyle": linestyles[(idx // color_count) % len(linestyles)],
            }
        )
    return styles


def _decimals_from_tick_step(step):
    """Infer decimals from a numeric step size using plain decimal notation."""
    if not np.isfinite(step) or step <= 0:
        return 2
    text = f"{float(step):.12f}".rstrip("0").rstrip(".")
    if "." not in text:
        return 0
    return min(len(text.split(".")[1]), 8)


def compute_decimals_from_uncertainty(yerr):
    """Choose decimals from the smallest positive uncertainty (1 significant digit rule)."""
    if yerr is None:
        return None

    arr = np.asarray(yerr, dtype=float)
    valid = arr[np.isfinite(arr) & (arr > 0.0)]
    if valid.size == 0:
        return None

    err_ref = float(np.min(valid))
    exponent = int(np.floor(np.log10(err_ref)))
    decimals = max(0, -exponent)
    return min(decimals, 8)


def _plain_number_formatter(decimals):
    """Return a formatter that never uses scientific notation."""

    def _format_tick(value, _pos):
        if not np.isfinite(value):
            return ""
        formatted = f"{value:.{decimals}f}"
        if formatted.startswith("-0"):
            as_float = float(formatted)
            if abs(as_float) < 10 ** (-(decimals + 1)):
                formatted = formatted[1:]
        return formatted

    return FuncFormatter(_format_tick)


def format_y_axis(ax, yerr=None, decimals=None):
    """Format Y ticks in plain decimal notation with precision tied to uncertainty."""
    if decimals is None:
        decimals = compute_decimals_from_uncertainty(yerr)

    if decimals is None:
        ticks = np.asarray(ax.get_yticks(), dtype=float)
        finite_ticks = ticks[np.isfinite(ticks)]
        if ax.get_yscale() == "log":
            positive_ticks = finite_ticks[finite_ticks > 0.0]
            if positive_ticks.size:
                min_tick = float(np.min(positive_ticks))
                decimals = max(0, int(np.ceil(-np.log10(min_tick)))) if min_tick < 1.0 else 0
            else:
                decimals = 2
        else:
            if finite_ticks.size >= 2:
                diffs = np.diff(np.unique(np.sort(finite_ticks)))
                diffs = diffs[diffs > 0.0]
                step = float(np.min(diffs)) if diffs.size else np.nan
                decimals = _decimals_from_tick_step(step)
            else:
                decimals = 2

    ax.yaxis.set_major_formatter(_plain_number_formatter(int(decimals)))
    ax.yaxis.offsetText.set_visible(False)


