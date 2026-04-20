"""Shared matplotlib style helpers for report-ready plots."""

import matplotlib.pyplot as plt


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


