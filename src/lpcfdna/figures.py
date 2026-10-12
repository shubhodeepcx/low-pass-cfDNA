"""Figures.

Three charts, each answering one question a reader will actually ask.

Palette: three categorical hues used in fixed order and never cycled, with
every series labelled directly on the line rather than in a legend box, so
identity never rests on colour alone. Reference lines are drawn in muted ink
rather than a series colour, because a reference is not a fourth series.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import (
    FixedFormatter,
    FixedLocator,
    NullFormatter,
    NullLocator,
)

from .analysis import limit_of_detection
from .junction import CapturePanel, JunctionModel, gigabases, required_coverage

# Headless by default: these figures are built in CI and on servers, where no
# display exists and the interactive backends fail at import time.
matplotlib.use("Agg")

__all__ = ["make_all", "fig_depth_plateau", "fig_power_curves",
           "fig_junction_cost"]

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
INK_MUTED = "#8a8984"
GRID = "#e4e3df"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")

LESION_LABELS = {
    "del(20q)": "del(20q), 16 Mb",
    "del(20q)_small": "del(20q), 5 Mb",
    "del(5q)": "del(5q), 87 Mb",
    "del(7q)": "del(7q), 55 Mb",
    "complex": "complex, 434 Mb",
    "minus7": "monosomy 7, 159 Mb",
    "plus8": "trisomy 8, 145 Mb",
}


def _style(ax: plt.Axes) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SOFT, labelsize=9, length=3)
    ax.grid(True, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def _log_ticks(ax: plt.Axes, which: str, values, labels) -> None:
    """Explicit ticks on a log axis, with the minor decade ticks suppressed."""
    axis = ax.xaxis if which == "x" else ax.yaxis
    axis.set_major_locator(FixedLocator(list(values)))
    axis.set_major_formatter(FixedFormatter(list(labels)))
    axis.set_minor_locator(NullLocator())
    axis.set_minor_formatter(NullFormatter())


def _titles(ax: plt.Axes, headline: str, subtitle: str) -> None:
    ax.set_title(headline, color=INK, fontsize=13, fontweight="bold",
                 loc="left", pad=16)
    ax.text(0.0, 1.018, subtitle, transform=ax.transAxes, color=INK_MUTED,
            fontsize=9.5)


def _save(fig: plt.Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return path


def fig_depth_plateau(sweep: pd.DataFrame, out_dir: Path,
                      lesions=("del(20q)", "del(5q)", "complex")) -> Path:
    """Detection limit against depth. The curves flatten, and that is the point."""
    offsets = {lesions[0]: 0, lesions[1]: 2, lesions[2]: -2}
    fig, ax = plt.subplots(figsize=(7.6, 4.7), facecolor=SURFACE)
    _style(ax)
    plotted: list[float] = []

    for colour, lesion in zip(SERIES, lesions, strict=True):
        group = sweep[sweep["lesion"] == lesion]
        points = sorted((float(depth), limit_of_detection(part))
                        for depth, part in group.groupby("depth"))
        xs = [p[0] for p in points if not np.isnan(p[1])]
        ys = [p[1] for p in points if not np.isnan(p[1])]
        if not xs:
            continue
        plotted.extend(ys)
        ax.plot(xs, ys, "-o", color=colour, linewidth=2.0, markersize=8,
                markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
        ax.annotate(LESION_LABELS.get(lesion, lesion), xy=(xs[-1], ys[-1]),
                    xytext=(10, offsets[lesion]), textcoords="offset points",
                    color=colour, fontsize=10, fontweight="bold",
                    va="center")

    ax.axvspan(0.5, 6.0, color=INK_MUTED, alpha=0.07, zorder=1, linewidth=0)
    ax.annotate("past about 0.5x the gains go small:\n"
                "bias limited, not depth limited",
                xy=(1.6, 0.30), color=INK_SOFT, fontsize=9.5, ha="center")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.008, 20)
    # Limits from the data rather than hard coded, so a point can never fall
    # off the chart without anyone noticing.
    low, high = min(plotted), max(plotted)
    ax.set_ylim(low * 0.8, high * 1.25)
    _log_ticks(ax, "x", [0.01, 0.05, 0.1, 0.5, 1.0, 5.0],
               ["0.01x", "0.05x", "0.1x", "0.5x", "1x", "5x"])
    ticks = [t for t in (0.03, 0.05, 0.075, 0.1, 0.15, 0.25, 0.4, 0.5)
             if low * 0.8 <= t <= high * 1.25]
    _log_ticks(ax, "y", ticks, [f"{t * 100:g}%" for t in ticks])
    ax.set_xlabel("sequencing depth", color=INK_SOFT, fontsize=10)
    ax.set_ylabel("smallest detectable tumour fraction", color=INK_SOFT,
                  fontsize=10)
    _titles(ax, "More sequencing stops helping",
            "limit of detection at 90% power, genome wide false positive rate 5%")
    return _save(fig, out_dir / "fig1_depth_plateau.png")


def fig_power_curves(sweep: pd.DataFrame, out_dir: Path, depth: float = 0.1,
                     lesions=("del(20q)_small", "del(5q)", "complex")) -> Path:
    """Detection probability against tumour fraction, with the copy neutral floor."""
    offsets = {lesions[0]: -1, lesions[1]: 17, lesions[2]: 0}
    fig, ax = plt.subplots(figsize=(7.6, 4.7), facecolor=SURFACE)
    _style(ax)

    at_depth = sweep[sweep["depth"] == depth]
    for colour, lesion in zip(SERIES, lesions, strict=True):
        group = at_depth[at_depth["lesion"] == lesion].sort_values(
            "tumour_fraction")
        if group.empty:
            continue
        ax.plot(group["tumour_fraction"], group["detected"], "-o",
                color=colour, linewidth=2.0, markersize=8,
                markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
        ax.annotate(LESION_LABELS.get(lesion, lesion),
                    xy=(group["tumour_fraction"].iloc[-1],
                        group["detected"].iloc[-1]),
                    xytext=(10, offsets[lesion]), textcoords="offset points",
                    color=colour, fontsize=10, fontweight="bold", va="center")

    floor = float(at_depth[at_depth["copy_neutral"]]["detected"].mean())
    ax.axhline(floor, color=INK_MUTED, linewidth=2.0, linestyle=(0, (4, 3)),
               zorder=2)
    ax.annotate(f"balanced translocations sit at {floor * 100:.0f}%, "
                "the false positive rate,\n"
                "at every tumour fraction and every depth",
                xy=(0.0, floor), xytext=(0, -14), textcoords="offset points",
                color=INK_SOFT, fontsize=9.5, va="top")

    ax.set_xlabel("tumour fraction", color=INK_SOFT, fontsize=10)
    ax.set_ylabel("probability of detection", color=INK_SOFT, fontsize=10)
    _titles(ax, f"Detection probability at {depth:g}x",
            f"{int(sweep['n_replicates'].iloc[0])} simulated libraries per point")
    ax.set_ylim(-0.26, 1.08)
    ax.set_xlim(-0.005, 0.56)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "25%", "50%", "75%", "100%"])
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_xticklabels(["0", "10%", "20%", "30%", "40%"])
    return _save(fig, out_dir / "fig2_power_curves.png")


def fig_junction_cost(out_dir: Path, model: JunctionModel | None = None,
                      panel: CapturePanel | None = None) -> Path:
    """Sequencing required for junction level evidence, targeted or not."""
    model = model or JunctionModel()
    panel = panel or CapturePanel()
    fractions = np.array([0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.25, 0.40])
    coverage = np.array([required_coverage(f, model=model) for f in fractions])

    fig, ax = plt.subplots(figsize=(7.6, 4.7), facecolor=SURFACE)
    _style(ax)
    for colour, values, label in (
        (SERIES[0], np.array([gigabases(c) for c in coverage]),
         "untargeted whole genome"),
        (SERIES[1], np.array([gigabases(c, panel.target_bp) for c in coverage]),
         f"{panel.target_bp / 1e6:.1f} Mb capture panel"),
    ):
        ax.plot(fractions, values, "-o", color=colour, linewidth=2.0,
                markersize=8, markeredgecolor=SURFACE, markeredgewidth=2,
                zorder=3)
        ax.annotate(label, xy=(fractions[-1], values[-1]), xytext=(10, 0),
                    textcoords="offset points", color=colour, fontsize=10,
                    fontweight="bold", va="center")

    ax.axhspan(1, 15, color=INK_MUTED, alpha=0.07, zorder=1, linewidth=0)
    ax.annotate("a routine low pass run\ndelivers roughly this much",
                xy=(0.011, 3.0), color=INK_SOFT, fontsize=9.5, va="center")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.008, 0.95)
    _log_ticks(ax, "x", [0.01, 0.02, 0.05, 0.1, 0.25, 0.4],
               ["1%", "2%", "5%", "10%", "25%", "40%"])
    _log_ticks(ax, "y", [0.01, 0.1, 1, 10, 100, 1000, 10000],
               ["0.01", "0.1", "1", "10", "100", "1,000", "10,000"])
    ax.set_xlabel("tumour fraction", color=INK_SOFT, fontsize=10)
    ax.set_ylabel("sequencing per sample (Gbp, log scale)", color=INK_SOFT,
                  fontsize=10)
    _titles(ax, "What a balanced translocation costs to see",
            f"coverage for {model.min_reads} junction reads at 90% power; "
            "the gap is evidence, not algorithm")
    return _save(fig, out_dir / "fig3_junction_cost.png")


def make_all(sweep: pd.DataFrame, out_dir: Path) -> list[Path]:
    """Build every figure and return the paths written."""
    return [fig_depth_plateau(sweep, out_dir),
            fig_power_curves(sweep, out_dir),
            fig_junction_cost(out_dir)]
