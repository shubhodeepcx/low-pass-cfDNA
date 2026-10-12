"""Turning the raw sweep into the numbers a reader actually wants.

The sweep records a detection probability per condition. What a clinician or a
study designer asks instead is "what is the smallest tumour fraction this can
see", which is a quantity the sweep only brackets. This module does the
interpolation, with the sampling grid's limits made explicit rather than
hidden.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["limit_of_detection", "lod_table", "lod_matrix",
           "depth_improvement", "copy_neutral_summary"]

DEFAULT_POWER = 0.90


def limit_of_detection(frame: pd.DataFrame, power: float = DEFAULT_POWER
                       ) -> float:
    """Smallest tumour fraction reaching ``power``, linearly interpolated.

    Returns ``nan`` when the grid never reaches the target, which is a real
    answer and not a missing value: it means the lesion was undetectable at
    every tumour fraction tested. Treating that as missing is how a sweep ends
    up quietly reporting only its successes.
    """
    ordered = frame.sort_values("tumour_fraction")
    x = ordered["tumour_fraction"].to_numpy(dtype=float)
    y = ordered["detected"].to_numpy(dtype=float)
    reached = np.flatnonzero(y >= power)
    if reached.size == 0:
        return float("nan")
    i = int(reached[0])
    if i == 0:
        return float(x[0])
    x0, x1, y0, y1 = x[i - 1], x[i], y[i - 1], y[i]
    if y1 == y0:
        return float(x1)
    return float(x0 + (power - y0) * (x1 - x0) / (y1 - y0))


def lod_table(sweep: pd.DataFrame, power: float = DEFAULT_POWER) -> pd.DataFrame:
    """One row per lesion and depth, carrying the interpolated limit."""
    records = []
    for (lesion, depth), group in sweep.groupby(["lesion", "depth"]):
        records.append(dict(
            lesion=lesion,
            lesion_mb=float(group["lesion_mb"].iloc[0]),
            copy_neutral=bool(group["copy_neutral"].iloc[0]),
            depth=float(depth),
            lod=limit_of_detection(group, power),
            max_detection=float(group["detected"].max()),
        ))
    return pd.DataFrame(records).sort_values(["lesion_mb", "depth"])


def lod_matrix(sweep: pd.DataFrame, power: float = DEFAULT_POWER) -> pd.DataFrame:
    """Lesion size by depth matrix of detection limits, copy neutral excluded."""
    table = lod_table(sweep, power)
    return (table[~table["copy_neutral"]]
            .pivot(index="lesion_mb", columns="depth", values="lod"))


def depth_improvement(sweep: pd.DataFrame, low: float, high: float,
                      power: float = DEFAULT_POWER) -> pd.Series:
    """Relative improvement in detection limit between two depths, in percent.

    This is the quantity behind the headline claim about the plateau: going
    from 0.01x to 0.1x buys a great deal, going from 0.5x to 5x buys almost
    nothing.
    """
    matrix = lod_matrix(sweep, power)
    for depth in (low, high):
        if depth not in matrix.columns:
            raise KeyError(f"depth {depth} is not in the sweep")
    return ((matrix[low] - matrix[high]) / matrix[low] * 100.0).round(1)


def copy_neutral_summary(sweep: pd.DataFrame) -> dict:
    """Detection statistics restricted to balanced rearrangements.

    If the detector is honest these should sit at the false positive rate, with
    no trend in either tumour fraction or depth. Any deviation is a bug in the
    simulation or the detector, not a discovery.
    """
    cn = sweep[sweep["copy_neutral"]]
    if cn.empty:
        raise ValueError("the sweep contains no copy neutral conditions")
    by_depth = cn.groupby("depth")["detected"].mean()
    by_fraction = cn.groupby("tumour_fraction")["detected"].mean()
    fractions = by_fraction.index.to_numpy(dtype=float)
    slope = float(np.polyfit(fractions, by_fraction.to_numpy(), 1)[0])
    return dict(
        n_conditions=int(len(cn)),
        mean_detection=float(cn["detected"].mean()),
        max_detection=float(cn["detected"].max()),
        by_depth=by_depth.round(4).to_dict(),
        trend_with_tumour_fraction=round(slope, 4),
    )
