"""Interpolation and summary statistics, including the end to end smoke test."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lpcfdna.analysis import (
    copy_neutral_summary,
    depth_improvement,
    limit_of_detection,
    lod_matrix,
    lod_table,
)
from lpcfdna.sweep import SweepConfig, run_sweep


def _frame(fractions, detected, **extra):
    base = dict(lesion="x", copy_neutral=False, lesion_mb=50.0, depth=0.1,
                n_replicates=100)
    base.update(extra)
    return pd.DataFrame([dict(tumour_fraction=f, detected=d, localised=d,
                              **base) for f, d in zip(fractions, detected, strict=True)])


def test_interpolation_lands_between_the_bracketing_points():
    frame = _frame([0.05, 0.10], [0.80, 1.00])
    assert limit_of_detection(frame, 0.90) == pytest.approx(0.075)


def test_target_already_met_at_the_first_point():
    frame = _frame([0.05, 0.10], [0.95, 1.00])
    assert limit_of_detection(frame, 0.90) == 0.05


def test_never_reached_is_nan_not_a_silent_maximum():
    frame = _frame([0.05, 0.10, 0.40], [0.05, 0.08, 0.20])
    assert np.isnan(limit_of_detection(frame, 0.90))


def test_flat_segment_does_not_divide_by_zero():
    frame = _frame([0.05, 0.10], [0.90, 0.90])
    assert limit_of_detection(frame, 0.90) == pytest.approx(0.05)


def test_unordered_input_is_sorted_first():
    frame = _frame([0.10, 0.05], [1.00, 0.80])
    assert limit_of_detection(frame, 0.90) == pytest.approx(0.075)


def test_lod_table_and_matrix_agree():
    rows = pd.concat([
        _frame([0.05, 0.10], [0.80, 1.00], lesion="a", lesion_mb=50.0, depth=0.1),
        _frame([0.05, 0.10], [0.95, 1.00], lesion="a", lesion_mb=50.0, depth=1.0),
    ])
    table = lod_table(rows)
    matrix = lod_matrix(rows)
    assert set(table["depth"]) == {0.1, 1.0}
    assert matrix.loc[50.0, 1.0] < matrix.loc[50.0, 0.1]


def test_depth_improvement_reports_percent_and_checks_its_arguments():
    rows = pd.concat([
        _frame([0.05, 0.10], [0.80, 1.00], lesion="a", depth=0.1),
        _frame([0.05, 0.10], [0.95, 1.00], lesion="a", depth=1.0),
    ])
    improvement = depth_improvement(rows, 0.1, 1.0)
    assert improvement.iloc[0] == pytest.approx(33.3, abs=0.2)
    with pytest.raises(KeyError):
        depth_improvement(rows, 0.1, 9.9)


def test_copy_neutral_summary_needs_copy_neutral_rows():
    with pytest.raises(ValueError):
        copy_neutral_summary(_frame([0.1], [0.05]))


@pytest.mark.slow
def test_quick_sweep_runs_end_to_end():
    """A miniature of the real study. Guards the whole pipeline, not a unit."""
    config = SweepConfig(
        bin_size=2_000_000, depths=(0.1,),
        tumour_fractions=(0.02, 0.10, 0.40),
        lesions=("complex", "t(9;22)"),
        n_panel=25, n_calibration=120, n_fpr_check=60, n_replicates=25,
    )
    rows, calibration = run_sweep(config, verbose=False)
    sweep = pd.DataFrame(rows)

    assert len(rows) == config.n_conditions
    assert len(calibration) == 1
    assert abs(calibration[0]["achieved_fpr"] - 0.05) < 0.10

    big = sweep[(sweep["lesion"] == "complex") & (sweep["tumour_fraction"] == 0.40)]
    assert float(big["detected"].iloc[0]) > 0.9

    balanced = sweep[sweep["copy_neutral"]]
    assert float(balanced["detected"].max()) < 0.25
    summary = copy_neutral_summary(sweep)
    assert summary["n_conditions"] == 3
