"""The detector, its calibration, and the properties a caller must satisfy."""

from __future__ import annotations

import numpy as np
import pytest

from lpcfdna.detect import (
    BinCountDetector,
    candidate_segments,
    gc_correct,
    is_localised,
    localisation_overlap,
)
from lpcfdna.genome import lesion_copy_number
from lpcfdna.simulate import simulate_library, simulate_panel


@pytest.fixture(scope="module")
def fitted(genome):
    """A detector fitted and calibrated on disjoint sets of normal libraries."""
    rng = np.random.default_rng(17)
    detector = BinCountDetector(genome, simulate_panel(genome, rng, n=60,
                                                       depth=0.1))
    detector.calibrate(simulate_panel(genome, rng, n=600, depth=0.1))
    return detector


def test_gc_correction_flattens_the_gc_trend(genome, rng):
    counts = simulate_library(genome, rng, depth=0.5)
    corrected = gc_correct(counts, genome.gc, genome.usable)
    usable = genome.usable

    def gc_slope(values):
        return abs(float(np.polyfit(genome.gc[usable], values[usable], 1)[0])
                   / float(values[usable].mean()))

    assert gc_slope(corrected) < 0.5 * gc_slope(counts)


def test_gc_correction_handles_an_empty_library(genome):
    out = gc_correct(np.zeros(genome.n_bins), genome.gc, genome.usable)
    assert (out == 0).all()


def test_gc_correction_rejects_mismatched_shapes(genome):
    with pytest.raises(ValueError):
        gc_correct(np.zeros(5), genome.gc, genome.usable)


def test_candidate_segments_are_well_formed(genome):
    segments = candidate_segments(genome)
    assert segments.ndim == 2 and segments.shape[1] == 2
    assert (segments[:, 1] > segments[:, 0]).all()
    assert segments[:, 1].max() <= genome.n_bins
    # Every segment lies inside one chromosome arm.
    for lo, hi in segments[::37]:
        assert len(np.unique(genome.chrom[lo:hi])) == 1
        assert len(np.unique(genome.arm[lo:hi])) == 1


def test_detector_requires_a_usable_panel(genome):
    rng = np.random.default_rng(1)
    with pytest.raises(ValueError):
        BinCountDetector(genome, simulate_panel(genome, rng, n=2, depth=0.1))
    with pytest.raises(ValueError):
        BinCountDetector(genome, np.zeros((5, 10)))


def test_calibration_switches_the_null(genome, fitted):
    assert fitted.calibrated
    assert fitted.n_segments > 1000
    assert fitted.threshold(0.05) > fitted.threshold(0.20)


def test_calibration_rejects_mismatched_libraries(genome, fitted):
    with pytest.raises(ValueError):
        fitted.calibrate(np.zeros((4, 7)))


def test_threshold_rejects_impossible_rates(fitted):
    for bad in (0.0, 1.0, -0.2):
        with pytest.raises(ValueError):
            fitted.threshold(bad)


def test_false_positive_rate_is_close_to_nominal(genome, fitted):
    """The central guarantee. Out of sample, on libraries nothing has seen."""
    rng = np.random.default_rng(999)
    threshold = fitted.threshold(0.05)
    calls = [fitted.call(simulate_library(genome, rng, depth=0.1),
                         threshold=threshold).detected for _ in range(400)]
    realised = float(np.mean(calls))
    # 400 draws at a nominal 5 percent gives a 95 percent interval of
    # roughly plus or minus 2.1 points; allow 3 for calibration noise.
    assert abs(realised - 0.05) < 0.03, f"realised false positive rate {realised}"


def test_an_unmistakable_lesion_is_detected_and_localised(genome, fitted):
    rng = np.random.default_rng(4)
    cn = lesion_copy_number(genome, "minus7")
    result = fitted.call(simulate_library(genome, rng, depth=0.1,
                                          tumour_copy_number=cn,
                                          tumour_fraction=0.6))
    assert result.detected
    assert result.n_called_segments > 0
    assert is_localised(result.called_bins, cn)
    assert localisation_overlap(result.called_bins, cn) > 0.5


def test_detection_increases_with_tumour_fraction(genome, fitted):
    """Monotonicity in the only direction that can be asserted cheaply."""
    rng = np.random.default_rng(21)
    cn = lesion_copy_number(genome, "del(5q)")
    threshold = fitted.threshold(0.05)

    def rate(fraction, n=60):
        return float(np.mean([
            fitted.call(simulate_library(genome, rng, depth=0.1,
                                         tumour_copy_number=cn,
                                         tumour_fraction=fraction),
                        threshold=threshold).detected for _ in range(n)]))

    assert rate(0.02) < rate(0.25)
    assert rate(0.40) > 0.8


def test_balanced_rearrangements_are_called_at_the_false_positive_rate(
        genome, fitted):
    """Not a bug to fix. A property of every depth based method."""
    rng = np.random.default_rng(77)
    cn = lesion_copy_number(genome, "t(15;17)")
    threshold = fitted.threshold(0.05)
    rate = float(np.mean([
        fitted.call(simulate_library(genome, rng, depth=0.1,
                                     tumour_copy_number=cn,
                                     tumour_fraction=0.40),
                    threshold=threshold).detected for _ in range(300)]))
    assert rate < 0.12


def test_localisation_overlap_handles_a_normal_truth(genome):
    diploid = np.full(genome.n_bins, 2.0)
    called = np.zeros(genome.n_bins, dtype=bool)
    assert localisation_overlap(called, diploid) == 0.0
    assert not is_localised(called, diploid)


def test_call_result_reports_what_it_called(genome, fitted):
    rng = np.random.default_rng(31)
    cn = lesion_copy_number(genome, "plus8")
    result = fitted.call(simulate_library(genome, rng, depth=0.1,
                                          tumour_copy_number=cn,
                                          tumour_fraction=0.5))
    assert result.called_mb == int(result.called_bins.sum())
    assert np.isfinite(result.max_statistic)
    assert result.max_statistic > result.threshold
