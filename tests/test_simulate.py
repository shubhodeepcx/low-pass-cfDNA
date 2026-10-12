"""The coverage model, especially the mixture arithmetic.

The test that matters most here is the closed form check on the mixture ratio.
A single copy loss at tumour fraction f depresses depth by f/2, not by f, and
an off by two there would halve every detection limit the study reports while
leaving all the plumbing apparently healthy.
"""

from __future__ import annotations

import numpy as np
import pytest

from lpcfdna.detect import gc_correct
from lpcfdna.genome import lesion_copy_number
from lpcfdna.simulate import (
    LibraryModel,
    mixture_ratio,
    simulate_library,
    simulate_panel,
)


def test_mixture_ratio_closed_form():
    cn = np.array([1.0, 2.0, 3.0])
    assert mixture_ratio(cn, 0.0) == pytest.approx([1.0, 1.0, 1.0])
    assert mixture_ratio(cn, 1.0) == pytest.approx([0.5, 1.0, 1.5])
    assert mixture_ratio(cn, 0.5) == pytest.approx([0.75, 1.0, 1.25])


def test_mixture_ratio_rejects_impossible_fractions():
    with pytest.raises(ValueError):
        mixture_ratio(np.array([2.0]), 1.5)
    with pytest.raises(ValueError):
        mixture_ratio(np.array([2.0]), -0.1)


def test_counts_are_non_negative_and_excluded_bins_are_empty(genome, rng):
    counts = simulate_library(genome, rng, depth=0.1)
    assert counts.shape == (genome.n_bins,)
    assert (counts >= 0).all()
    assert (counts[genome.excluded] == 0).all()


def test_depth_scales_the_read_count(genome, rng, ideal_library):
    """Checked on the ideal protocol, where the only randomness left is Poisson.

    On the default protocol this ratio carries the two libraries' GC curve
    draws as well, which is real variation rather than a bug, so it belongs in
    a test about bias and not in a test about depth.
    """
    low = simulate_library(genome, rng, depth=0.05, model=ideal_library).sum()
    high = simulate_library(genome, rng, depth=0.50, model=ideal_library).sum()
    assert high / low == pytest.approx(10.0, rel=0.02)


@pytest.mark.parametrize("fraction", [0.2, 0.5, 1.0])
def test_single_copy_loss_depresses_depth_by_half_the_fraction(
        genome, rng, ideal_library, fraction):
    """Checked at high depth, where counting noise is negligible.

    Chromosome 7 has its own GC and mappability profile, so even a tumour free
    library does not put its median exactly on the genome median. The test
    measures that offset on a diploid library first and then asserts the
    mixture term against it, which isolates the arithmetic under test from the
    bias structure around it.
    """
    cn = lesion_copy_number(genome, "minus7")

    def on_off_ratio(tumour_fraction: float, seed: int) -> float:
        counts = simulate_library(genome, np.random.default_rng(seed),
                                  depth=200.0, tumour_copy_number=cn,
                                  tumour_fraction=tumour_fraction,
                                  model=ideal_library)
        corrected = gc_correct(counts, genome.gc, genome.usable)
        on_seven = (genome.chrom == 7) & genome.usable
        off_seven = (genome.chrom != 7) & genome.usable
        return (float(np.median(corrected[on_seven]))
                / float(np.median(corrected[off_seven])))

    baseline = on_off_ratio(0.0, seed=101)
    observed = on_off_ratio(fraction, seed=101)
    assert observed / baseline == pytest.approx(1.0 - fraction / 2.0, abs=0.004)


def test_gain_raises_depth(genome, rng, ideal_library):
    cn = lesion_copy_number(genome, "plus8")
    counts = simulate_library(genome, rng, depth=200.0, tumour_copy_number=cn,
                              tumour_fraction=0.5, model=ideal_library)
    corrected = gc_correct(counts, genome.gc, genome.usable)
    on = float(np.median(corrected[(genome.chrom == 8) & genome.usable]))
    off = float(np.median(corrected[(genome.chrom != 8) & genome.usable]))
    assert on > off * 1.1


def test_balanced_rearrangement_is_indistinguishable_from_normal(genome):
    """The premise of the whole copy neutral argument, checked directly."""
    normal = simulate_library(genome, np.random.default_rng(5), depth=1.0,
                              lesion="normal", tumour_fraction=0.0)
    balanced = simulate_library(genome, np.random.default_rng(5), depth=1.0,
                                lesion="t(9;22)", tumour_fraction=0.9)
    assert np.array_equal(normal, balanced)


def test_zero_tumour_fraction_matches_a_normal_library(genome):
    with_lesion = simulate_library(genome, np.random.default_rng(3), depth=0.5,
                                   lesion="complex", tumour_fraction=0.0)
    plain = simulate_library(genome, np.random.default_rng(3), depth=0.5,
                             lesion="normal", tumour_fraction=0.0)
    assert np.array_equal(with_lesion, plain)


def test_per_library_bias_adds_variance_that_depth_cannot_remove(genome):
    """The floor exists. Without it the whole plateau finding would be false."""
    deep = dict(depth=50.0)
    noisy = LibraryModel(per_library_bias_sd=0.08, gc_strength_sd=0.0)
    quiet = LibraryModel(per_library_bias_sd=0.0, gc_strength_sd=0.0)

    def spread(model):
        rng = np.random.default_rng(11)
        panel = simulate_panel(genome, rng, n=12, model=model, **deep)
        corrected = np.vstack([gc_correct(r, genome.gc, genome.usable)
                               for r in panel])
        usable = genome.usable
        return float(np.mean(corrected[:, usable].std(axis=0, ddof=1)
                             / corrected[:, usable].mean(axis=0)))

    assert spread(noisy) > 3 * spread(quiet)


def test_panel_rows_are_independent(genome, rng):
    panel = simulate_panel(genome, rng, n=5, depth=0.1)
    assert panel.shape == (5, genome.n_bins)
    for i in range(1, 5):
        assert not np.array_equal(panel[0], panel[i])


def test_invalid_arguments_are_rejected(genome, rng):
    with pytest.raises(ValueError):
        simulate_library(genome, rng, depth=0.0)
    with pytest.raises(ValueError):
        simulate_library(genome, rng, tumour_fraction=2.0)
    with pytest.raises(ValueError):
        simulate_library(genome, rng, tumour_copy_number=np.ones(5))
    with pytest.raises(ValueError):
        simulate_panel(genome, rng, n=0)
    with pytest.raises(ValueError):
        LibraryModel(read_length_bp=0)
