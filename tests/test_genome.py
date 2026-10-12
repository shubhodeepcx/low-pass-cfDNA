"""Geometry, covariates and lesion definitions."""

from __future__ import annotations

import numpy as np
import pytest

from lpcfdna.genome import (
    CENTROMERES,
    CHROM_LENGTHS,
    COPY_NEUTRAL,
    LESIONS,
    BinnedGenome,
    Segment,
    gc_efficiency,
    lesion_copy_number,
    lesion_size_mb,
)


def test_bin_count_matches_chromosome_lengths(genome):
    expected = sum(length // genome.bin_size
                   for length in CHROM_LENGTHS.values())
    assert genome.n_bins == expected
    assert len(genome) == expected


def test_bins_do_not_overrun_their_chromosome(genome):
    for chrom, length in CHROM_LENGTHS.items():
        ends = genome.end[genome.chrom == chrom]
        assert ends.max() <= length


def test_every_bin_has_exactly_one_arm(genome):
    assert set(np.unique(genome.arm)) <= {0, 1}
    for chrom in CHROM_LENGTHS:
        p = genome.arm_mask(chrom, "p")
        q = genome.arm_mask(chrom, "q")
        assert not (p & q).any()
        assert (p | q).sum() == (genome.chrom == chrom).sum()


def test_arm_split_follows_the_centromere(genome):
    for chrom in CHROM_LENGTHS:
        p_starts = genome.start[genome.arm_mask(chrom, "p")]
        q_starts = genome.start[genome.arm_mask(chrom, "q")]
        if p_starts.size:
            assert p_starts.max() < CENTROMERES[chrom]
        if q_starts.size:
            assert q_starts.min() >= CENTROMERES[chrom]


def test_covariates_are_in_range(genome):
    assert 0.3 < float(genome.gc.mean()) < 0.5
    assert genome.gc.min() >= 0.25 and genome.gc.max() <= 0.65
    assert genome.mappability.min() >= 0.0
    assert genome.mappability.max() <= 1.0
    assert 0.0 < genome.usable.mean() <= 1.0


def test_genome_is_reproducible_from_its_seed():
    a = BinnedGenome(1_000_000, seed=7)
    b = BinnedGenome(1_000_000, seed=7)
    assert np.array_equal(a.gc, b.gc)
    assert np.array_equal(a.shared_bias, b.shared_bias)
    assert not np.array_equal(a.gc, BinnedGenome(1_000_000, seed=8).gc)


def test_bad_bin_size_is_rejected():
    with pytest.raises(ValueError):
        BinnedGenome(0)


def test_arm_label_and_mask_reject_nonsense(genome):
    with pytest.raises(ValueError):
        genome.arm_mask(5, "r")
    assert genome.arm_label(int(np.flatnonzero(genome.chrom == 5)[0])) == "5p"


def test_gc_efficiency_peaks_where_it_should():
    gc = np.array([0.30, 0.42, 0.55])
    eff = gc_efficiency(gc, strength=6.0)
    assert eff[1] == pytest.approx(1.0)
    assert eff[0] < eff[1] and eff[2] < eff[1]
    # A stronger curve must be at least as punishing away from the peak.
    assert gc_efficiency(gc, 12.0)[0] < eff[0]


@pytest.mark.parametrize("lesion", sorted(LESIONS))
def test_every_lesion_is_constructible(genome, lesion):
    cn = lesion_copy_number(genome, lesion)
    assert cn.shape == (genome.n_bins,)
    assert np.isfinite(cn).all()
    assert cn.min() >= 0


def test_unknown_lesion_names_the_alternatives(genome):
    with pytest.raises(KeyError, match="del\\(5q\\)"):
        lesion_copy_number(genome, "del(99q)")


def test_losses_and_gains_point_the_right_way(genome):
    loss = lesion_copy_number(genome, "del(5q)")
    assert (loss[genome.arm_mask(5, "q")] <= 2).all()
    assert (loss[genome.arm_mask(5, "q")] == 1).any()
    assert (loss[genome.chrom != 5] == 2).all()

    gain = lesion_copy_number(genome, "plus8")
    assert (gain[genome.chrom == 8] == 3).all()
    assert (gain[genome.chrom != 8] == 2).all()


@pytest.mark.parametrize("lesion", sorted(COPY_NEUTRAL | {"normal"}))
def test_copy_neutral_lesions_leave_every_bin_diploid(genome, lesion):
    assert (lesion_copy_number(genome, lesion) == 2).all()
    assert lesion_size_mb(genome, lesion) == 0.0


def test_lesion_sizes_are_physically_plausible(genome):
    assert 60 <= lesion_size_mb(genome, "del(5q)") <= 110
    whole_eight = CHROM_LENGTHS[8] / 1e6
    assert abs(lesion_size_mb(genome, "plus8") - whole_eight) < 1.0
    sizes = {name: lesion_size_mb(genome, name)
             for name in ("del(20q)_small", "del(5q)", "complex")}
    assert sizes["del(20q)_small"] < sizes["del(5q)"] < sizes["complex"]


def test_segment_validation_rejects_bad_geometry():
    with pytest.raises(ValueError):
        Segment(5, "r", 0.0, 1.0, 1)
    with pytest.raises(ValueError):
        Segment(5, "q", 0.8, 0.2, 1)
    with pytest.raises(ValueError):
        Segment(99, "q", 0.0, 1.0, 1)


def test_lesion_extent_is_stable_across_bin_sizes(genome, coarse_genome):
    fine = lesion_size_mb(genome, "del(5q)")
    coarse = lesion_size_mb(coarse_genome, "del(5q)")
    assert abs(fine - coarse) < 6.0
