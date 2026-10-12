"""The closed form junction evidence calculation.

Two properties are worth nailing down. The expected count must not depend on
genome length or fragment length, because the derivation says they cancel; and
the required coverage must scale as one over the tumour fraction, because that
scaling is the whole argument for why low pass sequencing cannot get there.
"""

from __future__ import annotations

import numpy as np
import pytest

from lpcfdna.genome import GENOME_HAPLOID_BP, SEQUENCED_GENOME_BP
from lpcfdna.junction import (
    CapturePanel,
    JunctionModel,
    detection_power,
    expected_junction_reads,
    gigabases,
    power_table,
    report,
    required_coverage,
)


def test_expected_reads_matches_the_closed_form():
    model = JunctionModel(fragment_bp=200.0, anchor_bp=0.0)
    assert model.usable_fraction == pytest.approx(1.0)
    assert expected_junction_reads(10.0, 0.20, model) == pytest.approx(1.0)


def test_expected_reads_is_linear_in_both_arguments():
    model = JunctionModel()
    base = expected_junction_reads(1.0, 0.1, model)
    assert expected_junction_reads(2.0, 0.1, model) == pytest.approx(2 * base)
    assert expected_junction_reads(1.0, 0.2, model) == pytest.approx(2 * base)


def test_fragment_length_cancels_once_the_anchor_is_zero():
    short = JunctionModel(fragment_bp=120.0, anchor_bp=0.0)
    long = JunctionModel(fragment_bp=320.0, anchor_bp=0.0)
    assert (expected_junction_reads(1.0, 0.1, short)
            == pytest.approx(expected_junction_reads(1.0, 0.1, long)))


def test_longer_anchors_cost_sensitivity():
    lenient = JunctionModel(anchor_bp=10.0)
    strict = JunctionModel(anchor_bp=60.0)
    assert (expected_junction_reads(1.0, 0.1, strict)
            < expected_junction_reads(1.0, 0.1, lenient))


def test_a_fragment_too_short_for_two_anchors_is_rejected():
    with pytest.raises(ValueError):
        JunctionModel(fragment_bp=30.0, anchor_bp=20.0)
    with pytest.raises(ValueError):
        JunctionModel(min_reads=0)


def test_power_is_monotone_and_bounded():
    model = JunctionModel()
    powers = [detection_power(c, 0.1, model) for c in (0.1, 1, 10, 100, 1000)]
    assert all(0.0 <= p <= 1.0 for p in powers)
    assert powers == sorted(powers)
    assert powers[0] < 1e-6
    assert powers[-1] > 0.999


def test_required_coverage_scales_inversely_with_tumour_fraction():
    a = required_coverage(0.05)
    b = required_coverage(0.10)
    assert a / b == pytest.approx(2.0, rel=0.02)


def test_required_coverage_round_trips_through_power():
    for fraction in (0.02, 0.1, 0.4):
        coverage = required_coverage(fraction, power=0.9)
        assert detection_power(coverage, fraction) == pytest.approx(0.9,
                                                                    abs=0.01)


def test_required_coverage_edge_cases():
    assert required_coverage(0.0) == float("inf")
    with pytest.raises(ValueError):
        required_coverage(0.1, power=1.0)


def test_low_pass_yields_essentially_no_junction_evidence():
    """The headline number: four thousandths of a read, not a probability."""
    reads = expected_junction_reads(0.1, 0.10)
    assert reads < 0.01
    assert detection_power(0.1, 0.10) < 1e-6


def test_capture_enrichment_is_large_and_consistent():
    panel = CapturePanel(n_fusions=10, bp_per_region=150_000)
    assert panel.target_bp == pytest.approx(1.5e6)
    assert panel.genome_fraction < 0.001
    assert panel.enrichment > 1000
    assert (gigabases(100.0, panel.target_bp) * panel.enrichment
            == pytest.approx(gigabases(100.0), rel=1e-9))


def test_power_table_is_self_consistent():
    rows = power_table([0.05, 0.10])
    assert len(rows) == 2
    assert rows[0]["required_coverage"] > rows[1]["required_coverage"]
    for row in rows:
        assert row["wgs_gbp"] > row["capture_gbp"]
        assert np.isfinite(row["reads_at_0_1x"])


def test_report_mentions_the_things_a_reader_needs():
    text = report()
    for token in ("fragment length", "anchor", "enrichment", "capture"):
        assert token in text
    # The project writes no long dashes; chr() keeps this file clean too.
    assert chr(0x2014) not in text


def test_sequencing_cost_is_charged_against_the_whole_genome():
    """The bin grid is autosomal; a sequencing run is not.

    Quoting untargeted cost against the autosomes alone would understate a
    real run by the length of the sex chromosomes, so ``gigabases`` defaults
    to the sequenced genome rather than the binned one.
    """
    assert SEQUENCED_GENOME_BP > GENOME_HAPLOID_BP
    assert gigabases(1.0) == pytest.approx(SEQUENCED_GENOME_BP / 1e9)
    assert CapturePanel().enrichment == pytest.approx(
        SEQUENCED_GENOME_BP / CapturePanel().target_bp)
