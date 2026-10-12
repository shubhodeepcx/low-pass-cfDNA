"""Detectability of cytogenetic lesions in low pass cell free DNA.

A simulation study of what a bin count copy number caller can and cannot see
when you try to read a haematological karyotype from a blood draw, and a
closed form answer for the lesions it can never see at all.

The short version of the findings:

1. Past roughly 0.5x coverage, extra sequencing buys almost nothing, because
   library to library bias rather than counting noise sets the floor.
2. Lesions below about 30 megabases need tumour fractions high enough that the
   diagnosis was never in doubt.
3. Balanced translocations are detected at exactly the false positive rate, and
   the evidence needed to see them is three orders of magnitude beyond what a
   low pass run delivers. That is a coverage allocation problem, not an
   algorithm problem.

Typical use::

    from lpcfdna import BinnedGenome, simulate_panel, BinCountDetector
    import numpy as np

    genome = BinnedGenome(500_000, seed=1)
    rng = np.random.default_rng(0)
    detector = BinCountDetector(genome, simulate_panel(genome, rng, 100))
    detector.calibrate(simulate_panel(genome, rng, 1000))
"""

from __future__ import annotations

from .analysis import limit_of_detection, lod_matrix, lod_table
from .detect import (
    BinCountDetector,
    CallResult,
    is_localised,
    localisation_precision,
)
from .genome import (
    COPY_NEUTRAL,
    GENOME_HAPLOID_BP,
    LESIONS,
    SEQUENCED_GENOME_BP,
    BinnedGenome,
    lesion_copy_number,
    lesion_size_mb,
)
from .junction import (
    CapturePanel,
    JunctionModel,
    detection_power,
    expected_junction_reads,
    gigabases,
    required_coverage,
)
from .simulate import LibraryModel, simulate_library, simulate_panel
from .sweep import SweepConfig, run_sweep

__version__ = "1.0.0"

__all__ = [
    "__version__",
    "BinnedGenome", "LESIONS", "COPY_NEUTRAL", "lesion_copy_number",
    "lesion_size_mb", "GENOME_HAPLOID_BP", "SEQUENCED_GENOME_BP",
    "LibraryModel", "simulate_library", "simulate_panel",
    "BinCountDetector", "CallResult", "is_localised",
    "localisation_precision",
    "JunctionModel", "CapturePanel", "expected_junction_reads",
    "detection_power", "required_coverage", "gigabases",
    "SweepConfig", "run_sweep",
    "limit_of_detection", "lod_table", "lod_matrix",
]
