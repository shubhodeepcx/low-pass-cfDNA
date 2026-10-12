"""Shared fixtures.

The genome is session scoped because building it is the slowest part of a test
run and because every test must share one genome: a detector fitted on one
genome's bias structure is meaningless against another's.
"""

from __future__ import annotations

import numpy as np
import pytest

from lpcfdna.genome import BinnedGenome
from lpcfdna.simulate import LibraryModel


@pytest.fixture(scope="session")
def genome() -> BinnedGenome:
    return BinnedGenome(bin_size=500_000, seed=1)


@pytest.fixture(scope="session")
def coarse_genome() -> BinnedGenome:
    """Two megabase bins, for tests that only need the structure to be right."""
    return BinnedGenome(bin_size=2_000_000, seed=1)


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(12345)


@pytest.fixture(scope="session")
def ideal_library() -> LibraryModel:
    """No library to library bias, so closed form checks are exact."""
    return LibraryModel(per_library_bias_sd=0.0, gc_strength_sd=0.0,
                        name="ideal")
