"""Low pass cell free DNA library simulation.

The coverage model is a product of five terms, and the split between them is
the only part of this package that really matters for the conclusions:

===========================  =================================  ==============
term                         behaviour across libraries         removable?
===========================  =================================  ==============
shared systematic bias       identical in every library         yes, by a panel
per library bias             redrawn for every library          no
GC efficiency curve          same shape, per library strength   partly
mappability                  fixed property of the reference    yes, by masking
Poisson counting noise       shrinks as depth grows             yes, with depth
===========================  =================================  ==============

The per library term is the one that sets a floor. Without it a simulation
reports detection limits that keep improving as coverage grows, which is
flattering and wrong: real callers plateau, and they plateau because library to
library variation does not average away no matter how many reads you buy.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .genome import BinnedGenome, gc_efficiency, lesion_copy_number

__all__ = ["LibraryModel", "simulate_library", "simulate_panel",
           "mixture_ratio"]


@dataclass(frozen=True)
class LibraryModel:
    """Nuisance parameters of the sequencing protocol.

    Attributes
    ----------
    read_length_bp:
        Used only to turn a coverage figure into a read count per bin.
    per_library_bias_sd:
        Standard deviation of the per bin log bias that is redrawn for every
        library. This single number sets the detection floor, so it is the
        first thing to fit against real data before trusting any absolute
        number this package produces.
    gc_strength_mean, gc_strength_sd:
        Mean and spread of the per library GC curve sharpness.
    """

    read_length_bp: float = 150.0
    per_library_bias_sd: float = 0.045
    gc_strength_mean: float = 6.0
    gc_strength_sd: float = 1.2
    name: str = "default"
    notes: str = field(default="", repr=False)

    def __post_init__(self) -> None:
        if self.read_length_bp <= 0:
            raise ValueError("read_length_bp must be positive")
        if self.per_library_bias_sd < 0:
            raise ValueError("per_library_bias_sd must be non negative")
        if self.gc_strength_sd < 0:
            raise ValueError("gc_strength_sd must be non negative")


#: A protocol with no library to library bias at all. Not realistic, but it is
#: the right control for isolating counting noise, and the tests use it to
#: check the mixture arithmetic in closed form.
IDEAL_LIBRARY = LibraryModel(per_library_bias_sd=0.0, gc_strength_sd=0.0,
                             name="ideal",
                             notes="no per library bias; for tests only")


def mixture_ratio(tumour_copy_number: np.ndarray,
                  tumour_fraction: float) -> np.ndarray:
    """Expected depth relative to a diploid library.

    A fraction ``f`` of molecules come from a clone with copy number ``c`` and
    the rest from diploid cells, so the expected depth scales as
    ``(f * c + (1 - f) * 2) / 2``. A single copy loss at fraction ``f``
    therefore depresses depth by ``f / 2``, not by ``f``, which is the factor
    of two most easily got wrong.
    """
    if not 0.0 <= tumour_fraction <= 1.0:
        raise ValueError("tumour_fraction must lie in [0, 1]")
    return (tumour_fraction * tumour_copy_number
            + (1.0 - tumour_fraction) * 2.0) / 2.0


def simulate_library(
    genome: BinnedGenome,
    rng: np.random.Generator,
    depth: float = 0.1,
    lesion: str = "normal",
    tumour_fraction: float = 0.0,
    model: LibraryModel | None = None,
    tumour_copy_number: np.ndarray | None = None,
) -> np.ndarray:
    """Simulate per bin read counts for one cell free DNA library.

    Parameters
    ----------
    depth:
        Genome wide mean coverage, so 0.1 means a 0.1x low pass run.
    lesion:
        Name from :data:`lpcfdna.genome.LESIONS`. Ignored when
        ``tumour_copy_number`` is supplied.
    tumour_fraction:
        Fraction of molecules originating from the abnormal clone.
    tumour_copy_number:
        Precomputed copy number vector. Passing this skips a lookup that would
        otherwise repeat tens of thousands of times inside a sweep.

    Returns
    -------
    numpy.ndarray
        Read counts per bin, float valued, zero in excluded bins.
    """
    if depth <= 0:
        raise ValueError("depth must be positive")
    model = model or LibraryModel()

    cn = (lesion_copy_number(genome, lesion)
          if tumour_copy_number is None else tumour_copy_number)
    if cn.shape != (genome.n_bins,):
        raise ValueError("tumour_copy_number has the wrong shape for this genome")

    reads_per_bin = depth * genome.bin_size / model.read_length_bp
    strength = rng.normal(model.gc_strength_mean, model.gc_strength_sd)
    strength = max(strength, 0.0)
    per_library = np.exp(rng.normal(0.0, model.per_library_bias_sd,
                                    genome.n_bins))

    expected = (reads_per_bin
                * genome.shared_bias
                * per_library
                * gc_efficiency(genome.gc, strength)
                * genome.mappability
                * mixture_ratio(cn, tumour_fraction))
    expected[genome.excluded] = 0.0
    return rng.poisson(expected).astype(float)


def simulate_panel(
    genome: BinnedGenome,
    rng: np.random.Generator,
    n: int = 100,
    depth: float = 0.1,
    model: LibraryModel | None = None,
) -> np.ndarray:
    """Simulate ``n`` normal libraries, shape ``(n, n_bins)``.

    Used both for the reference panel a detector is fitted on and for the
    held out normals it is calibrated against. Those two sets must be
    disjoint, which is why this returns a plain array rather than caching
    anything.
    """
    if n <= 0:
        raise ValueError("n must be positive")
    model = model or LibraryModel()
    out = np.empty((n, genome.n_bins))
    diploid = np.full(genome.n_bins, 2.0)
    for i in range(n):
        out[i] = simulate_library(genome, rng, depth=depth, model=model,
                                  tumour_fraction=0.0,
                                  tumour_copy_number=diploid)
    return out
