"""Binned reference genome and the karyotype lesions we simulate on it.

Everything downstream works in bin space, so this module owns the bin grid, the
per-bin covariates that drive coverage bias (GC content and mappability), and
the mapping from a named cytogenetic lesion to a per-bin tumour copy number.

GC content and mappability are modelled rather than read from a reference
FASTA. That keeps the package dependency free and lets the whole study run on a
laptop. A production pipeline swaps in measured per-bin values at exactly two
attribute assignments, :attr:`BinnedGenome.gc` and
:attr:`BinnedGenome.mappability`, and nothing else in the package changes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "CHROM_LENGTHS",
    "CENTROMERES",
    "BinnedGenome",
    "LESIONS",
    "COPY_NEUTRAL",
    "lesion_copy_number",
    "lesion_size_mb",
    "gc_efficiency",
    "GENOME_HAPLOID_BP",
    "SEQUENCED_GENOME_BP",
]

#: GRCh38 primary assembly autosome lengths, in base pairs.
CHROM_LENGTHS: dict[int, int] = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
    6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
    11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
    16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468,
}

#: Approximate GRCh38 centromere midpoints, used only to split p from q.
CENTROMERES: dict[int, int] = {
    1: 123400000, 2: 93900000, 3: 90900000, 4: 50000000, 5: 48800000,
    6: 59800000, 7: 60100000, 8: 45200000, 9: 43000000, 10: 39800000,
    11: 53400000, 12: 35500000, 13: 17700000, 14: 17200000, 15: 19000000,
    16: 36800000, 17: 25100000, 18: 18500000, 19: 26200000, 20: 28100000,
    21: 12000000, 22: 15000000,
}

#: Total autosomal length, which is the span the bin grid covers.
GENOME_HAPLOID_BP: int = sum(CHROM_LENGTHS.values())

#: What a whole genome run actually sequences: the autosomes plus the sex
#: chromosomes and the unplaced contigs that come along with them. Sequencing
#: cost is charged against this, not against the autosomes alone, so the
#: junction calculation uses it wherever gigabases enter the arithmetic.
SEQUENCED_GENOME_BP: int = 3_100_000_000


def gc_efficiency(gc: np.ndarray, strength: float, peak: float = 0.42,
                  width: float = 0.10) -> np.ndarray:
    """Relative capture efficiency as a function of bin GC content.

    A unimodal curve peaking near the genome wide mean GC, which is the shape
    library preparation and PCR actually produce. ``strength`` scales how
    sharply efficiency falls away from the peak; it is drawn per library in
    :func:`lpcfdna.simulate.simulate_library` because it varies between
    libraries even on one protocol.
    """
    return np.exp(-strength * ((gc - peak) / width) ** 2)


@dataclass(frozen=True)
class Segment:
    """One contiguous copy number change, expressed as a fraction of an arm."""

    chrom: int
    arm: str
    frac_start: float
    frac_end: float
    copy_number: int

    def __post_init__(self) -> None:
        if self.arm not in ("p", "q"):
            raise ValueError(f"arm must be 'p' or 'q', got {self.arm!r}")
        if not 0.0 <= self.frac_start < self.frac_end <= 1.0:
            raise ValueError("require 0 <= frac_start < frac_end <= 1")
        if self.chrom not in CHROM_LENGTHS:
            raise ValueError(f"unknown autosome {self.chrom}")


def _whole(chrom: int, copy_number: int) -> tuple[Segment, Segment]:
    """A whole chromosome change, as its two arms."""
    return (Segment(chrom, "p", 0.0, 1.0, copy_number),
            Segment(chrom, "q", 0.0, 1.0, copy_number))


#: Named lesions. Sizes deliberately span two orders of magnitude, because
#: lesion size is one of the axes the detectability study sweeps over. The
#: copy neutral entries carry no segments at all, which is the point: a
#: balanced rearrangement changes no bin's copy number.
LESIONS: dict[str, tuple[Segment, ...]] = {
    "del(5q)": (Segment(5, "q", 0.35, 1.00, 1),),
    "del(7q)": (Segment(7, "q", 0.45, 1.00, 1),),
    "del(13q)": (Segment(13, "q", 0.30, 0.65, 1),),
    "del(20q)": (Segment(20, "q", 0.30, 0.75, 1),),
    "del(20q)_small": (Segment(20, "q", 0.45, 0.60, 1),),
    "minus7": _whole(7, 1),
    "plus8": _whole(8, 3),
    "complex": (
        Segment(5, "q", 0.35, 1.00, 1),
        *_whole(7, 1),
        *_whole(8, 3),
        Segment(17, "p", 0.00, 0.90, 1),
        Segment(12, "p", 0.00, 0.60, 1),
    ),
    "t(9;22)": (),
    "t(15;17)": (),
    "t(8;21)": (),
    "inv(16)": (),
    "normal": (),
}

#: Lesions that are balanced, so invisible to any depth based method.
COPY_NEUTRAL: frozenset[str] = frozenset(
    {"t(9;22)", "t(15;17)", "t(8;21)", "inv(16)"})


class BinnedGenome:
    """Fixed width autosomal bins plus the covariates a coverage model needs.

    Parameters
    ----------
    bin_size:
        Bin width in base pairs. 500 kb is the default because it is the scale
        at which low pass callers operate: small enough to localise an arm
        level lesion, large enough that a 0.1x library still puts hundreds of
        reads in each bin.
    seed:
        Seed for the covariates. The same seed gives the same genome, which
        matters because the reference panel and the test libraries must share
        one genome for the comparison to mean anything.
    mappability_floor:
        Bins below this mappability are excluded outright, standing in for
        centromeric and other unalignable regions that every real pipeline
        drops.
    """

    def __init__(self, bin_size: int = 500_000, seed: int = 0,
                 mappability_floor: float = 0.55) -> None:
        if bin_size <= 0:
            raise ValueError("bin_size must be positive")
        self.bin_size = int(bin_size)
        self.seed = int(seed)
        self.mappability_floor = float(mappability_floor)
        rng = np.random.default_rng(seed)

        chrom, start, arm = [], [], []
        for c, length in CHROM_LENGTHS.items():
            n = length // self.bin_size
            if n == 0:
                continue
            starts = np.arange(n) * self.bin_size
            chrom.append(np.full(n, c, dtype=np.int16))
            start.append(starts)
            arm.append((starts >= CENTROMERES[c]).astype(np.int8))

        self.chrom = np.concatenate(chrom)
        self.start = np.concatenate(start)
        self.end = self.start + self.bin_size
        self.arm = np.concatenate(arm)
        self.n_bins = int(self.chrom.size)

        self.gc = self._simulate_gc(rng)
        self.mappability = self._simulate_mappability(rng)
        self.excluded = self.mappability < self.mappability_floor

        # Systematic per bin bias shared by every library on this protocol.
        # A reference panel can learn and remove it, which is the whole reason
        # reference panels exist.
        self.shared_bias = np.exp(rng.normal(0.0, 0.14, self.n_bins))

    # -- construction helpers ---------------------------------------------

    def _simulate_gc(self, rng: np.random.Generator) -> np.ndarray:
        """Large scale isochore waves plus local variation, centred on 0.41."""
        pos_mb = self.start / 1e6
        waves = (0.020 * np.sin(2 * np.pi * pos_mb / 12.0 + self.chrom)
                 + 0.012 * np.sin(2 * np.pi * pos_mb / 3.1 + 2 * self.chrom))
        return np.clip(0.41 + waves + rng.normal(0, 0.025, self.n_bins),
                       0.25, 0.65)

    def _simulate_mappability(self, rng: np.random.Generator) -> np.ndarray:
        """Mostly high, degrading towards centromeres."""
        base = np.clip(rng.beta(18.0, 1.6, self.n_bins), 0.0, 1.0)
        centromere = np.array([CENTROMERES[c] for c in self.chrom])
        dist_mb = np.abs(self.start - centromere) / 1e6
        return base * (1.0 - 0.55 * np.exp(-dist_mb / 2.5))

    # -- lookups -----------------------------------------------------------

    @property
    def usable(self) -> np.ndarray:
        """Boolean mask of bins any caller is allowed to look at."""
        return ~self.excluded

    def arm_mask(self, chrom: int, arm: str) -> np.ndarray:
        """Boolean mask selecting one chromosome arm."""
        if arm not in ("p", "q"):
            raise ValueError(f"arm must be 'p' or 'q', got {arm!r}")
        return (self.chrom == chrom) & (self.arm == (0 if arm == "p" else 1))

    def arm_label(self, index: int) -> str:
        """Human readable label for the arm a bin belongs to, such as ``5q``."""
        return f"{self.chrom[index]}{'p' if self.arm[index] == 0 else 'q'}"

    def __len__(self) -> int:
        return self.n_bins

    def __repr__(self) -> str:
        return (f"BinnedGenome(bin_size={self.bin_size}, seed={self.seed}, "
                f"n_bins={self.n_bins}, usable={int(self.usable.sum())})")


def lesion_copy_number(genome: BinnedGenome, lesion: str) -> np.ndarray:
    """Per bin copy number in the abnormal clone. Normal diploid is 2.

    Raises
    ------
    KeyError
        If the lesion name is not in :data:`LESIONS`.
    """
    try:
        segments = LESIONS[lesion]
    except KeyError as exc:
        known = ", ".join(sorted(LESIONS))
        raise KeyError(f"unknown lesion {lesion!r}; known: {known}") from exc

    cn = np.full(genome.n_bins, 2.0)
    for seg in segments:
        idx = np.flatnonzero(genome.arm_mask(seg.chrom, seg.arm))
        if idx.size == 0:
            continue
        lo = int(round(idx.size * seg.frac_start))
        hi = int(round(idx.size * seg.frac_end))
        cn[idx[lo:hi]] = seg.copy_number
    return cn


def lesion_size_mb(genome: BinnedGenome, lesion: str) -> float:
    """Total megabases whose copy number differs from diploid."""
    cn = lesion_copy_number(genome, lesion)
    return float(np.count_nonzero(cn != 2) * genome.bin_size / 1e6)
