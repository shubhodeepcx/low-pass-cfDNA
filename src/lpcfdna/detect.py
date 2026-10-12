"""A bin count copy number detector of the WisecondorX and ichorCNA family.

This is a faithful reimplementation of the *approach*, not of either tool:
within library normalisation, GC correction, per bin z scores against a
reference panel, then a segment level statistic over chromosome arms and
sliding windows.

It exists so that the detectability study measures a method of the right class
with its false positive rate pinned down, rather than asserting a sensitivity
figure at unknown specificity. Threshold calibration is the load bearing idea
here. Roughly two and a half thousand candidate segments are tested per
library, so a per segment threshold would call something in nearly every normal
sample. The threshold used is instead the upper quantile of the *maximum*
segment statistic across held out normal libraries, which controls the genome
wide false positive rate directly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .genome import BinnedGenome

__all__ = ["gc_correct", "candidate_segments", "CallResult",
           "BinCountDetector", "localisation_overlap", "is_localised"]

DEFAULT_WINDOWS: tuple[int, ...] = (8, 16, 32, 64, 128, 256)


def gc_correct(counts: np.ndarray, gc: np.ndarray, usable: np.ndarray,
               n_knots: int = 24, min_bins_per_knot: int = 12) -> np.ndarray:
    """Normalise a library to its own total, then divide out its GC curve.

    The GC curve is estimated as a median within GC quantile bins and then
    interpolated, which is a robust stand in for the LOESS fit the published
    tools use. Medians rather than means matter here: a real copy number
    lesion must not drag the correction towards itself.
    """
    out = np.asarray(counts, dtype=float).copy()
    if out.shape != usable.shape:
        raise ValueError("counts and usable must have the same shape")

    total = out[usable].sum()
    if total <= 0:
        return np.zeros_like(out)
    out = out / total * usable.sum()

    edges = np.quantile(gc[usable], np.linspace(0.0, 1.0, n_knots + 1))
    centres, medians = [], []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        sel = usable & (gc >= lo) & (gc <= hi)
        if np.count_nonzero(sel) >= min_bins_per_knot:
            centres.append(0.5 * (lo + hi))
            medians.append(float(np.median(out[sel])))
    if len(centres) < 4:
        return out

    expected = np.interp(gc, np.asarray(centres), np.asarray(medians))
    overall = float(np.median(medians))
    with np.errstate(invalid="ignore", divide="ignore"):
        corrected = np.where(expected > 0, out / expected * overall, 0.0)
    corrected[~np.isfinite(corrected)] = 0.0
    return corrected


def candidate_segments(genome: BinnedGenome,
                       windows: tuple[int, ...] = DEFAULT_WINDOWS
                       ) -> np.ndarray:
    """Half open ``(lo, hi)`` bin ranges to test: whole arms and sliding windows.

    Returned as an ``(n, 2)`` integer array so segment statistics can be
    computed with prefix sums rather than a Python loop, which is what makes a
    hundred thousand library calls finish in minutes rather than hours.
    """
    segments: list[tuple[int, int]] = []
    for chrom in np.unique(genome.chrom):
        for arm in (0, 1):
            idx = np.flatnonzero((genome.chrom == chrom) & (genome.arm == arm))
            if idx.size < 8:
                continue
            lo0, hi0 = int(idx[0]), int(idx[-1]) + 1
            segments.append((lo0, hi0))
            for width in windows:
                if width >= idx.size:
                    continue
                step = max(1, width // 2)
                for offset in range(0, idx.size - width + 1, step):
                    segments.append((lo0 + offset, lo0 + offset + width))
    return np.asarray(segments, dtype=np.int64)


@dataclass(frozen=True)
class CallResult:
    """What the detector concluded about one library."""

    detected: bool
    max_statistic: float
    threshold: float
    called_bins: np.ndarray
    n_called_segments: int

    @property
    def called_mb(self) -> int:
        """Number of bins called, which the caller can scale to megabases."""
        return int(np.count_nonzero(self.called_bins))


class BinCountDetector:
    """Fit on a panel of normal libraries, calibrate on held out normals.

    Parameters
    ----------
    genome:
        The same :class:`~lpcfdna.genome.BinnedGenome` the libraries came from.
    panel:
        Normal libraries used to learn per bin means and standard deviations.
    windows:
        Sliding window widths in bins. Whole arms are always included.
    min_segment_bins:
        Segments with fewer usable bins than this are not tested at all.
    """

    def __init__(self, genome: BinnedGenome, panel: np.ndarray,
                 windows: tuple[int, ...] = DEFAULT_WINDOWS,
                 min_segment_bins: int = 8,
                 min_usable_fraction: float = 0.5) -> None:
        panel = np.atleast_2d(np.asarray(panel, dtype=float))
        if panel.shape[1] != genome.n_bins:
            raise ValueError("panel columns do not match the genome's bins")
        if panel.shape[0] < 3:
            raise ValueError("need at least 3 panel libraries to estimate a sd")

        self.genome = genome
        self.min_segment_bins = int(min_segment_bins)
        self.n_panel = int(panel.shape[0])

        corrected = np.vstack([gc_correct(row, genome.gc, genome.usable)
                               for row in panel])
        self.panel_mean = corrected.mean(axis=0)
        self.panel_sd = corrected.std(axis=0, ddof=1)

        # Drop bins the panel itself cannot characterise. A bin with zero
        # variance across normals is not informative, it is broken.
        unusable = (self.panel_mean <= 0) | (self.panel_sd <= 0)
        self.usable = genome.usable & ~unusable
        if self.usable.mean() < min_usable_fraction:
            raise RuntimeError(
                f"only {self.usable.mean():.1%} of bins are usable; "
                "check the genome or the panel depth")

        self.segments = candidate_segments(genome, windows)
        self._null_max = self._leave_one_out_null(corrected)
        self._calibrated = False

    # -- internals ---------------------------------------------------------

    def _zscores(self, corrected: np.ndarray, mean: np.ndarray | None = None,
                 sd: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        mean = self.panel_mean if mean is None else mean
        sd = self.panel_sd if sd is None else sd
        ok = self.usable & (sd > 0)
        z = np.zeros(self.genome.n_bins)
        z[ok] = (corrected[ok] - mean[ok]) / sd[ok]
        return z, ok

    def _segment_statistics(self, z: np.ndarray, ok: np.ndarray) -> np.ndarray:
        """Mean z over each segment's usable bins, scaled by the square root of
        their count. Computed with prefix sums, so cost is independent of how
        many segments overlap."""
        cum_z = np.concatenate(([0.0], np.cumsum(np.where(ok, z, 0.0))))
        cum_n = np.concatenate(([0], np.cumsum(ok.astype(np.int64))))
        lo, hi = self.segments[:, 0], self.segments[:, 1]
        n = cum_n[hi] - cum_n[lo]
        total = cum_z[hi] - cum_z[lo]
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(n >= self.min_segment_bins,
                            total / np.sqrt(np.maximum(n, 1)), np.nan)

    def _max_statistic(self, corrected: np.ndarray, mean=None, sd=None) -> float:
        z, ok = self._zscores(corrected, mean, sd)
        stats = self._segment_statistics(z, ok)
        return float(np.nanmax(np.abs(stats))) if np.isfinite(stats).any() else np.nan

    def _leave_one_out_null(self, corrected: np.ndarray) -> np.ndarray:
        """Provisional null from the fitting panel, used until calibration.

        Slightly optimistic, because each library is scored against a panel it
        helped build. :meth:`calibrate` replaces it, and
        :attr:`calibrated` says which one is in force.
        """
        values = np.empty(corrected.shape[0])
        for i in range(corrected.shape[0]):
            others = np.delete(corrected, i, axis=0)
            values[i] = self._max_statistic(corrected[i], others.mean(axis=0),
                                            others.std(axis=0, ddof=1))
        return values

    # -- public API --------------------------------------------------------

    @property
    def calibrated(self) -> bool:
        """True once :meth:`calibrate` has supplied a held out null."""
        return self._calibrated

    @property
    def n_segments(self) -> int:
        return int(self.segments.shape[0])

    def calibrate(self, normal_libraries: np.ndarray) -> None:
        """Set the null distribution from libraries the panel has never seen.

        These must be normal and must not appear in the fitting panel. A
        thousand of them is not excessive: at three hundred the 95th percentile
        is noisy enough to inflate the realised false positive rate by half
        again.
        """
        normal_libraries = np.atleast_2d(np.asarray(normal_libraries,
                                                    dtype=float))
        if normal_libraries.shape[1] != self.genome.n_bins:
            raise ValueError("calibration libraries do not match the genome")
        values = np.empty(normal_libraries.shape[0])
        for i, row in enumerate(normal_libraries):
            values[i] = self._max_statistic(
                gc_correct(row, self.genome.gc, self.genome.usable))
        self._null_max = values
        self._calibrated = True

    def threshold(self, false_positive_rate: float = 0.05) -> float:
        """Genome wide threshold at the requested false positive rate."""
        if not 0.0 < false_positive_rate < 1.0:
            raise ValueError("false_positive_rate must lie in (0, 1)")
        return float(np.quantile(self._null_max, 1.0 - false_positive_rate))

    def call(self, counts: np.ndarray, false_positive_rate: float = 0.05,
             threshold: float | None = None) -> CallResult:
        """Score one library and decide whether anything is called."""
        corrected = gc_correct(counts, self.genome.gc, self.genome.usable)
        z, ok = self._zscores(corrected)
        stats = self._segment_statistics(z, ok)
        thr = (self.threshold(false_positive_rate) if threshold is None
               else float(threshold))

        hit = np.abs(stats) > thr
        called = np.zeros(self.genome.n_bins, dtype=bool)
        for lo, hi in self.segments[hit]:
            called[lo:hi] = True

        finite = np.isfinite(stats)
        return CallResult(
            detected=bool(hit.any()),
            max_statistic=float(np.nanmax(np.abs(stats))) if finite.any() else float("nan"),
            threshold=thr,
            called_bins=called,
            n_called_segments=int(np.count_nonzero(hit)),
        )


def localisation_overlap(called_bins: np.ndarray,
                         tumour_copy_number: np.ndarray) -> float:
    """Fraction of the true lesion that the calls actually cover."""
    truth = tumour_copy_number != 2
    if not truth.any():
        return 0.0
    return float(np.count_nonzero(called_bins & truth) / np.count_nonzero(truth))


def is_localised(called_bins: np.ndarray, tumour_copy_number: np.ndarray,
                 min_overlap: float = 0.5) -> bool:
    """Whether the calls recovered at least ``min_overlap`` of the lesion."""
    return localisation_overlap(called_bins, tumour_copy_number) >= min_overlap
