"""The detectability sweep.

For every combination of lesion, tumour fraction and sequencing depth, estimate
the probability that the bin count detector calls it, with the genome wide
false positive rate held fixed by calibration.

Two details are deliberate and easy to get wrong. The reference panel and the
calibration set are regenerated at every depth, because a panel built at one
depth does not transfer to another. And the calibration libraries, the panel
libraries and the normal libraries used to measure the achieved false positive
rate are three disjoint sets, so the reported rate is an out of sample number
rather than a restatement of the threshold.
"""

from __future__ import annotations

import csv
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .detect import BinCountDetector, is_localised, localisation_precision
from .genome import COPY_NEUTRAL, BinnedGenome, lesion_copy_number, lesion_size_mb
from .simulate import LibraryModel, simulate_library, simulate_panel

__all__ = ["SweepConfig", "run_sweep", "write_csv"]

DEFAULT_LESIONS = ("plus8", "minus7", "del(5q)", "del(7q)", "del(13q)",
                   "del(20q)", "del(20q)_small", "complex",
                   "t(9;22)", "t(15;17)", "t(8;21)", "inv(16)")


@dataclass
class SweepConfig:
    """Everything the sweep needs, in one object so a run is reproducible."""

    bin_size: int = 500_000
    genome_seed: int = 1
    depths: tuple[float, ...] = (0.01, 0.05, 0.1, 0.5, 1.0, 5.0)
    tumour_fractions: tuple[float, ...] = (0.01, 0.02, 0.03, 0.05, 0.075,
                                           0.10, 0.15, 0.25, 0.40)
    lesions: tuple[str, ...] = DEFAULT_LESIONS
    n_panel: int = 100
    n_calibration: int = 2000
    n_fpr_check: int = 1000
    n_replicates: int = 200
    target_fpr: float = 0.05
    seed: int = 20261012
    library: LibraryModel = field(default_factory=LibraryModel)

    @property
    def n_conditions(self) -> int:
        return len(self.depths) * len(self.lesions) * len(self.tumour_fractions)

    @property
    def n_calls(self) -> int:
        return self.n_conditions * self.n_replicates


def run_sweep(config: SweepConfig | None = None, verbose: bool = True
              ) -> tuple[list[dict], list[dict]]:
    """Run the sweep. Returns ``(detection_rows, calibration_rows)``."""
    config = config or SweepConfig()
    genome = BinnedGenome(config.bin_size, seed=config.genome_seed)
    copy_numbers = {name: lesion_copy_number(genome, name)
                    for name in config.lesions}
    sizes = {name: lesion_size_mb(genome, name) for name in config.lesions}
    diploid = np.full(genome.n_bins, 2.0)

    rows: list[dict] = []
    calibration: list[dict] = []
    started = time.time()
    master = np.random.SeedSequence(config.seed)
    streams = master.spawn(len(config.depths))

    for depth, stream in zip(config.depths, streams, strict=True):
        rng = np.random.default_rng(stream)
        panel = simulate_panel(genome, rng, n=config.n_panel, depth=depth,
                               model=config.library)
        detector = BinCountDetector(genome, panel)
        detector.calibrate(simulate_panel(genome, rng, n=config.n_calibration,
                                          depth=depth, model=config.library))
        thr = detector.threshold(config.target_fpr)

        # Fresh normals, seen by neither the panel nor the calibration set.
        achieved = float(np.mean([
            detector.call(simulate_library(genome, rng, depth=depth,
                                           model=config.library,
                                           tumour_copy_number=diploid),
                          threshold=thr).detected
            for _ in range(config.n_fpr_check)]))
        calibration.append(dict(
            depth=depth, threshold=round(thr, 4),
            target_fpr=config.target_fpr, achieved_fpr=round(achieved, 4),
            n_fresh_normals=config.n_fpr_check,
            n_panel=config.n_panel, n_calibration=config.n_calibration,
            mean_reads_per_bin=round(float(panel.mean()), 1),
            n_segments=detector.n_segments,
            usable_bins=int(detector.usable.sum()),
        ))
        if verbose:
            print(f"depth {depth:>6}x  threshold {thr:6.2f}  "
                  f"achieved FPR {achieved:.3f}", flush=True)

        for lesion in config.lesions:
            cn = copy_numbers[lesion]
            for fraction in config.tumour_fractions:
                detected = localised = 0
                precisions = []
                for _ in range(config.n_replicates):
                    counts = simulate_library(
                        genome, rng, depth=depth, tumour_fraction=fraction,
                        model=config.library, tumour_copy_number=cn)
                    result = detector.call(counts, threshold=thr)
                    detected += result.detected
                    localised += is_localised(result.called_bins, cn)
                    if result.detected:
                        precisions.append(
                            localisation_precision(result.called_bins, cn))
                rows.append(dict(
                    lesion=lesion,
                    copy_neutral=lesion in COPY_NEUTRAL,
                    lesion_mb=round(sizes[lesion], 1),
                    depth=depth,
                    tumour_fraction=fraction,
                    n_replicates=config.n_replicates,
                    detected=round(detected / config.n_replicates, 4),
                    localised=round(localised / config.n_replicates, 4),
                    precision=(round(float(np.mean(precisions)), 4)
                               if precisions else ""),
                ))
            if verbose:
                print(f"  {lesion:<16} done  "
                      f"({time.time() - started:6.0f}s)", flush=True)

    if verbose:
        print(f"{len(rows)} conditions, {config.n_calls} calls, "
              f"{time.time() - started:.0f}s", flush=True)
    return rows, calibration


def write_csv(rows: list[dict], path: Path) -> None:
    """Write a list of uniform dictionaries as a CSV file."""
    if not rows:
        raise ValueError("nothing to write")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
