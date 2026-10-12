"""What junction level evidence costs, for the lesions depth cannot see.

A balanced rearrangement changes no bin's copy number, so the power of any
depth based caller equals its false positive rate exactly. That is not a
limitation of one implementation, it is a property of the whole class. The
question worth answering is therefore not how to improve the caller but how
much evidence the alternative needs.

For a heterozygous balanced rearrangement at tumour fraction ``f``, the
expected number of library molecules spanning the junction is

.. math::

    E = c_\\text{phys} \\cdot \\frac{f}{2}

with ``c_phys`` the physical, that is fragment level, coverage. Genome length
cancels, and so does fragment length. One of the two homologues carries the
derivative chromosome, which is where the factor of two comes from.

For cell free DNA with 150 base paired end reads on roughly 167 base
fragments, the read pair covers the whole fragment, so split reads are the only
junction evidence on offer. Discordant pair evidence needs an insert longer
than the read pair, which cell free DNA does not provide.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import poisson

from .genome import GENOME_HAPLOID_BP

__all__ = ["JunctionModel", "CapturePanel", "expected_junction_reads",
           "detection_power", "required_coverage", "gigabases",
           "power_table"]

#: Modal cell free DNA fragment length in bases, the mononucleosomal peak.
CFDNA_FRAGMENT_BP: float = 167.0


@dataclass(frozen=True)
class JunctionModel:
    """Assumptions behind the junction evidence calculation.

    Attributes
    ----------
    fragment_bp:
        Modal fragment length.
    anchor_bp:
        Unique sequence required on each side of the breakpoint for a split
        read to be mappable with confidence.
    min_reads:
        Supporting reads required before a junction is believed. Two is
        optimistic in the presence of mapping artefacts; three is the usual
        working minimum.
    allele_copies:
        Derivative chromosome copies per diploid genome in the clone. One, for
        a heterozygous balanced rearrangement.
    """

    fragment_bp: float = CFDNA_FRAGMENT_BP
    anchor_bp: float = 20.0
    min_reads: int = 3
    allele_copies: float = 1.0

    def __post_init__(self) -> None:
        if self.fragment_bp <= 2 * self.anchor_bp:
            raise ValueError("fragment is too short to carry two anchors")
        if self.min_reads < 1:
            raise ValueError("min_reads must be at least 1")

    @property
    def usable_fraction(self) -> float:
        """Share of spanning fragments that yield a usable split read."""
        return (self.fragment_bp - 2 * self.anchor_bp) / self.fragment_bp


@dataclass(frozen=True)
class CapturePanel:
    """A hybrid capture panel tiling known breakpoint regions.

    The recurrent haematological fusions have breakpoints confined to known
    introns, from a few kilobases up to around 140 kilobases each. Tiling them
    is a small target, and coverage spent there is coverage not wasted on the
    99.95 percent of the genome that carries no junction.
    """

    n_fusions: int = 10
    bp_per_region: float = 150_000.0
    name: str = "recurrent haematological fusions"

    @property
    def target_bp(self) -> float:
        return self.n_fusions * self.bp_per_region

    @property
    def genome_fraction(self) -> float:
        return self.target_bp / GENOME_HAPLOID_BP

    @property
    def enrichment(self) -> float:
        """How much more coverage a fixed budget buys inside the target."""
        return GENOME_HAPLOID_BP / self.target_bp


def expected_junction_reads(coverage: float, tumour_fraction: float,
                            model: JunctionModel | None = None) -> float:
    """Expected usable junction spanning reads for one rearrangement."""
    model = model or JunctionModel()
    if coverage < 0:
        raise ValueError("coverage must be non negative")
    if not 0.0 <= tumour_fraction <= 1.0:
        raise ValueError("tumour_fraction must lie in [0, 1]")
    return (coverage * tumour_fraction * model.allele_copies / 2.0
            * model.usable_fraction)


def detection_power(coverage: float, tumour_fraction: float,
                    model: JunctionModel | None = None) -> float:
    """Probability of seeing at least ``min_reads`` junction reads.

    Poisson, because junction spanning fragments are a thinned sample of an
    already random fragmentation process. The thinning probability is tiny and
    the fragment count is huge, which is exactly the Poisson limit.
    """
    model = model or JunctionModel()
    lam = expected_junction_reads(coverage, tumour_fraction, model)
    return float(poisson.sf(model.min_reads - 1, lam))


def required_coverage(tumour_fraction: float, power: float = 0.90,
                      model: JunctionModel | None = None,
                      max_coverage: float = 1e7) -> float:
    """Physical coverage needed to reach the requested power.

    Solved by bisection on a log scale. Power is monotone in coverage, so the
    bracket is safe.
    """
    model = model or JunctionModel()
    if not 0.0 < power < 1.0:
        raise ValueError("power must lie in (0, 1)")
    if tumour_fraction <= 0:
        return float("inf")
    if detection_power(max_coverage, tumour_fraction, model) < power:
        return float("inf")

    lo, hi = 1e-6, max_coverage
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if detection_power(mid, tumour_fraction, model) < power:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi))


def gigabases(coverage: float, target_bp: float = GENOME_HAPLOID_BP) -> float:
    """Sequenced gigabases implied by a coverage over a target of this size."""
    return coverage * target_bp / 1e9


def power_table(tumour_fractions, power: float = 0.90,
                model: JunctionModel | None = None,
                panel: CapturePanel | None = None) -> list[dict]:
    """Rows of coverage and sequencing cost per tumour fraction."""
    model = model or JunctionModel()
    panel = panel or CapturePanel()
    rows = []
    for f in tumour_fractions:
        cov = required_coverage(f, power, model)
        rows.append(dict(
            tumour_fraction=float(f),
            required_coverage=cov,
            wgs_gbp=gigabases(cov),
            capture_gbp=gigabases(cov, panel.target_bp),
            reads_at_0_1x=expected_junction_reads(0.1, f, model),
        ))
    return rows


def report(model: JunctionModel | None = None,
           panel: CapturePanel | None = None) -> str:
    """A plain text summary, used by the command line interface."""
    model = model or JunctionModel()
    panel = panel or CapturePanel()
    fractions = [0.01, 0.02, 0.05, 0.10, 0.25, 0.50]
    rows = power_table(fractions, model=model, panel=panel)

    out = ["COPY NEUTRAL DETECTION: WHAT THE EVIDENCE COSTS", "",
           f"fragment length            {model.fragment_bp:.0f} bp",
           f"anchor each side           {model.anchor_bp:.0f} bp",
           f"usable spanning fragments  {model.usable_fraction:.2f}",
           f"reads required             {model.min_reads}", "",
           "expected junction reads = coverage * tumour_fraction / 2 * "
           f"{model.usable_fraction:.2f}", "",
           f"{'fraction':>9} | {'at 0.1x':>10} | {'coverage':>10} | "
           f"{'WGS Gbp':>9} | {'capture Gbp':>12}",
           "-" * 62]
    for r in rows:
        out.append(f"{r['tumour_fraction']:>9.2f} | {r['reads_at_0_1x']:>10.4f} | "
                   f"{r['required_coverage']:>9.0f}x | {r['wgs_gbp']:>9.0f} | "
                   f"{r['capture_gbp']:>12.2f}")
    out += ["",
            f"panel target               {panel.target_bp/1e6:.1f} Mb "
            f"({panel.genome_fraction*100:.3f}% of the genome)",
            f"enrichment                 {panel.enrichment:.0f}x"]
    return "\n".join(out)
