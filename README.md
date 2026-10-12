# low-pass-cfDNA

**What a bin count copy number caller can and cannot see when you try to read a
haematological karyotype from a blood draw.** Measured, not asserted.

[![tests](https://github.com/shubhodeepcx/low-pass-cfDNA/actions/workflows/ci.yml/badge.svg)](https://github.com/shubhodeepcx/low-pass-cfDNA/actions/workflows/ci.yml)

---

## Why this exists

The obvious project here is to build a better low pass CNV caller. That is a
bad idea. WisecondorX and ichorCNA have had years of tuning on within sample
normalisation, GC and mappability correction and reference panel construction,
and a new entrant competes against the incumbent's strongest axis.

Two things are genuinely missing from the literature, and this repository
addresses both.

**An open sensitivity map for myeloid lesions.** Published detection limits for
low pass cell free DNA are mostly derived for prenatal whole chromosome
trisomies. The lesions that matter in marrow disease, del(5q), del(7q),
del(20q), trisomy 8 and complex karyotype, span two orders of magnitude in
size, and nobody has published the detectability surface across tumour
fraction, depth and lesion size with the false positive rate held fixed.

**The copy neutral blind spot, quantified.** Balanced translocations change no
bin's depth, so no depth based method can ever call them. This repository shows
that empirically and then answers the question that follows: what would
junction level evidence actually cost?

## The three findings

### 1. Past about 0.5x, more sequencing buys almost nothing

![limit of detection against depth](figures/fig1_depth_plateau.png)

Going from 0.01x to 0.1x improves the limit of detection by roughly half.
Going from 0.5x to 5x, ten times the sequencing and ten times the cost,
improves it by a few percent. Past the half coverage mark the detector is
limited by library to library bias rather than counting noise, so a lab running
1x or 5x cell free DNA for cytogenetics is paying for coverage it cannot use.

### 2. Small lesions are out of reach

![power curves](figures/fig2_power_curves.png)

At a realistic 0.1x, a 434 Mb complex karyotype is detectable at about 7 percent
tumour fraction. A 5 Mb del(20q) needs nearly 40 percent, which is a tumour
fraction at which nobody needed the test. Detection scales with roughly the
square root of lesion size, so halving the detectable tumour fraction takes a
lesion four times larger.

### 3. Balanced translocations are detected at exactly the false positive rate

![junction cost](figures/fig3_junction_cost.png)

For a heterozygous balanced rearrangement at tumour fraction `f`, the expected
number of library molecules spanning the junction is

```
E = physical_coverage * f / 2
```

Genome length cancels. Fragment length cancels. At 0.1x and 10 percent tumour
fraction that is **0.004 reads**, four thousandths of one read, not a
probability. Three junction reads at 90 percent power needs about 140x
coverage, roughly 434 Gbp per sample, against 0.21 Gbp from a 1.5 Mb capture
panel tiling the recurrent breakpoint regions.

The conclusion is unwelcome but clean: **the copy neutral problem is a coverage
allocation gap, not an algorithm gap.** No amount of fragmentomics or machine
learning retrieves information from four thousandths of a read. The fix is to
point the sequencing at the junctions.

Full write up with every number: [docs/RESULTS.md](docs/RESULTS.md). The paper
in IEEE format: [paper/main.pdf](paper/main.pdf).

## Install and run

```bash
git clone https://github.com/shubhodeepcx/low-pass-cfDNA
cd low-pass-cfDNA
pip install -e ".[dev]"

pytest -q                # the test suite
lpcfdna sweep --quick    # a few minutes, indicative numbers
lpcfdna sweep            # the real thing, about five minutes
lpcfdna figures          # charts from the CSVs
lpcfdna junction         # the copy neutral evidence calculation
lpcfdna summary          # headline numbers, printed
```

No reference genome, no BAM files, no GPU. The whole study runs on a laptop.

## How it works

### The coverage model

Reads are simulated per 500 kb bin across the 22 autosomes. The expected count
in a bin is the product of depth, a shared systematic bias, a per library bias,
a GC efficiency curve, mappability and the copy number ratio implied by the
lesion and the tumour fraction. Counts are Poisson.

The split between the bias terms is the part that matters:

| Term | Behaviour across libraries | Removable? |
|---|---|---|
| shared systematic bias | identical in every library | yes, by a reference panel |
| per library bias | redrawn every library | **no, this sets the floor** |
| GC efficiency curve | same shape, per library strength | partly |
| mappability | fixed property of the reference | yes, by masking |
| Poisson counting noise | shrinks with depth | yes, with depth |

That third row is why detection limits plateau instead of improving forever.
Any simulation without a per library bias term reports limits that keep getting
better with coverage, which is flattering and wrong.

GC content and mappability are modelled rather than read from a reference
FASTA. A production pipeline substitutes measured per bin values at exactly two
attribute assignments and nothing else changes.

### The detector

Within library normalisation, binned median GC correction, per bin z scores
against a 100 library reference panel, then a Stouffer style statistic over
about 2,500 candidate segments: every chromosome arm plus sliding windows of 8
to 256 bins.

**The threshold is calibrated, not chosen.** It is the 95th percentile of the
maximum segment statistic across 1,000 held out normal libraries, which
controls the genome wide false positive rate at 5 percent. This is the single
most important methodological point in the repository. A per segment threshold
applied across 2,500 segments calls something in nearly every normal library,
and a sensitivity figure reported at uncontrolled specificity means nothing.
`results/calibration.csv` records the achieved rate at every depth so the claim
can be checked rather than believed.

One incidental finding worth repeating: calibrating a genome wide threshold at
the 5 percent level takes on the order of 1,000 normal libraries. With 300 the
quantile estimate is noisy enough to deliver 8 percent instead of 5. For a
small lab without a large normal cohort that is a real obstacle, and it is
separate from the better known requirement that the panel match the assay
chemistry.

## Layout

```
src/lpcfdna/
  genome.py     binned genome, covariates, lesion definitions
  simulate.py   the coverage model
  detect.py     the bin count detector and its calibration
  junction.py   closed form junction evidence calculation
  sweep.py      the detectability sweep
  analysis.py   limits of detection and summary statistics
  figures.py    the three charts
  cli.py        command line interface
tests/          80 tests, property based where it matters
results/        sweep.csv, calibration.csv, limit_of_detection.csv
figures/        three PNG figures
docs/           RESULTS.md, METHODS.md, LIMITATIONS.md
paper/          IEEE format manuscript, LaTeX source and build
```

## What this is not

Not a clinical tool, and not a claim about any real assay. It is a simulation
study whose purpose is to bound what is possible before anyone spends money on
sequencing. The numbers are only as good as the bias model, and the bias model
is a reasoned guess rather than a measurement.

The honest reading of every absolute figure here is "order of magnitude, under
a stated model". The qualitative findings, the plateau, the size dependence and
the copy neutral floor, do not depend on the bias parameter. Where the plateau
sits does. [docs/LIMITATIONS.md](docs/LIMITATIONS.md) is the full list, written
before anyone had to ask.

## What would make these numbers real

1. Fit the bias model to real low pass cell free DNA libraries. Highest value
   single step, and it turns every number here from an estimate into a
   prediction.
2. Run the actual WisecondorX and ichorCNA binaries on the simulated libraries.
   The detector here is of the right class but is not either tool.
3. Replace modelled GC and mappability with measured per bin values from
   GRCh38.
4. Validate against cases with a known karyotype and a measured tumour
   fraction.

## Licence

MIT. See [LICENSE](LICENSE).
