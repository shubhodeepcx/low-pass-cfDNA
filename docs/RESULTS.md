# Results

Three findings, ordered by how much they should change what someone does.

Everything below comes from `results/sweep.csv`: 540 conditions, 200 simulated
libraries each, 108,000 detector calls, with the genome wide false positive
rate held at 5 percent by calibration against 1,000 held out normal libraries
at every depth. Reproduce with `lpcfdna sweep && lpcfdna summary`.

## 0. The false positive rate really is controlled

Before any sensitivity claim means anything, the specificity has to be pinned
down and shown.

| Depth | Threshold | Target FPR | Achieved FPR |
|---|---|---|---|
| 0.01x | 4.87 | 5% | 7.4% |
| 0.05x | 5.07 | 5% | 5.4% |
| 0.1x | 5.73 | 5% | 5.8% |
| 0.5x | 6.85 | 5% | 3.6% |
| 1x | 6.88 | 5% | 3.2% |
| 5x | 7.33 | 5% | 4.8% |

Achieved rates are measured on 500 normal libraries that neither the panel nor
the calibration set has seen. They average 5.0 percent and scatter from 3.2 to
7.4 percent, which is quantile estimation noise rather than bias: the standard
error on a 5 percent rate from 500 draws is about 1 percentage point, and the
threshold itself is a 95th percentile estimated from 1,000 samples.

That scatter is itself worth reporting. **Calibrating a genome wide threshold
to 5 percent is harder than it looks.** With 300 calibration libraries instead
of 1,000, the achieved rate lands near 8 percent, a 60 percent inflation from
calibration noise alone with every other part of the method correct. For a lab
without a large normal cohort this is a real obstacle, separate from the better
known requirement that the panel match the assay chemistry.

## 1. Past about 0.5x, more sequencing buys almost nothing

![limit of detection against depth](../figures/fig1_depth_plateau.png)

Smallest tumour fraction detected with 90 percent probability:

| Lesion | 0.01x | 0.05x | 0.1x | 0.5x | 1x | 5x |
|---|---|---|---|---|---|---|
| del(20q), 5.5 Mb | not reached | not reached | 39.5% | 37.8% | 37.2% | 37.0% |
| del(20q), 16 Mb | not reached | 35.3% | 27.0% | 23.5% | 23.3% | 23.4% |
| del(13q), 33.5 Mb | 38.3% | 23.0% | 21.1% | 14.3% | 13.9% | 13.9% |
| del(7q), 54 Mb | 35.9% | 22.1% | 16.3% | 13.7% | 13.4% | 12.8% |
| del(5q), 86 Mb | 24.9% | 14.0% | 13.6% | 9.7% | 9.6% | 9.4% |
| trisomy 8, 145 Mb | 23.7% | 13.1% | 10.0% | 8.1% | 7.3% | 7.4% |
| monosomy 7, 159 Mb | 23.2% | 13.2% | 10.5% | 9.2% | 9.0% | 8.1% |
| complex, 434 Mb | 15.9% | 8.9% | 7.2% | 6.8% | 5.3% | 5.0% |

"Not reached" means the lesion was never detected at 90 percent power at any
tumour fraction up to 40 percent. It is a result, not a missing value.

Going from 0.01x to 0.1x improves the detection limit by 45 to 58 percent,
median 55. Going from 0.5x to 5x, **ten times the sequencing and ten times the
cost**, improves it by 0.3 to 27 percent, median 4.9. Six of the eight lesions
gain under 9 percent. The one large gain, 27 percent for the complex karyotype,
is the easiest case, where the limit is already so low that the tumour fraction
grid is coarse relative to the remaining headroom.

This is the per library bias term doing its work. At 0.01x the detector is
counting noise limited and depth helps enormously. By 0.5x the counting noise
has fallen below the library to library bias, and depth stops mattering. **A
lab running 1x or 5x cell free DNA for cytogenetics is paying for coverage it
cannot use.** The route past the plateau is a better reference panel or a
different kind of evidence, not more reads.

The caveat that matters, stated once and in full: the position of the plateau
is set by the per library bias parameter, which here is a reasoned guess rather
than a measurement. The *existence* of a plateau does not depend on that
parameter, because any non zero per library bias produces one. Where it sits
does. See [LIMITATIONS.md](LIMITATIONS.md).

## 2. Small lesions are out of reach, and small lesions are most of the list

![power curves at 0.1x](../figures/fig2_power_curves.png)

At a realistic 0.1x, the 434 Mb complex karyotype is detectable at 7.2 percent
tumour fraction. The 5.5 Mb del(20q) needs 39.5 percent, a tumour fraction at
which the diagnosis was never in doubt.

Fitting the detection limit against lesion size on log axes gives a slope of
-0.40 at 0.1x and -0.46 at 5x, close to the -0.5 the segment statistic implies:
**detection scales with the inverse square root of lesion size.** Halving the
detectable tumour fraction takes a lesion four times larger. That is not a
tuning problem, it is what averaging over bins does.

Localisation is not the bottleneck. In conditions at 0.1x where detection
reached 90 percent, the calls covered at least half the true lesion 94 percent
of the time. When this detector fires, it fires in the right place. It simply
does not fire often enough.

## 3. Balanced translocations are detected at exactly the false positive rate

Across all 108 copy neutral conditions, four rearrangements by six depths by
nine tumour fractions up to 40 percent, mean detection was **5.3 percent**
against a 5 percent false positive rate. The highest single value was 11.5
percent. The trend with tumour fraction is -0.007 per unit fraction, which is
zero to within noise.

There is no trend because there is nothing to detect. The simulator confirms
the premise directly: a balanced rearrangement at 90 percent tumour fraction
produces a library byte identical to a normal one at the same seed
(`test_balanced_rearrangement_is_indistinguishable_from_normal`).

This is not a weakness of this implementation. It is a property of the entire
class of depth based methods, and it is why t(9;22), t(15;17), t(8;21) and
inv(16), the lesions that define the diseases, are invisible to WisecondorX,
ichorCNA and every successor built on bin counts.

### What junction level evidence would cost

![junction cost](../figures/fig3_junction_cost.png)

For a heterozygous balanced rearrangement at tumour fraction `f`, the expected
number of library molecules spanning the junction is

```
E = physical_coverage * f / 2
```

Genome length cancels. Fragment length cancels. The derivation is in
[METHODS.md](METHODS.md) section 8. With a 20 base anchor required each side on
a 167 base cell free DNA fragment, 76 percent of spanning molecules yield a
usable split read, so expected usable evidence is `0.38 * c * f`.

At 0.1x and 10 percent tumour fraction that is **0.004 reads**. Not a 0.004
probability of detection. Four thousandths of one read.

Detection needing three junction reads at 90 percent power requires:

| Tumour fraction | Physical coverage | Untargeted WGS | 1.5 Mb capture panel |
|---|---|---|---|
| 1% | 1,400x | 4,339 Gbp | 2.10 Gbp |
| 2% | 700x | 2,170 Gbp | 1.05 Gbp |
| 5% | 280x | 868 Gbp | 0.42 Gbp |
| 10% | 140x | 434 Gbp | 0.21 Gbp |
| 25% | 56x | 174 Gbp | 0.08 Gbp |
| 50% | 28x | 87 Gbp | 0.04 Gbp |

A routine low pass cell free DNA run delivers something between 1 and 15 Gbp.
Untargeted whole genome sequencing is off by two to three orders of magnitude.
A capture panel over the recurrent breakpoint regions is comfortably inside
budget, because those regions total about 1.5 Mb, 0.05 percent of the genome,
for an enrichment of roughly 2,000 fold.

**The copy neutral problem is a coverage allocation gap, not an algorithm
gap.** No amount of fragmentomics, machine learning or clever normalisation
retrieves information from four thousandths of a read. The fix is to point the
sequencing at the junctions.

## What this does to the research direction

The speculative half of this project, detecting balanced rearrangements from
untargeted low pass cell free DNA, is closed by the arithmetic above. That is
worth more than another month of trying.

What is left open, and narrower for it:

1. **Minimal breakpoint panel design.** Which bases must a capture panel tile
   to cover the recurrent haematological fusions, and how few can it be? This
   is a pure optimisation problem on public breakpoint coordinates and needs no
   new data at all.
2. **Junction calling at single digit read support.** Three reads is the
   working minimum because one or two are indistinguishable from mapping
   artefacts. Lowering that threshold defensibly would cut the required
   coverage proportionally.
3. **Aggregate fragmentomic signal.** Whether fragment length and nucleosome
   positioning can flag *that* an abnormal clone is present without saying
   *which* rearrangement it carries. A different question from the one closed
   above, and still open.

## What would make these numbers real

1. Fit the bias model to real low pass cell free DNA libraries. Highest value
   single step, and the one that converts every figure here from an estimate
   into a prediction.
2. Run the actual WisecondorX and ichorCNA binaries on these simulated
   libraries. The detector here is of the right class but is not either tool,
   and this is the comparison a reviewer asks for first.
3. Replace modelled GC and mappability with measured per bin values from
   GRCh38.
4. Validate against cases with a known karyotype and a measured tumour
   fraction.
