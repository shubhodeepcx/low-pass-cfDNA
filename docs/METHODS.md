# Methods

Enough detail to reimplement the study from scratch, and enough honesty about
the choices to argue with it.

## 1. The binned genome

The 22 autosomes are divided into fixed 500 kb bins, giving 5,738 bins in
total. 500 kb is the scale low pass callers actually work at: small enough to
localise an arm level lesion, large enough that even a 0.1x library deposits
around 200 reads per bin.

Each bin carries three covariates.

**GC content** is modelled as two superposed sinusoids standing in for isochore
structure, plus local Gaussian variation, centred on 0.41 and clipped to
[0.25, 0.65]. The amplitudes are chosen so the spread of per bin GC matches
what a real 500 kb binning of GRCh38 produces.

**Mappability** is drawn from a Beta(18, 1.6) distribution, which is mostly
near one with a tail of poor bins, then multiplied by a factor that decays
towards each centromere. Bins below 0.55 are excluded outright, which stands in
for the centromeric and repeat regions every real pipeline drops.

**Shared systematic bias** is a log normal multiplier with standard deviation
0.14, identical in every library simulated on the same genome. This is the
reproducible part of coverage bias that a reference panel exists to remove.

All three are modelled rather than read from a reference FASTA, which is what
keeps the package free of large data dependencies. A production pipeline
substitutes measured per bin values at `BinnedGenome.gc` and
`BinnedGenome.mappability` and nothing downstream changes.

## 2. The coverage model

For bin `i` in a library at depth `d`, the expected read count is

```
E[count_i] = (d * bin_size / read_length)
             * shared_bias_i
             * per_library_bias_i
             * gc_efficiency(gc_i, strength)
             * mappability_i
             * ratio_i
```

and the observed count is Poisson about that expectation.

`gc_efficiency` is a Gaussian in GC peaking at 0.42, with a per library
`strength` drawn from N(6.0, 1.2). Drawing the strength per library rather than
fixing it reflects the real situation: the shape of the GC bias is a protocol
property, its sharpness varies between libraries on the same protocol.

`per_library_bias_i` is log normal with standard deviation 0.045, **redrawn for
every library**. This is the single most consequential parameter in the study.
It is the part of coverage bias that a reference panel cannot learn, because it
is different in every sample, and it is therefore what stops detection limits
improving indefinitely with depth.

`ratio_i` encodes the lesion and the tumour fraction:

```
ratio_i = (f * copy_number_i + (1 - f) * 2) / 2
```

A single copy loss at tumour fraction `f` depresses expected depth by `f / 2`,
not by `f`. That factor of two is the easiest thing in this whole study to get
wrong, and it is checked directly in `tests/test_simulate.py` by measuring the
chromosome specific bias offset on a diploid library and then asserting the
mixture term against it.

## 3. Lesions

Nine copy number lesions spanning 5 Mb to 434 Mb, plus four balanced
rearrangements carrying no copy number change at all.

| Name | Size | What it stands for |
|---|---|---|
| del(20q) small | 5.5 Mb | the small end of what gets reported |
| del(20q) | 16 Mb | a common myeloid deletion |
| del(13q) | 33.5 Mb | mid sized |
| del(7q) | 54.5 Mb | adverse risk in myeloid disease |
| del(5q) | 86.5 Mb | the defining lesion of 5q minus syndrome |
| trisomy 8 | 145 Mb | the commonest myeloid trisomy |
| monosomy 7 | 159 Mb | whole chromosome loss |
| complex | 434 Mb | five simultaneous lesions |
| t(9;22), t(15;17), t(8;21), inv(16) | zero | copy neutral by construction |

The copy neutral entries are not a formality. They are the control that proves
the detector is honest: anything other than detection at the false positive
rate would mean the simulation or the caller is leaking information that does
not exist.

## 4. The detector

A reimplementation of the WisecondorX and ichorCNA *approach*, not of either
tool.

1. **Within library normalisation.** Divide each bin by the library total.
2. **GC correction.** Split usable bins into 24 GC quantile groups, take the
   median of each, interpolate between group centres, divide out. Medians
   rather than means, so a real lesion cannot drag the correction towards
   itself.
3. **Reference panel z scores.** Per bin mean and standard deviation across a
   panel of 100 normal libraries. Bins with zero panel variance are dropped as
   broken rather than treated as informative.
4. **Segment statistics.** For each candidate segment, the mean z over its
   usable bins scaled by the square root of their count, which is the Stouffer
   combination under independence. Candidate segments are every chromosome arm
   plus sliding windows of 8, 16, 32, 64, 128 and 256 bins at half window
   steps, around 2,500 in total.
5. **Call.** Any segment whose absolute statistic exceeds the threshold.

Segment statistics are computed with prefix sums, so the cost does not grow
with how many segments overlap. That is what makes a hundred thousand library
calls finish in minutes.

## 5. Threshold calibration

This is the load bearing step.

With around 2,500 candidate segments tested per library, a threshold set per
segment would call something in nearly every normal sample. The threshold here
is instead the upper quantile of the **maximum** segment statistic across a set
of held out normal libraries, which controls the genome wide false positive
rate directly.

Three sets of libraries are kept disjoint:

- **Panel**, 100 libraries, used to learn per bin means and standard
  deviations.
- **Calibration**, 1,000 libraries, used to build the null distribution of the
  maximum statistic.
- **False positive check**, 500 libraries, used to measure the achieved rate.

Only the third gives an out of sample number, and it is the one reported in
`results/calibration.csv`.

The class also carries a leave one out null built from the panel, used before
`calibrate` is called. It is slightly optimistic, because each library is
scored against a panel it helped build, and `BinCountDetector.calibrated`
reports which null is in force so the distinction is never silent.

**On panel size.** A thousand calibration libraries is not excessive. At 300
the 95th percentile estimate is noisy enough that the achieved false positive
rate lands near 8 percent instead of 5. That is a 60 percent inflation from
calibration noise alone, with every other part of the method correct.

## 6. The sweep

Six depths from 0.01x to 5x, nine tumour fractions from 1 to 40 percent, ten
lesions, 200 replicate libraries per condition. 540 conditions and 108,000
detector calls.

The panel, the calibration set and the threshold are rebuilt at every depth,
because a panel built at one depth does not transfer to another. Random streams
are spawned from one seed sequence per depth, so the whole sweep is
reproducible and the depths are independent of each other.

Two quantities are recorded per condition:

- **detected**, the fraction of replicates where any segment was called.
- **localised**, the fraction where the calls covered at least half the true
  lesion.

The gap between them answers a question readers ask: when the detector fires,
is it firing in the right place? It is. At 0.1x, in conditions where detection
reached 90 percent, localisation succeeded 95 percent of the time. Detection is
the bottleneck, not localisation.

## 7. Limits of detection

The sweep records probabilities on a grid. The limit of detection is the
smallest tumour fraction reaching 90 percent detection, linearly interpolated
between the two bracketing grid points.

When the grid never reaches 90 percent the limit is recorded as missing, and
that is a real answer rather than a gap: it means the lesion was undetectable
at every tumour fraction tested, up to 40 percent. Dropping those rows is how
a sweep ends up quietly reporting only its successes.

## 8. Junction level evidence

For the lesions depth cannot see, the question becomes how much evidence the
alternative needs.

Consider a heterozygous balanced rearrangement in a clone making up fraction
`f` of the cell free DNA. A library molecule carries the junction only if it
came from the derivative chromosome, which is one of the two homologues, and
only if it physically spans the breakpoint.

Let `N` be the number of molecules sequenced, `L` the fragment length and `G`
the haploid genome length. The probability that a given molecule spans a
particular point is `L / G`. The fraction of molecules from the derivative
allele is `f / 2`. So

```
E[spanning molecules] = N * (L / G) * (f / 2)
```

and since physical coverage is `c = N * L / G`, this collapses to

```
E = c * f / 2
```

Genome length cancels. Fragment length cancels. The result depends only on
coverage and tumour fraction, which is why it is worth stating as a formula
rather than a simulation.

Not every spanning molecule yields usable evidence. A split read needs enough
unique sequence on each side of the breakpoint to map confidently, taken here
as 20 bases, so on a 167 base fragment 76 percent of spanning molecules
qualify. Discordant pair evidence would need the junction to fall between the
two reads of a pair, which for 150 base paired end reads on 167 base fragments
never happens: the pair covers the whole fragment. For cell free DNA, split
reads are the only junction evidence on offer.

Detection is then a Poisson tail probability: the chance of observing at least
`k` supporting reads when the expectation is `0.38 * c * f`. Three is the
working minimum for `k`, because one or two supporting reads are
indistinguishable from mapping artefacts at scale.

The Poisson assumption is the right one here. Junction spanning molecules are a
heavily thinned sample of an enormous fragment population, which is exactly the
regime where the binomial becomes Poisson.
