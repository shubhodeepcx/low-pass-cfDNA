# Limitations

Written before anyone had to ask. Every item here is a real reason to doubt a
number in this repository, ordered by how much damage it does.

## 1. The bias model is a reasoned guess, not a measurement

The per library bias standard deviation of 0.045 is the most consequential
number in the study, and it was chosen to be plausible rather than fitted to
data. It sets where the detection plateau sits.

What this does and does not threaten:

- **Does not threaten:** that a plateau exists. Any non zero per library bias
  produces one, because that component does not shrink with depth while
  counting noise does. The qualitative finding is structural.
- **Does threaten:** the claim that the plateau sits at about 0.5x. A smaller
  bias moves it deeper, a larger bias moves it shallower. Take the specific
  depth as order of magnitude.

The fix is to estimate the parameter from a set of real low pass cell free DNA
libraries by measuring per bin variance across libraries after GC correction
and panel normalisation, then subtracting the Poisson component. That is a
morning of work given access to the data, and it would convert every absolute
figure here from an estimate into a prediction.

## 2. The detector is of the right class but is not WisecondorX

The comparison a reviewer will ask for first is against the actual binaries.
This repository reimplements the approach, faithfully as far as the published
descriptions go, but WisecondorX carries years of refinements that a
reimplementation does not.

The direction of the error is not obvious. The reimplementation may be weaker
than the real tool, in which case the detection limits here are pessimistic.
Against that, this detector tests around 2,500 candidate segments with a
genome wide calibrated threshold, which is a more careful multiple testing
treatment than some published pipelines apply, so on specificity it may be
stricter.

Either way, the copy neutral result is unaffected: no bin count method of any
quality can see a lesion that changes no bin count.

## 3. GC and mappability are modelled, not measured

The GC covariate is two sinusoids and some noise. Real GC structure is
correlated with gene density, replication timing and repeat content in ways
this does not capture, and real mappability has sharp structure around
segmental duplications rather than the smooth centromeric decay used here.

The consequence is that the bias structure is right in character and wrong in
detail. Since the detector removes the shared component by panel normalisation
anyway, the effect on detection limits should be small, but "should be" is
doing work in that sentence and it has not been verified against measured
tracks.

## 4. Lesions are clean, real karyotypes are not

Every simulated lesion has sharp boundaries at exact copy number one or three.
Real samples bring subclonal structure, where part of the clone carries the
lesion and part does not, which is effectively a lower tumour fraction and is
therefore partly covered by the sweep. They also bring whole genome doubling,
which changes the baseline ploidy and would need a separate treatment, and
chromothripsis, where copy number oscillates over short distances in ways a
segment statistic handles badly.

The complex karyotype case is five simultaneous clean lesions, which is a
friendlier object than a real complex karyotype.

## 5. One tumour fraction per library

The simulation mixes exactly two populations: one abnormal clone and normal
diploid background. Real cell free DNA from a marrow disorder carries multiple
subclones at different fractions, and the clone that carries del(5q) may not be
the clone that carries trisomy 8.

For detection of any single lesion this is a reasonable approximation. For
anything about clonal architecture it is not, and this repository makes no
claims there.

## 6. The junction calculation assumes a known breakpoint region

The closed form says what coverage is needed to observe junction spanning
molecules. It assumes you know where to look, which is true for the recurrent
fusions and false for novel rearrangements. A genome wide search for unknown
junctions carries a multiple testing burden the calculation does not model, so
the stated coverage requirements are optimistic for discovery and about right
for targeted confirmation.

It also assumes one derivative chromosome per diploid genome. A balanced
translocation present in both homologues, or amplified as part of a complex
rearrangement, would be easier to see by the corresponding factor.

## 7. Mapping artefacts are not modelled

The requirement for three supporting reads is a stand in for the real problem,
which is that one or two reads spanning an apparent junction are frequently
mapping errors, especially near repeats. A proper treatment would model the
artefact rate and derive the read threshold from a target false discovery rate
rather than assuming three.

Since the required coverage scales linearly with the read threshold, and the
threshold could plausibly be anywhere from two to five, the junction cost
numbers carry roughly a factor of two of uncertainty from this alone. That does
not change the conclusion, because the gap being argued about is three orders
of magnitude.

## 8. The false positive rate is controlled, not the false discovery rate

A 5 percent genome wide false positive rate means 5 percent of normal libraries
produce at least one call. It says nothing about what fraction of calls in a
real cohort would be wrong, which depends on how many abnormal samples the
cohort contains. For a screening application, where most samples are normal,
the relevant quantity is the false discovery rate and it would be considerably
less flattering.

## 9. No validation against real data

Nothing here has been checked against a sample with a known karyotype and a
measured tumour fraction. Every number is internal to the model. The study is
useful for bounding what is possible before spending money on sequencing, and
it is not evidence about any particular assay's performance.

## 10. Simulation studies flatter their own assumptions

The general risk. A simulation can only be wrong in the ways its author
imagined, and the detector was written by the same person who wrote the
simulator, against the same mental model of what matters.

The partial defence is that the main findings are not fitted quantities. The
copy neutral result is a theorem about depth based methods rather than an
empirical measurement. The plateau follows from the structure of the bias
model rather than its parameter values. The square root scaling of detection
with lesion size follows from the segment statistic. The numbers that are
genuinely at risk are the specific detection limits, and those are stated with
the model they depend on.
