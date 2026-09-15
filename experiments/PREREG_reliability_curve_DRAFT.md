# How many contexts does a direction-instability ranking need?

**Status:** DRAFT — not frozen. To be frozen with `prereg freeze` before any
analysis is run on the datasets below.

## Foreknowledge

Direction instability (D = 1 − mean pairwise cosine of a perturbation's
signatures across contexts) has been computed and published on five datasets
across two companion manuscripts. Held-out context prediction has been computed
on exactly one.

| dataset | D computed | held-out prediction | role here |
|---|---|---|---|
| LINCS L1000 Phase I (GSE92742) | companion, atlas, paper 12 | **yes — 66 folds** | diagnostic only, cannot contribute confirmatory evidence |
| LINCS L1000 Phase II (GSE70138) | no | no | **primary confirmatory** |
| Tahoe-100M | atlas | no | replication; marginal D distribution previously seen |
| JUMP-CP CRISPR / ORF | atlas | no | technical mechanism test only |

The author has seen, on Phase I: mean raw Spearman ρ = −0.602 against held-out
cosine, negative in all 66 folds; the transport-stable variant at +0.497; and
the full λ sensitivity curve. That comparison was scored on signed correlations
between two oppositely-oriented statistics, so the reported 66/66 advantage
reflects orientation rather than prediction. No reliability or rank-stability
analysis has been run on any dataset.

Nothing below has been computed anywhere.

## Question

Ranking perturbations by a consistency score assumes the ranking is stable. A
drug scored from 9 cell lines and one scored from 60 receive numbers on the
same scale and are ranked against each other. This registration asks at what
number of contexts that ranking becomes stable enough to act on, and whether
context diversity matters beyond context count.

D is the complement of pairwise phase consistency `@vinck2010ppc`, which
averages over distinct pairs and is therefore unbiased for its population
quantity. Its *variance*, and the composition of contexts it averages over, are
not controlled. Low between-unit reliability caps any correlation a ranking can
show `@hedge2018reliability`; the analogous audit in neural population geometry
found norm-based metrics drifting with sample size and cross-session ICC below
0.15 for 12 of 16 metrics `@tower2026neurogeometry`.

## Hypotheses

**H1.** Within-perturbation subsampling holds E[D] approximately constant as the
number of sampled contexts K varies, in contrast to the √n drift reported for
norm-based neural geometry metrics.

**H2.** Top-N rank overlap between disjoint half-samples increases monotonically
with K, for N ∈ {10, 50, 100}.

**H3.** At matched K, perturbations whose sampled contexts are more diverse show
higher rank stability than perturbations whose contexts are narrow, so context
count alone does not determine reliability.

**H4.** A per-perturbation delete-one-context jackknife standard error is
negatively associated with that perturbation's rank displacement between
disjoint half-samples, making it usable as a shipped reliability flag.

**H5.** The H2 relationship holds in at least two of the three datasets.

H2 carries the design. If rank stability does not vary with K, there is no
reliability problem to characterize and H3 and H4 describe nothing.

H1 is registered as a result that would strengthen D rather than embarrass it:
confirmation establishes that D is free of the sample-size bias that invalidated
the neural geometry metrics, and that only its variance and coverage require
handling.

## Method

For each perturbation with at least 2m eligible contexts, draw two disjoint
sets of m contexts, compute D independently in each, and rank all eligible
perturbations within each half. Repeat over 200 balanced partitions.
m ∈ {3, 4, 5, 6, 8, 10, 15, 20}, subject to eligibility.

Reliability at each m is reported as a curve, not a test: top-N rank overlap
for N ∈ {10, 50, 100}, and Spearman ρ between the two half-sample rankings. No
reliability threshold is set. A reader selects an operating point from the curve
according to what a false follow-up costs them.

Context diversity for H3 is the mean pairwise distance among a perturbation's
sampled contexts in a reference space defined before analysis: baseline
expression profiles for the transcriptomic datasets, untreated-well morphology
for JUMP-CP. Perturbations are split at the median diversity within each m
stratum.

Uncertainty for H4 is the delete-one-context jackknife standard error of D.

## Datasets

**LINCS L1000 Phase II (GSE70138).** Primary. Chemical perturbations, bulk
transcriptomic, cell lines as contexts. Dose and time held fixed where metadata
permit; unmatched conditions excluded rather than pooled.

**Tahoe-100M.** Replication. Chemical perturbations, single-cell transcriptomic
aggregated to one signature per drug-dose-cell-line condition, cell lines as
contexts.

**JUMP-CP CRISPR and ORF.** Mechanism test. Genetic perturbations, morphological
profiles, **plates as contexts**. `Metadata_Source` is single-valued within
perturbation, so the only available context axis is technical.

## Inference criteria

| hypothesis | holds when |
|---|---|
| H1 | median \|Spearman ρ(D, m)\| across perturbations < 0.15 in every dataset, against ρ = +0.989 reported for the uncorrected neural metric |
| H2 | top-N overlap is monotonically non-decreasing in m for all three N in the primary dataset, and strictly increasing from the smallest to the largest m |
| H3 | top-50 overlap for high-diversity perturbations exceeds that for low-diversity perturbations at matched m, in a majority of m strata |
| H4 | Spearman ρ between jackknife SE and half-sample rank displacement is positive in the primary dataset, with a context-level bootstrap CI excluding zero |
| H5 | H2 holds in at least two of three datasets |

Uncertainty is estimated by hierarchical bootstrap over contexts and then
perturbations within context, 10,000 replicates. Half-sample partitions are not
independent replicates and no test treats them as such.

**All hypotheses are void if H2 fails.**

**A dataset is void if fewer than 200 perturbations are eligible at m = 5.**

**H3 is void if context diversity does not vary within m strata, quantified as
an interquartile range below 10% of the median.**

**JUMP-CP results cannot support any claim about biological generalization.**
Plate is technical variation. A positive JUMP result establishes the statistical
mechanism and nothing about transfer across cell types, whatever it returns.

## Stopping

If the primary dataset returns a null on H2, the analysis stops and the null is
reported. Datasets are run against this frozen plan in the order LINCS Phase II,
Tahoe, JUMP-CP; results from any dataset do not modify the hypotheses or the
criteria above.

## Sample size

Eligibility at m = 5 requires 10 contexts per perturbation. Phase I carried
8,949 drugs with ≥ 5 cell lines, of which a minority reach 10, so the eligible
set is expected in the hundreds to low thousands rather than the full panel.
The design can resolve whether rank stability varies with K across the tested
range. It cannot resolve the shape of that curve below m = 3, and it cannot
separate cell-line identity from cell-line count for perturbations profiled in
a fixed panel.
