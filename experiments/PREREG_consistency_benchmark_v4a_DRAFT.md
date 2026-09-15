# How much context coverage do perturbation rankings need, and does it depend on which transformations the score treats as irrelevant?

**Status:** DRAFT v4a — not frozen. To be frozen with `prereg freeze` after the
algebraic redundancy audit, the metadata-only feasibility audit, and the control
lineage audit described below, and before any consistency score is computed.

## Foreknowledge

Direction instability, D = 1 − mean pairwise cosine of a perturbation's
signatures across contexts, has been computed and published on five datasets
across two companion manuscripts. Held-out context prediction has been computed
on exactly one. No reliability, rank-stability or metric-comparison analysis has
been run on any of them.

| dataset | D computed | held-out prediction | role here |
|---|---|---|---|
| LINCS L1000 Phase I (GSE92742) | companion, atlas, paper 12 | **yes — 66 folds** | pilot and diagnostic only; cannot contribute confirmatory evidence |
| LINCS L1000 Phase II (GSE70138) | no | no | **primary confirmatory** |
| Tahoe-100M | atlas | no | replication; marginal D distribution previously seen |
| JUMP-CP CRISPR and ORF | atlas | no | two related technical replications, one dataset family |

Seen on Phase I: mean raw Spearman ρ = −0.602 against held-out cosine, negative
in all 66 folds; a variance-penalized variant at +0.497; the full penalty-weight
sensitivity curve. That comparison was scored on signed correlations between two
oppositely-oriented statistics, so its reported advantage reflects orientation
rather than prediction. This registration's orientation protocol exists because
of that failure.

## Question

A perturbation profiled in 9 contexts and one profiled in 60 receive scores on
the same scale and are ranked against each other, and the top of that ranking
decides what gets followed up. This registration measures how many contexts a
consistency ranking needs before it reproduces, and whether that requirement
depends on which transformations of the response the score treats as irrelevant.

**Estimand.** Reliability of rankings over the finite available context panel of
each dataset, not over a superpopulation of possible cell types. D is an
order-two U-statistic with a bounded cosine kernel; under a non-degenerate
first-order projection its variance is order 1/K.

**Closest prior work.** `@dewolf` varied the number and identity of LINCS cell
lines and time points, recommending at least three cell lines with marginal
benefit beyond roughly six or seven, for transcriptional characterization rather
than ranking stability across perturbations. `@limpavlidis` report low
consistency between CMap 1 and CMap 2 prioritization results. CMap Touchstone's
nine-cell-line core is a resource design, not a validated reliability threshold.
No published work constructs independent rankings from disjoint context subsets
over a fixed perturbation cohort and measures how ranking agreement changes with
context count, and none compares consistency metrics on that axis.

## Metric set

Entries are chosen one per **invariance class** — the transformation group under
which each score is unchanged, in the sense of Stevens' scale taxonomy
generalized to group actions. Invariance is the design principle; it does not by
itself prevent redundancy, because rank stability observes only orderings and is
blind to any monotone relationship between two scores. The redundancy audit
below is what verifies non-duplication.

| # | class | a perturbation counts as consistent even if… | metric |
|---|---|---|---|
| 1 | independent per-context positive rescaling | it is ten times stronger in one context | C_cos = 1 − D |
| 2 | independent per-context strictly monotone | effect sizes are ordered the same but spaced differently | mean pairwise Spearman |
| 3 | everything but signed support membership | only which features move, and in which direction | mean signed top-k overlap |
| 4 | common affine only, not context-specific | — location and scale disagreement both count against it | mean pairwise Lin concordance |
| 5 | *relative*, background-dependent | its response is identifiable among distractors | per-perturbation average precision |

Class 5 is a **relative identifiability criterion**, not an intrinsic geometry:
average precision depends on the distractor corpus and inherits the invariances
of whatever similarity produced its ranking. It is reported as its own class.
Its background is frozen once from the full panel and held identical across all
m; otherwise m changes the null.

Mean pairwise Pearson is admitted only if the redundancy audit separates it from
class 1. The frozen preprocessing does not feature-center signatures, so Pearson
is not algebraically identical to cosine, but the two may be empirically
near-duplicate.

Energy distance and maximum mean discrepancy are a declared Tahoe-only
extension. They are excluded from H3, whose coverage measure is itself an energy
distance.

**Contamination robustness** is a recognized missing class: every pairwise mean
can be moved by one aberrant context, which enters m − 1 of the pairs. A robust
consensus deviation — median angular deviation from a jointly estimated
spherical median — is registered as **exploratory**, excluded from the headline
metric count, with its orientation frozen before any score is inspected.

## The control, which is a unit test

Including self-comparisons in the average gives, exactly and unconditionally,

    V_m = (1/m²) ΣᵢΣⱼ cos(sᵢ, sⱼ) = 1/m + ((m−1)/m) · c̄

where c̄ is the mean over distinct pairs, so V_m − c̄ = (1 − c̄)/m. The bias is
(1 − θ)/m for population consistency θ, and the m = 3 → 20 drift is
(17/60)(1 − θ) ≈ 0.283(1 − θ) — a coefficient, not a constant. At the θ ≈ 0.81
implied by published pan-HDAC values the predicted drift is about 0.05; at
θ ≈ 0 it approaches 0.28. Removing this bias is what yields pairwise phase
consistency `@vinck2010ppc`; the derivation for the squared phase-locking value
is `@aydore2013plv`.

**This identity holds exactly on the same sample, so it is a software unit test
rather than an empirical hypothesis.** Deviation beyond floating-point tolerance
indicates a pipeline error, not a finding. The registered prediction is the
conditional quantity Δ̂ₚ = (17/60)(1 − θ̂ₚ), with θ̂ₚ estimated from the
full-panel off-diagonal mean and never from the V_m realization being tested.

**V_m is rank-degenerate with C_cos.** Being affine in c̄ at fixed m, it induces
an identical ordering and therefore an identical stability curve. It appears on
the bias curve only, is never counted as a metric, and contributes nothing to
any stability comparison.

Kendall's W is excluded: for the complete rankings used here, mean pairwise
Spearman = (mW − 1)/(m − 1), an affine transform at fixed m, so W induces
identical rankings and identical stability curves. Tie corrections are
negligible at 978 features and do not rescue it.

## Common-cohort constraint

Eligibility at subset size m requires 2m contexts, so a naive curve across m
compares different perturbations at every point and can manufacture a trend from
population change alone. Every comparison across m is made on a **fixed cohort**.

*Feasibility audit, metadata only.* Before freezing, context counts per
perturbation are tabulated from metadata; no signature is loaded and no score is
computed. The primary grid is the largest grid whose common cohort — the
perturbations eligible at its largest m — holds at least 650 perturbations.
Candidate grid m ∈ {3, 4, 5, 6, 8, 10, 15, 20}, truncated from the top until met.
The realized grid is recorded in the frozen file.

Pairwise secondary: for each pair m_a < m_b, stability at both sizes is computed
on the cohort eligible at m_b. Each pair defines a **pair-specific conditional
estimand**; points are reported separately and never stitched into one curve.

## Hypotheses

**H1.** V_m − c̄_m = (1 − c̄_m)/m reproduces to floating-point tolerance at every
m, and the bias curve matches Δ̂ₚ = (17/60)(1 − θ̂ₚ) within its 90% interval.
*Registered as a pipeline unit test; deviation indicates a bug, not a result.*

**H2.** Chance-corrected top-decile ranking stability for C_cos increases with
the number of sampled contexts.

**H3.** Within a perturbation and at fixed m, context subsets providing better
baseline-covariate coverage of a candidate pool yield more accurate prediction in
evaluation contexts excluded from coverage estimation, subset selection and score
construction, than equally sized worse-covering subsets.

**H4.** Delete-one-context jackknife standard error is positively associated with
half-sample percentile-rank displacement.

**H5.** H2's contrast holds in at least two of three dataset families.

**H6.** The number of contexts required to reach 0.5 chance-corrected top-decile
stability differs across invariance classes, so which transformations a score
treats as irrelevant changes how much coverage it needs.

**H7.** At matched m, chance-corrected top-decile stability is lower among
sparsely-profiled perturbations than among well-profiled ones, so a coverage
requirement estimated on the well-covered cohort understates what the sparse
regime needs.

H6 carries the design. H2 establishes that the curve exists; H6 is the claim that
the benchmark is for, and it is the outcome a practitioner acts on.

**Failure of H2 for C_cos neither terminates nor voids the comparator analyses.**
If one metric's curve is flat while another's rises, that contrast is a result of
this study and is reported as such.

## Method

Within each dataset, for each perturbation in the cohort, draw two disjoint sets
of m contexts, score independently in each under every metric, and rank the
cohort within each half. Repeat over 200 balanced partitions, seeds fixed;
partitions increase if the Monte Carlo standard error of any primary curve point
exceeds 0.01.

**Primary statistic.** Chance-corrected Jaccard overlap of the top decile between
half-sample rankings, (J_obs − J_exp)/(1 − J_exp). Because all items tied at a
decile boundary are included, lists need not hold exactly N = qM items, so J_exp
uses **realized** list sizes a and b via E[|A ∩ B|] = ab/M. Top lists contain the
most-consistent perturbations under the frozen orientation.

**Headline quantity.** m₀.₅ⱼ = min{m in the realized grid : Ŝⱼ(m) ≥ 0.5},
evaluated at observed grid points with no monotone smoothing in the confirmatory
result. A metric that never reaches 0.5 is reported as m₀.₅ > m_max, right-
censored, never dropped. Where realized grids differ across datasets, the claim
is stated per family and never pooled.

**Secondary.** Spearman ρ between half-sample rankings over the full cohort. A
label-permutation null run through the complete ranking pipeline, preserving the
actual control structure, is registered as a sensitivity analysis on the analytic
correction and as the leakage floor.

**External utility.** Held-out-context cosine against the training consensus and
held-out top-feature Jaccard, as functions of m, for every metric.

**Coverage for H3.** Each perturbation's contexts are split once by seed into a
candidate pool and an evaluation set. Coverage estimation, subset selection and
score construction use the candidate pool only; evaluation contexts are seen by
none of the three. Coverage is energy distance between subset and candidate pool
on baseline covariates only. The design is paired within perturbation: matched-
size subsets from the best and worst coverage deciles, compared on the same
evaluation contexts. |D_m − D_full| is a secondary reference-reconstruction
diagnostic only — better-covering subsets are selected for resembling the pool
that defines D_full, so their closeness to it partly follows from the selection
rule.

**Compute.** Per-metric wall-clock and peak memory are reported as a secondary
axis.

## Orientation and redundancy protocol

**Orientation** is frozen in a static registry before any score is computed:
every metric is defined so that higher means more consistent. No data-dependent
sign flip is permitted at any stage.

**Redundancy is audited algebraically first.** Known identities — Kendall's W
with mean pairwise Spearman, squared Euclidean distance on unit vectors with
cosine, Euclidean centroid dispersion with the pairwise mean, resultant length
with V_m — are resolved on paper and the redundant entries excluded with the
identity stated. Only then is the empirical audit run: a full-panel Spearman
matrix among all admitted metrics, published in full regardless of outcome. Pairs
at |ρ| ≥ 0.95 are declared one effective entry; pairs in 0.80–0.95 are retained
and flagged in the multiplicity accounting.

**Near-duplicates are never averaged.** A canonical representative is chosen and
the other reported as a supplement. Averaging would create an unvalidated
composite that silently double-weights one construct.

## Inference criteria

| hypothesis | holds when |
|---|---|
| H1 | identity reproduces within floating-point tolerance and the bias curve lies inside the 90% interval of Δ̂ₚ at every m |
| H2 | the prespecified linear contrast of chance-corrected top-decile stability against log m is positive for C_cos, context-block bootstrap 95% CI excluding zero, on the fixed cohort in Phase II |
| H3 | within-perturbation paired difference in held-out prediction favors better coverage, crossed-bootstrap 95% CI excluding zero, in the primary dataset |
| H4 | Spearman ρ between jackknife SE and percentile-rank displacement is positive, CI excluding zero, in the primary dataset |
| H5 | H2's contrast is positive with CI excluding zero in at least two of three dataset families |
| H6 | the range of m₀.₅ across admitted invariance classes spans at least one grid step, with non-overlapping simultaneous bands for the extreme pair, in the primary dataset |
| H7 | stability at matched m is lower in the sparsest band than in the best-covered band, context-block bootstrap 95% CI on the paired difference excluding zero, at every shared m in the primary dataset |

C_cos is the confirmatory metric and carries H2's single formal test. Comparisons
across the metric set are reported as **simultaneous bands** from a
Westfall–Young permutation max-statistic resampling perturbations as primary
units, not as per-metric p-values. H3–H5 are secondary and Holm-corrected within
that family. Half-sample partitions are Monte Carlo replicates and are never
treated as independent observations.

**A dataset's top-decile analysis is void if its common cohort holds fewer than
650 perturbations.** For two random top-N sets from M items, E[|A ∩ B|] = N²/M
and SE(X/N) ≈ (1 − q)/√M at q = N/M; detecting a 0.10 increase over chance at 80%
power and two-sided α = 0.05 requires M ≳ [(1.96 + 0.84)(1 − 0.10)/0.10]² ≈ 635.
Full-cohort Spearman remains descriptive when this fails. Repeated partitions
reduce Monte Carlo error and cannot substitute for cohort size.

**H2's estimand is the common cohort, not the dataset.** The fixed cohort holds
the best-sampled perturbations, so the primary curve describes the well-covered
regime. H7 tests whether it transfers to the sparse regime rather than assuming
it does.

**Panel-size stratification for H7.** Perturbations are assigned to bands by
total available contexts K, with band edges set by the metadata-only feasibility
audit so that each band clears the 650 floor and the bands partition the dataset.
Within each band the stability curve is computed at every m the band supports.
Bands are compared **only at m values they share**, so the contrast is at matched
context count and never across different subset sizes.

Bands differ in more than context count: a perturbation profiled in few contexts
is typically less studied and may have weaker or noisier response. H7 therefore
separates context count from everything else that accompanies being sparse, and
a difference at matched m is attributable to the latter. That is the informative
outcome, not a confound to be regretted — it would mean a count-based
recommendation is insufficient on its own.

**JUMP-CP results cannot support any claim about biological generalization.**
Plate is the only available context axis. A positive JUMP result establishes the
statistical mechanism on technical contexts and nothing about transfer across
cell types, whatever it returns.

**Half-sample independence is an empirical claim about provenance, not an
assumption.** If both halves compute signatures against a shared control
population, measured stability is partly common-control noise. Control lineage is
audited from metadata before freezing; where disjoint-control construction is
feasible it is the primary, otherwise a shared-versus-disjoint sensitivity
analysis is registered. Control units are never split below the well or plate
level.

## Continuation

Continuation to Tahoe and JUMP depends on feasibility, never on Phase II's
direction or significance: the estimands must be computable, cohorts must clear
the size floor, curves must be non-degenerate. A null on H2 in Phase II is
reported and the replications still run, because the question is about generality
and answering it conditionally on a positive first result would not answer it.

## Frozen implementation choices

Replicate aggregation, dose and time selection, feature filtering, normalization,
sphering, zero-norm handling and missing-context rules are fixed per dataset in
`config/` before freezing. The signed top-k construction — whether the positive
set is the k largest signed values or the k largest absolute values — is frozen
explicitly, with k prespecified at more than one value. Ties receive midranks and
all items tied at a decile boundary are included. D_full is computed on the
candidate pool and its dependence on the subsets evaluated against it is stated
wherever reported. JUMP CRISPR and ORF are two related technical replications
within one dataset family.

## Interpretation

| observed | conclusion |
|---|---|
| unit test passes, stability rises with m, m₀.₅ differs across classes | context requirement is metric-dependent; the per-class m₀.₅ table is the deliverable |
| unit test passes, stability rises, m₀.₅ identical across classes | context requirement is a property of the data, not of the invariance choice |
| stability flat for all metrics across the realized grid | rankings are stable throughout the tested range and no coverage guidance is warranted |
| H7 confirmed — sparse band less stable at matched m | context count is necessary but not sufficient; the recommendation must be a per-perturbation reliability estimate, not a threshold on K |
| H7 null — bands coincide at matched m | the coverage requirement transfers, and the headline m₀.₅ applies across the dataset rather than only to well-covered perturbations |
| unit test fails | pipeline error; no scientific conclusion is drawn until it is resolved |

## Sample size

The realized cohort is set by the feasibility audit and recorded before freezing.
The design can resolve whether ranking stability varies with context count across
the realized grid on a fixed cohort, whether m₀.₅ differs across invariance
classes by at least one grid step, and whether stability at matched m differs
between the sparsest and best-covered panel-size bands. It cannot resolve behavior below the smallest
realized m, cannot separate context identity from context count for perturbations
on a fixed panel, and cannot establish that any metric is correct — only how much
coverage each one needs before its ranking reproduces.

## Citations to resolve before freezing

`@dewolf`, `@limpavlidis`, `@vinck2010ppc`, `@aydore2013plv`. Resolve with
`citations resolve` and pin each record before `prereg freeze`.
