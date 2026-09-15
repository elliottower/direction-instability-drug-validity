# How do context count and coverage affect the reliability of direction-instability rankings?

**Status:** DRAFT v2 — not frozen. To be frozen with `prereg freeze` after the
metadata-only feasibility audit below and before any score is computed.

## Foreknowledge

Direction instability (D = 1 − mean pairwise cosine of a perturbation's
signatures across contexts; lower D means greater consistency) has been
computed and published on five datasets across two companion manuscripts.
Held-out context prediction has been computed on exactly one.

| dataset | D computed | held-out prediction | role here |
|---|---|---|---|
| LINCS L1000 Phase I (GSE92742) | companion, atlas, paper 12 | **yes — 66 folds** | pilot and diagnostic only; cannot contribute confirmatory evidence |
| LINCS L1000 Phase II (GSE70138) | no | no | **primary confirmatory** |
| Tahoe-100M | atlas | no | replication; marginal D distribution previously seen |
| JUMP-CP CRISPR and ORF | atlas | no | two related technical replications, not two modalities |

Seen on Phase I: mean raw Spearman ρ = −0.602 against held-out cosine, negative
in all 66 folds; the transport-stable variant at +0.497; the full λ sensitivity
curve. That comparison was scored on signed correlations between two
oppositely-oriented statistics, so the reported 66/66 advantage reflects
orientation, not prediction. No reliability or rank-stability analysis has been
run on any dataset.

## Question

A drug scored from 9 cell lines and one scored from 60 receive numbers on the
same scale and are ranked against each other. D is the complement of pairwise
phase consistency `@vinck2010ppc` and is unbiased for its population quantity,
but its variance and the composition of contexts it averages over are
uncontrolled. Low between-unit reliability caps any correlation a ranking can
support `@hedge2018reliability`. This registration asks how ranking stability
varies with the number of sampled contexts, and whether which contexts were
sampled matters beyond how many.

**Estimand.** Reliability of rankings over the finite available context panel
of each dataset, not over a superpopulation of possible cell types.

## The common-cohort constraint

Eligibility at subset size m requires 2m contexts, so a naive curve across m
compares different perturbations at every point and can manufacture a trend
from population change alone. Every comparison across m is therefore made on a
**fixed cohort**.

*Feasibility audit, metadata only.* Before freezing, context counts per
perturbation are tabulated from metadata. No signature is loaded and no D is
computed. The primary grid is the largest m grid whose common cohort — the
perturbations eligible at the largest m in the grid — contains at least 650
perturbations (justified below). Candidate grid m ∈ {3, 4, 5, 6, 8, 10, 15, 20},
truncated from the top until the constraint is met. The realized grid is
recorded in the frozen file.

Secondary pairwise analysis: for each pair m_a < m_b, stability at both sizes is
computed on the cohort eligible at m_b, preserving perturbations the primary
cohort excludes.

## Hypotheses

**H1.** For each m, the mean within-perturbation difference between subsampled
D_m and full-panel D_full is reported as a bias curve with 90% confidence
intervals, against ρ = +0.989 for the uncorrected neural-geometry metric as a
positive control. *Registered as a model check; no pass/fail criterion.*

**H2.** Chance-corrected top-decile ranking stability increases with the number
of sampled contexts.

**H3.** Within a perturbation and at fixed m, context subsets that better cover
that perturbation's full available panel produce smaller |D_m − D_full| and more
accurate held-out-context prediction than worse-covering subsets of the same
size.

**H4.** Delete-one-context jackknife standard error is **positively** associated
with a perturbation's half-sample percentile-rank displacement.

**H5.** H2 holds in at least two of three dataset families (LINCS Phase II,
Tahoe, JUMP).

H2 carries the design. If ranking stability does not vary with context count on
a fixed cohort, H3 and H4 describe nothing.

H3 is registered as an outcome that can strengthen D rather than only embarrass
it: if coverage rather than count drives error, a practitioner can act on which
contexts to add rather than being told only to collect more.

## Method

Within each dataset, for each perturbation in the cohort, draw two disjoint
sets of m contexts, compute D independently in each, and rank the cohort within
each half. Repeat over 200 balanced partitions, seeds fixed; partitions are
increased if the Monte Carlo standard error of any primary curve point exceeds
0.01.

**Primary statistic.** Chance-corrected Jaccard overlap of the top decile
between half-sample rankings: (J_obs − J_exp)/(1 − J_exp), with J_exp derived
from the hypergeometric expectation N²/M for a cohort of size M. Top lists
contain the **lowest** D values.

**Secondary.** Spearman ρ between the two half-sample rankings over the full
cohort. Top-N overlap for N ∈ {10, 50, 100} is reported as a user-facing display
only and carries no hypothesis. Rank-biased overlap is a sensitivity analysis.

**External utility.** Held-out-context cosine against the training consensus,
and held-out top-feature Jaccard, both as functions of m. Internal rank
reproducibility alone does not establish that a stable ranking is useful.

**Coverage for H3.** Representativeness of a subset relative to the
perturbation's full panel, computed on baseline covariates only — untreated
expression profiles, and for imaging untreated wells normalized without
reference to treated profiles. Primary measure is energy distance between the
sampled subset and the full panel; secondary is maximum nearest-neighbor
distance from panel to subset. Feature set, distance metric and scaling are
frozen before analysis; no dimensionality reduction. The design is paired within
perturbation: for each drug and each m, generate matched-size subsets in the
best and worst coverage deciles and compare them to each other, holding drug
identity, m and available panel fixed.

**Uncertainty for H4.** Delete-one-context jackknife standard error computed
from the observed m-context subset alone, so it is available to a practitioner
who has only that subset. Displacement is the absolute difference in percentile
rank between halves.

## Inference criteria

| hypothesis | holds when |
|---|---|
| H1 | *no criterion; bias curve with 90% CI reported for every m* |
| H2 | the prespecified linear contrast of chance-corrected top-decile stability against log m is positive, with a context-block bootstrap 95% CI excluding zero, on the fixed cohort in Phase II |
| H3 | within-perturbation paired difference in \|D_m − D_full\| between best- and worst-coverage subsets is positive, CI excluding zero, in the primary dataset |
| H4 | Spearman ρ between jackknife SE and percentile-rank displacement is positive, CI excluding zero, in the primary dataset |
| H5 | the H2 contrast is positive with CI excluding zero in at least two of three dataset families |

H2 is primary. H3–H5 are secondary and Holm-corrected within that family.
Uncertainty is estimated by a two-way crossed bootstrap resampling
perturbations as primary units, with context partitioning rerun inside each
replicate, 10,000 replicates. Half-sample partitions are Monte Carlo replicates
and are never treated as independent observations.

**All hypotheses are void if H2 fails.**

**A dataset's top-decile analysis is void if its common cohort holds fewer than
650 perturbations.** For two random top-N sets drawn from M items,
E[|A ∩ B|] = N²/M and SE(X/N) ≈ (1 − q)/√M at q = N/M. Detecting a 0.10
increase in overlap proportion over chance at 80% power, two-sided α = 0.05,
requires M ≳ [(1.96 + 0.84)(1 − 0.10)/0.10]² ≈ 635. Full-cohort Spearman results
remain descriptive when this fails. Repeated partitions reduce Monte Carlo error
and cannot substitute for cohort size.

**JUMP-CP results cannot support any claim about biological generalization.**
Plate is the only available context axis; `Metadata_Source` is single-valued
within perturbation. A positive JUMP result establishes the statistical
mechanism on technical contexts and nothing about transfer across cell types,
whatever it returns.

## Continuation

Continuation to Tahoe and JUMP depends on feasibility, never on Phase II's
direction or significance: the registered estimands must be computable, the
cohorts must clear the size floor, and the curves must be non-degenerate. A null
on H2 in Phase II is reported and the replications still run, because the
question this registration asks is about generality and answering it
conditionally on a positive first result would not answer it.

## Frozen implementation choices

Replicate aggregation, dose and time selection, feature filtering,
normalization, sphering, zero-norm handling and missing-context rules are fixed
per dataset in the accompanying `config/` before freezing. Ties receive
midranks, and all items tied at a top-decile boundary are included. D_full is
computed on a perturbation's entire available panel and its dependence on the
subsets evaluated against it is stated wherever reported. JUMP CRISPR and ORF
count as two related technical replications within one dataset family.

## Interpretation

| observed | conclusion |
|---|---|
| bias flat, stability rises with m | variance-limited; a per-perturbation uncertainty flag is the deliverable |
| bias drifts with m, stability flat | D inherits a sample-size artifact and needs normalization, not shrinkage |
| both | count affects the estimator and the ranking, and neither correction alone suffices |
| neither | D rankings are stable across the tested range and no reliability adjustment is warranted |

## Sample size

The realized cohort is set by the feasibility audit and recorded before
freezing. The design can resolve whether ranking stability varies with context
count across the realized grid, on a cohort held constant. It cannot resolve
behavior below the smallest realized m, and it cannot separate context identity
from context count for perturbations profiled on a fixed panel.
