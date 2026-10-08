# What could the unique-target permutation have detected, and what holds its null up at 0.37?

**Date:** 2026-10-08
**Status:** DRAFT. Not frozen. Revises
`PREREG_H3_S2T_POWER_ENVELOPE_DRAFT.md` after a code review of its designs. H1
becomes three geometry-preserving arms, H2 loses the cross-paper question it was
built on, and H3 becomes a repeated-cohort power estimate.
**Commit SHA:** pending
**Kind:** post-hoc diagnostic. No registered statistic is recomputed, no
registered criterion is changed, and no outcome here can move the gate at
`results/03c_h3_sensitivity/h3_sensitivity_results.json`, which is closed:
`h3_interpretation_survives: false`.
**Relates to:** `PREREG_H3_S1S3_CORRECTED_BASELINE.md`, frozen at 7f57136, whose
S2-T failed; `results/03j_reference_axis_geometry/reference_axis_geometry.json`,
the descriptive geometry of the reference;
`results/03k_partial_spearman_estimators/estimator_comparison.json`, which
settled what H2 was originally for.

## The questions

S2-T reassigns target directions among targets and found the association
unchanged, p_perm = 0.4521. Three readings of that are consistent with the
recorded numbers and they differ in what the paper may say.

1. The association does not depend on which target a drug carries.
2. The association is a property of projecting signatures onto any axis, so no
   reference could have carried target information through it.
3. The shRNA directions are too mutually similar for reassignment to change the
   axis enough, and a target-specific association would have survived the
   permutation regardless.

## Foreknowledge

Everything in the registered run is known, including every verdict. So is the
reference geometry: the 258 unique directions have median pairwise absolute
cosine 0.239 against an isotropic expectation of 0.026, a Gram participation
ratio of 13.0, a leading eigenvalue holding 26% of the trace, and 140 components
carrying 90% of it. The set has a strong common component and a long tail; it is
not thirteen-dimensional. Target cluster sizes are known: 131 of 258 targets
carry one drug, the largest carries 37, and the Kish effective cluster count is
82. That count warns against treating 258 targets as independent and is not a
plug-in effective sample size for this permutation distribution.

The recorded excess is known. S2-T refuses at p < 0.01, which on the recorded
null is an observed partial rho of 0.4599 or above, 0.0928 over the null median.
The observed 0.3722 sits 0.0877 below it. That figure describes the recorded
draws and is not a power calculation.

Rank residualization with an explicit intercept, the procedure registered at
f288507, is unbiased under conditional independence. The atlas project's
regression test reports that rank-then-residualize yields spurious partial
correlations of 0.05 to 0.25 and its comparator omits the intercept; measured on
that test's own model at n = 5000, the registered estimator has mean -0.0029
against the test's own unbiased band of 0.02, the no-intercept comparator has
mean +0.2646, and Pearson residualization gives -0.0022 with or without the
intercept. H2 was written to settle a contradiction between two of these papers.
There is no contradiction, so what remains of H2 is narrower and is stated as
such below.

## Hypotheses

Three substitutions replace the real reference, ordered by how much of its
geometry they keep, so that a difference from the recorded floor of 0.3671 is
attributable to the property each one destroys and to nothing else.

**H1a.** Replacing the 258 directions with independent isotropic unit directions
leaves the floor at 0.3671. This destroys pairwise similarity, spectrum, subspace
and orientation together, so it can only say whether a generic axis construction
is sufficient.

**H1b.** Replacing them with a set that has the same Gram matrix but a
Haar-random orientation in the ambient space leaves the floor at 0.3671. Every
pairwise cosine and every singular value survive; the orientation relative to the
drug-signature space does not.

**H1c.** Replacing them with a Haar rotation inside their own empirical row span
leaves the floor at 0.3671. The Gram matrix and the shared subspace both survive;
the identity of each direction does not.

**H2.** The registered estimator is centered on zero under a conditional
independence model fitted to this cohort's own covariate relationships at
n = 795.

**H3.** S2-T detects a planted target-specific contribution with probability at
least 0.80, at a contribution of 0.10 partial rho on the analysis scale, over
independently generated cohorts carrying the real directions, the real cluster
sizes and the real context counts.

H3 carries the diagnostic. H1a to H1c say what holds the floor up; only H3 says
which alternatives the test was capable of refusing, and a failing S2-T is
evidence about the data only for alternatives it could have refused.

**H3 is the outcome that strengthens the registered finding.** If detection
probability reaches 0.80 at or below the recorded excess of 0.0928, then reading
3 above is refuted and the failing S2-T means what it appears to mean.

## Inference criteria

Every median and every probability is reported with a Monte Carlo interval, and
every criterion is read against that interval rather than against a point.

| hypothesis | holds when |
|---|---|
| H1a | the 95% Monte Carlo interval for the arm's median contains 0.3671 |
| H1b | the same, for the Gram-preserving ambient rotation |
| H1c | the same, for the in-span rotation |
| H2 | the 95% Monte Carlo interval for the median registered partial rho contains 0 under both nulls |
| H3 | the estimated detection probability at a planted contribution of 0.10 is at least 0.80 and the lower bound of its 95% Wilson interval exceeds 0.70 |

No equivalence band is set on the H1 arms. A band of fixed width would be chosen
knowing the value it is tested against, and the question is which arms differ
from the recorded floor rather than whether a difference clears a threshold.

H3 is also run at planted contributions of 0.05 and 0.20 and the three points are
reported as a curve. The 0.10 level carries the criterion because it brackets the
recorded excess of 0.0928 from above; failing at 0.05 and passing at 0.20 is a
power curve and is reported as one, not as a verdict.

**Reading the H1 arms.** The arms are read jointly, through the property each one
destroys, and no single arm supports a claim on its own.

| pattern | reading |
|---|---|
| H1a holds | any arbitrary common-axis construction produces the floor |
| H1a fails, H1b holds | the inter-axis geometry alone is sufficient |
| H1b fails, H1c holds | overlap between the reference subspace and the drug-signature structure is sufficient |
| all three fail | the floor involves the exact orientation of the reference, which is implicated and not shown to be causal |

**All hypotheses are void if the identity permutation does not reproduce the
recorded observed statistic of 0.3721909618253966**, which is the check that the
diagnostic runs on the registered cohort.

**H3 is void if its simulation-validity gate fails**, below.

## Implementation, fixed before it runs

Every arm imports `partial_spearman`, `permuted_statistic`, `load_bundle` and
`target_representatives` from `experiments/03c_h3_sensitivity.py` at the sealed
hash 6290281e, and `unique_target_permutations` from `geometry.inference`, so the
diagnostic exercises the code the run used rather than a restatement of it. The
covariates stay the registered pair, mean pairwise difference norm `M_delta` and
context count `K`, and are recomputed from the signatures in every arm rather
than carried. Seeds are declared in the script and recorded in its output.

### H1a to H1c

Write the real direction matrix as `U = A S Vt`, where `U` is 258 by 978 with
unit rows and `r` is its numerical rank.

- **H1a** draws 258 independent standard normal vectors in 978 dimensions and
  normalizes each to unit length.
- **H1b** replaces `Vt` with a Haar-distributed `r` by 978 orthonormal frame, so
  `U* U*^T = U U^T` exactly and the rows stay unit length.
- **H1c** takes an orthonormal basis `Q` of the row span of `U`, draws a Haar
  rotation `R` of size `r`, and sets `U* = U Q R Q^T`, which preserves both
  `U U^T` and the span.

Each arm runs 2,000 replicate direction sets. Within a replicate, each target
receives its substituted direction, every drug of that target receives its
target's direction, and `P`, `E` and the registered partial correlation are
recomputed from the signatures. The arm's statistic is the median over
replicates, with a 95% interval from 2,000 bootstrap resamples of the replicate
medians. Haar matrices come from the QR decomposition of a standard normal
matrix with the sign correction that makes the distribution Haar.

### H2

Two nulls, both reported, neither claimed to be general.

The **parametric** null standardizes `M_delta` and `K`, regresses the recorded
`P` and `E` on them by ordinary least squares in the raw scale, and draws `P*`
and `E*` from the fitted means with independent Gaussian noise at the fitted
residual scales. It assumes linear conditional means, homoscedastic Gaussian
residuals and no target-cluster dependence, and it can place values outside the
natural ranges of `P` and `E`. The claim it licenses is that the estimator's
Monte Carlo median took a stated value under this linear-Gaussian null fitted to
the cohort's first-order covariate relationships, and not that the estimator is
unbiased under independence.

The **semiparametric** null is Freedman-Lane in form: fit `E` on the covariates
in the raw scale, permute the raw residuals across drugs, and add them back to
the per-drug fitted values. It preserves the fitted linear covariate relationship
and the empirical residual multiset. It does not preserve the marginal
distribution of `E`, and its validity rests on residual exchangeability, which
fails under nonlinear conditional means, heteroscedasticity or target-cluster
dependence. Three diagnostics are fixed here and reported with the arm: residual
variance against fitted value and against each covariate, the target-level
intraclass correlation of the residuals, and the marginal quantiles of `E`
before and after permutation. If the residual intraclass correlation exceeds
0.10, the permutation is applied to target-level residual blocks instead of to
drugs, and the switch is reported.

Each replicate is scored by three estimators on identical data: the registered
rank residualization with an intercept; Pearson residualization with an
intercept followed by ranking; and the no-intercept rank residualization, which
is present only to reproduce the software defect in the atlas project's
regression test and is labeled as a software-diagnostic arm wherever it appears.
2,000 replicates per null.

### H3

Each generated cohort reuses the real target labels, the real per-drug context
counts `K_c`, and the real unique direction `u_t` of each target. For a drug `d`
of target `t` with `K_d` contexts, in 978 dimensions:

```
q_d    ~ Uniform(0.3, 3.0) * a        the drug's quality, a the planted scale
j_dk   ~ Normal(0, 1)                 one scalar per context
eps_dk ~ Normal(0, 1) per coordinate  isotropic, unit variance per coordinate
S_dk   = q_d * u_t + j_dk * q_d * u_t + eps_dk
```

So every coordinate of the noise has unit variance and a noise vector has
expected norm `sqrt(978)`; the planted component lives on one axis and the same
`q_d` controls both mean alignment and context dispersion. Nothing is centered or
rescaled after construction. `P`, `E`, `M_delta` and `K` are then computed from
`S` exactly as the analysis computes them.

The planted scale `a` is not an effect size. Three levels are frozen by a pilot
that calibrates `a` to a target contribution on the analysis scale, where the
contribution is the difference between the cohort's observed partial correlation
and the median of its own unique-target permutation null. The pilot runs on seeds
drawn from a separate stream that the registered runs never use, and its only
output is the three values of `a` that reach contributions of approximately 0.05,
0.10 and 0.20. The pilot's detection outcomes are not inspected.

At each level, 200 cohorts are generated independently and each is tested with
1,999 unique-target permutations, giving a p-value resolution of 1/2000. Reported
per level: the proportion of cohorts with `p_perm < 0.01`, the proportion with a
contribution above 0.0928, the proportion with both, a 95% Wilson interval for
each proportion, and the achieved contribution distribution.

**Simulation-validity gate, fixed before any detection outcome is read.** At each
level, across the 200 cohorts, report the Spearman correlation between `q_d` and
`M_delta`, and the proportion of the unadjusted `P`-`E` association removed by
adjusting for `M_delta` and `K`. The level is void if the median absolute
Spearman correlation between `q_d` and `M_delta` exceeds 0.30, or if adjustment
removes more than half of the unadjusted association, because the generator would
then be planting a magnitude effect that the covariates absorb rather than a
target-specific one. Neither threshold is revised after a detection outcome is
seen, and the generator is not retuned to pass the gate.

Each arm runs once at each level. No tolerance is adjusted, no covariate is
added, and no statistic is chosen after seeing a result.

## Maximum claim under this registration

The paper may say which geometric property of the reference holds the S2-T null
up, what the registered estimator's median was under two stated conditional
independence models at this sample size, and with what probability S2-T would
have detected a planted target-specific contribution of a stated size under one
stated generator.

It may not say that the association contains no target-specific component. A
failing test is absence of support, and a power diagnostic bounds which
alternatives were refusable, so the strongest licensed form is: under the
prespecified generator, S2-T detected contributions of at least X partial rho
with Y% probability, and the observed statistic did not exceed its reassignment
null, so the data give no evidence for a target-specific contribution of that
detectable form and scale.

Nothing here licenses rerunning S2-T, changing its criterion, or reporting any of
these results as confirmatory. S1 and S1-T stand as registered: each satisfied
its prespecified criterion, and neither is discriminating evidence for target
specificity, because the observed association was typical of the unique-target
reassignment null.
