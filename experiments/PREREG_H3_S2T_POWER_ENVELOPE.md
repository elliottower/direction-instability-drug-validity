# What could the unique-target permutation have detected, and what holds its null up at 0.37?

**Date:** 2026-10-08
**Status:** FROZEN before the diagnostic ran. No H1, H2 or H3 statistic has been
computed, the pilot has not been run, and no rotation, null or cohort exists on
disk. Supersedes the drafts `PREREG_H3_S2T_POWER_ENVELOPE_DRAFT.md`,
`_v1a_DRAFT.md` and `_v1b_DRAFT.md`, which are kept unaltered: v1b's H3 criterion
was unattainable by construction, because the planted level was calibrated so the
median cohort contribution sat at the detection threshold and detection then
required a cohort to exceed that same threshold, which at most half of them can.
Detection here is the permutation rejection alone and the contribution is a
separate calibration gate.
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

Reading 2 is named here as shared-axis and reference-signature coupling, because
several mechanisms sit inside it: direct coupling through the shared substituted
axis, anisotropic signature covariance, overlap between the reference subspace
and shared drug-response programs, correlated geometry among the directions,
cluster imbalance interacting with a nonlinear statistic, and residual covariate
misspecification. H1a to H1c begin to separate them and do not settle them.

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
geometry they keep. Contrasts are interpreted according to the properties
preserved and randomized by each prespecified ensemble; each arm changes a bundle
of properties whose effects may interact, and no contrast need identify a unique
causal property. Each is an equivalence hypothesis against a margin fixed before
any arm runs, not a failure-to-reject rule.

**H1a.** The floor under 258 independent isotropic unit directions is
practically equivalent to the recorded floor. This destroys pairwise similarity,
spectrum, subspace and ambient orientation together.

**H1b.** The floor under a direction set with the same Gram matrix and a
Haar-random ambient orientation is practically equivalent to the recorded floor.
Every pairwise cosine and every singular value survive; the empirical row
subspace is replaced.

**H1c.** The floor under a Haar rotation inside the reference's own numerical row
span is practically equivalent to the recorded floor. The Gram matrix and the
subspace both survive; the coordinate alignment of each labeled direction with
the drug-signature space does not.

**H2a.** The registered estimator is practically equivalent to zero under a
linear-Gaussian conditional independence model fitted to this cohort's covariate
relationships at n = 795.

**H2b.** The same, under a drug-level Freedman-Lane residual permutation.

**H2c.** The same, under a target-level cluster-wild residual bootstrap.

H2a to H2c are separate hypotheses. One holding does not compensate for another
failing, and disagreement among them is a result about residual exchangeability
rather than a tie to be broken.

**H3.** At a generator scale whose median achieved contribution is calibrated to
the recorded excess of 0.09282365, S2-T rejects at `p_perm < 0.01` in at least
80% of independently generated cohorts carrying the real directions, the real
cluster sizes and the real context counts.

H3 carries the diagnostic. H1a to H1c say what holds the floor up; only H3 says
which alternatives the test was capable of refusing, and a failing S2-T is
evidence about the data only for alternatives it could have refused.

**H3 is the outcome that strengthens the registered finding.** A passing H3
disfavors reading 3, the masking explanation, for target-specific alternatives of
the prespecified generated form and scale. It does not establish power against
all biologically plausible target-specific alternatives, and it does not refute
every version of reading 3.

## Inference criteria

Equivalence is judged against a margin of 0.02, carried over unchanged from the
draft that proposed it before any arm existed. For an arm with median `m` and a
comparison value `v`, define `delta = m - v` and take a 90% interval for `delta`.

| verdict | when |
|---|---|
| equivalent | the whole 90% interval for `delta` lies inside [-0.02, +0.02] |
| different | the whole 90% interval lies outside [-0.02, +0.02], on one side |
| inconclusive | neither, which is reported as inconclusive and not as either verdict |

| hypothesis | comparison value | holds when |
|---|---|---|
| H1a, H1b, H1c | the recorded unique-target null median, 0.36711689048130913 | the arm's verdict is equivalent |
| H2a, H2b, H2c | zero | the arm's verdict is equivalent |
| H3 | — | at the calibrated level, the estimated probability of `p_perm < 0.01` is at least 0.80 |

The two intervals are built differently, because H1 compares against an
estimated median and H2 against a fixed constant.

For **H1**, 10,000 bootstrap iterations, each resampling the arm's 2,000
replicate statistics with replacement and, independently, the 10,000 recorded
null draws with replacement, taking the difference of the two medians, and
reading the 5th and 95th percentiles. The arms are not paired with the recorded
permutations and are not treated as paired.

For **H2**, 10,000 bootstrap iterations resampling only that arm's 2,000
simulated estimator statistics, giving an interval for their median, which is
then compared against fixed zero. The recorded null draws play no part.

Whichever interval results controls the verdict. No arm gains replicates after
its interval is seen; an interval too wide to decide returns inconclusive, which
is a registered outcome.

For **H3** the primary detection event is `p_perm < 0.01` alone. Its 95% Wilson
interval is reported, and the result states both that the lower bound rules out
power below its own value and that it does not establish power above 0.80. The
proportion of cohorts whose achieved contribution exceeds 0.09282365 is reported
descriptively and is not part of the detection event: the level is calibrated so
that proportion sits near one half, which is a property of the calibration and
not a measure of sensitivity.

H3 is also run at planted contributions of 0.05 and 0.20 and the three points are
reported as a curve. Only the 0.09282365 level carries the criterion.

**Reading the H1 arms.** Read jointly, each limited to the property it destroys.
The arms sample different ensembles and may behave non-monotonically, so the
patterns below are not exhaustive and any other pattern is reported as it falls.

| pattern | licensed reading |
|---|---|
| H1a equivalent | the floor is reproducible under the prespecified independent-isotropic axis ensemble |
| H1a different, H1b equivalent | preserved inter-axis Gram geometry is sufficient under a random ambient orientation |
| H1b different, H1c equivalent | preserving the empirical reference row subspace is sufficient; exact within-subspace orientation is not required |
| H1c different | the observed orientation within the empirical subspace is atypical relative to random in-span rotations |
| any arm inconclusive | that contrast is reported as inconclusive and supports no reading |

No arm licenses a universal or causal claim. The phrases "any axis
construction", "geometry alone" and "orientation causes the floor" are not used.

**All hypotheses are void if the identity permutation does not reproduce the
recorded observed statistic of 0.3721909618253966**, which is the check that the
diagnostic runs on the registered cohort.

## Implementation, fixed before it runs

Every arm imports `partial_spearman`, `permuted_statistic`, `load_bundle` and
`target_representatives` from `experiments/03c_h3_sensitivity.py` at the sealed
hash 6290281e, and `unique_target_permutations` from `geometry.inference`, so the
diagnostic exercises the code the run used rather than a restatement of it. The
covariates stay the registered pair, mean pairwise difference norm `M_delta` and
context count `K`, and are recomputed from the signatures in every arm rather
than carried.

Random numbers come from two disjoint streams. `numpy.random.SeedSequence`
entropy 20261008002 spawns, in order, the seeds for H1a, H1b, H1c, H2a, H2b, H2c
and the H3 registered runs. Entropy 20261008001 spawns the H3 pilot's seeds and
is used nowhere else. Every spawned seed is recorded in the output.

### H1a to H1c

Write the real direction matrix as `U = A S Vt` by singular value decomposition,
`U` being 258 by 978 with unit rows. The numerical rank `r` is the count of
singular values exceeding `S.max() * 1e-10`, a tolerance fixed here rather than
left to a library default. `r` and all 258 singular values are recorded in the
output. The 140-component figure that carries 90% of the Gram trace is a variance
summary and is not used as a rank.

- **H1a** draws 258 independent standard normal vectors in 978 dimensions and
  normalizes each to unit length.
- **H1b** replaces `Vt` with a Haar-distributed `r` by 978 orthonormal frame `W`,
  giving `U* = A S W`. Since `W W^T = I`, `U* U*^T = A S^2 A^T = U U^T`, so every
  pairwise cosine, every row norm and the whole spectrum survive.
- **H1c** takes an orthonormal basis `Q`, 978 by `r`, of the full numerical row
  span of `U`, draws an `r` by `r` Haar orthogonal `R`, and sets
  `U* = U Q R Q^T`, which stays in that span and preserves `U U^T`.

Haar matrices come from the QR decomposition of a standard normal matrix with the
sign correction on the diagonal of `R` that makes the distribution Haar.

Three computational integrity assertions run on every H1b and H1c replicate, with
tolerances frozen here at 1e-9 and no outcome attached: the maximum absolute
element difference between `U* U*^T` and `U U^T`; the maximum absolute deviation
of the row norms of `U*` from one; and, for H1c only, the maximum absolute
residual of `U*` after projection onto the row span of `U`. A replicate that
exceeds any tolerance refuses rather than reporting.

Each arm runs 2,000 replicate direction sets. Within a replicate, each target receives its
substituted direction, every drug of that target receives its target's direction,
and `P`, `E` and the registered partial correlation are recomputed from the
signatures. The arm's statistic is the median over its 2,000 replicates.

### H2a to H2c

All three nulls run, all three are reported, and none is selected on the basis of
a diagnostic computed from the data.

**H2a, parametric.** Standardize `M_delta` and `K`, regress the recorded `P` and
`E` on them by ordinary least squares in the raw scale, and draw `P*` and `E*`
from the fitted means with independent Gaussian noise at the fitted residual
scales. It assumes linear conditional means, homoscedastic Gaussian residuals and
no target-cluster dependence, and can place values outside the natural ranges of
`P` and `E`. It licenses a statement about the estimator's median under this
linear-Gaussian null fitted to the cohort's first-order covariate relationships,
and not about unbiasedness under independence.

**H2b, Freedman-Lane at the drug level.** Fit `E` on the covariates in the raw
scale, permute the raw residuals across drugs, add them back to the per-drug
fitted values. It preserves the fitted linear covariate relationship and the
empirical residual multiset. It does not preserve the marginal distribution of
`E`, and its validity rests on residual exchangeability across drugs.

**H2c, cluster-wild at the target level.** Fit the same nuisance model, keep each
target's residual vector whole, multiply every residual in a target by one shared
Rademacher weight drawn independently per target, and add the weighted residuals
to the fitted values. This retains unequal cluster sizes, which range from one to
37 drugs, and within-target residual dependence. It is a cluster wild bootstrap
and is named as one, not as a permutation.

Each replicate under each null is scored by three estimators on identical data:
the registered rank residualization with an intercept; Pearson residualization
with an intercept followed by ranking; and the no-intercept rank residualization,
present only to reproduce the software defect in the atlas project's regression
test and labeled as a software-diagnostic arm wherever it appears. 2,000
replicates per null.

Four quantities are reported descriptively with the arms, with no threshold and
no branch taken on any of them.

The Spearman correlation between the squared nuisance residual and the fitted
value. The target-level intraclass correlation of the nuisance residuals,
computed as the one-way random-effects ICC(1,1) by the ANOVA moment estimator for
unbalanced groups, over all 258 targets with the standard `n_0` correction; the
131 singleton targets contribute to the between-group sum of squares and carry no
within-group degrees of freedom, and their count travels with the value. The
marginal deciles of `E` before and after each resampling, which for H2c describe
the resampling distribution rather than a preserved empirical distribution, since
Rademacher weights symmetrize the residuals. The skewness of the nuisance
residuals before and after, for the same reason.

### H3

Each generated cohort reuses the real target labels, the real per-drug context
counts `K_d`, and the real unique direction `u_t` of each target. For a drug `d`
of target `t`, in 978 dimensions:

```
q_d    ~ Uniform(0.3, 3.0) * a        the drug's quality, a the planted scale
j_dk   ~ Normal(0, 1)                 one scalar per context k
eps_dk ~ Normal(0, 1) per coordinate  isotropic, unit variance per coordinate
S_dk   = q_d * u_t + j_dk * q_d * u_t + eps_dk
```

Every coordinate of the noise has unit variance, so a noise vector has expected
norm `sqrt(978)` while the planted component lives on one axis; calibration may
therefore need large `a` and may produce signature magnitudes unlike the real
cohort. Nothing is centered or rescaled after construction. The same `q_d`
controls mean alignment and context dispersion, which is what makes the
alternative target-specific. `P`, `E`, `M_delta` and `K` are computed from `S`
exactly as the analysis computes them.

The achieved contribution of a cohort is its observed registered partial
correlation minus the median of its own unique-target permutation null. The
planted scale `a` is not an effect size and is calibrated to the contribution.

**Pilot calibration, frozen here as an algorithm.** Candidate scales are twelve
points log-spaced from 0.1 to 10. Each candidate runs 20 pilot cohorts at 499
permutations each, and its achieved value is the median contribution over those
20 cohorts. For each target `c` in 0.05, 0.09282365 and 0.20:

```
  errors    = achieved[i] - c  for each candidate i, in ascending scale order
  brackets  = adjacent pairs (i, i+1) whose errors have opposite signs,
              or where either error is exactly zero
  if brackets is empty:
      take the candidate with the smallest |error|, ties to the lower scale;
      perform no extrapolation; record no-bracket for this target and stop
  pick      = the bracket of smallest width in log scale,
              ties to the one with the lower scale
  repeat exactly 3 times:
      a_new        = linear interpolation of the bracket endpoints in log scale
                     against their errors
      achieved_new = median contribution of 20 cohorts at a_new,
                     499 permutations each
      replace whichever bracket endpoint lies on the same side of c as
                     achieved_new with (a_new, achieved_new)
  selected  = the attempted scale, grid or confirmation, with the smallest
              |achieved - c|, ties to the lower scale
```

Nonmonotonic achieved medians need no special case: a bracket is a sign change in
the error, several brackets resolve by log width then by lower scale, and the
final selection is on smallest absolute error over every attempt rather than on
the last confirmation, so an overshoot or undershoot cannot be selected over a
better earlier attempt. The tolerance of 0.01 does not stop the search; it is the
calibration gate below, applied to the registered run. Every attempted scale and
its achieved contribution is published, grid and confirmation alike, not only the
selected ones. The pilot computes permutations only to obtain each cohort's null
median for the contribution; pilot p-values and pilot detection indicators are
not used to select a scale.

**Registered runs.** The pipeline executes once; within that execution it
generates exactly 200 cohorts at each of the three frozen scales, each tested
with 1,999 unique-target permutations, giving a p-value resolution of 1/2000.
Reported per level: the proportion of cohorts rejecting at `p_perm < 0.01`,
which is the detection event, with its 95% Wilson interval; the proportion whose
achieved contribution exceeds 0.09282365, descriptively and with its own
interval; and the distribution of achieved contributions.

**Calibration gate.** A level supports a power claim only if the median achieved
contribution across its 200 registered cohorts lies within 0.01 of its target. A
level outside that tolerance keeps its full registered output, including its
rejection proportion and Wilson interval, is labeled miscalibrated, and yields no
H3 power verdict. A miscalibrated level is never dropped from the report and its
nearest scale is never treated as representing the target: either would hide a
failure of the simulation design. This is the only gate, and it is on calibration
rather than on any property of the detection outcome.

Reported descriptively at every level, with no threshold attached: the Spearman
correlation between `q_d` and `M_delta`, the unadjusted and adjusted planted
association and the proportion of the unadjusted association removed by
adjustment, and the simulated distributions of `P`, `E` and `M_delta` against
their real-data distributions.

No tolerance is adjusted, no covariate is added, and no statistic is chosen after
seeing a result.

## Maximum claim under this registration

The paper may say which geometric property of the reference holds the S2-T null
up, what the registered estimator's median was under three stated conditional
independence models at this sample size, and with what probability S2-T rejects
at `p_perm < 0.01` under one stated generator whose median achieved contribution
is calibrated to 0.09282365.

It may not say that the association contains no target-specific component. A
failing test is absence of support, and a power diagnostic bounds which
alternatives were refusable, so the strongest licensed form is: under the
prespecified generator calibrated to a median contribution of X partial rho, S2-T
rejected in Y% of cohorts, and the observed statistic did not exceed its
reassignment null, so the data give no evidence for a target-specific
contribution of that detectable form and scale.

The floor is called shared-axis and reference-signature coupling, which names a
family of mechanisms rather than a proved algebraic cause, until H1a to H1c
report.

Nothing here licenses rerunning S2-T, changing its criterion, or reporting any of
these results as confirmatory. S1 and S1-T stand as registered: each satisfied
its prespecified criterion, and neither is discriminating evidence for target
specificity, because the observed association was typical of the unique-target
reassignment null.
