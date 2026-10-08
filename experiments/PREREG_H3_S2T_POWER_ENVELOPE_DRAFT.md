# What could the unique-target permutation have detected, and what holds its null up at 0.37?

**Date:** 2026-10-08
**Status:** DRAFT. Not frozen. A code review of H1 and H2 is outstanding.
**Commit SHA:** pending
**Kind:** post-hoc diagnostic. No registered statistic is recomputed, no
registered criterion is changed, and no outcome here can move the gate at
`results/03c_h3_sensitivity/h3_sensitivity_results.json`, which is closed:
`h3_interpretation_survives: false`.
**Relates to:** `PREREG_H3_S1S3_CORRECTED_BASELINE.md`, frozen at 7f57136, whose
S2-T failed; `results/03j_reference_axis_geometry/reference_axis_geometry.json`,
the descriptive geometry of the reference.

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
ratio of 13.0, and a leading eigenvalue holding 26% of the trace. Target cluster
sizes are known: 131 of 258 targets carry one drug, the largest carries 37, and
the Kish effective cluster count is 82.

The detection threshold is known. S2-T refuses at p < 0.01, which on the recorded
null is an observed partial rho of 0.4599 or above, a target-specific
contribution of 0.0928. The observed 0.3722 sits 0.0877 below it.

One anticipation is on record before H2 runs: the atlas analysis
(`direction-instability-atlas`, `paper/atlas_paper_v7_plosone.tex` line 411)
residualizes with Pearson and then ranks, and states that rank-then-residualize
produces spurious partial correlations of 0.05 to 0.25 under the null at large n.
The registration frozen at f288507 registered rank-then-residualize. The two
cannot both be right about this data.

## Hypotheses

**H1.** Replacing each target's direction with an isotropic random direction
leaves the null median where the unique-target permutation put it, 0.3671.

**H2.** Under independence at n = 795 with the cohort's own covariates, the
registered partial correlation is centered on zero.

**H3.** A cohort carrying the real directions, the real cluster sizes and the real
context counts, in which the association is target-specific by construction,
attenuates by more than 0.0928 under the unique-target permutation.

H3 carries the diagnostic. H1 and H2 characterize the floor; only H3 says whether
the test could have fired, and a failing S2-T is evidence about the biology only
if it could.

**H3 is the outcome that strengthens the registered finding.** If a constructed
target-specific association attenuates well past the threshold on this reference,
then reading 3 is refuted and the failing S2-T means what it appears to mean.

## Inference criteria

| hypothesis | holds when |
|---|---|
| H1 | the random-axis null median lies within 0.02 of 0.3671 |
| H2 | the median registered partial rho lies within 0.02 of 0 under both nulls, over 2,000 replicates each |
| H3 | p_perm < 0.01 and the attenuation exceeds 0.0928, at every planted strength |

H1 failing means the shRNA directions, not axis projection in general, hold the
floor up. H2 failing means S1's registered floor of 0.20 was measured against a
shifted null; S1's verdict stands as registered and the shift is reported as a
limitation of that criterion, not as a correction to it. H3 failing means the
paper states that this reference cannot discriminate targets well enough to test
target specificity, and makes no claim either way about whether the association
is target-specific.

**All three are void if the identity permutation does not reproduce the recorded
observed statistic**, which is the check that the diagnostic is running on the
registered cohort.

## Implementation, fixed before it runs

Every arm imports `partial_spearman`, `permuted_statistic`, `load_bundle` and
`target_representatives` from `experiments/03c_h3_sensitivity.py` at the sealed
hash 6290281e, and `unique_target_permutations` from `geometry.inference`, so the
diagnostic exercises the code the run used rather than a restatement of it. The
covariates stay the registered pair, mean pairwise difference norm and context
count. Permutation counts are 10,000 for H1 and H3 and 2,000 replicates for each of H2's two nulls.
Seeds are declared in the script and recorded in its output.

H1 draws, per replicate, 258 isotropic unit directions in the 978 landmark
dimensions, assigns each to its target, propagates to that target's drugs, and
recomputes `P`, `E` and the partial correlation. Covariates are recomputed, not
held, because a random axis changes nothing they depend on and recomputing them
removes the question.

H2 needs a null in which `P` and `E` are conditionally independent given the
covariates **while each keeps its dependence on them**, because that dependence
is what the bias feeds on: a free permutation of `E` across drugs would destroy
it and would under-detect the bias by construction. Two nulls run, both reported.

The parametric null follows the atlas test's generating model on this cohort's
own covariate matrix. `M_delta` and `K` are standardized and the recorded `P` and
`E` are each regressed on them by ordinary least squares in the raw scale; the
fitted coefficients and residual scales are then used to draw `P*` and `E*` with
independent Gaussian noise, so the only association between them runs through the
covariates. The semiparametric null keeps the real marginals instead: the raw
residuals of `E` on the covariates are permuted across drugs and added back to
the fitted values, leaving `E`'s covariate dependence intact and its link to `P`
broken.

Each replicate is scored twice, once by the registered rank-then-residualize
procedure and once by Pearson residualization followed by ranking, so the two
procedures are compared on one set of draws rather than on separate simulations.

H3 builds signatures rather than resampling them. Each drug of target `t` is
given a quality drawn uniformly on (0.3, 3.0) times a planted strength; its mean
sits that far along `u_t` and its contexts jitter along `u_t` by the same factor,
with isotropic noise of unit scale so the mean pairwise difference norm stays
effectively constant and the magnitude covariate cannot absorb the signal. The
real cluster sizes and the real per-drug context counts are reused. Planted
strengths are 0.5, 1.0 and 2.0, all three stated here, all three reported, and
no other strength is tried.

Three arms run once each. No tolerance is adjusted, no covariate is added, and no
statistic is chosen after seeing a result.

## Maximum claim under this registration

The paper may say what holds the S2-T null up, whether the registered estimator
is centered under independence at this sample size, and whether a target-specific
association of any planted strength would have been detected on this reference.
It may not say that the association is or is not target-specific in the data:
that question was settled by S2-T and settled against the claim. Nothing here
licenses rerunning S2-T, changing its criterion, or reporting any of these three
results as confirmatory.
