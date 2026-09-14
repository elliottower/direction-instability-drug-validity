# Pre-registration: Does the phenotype-projection result survive magnitude, coverage and shared-axis controls?

**Date:** 2026-09-13
**Status:** PRE-REGISTERED (frozen before the LINCS matrices are rebuilt or any adjusted statistic is computed)
**Commit SHA:** f288507
**Amends:** H3 in `PREREGISTRATION.md` (commit `249abaf`) and Entry 1 of
`DEVIATION_LOG.md`. This registers three sensitivity analyses of an existing
registered hypothesis. It registers no new confirmatory hypothesis and cannot
raise H3 above the evidential status its original registration supports. The
sensitivity analyses are labeled S1–S3 to keep them distinct from the paper's
registered H1–H3.

## Rationale

H3 is the one hypothesis in the paper that meets its registered criterion and
survives multiplicity correction. Its statistic is computed by
`geometry/direction_instability.py::phenotype_projected_instability`, which for
each pair of contexts projects the *difference* between two signatures onto the
unit on-target direction, takes the absolute value, and averages over pairs:

    P_c = C(K,2)^-1 * sum_{i<j} |(s_ci - s_cj) . u_c|

Three properties of that computation are not addressed by the original
registration. Each term scales with `||s_ci - s_cj||`, so `P` carries the
magnitude of the pairwise differences, while the enrichment it is correlated
against, `E_c = cos^2(mean signature, u_c)`, is magnitude-free. `E` is a
signless drug-target axis alignment, historically labeled on-target enrichment:
a squared cosine cannot separate alignment from anti-alignment. Both `P` and `E`
are constructed from the same shRNA consensus direction `u_c`, so an association
is available to them algebraically whether or not each drug is matched to its
own target. And cell-line coverage is already known to correlate with both.

## Foreknowledge of data or evidence

Recomputed from `results/03_phenotype_projection/phenotype_projection_results.json`
on 2026-09-13, recorded in `results/03b_h3_sensitivity/artifact_level_checks.json`:

| quantity | value |
|---|---|
| Spearman rho(P, E) | +0.3756 (p = 4.8e-28) |
| Spearman rho(D, E), raw direction instability | -0.0433 (p = 0.22) |
| Spearman rho(P, K), K = number of cell lines | +0.2555 |
| Spearman rho(E, K) | +0.2149 |
| partial rho(P, E) given K | +0.3397 |
| Spearman rho(P, D) | -0.2475 |
| range of P | 0.42 to 22.50 |

The registered H3 criteria were |rho| > 0.3 for the projected score and
|rho| < 0.15 for the raw score; both were met.

Also known were the reported removal sensitivities (HDAC-class removal
rho = 0.3759 at n = 775; removing the 20 lowest-raw-D drugs rho = 0.387;
removing the 20 highest-P drugs rho = 0.350), the bootstrap 95% interval on the
full set [0.315, 0.433], and the failed CRISPRi convergent-validity check
(rho = -0.12, p = 0.16, n = 131). None of those used signature-magnitude
adjustment or a target-assignment permutation.

Never computed: any signature magnitude for these drugs, any magnitude-adjusted
association, and any permutation that reassigns target directions. The LINCS
signature matrices were not retained and will be rebuilt from GEO GSE92742.

## Research questions or hypotheses

**S1.** The association between target-axis projected dispersion and signless
drug-target axis alignment persists after rank-based adjustment for
pairwise-difference magnitude and cell-line coverage. That quantity is
historically labeled on-target enrichment, and the deposited field keeps that
name.

**S2.** That adjusted association exceeds the one produced when target
assignments are permuted across drugs while both quantities continue to be
constructed from the same target-direction axis.

**S3.** The raw score remains practically equivalent to zero under the same
adjustment, preserving the contrast the original registration drew between the
projected and raw scores.

S2 addresses the central shared-axis concern. Without it, S1 is consistent with an association that
sharing the axis `u` produces for any drug set, whether or not each drug is
paired with its own target. S3 is the component that can return a result in the
paper's favor; each of the three can also weaken H3, and none is an
only-strengthening test.

## Inference criteria

Let `M_delta,c = C(K,2)^-1 * sum_{i<j} ||s_ci - s_cj||_2`, the mean pairwise
signature-difference norm, matched to the scale `P` actually carries. Let `K_c`
be the number of cell lines. Partial Spearman correlation is computed by one
frozen procedure, not by a library's `partial_corr`:

1. Convert each of `P`, `D`, `E`, `M_delta`, `K` to midranks over the complete-case cohort.
2. Regress ranked `P` on an intercept, ranked `M_delta` and ranked `K` by ordinary least squares.
3. Regress ranked `E` on the same covariates.
4. Take the Pearson correlation of the two residual vectors.

All ranks use average ranks for ties. The intercept is included explicitly. One
fixed complete-case cohort and covariate matrix serve S1, S2 and S3.

For S3, substitute ranked `D` for ranked `P`.

| hypothesis | holds when |
|---|---|
| S1 | partial rho(P, E given M_delta, K) >= 0.20 **and** the lower bound of a 95% drug-level bootstrap interval exceeds 0 |
| S2 | p_perm < 0.01, where p_perm = (1 + #{rho_b >= rho_obs}) / 10001 over 10,000 target permutations |
| S3 | the 90% bootstrap interval for partial rho(D, E given M_delta, K) lies entirely within (-0.15, 0.15) |

0.20 is the smallest effect treated as substantively meaningful here. It is not
derived from a power calculation.

Bootstrap intervals are percentile intervals. Within each of 10,000 drug-level
resamples, every variable is reranked using average ranks, both residual
regressions are refitted with `numpy.linalg.lstsq`, and the residual correlation
is recomputed. S1 uses the 2.5th and 97.5th percentiles; S3 uses the 5th and
95th. Any nonfinite replicate invalidates the run and is investigated, never
silently discarded or redrawn.

Seeds: bootstrap 20260913, permutation 20260914.

Permutation: permute the complete target-direction vectors `u_c` across drug
records as a multiset, so duplicate targets stay duplicated; recompute both
`P_c` and `E_c` from the permuted direction; hold `D_c`, `M_delta,c` and `K_c`
fixed; rerank the permuted `P` and `E` before residualizing, and recompute the
adjusted partial correlation. The observed statistic is the same one S1 uses,
partial rho(P, E given M_delta, K). Fixed points are allowed and no derangement
is imposed; the permutation acts on the target-direction vector attached to each
drug record, not independently on labels and direction matrices. Report the
empirical p-value and the 2.5th, 50th and 97.5th percentiles of the null.

Mean signature norm is a declared secondary sensitivity, reported without a
criterion. `P/M_delta` is descriptive, carries no criterion, and is missing where `M_delta` is zero, with no imputation. No alternative
magnitude definition is tested and selected afterwards.

**S1, S2 and S3 are required jointly: H3's current interpretation survives
these sensitivity analyses only if all three hold.** The original registered
numerical criterion has already been met and does not retroactively become
unmet. Component outcomes are reported separately, but no component compensates
for the failure of another.
This is an intersection-union gate, so no multiplicity correction is applied
across them.

**All three are void unless the rebuild validates (below).**

## Rebuild validation

An aggregate correlation can reproduce while individual records are mismatched,
so validation is per drug. Let `C` be the intersection of deposited and rebuilt
drug identifiers. Numerical comparison happens only on `C`:

    no rebuilt identifier absent from the deposited artifact
    at least 700 common records
    one record per unique drug identifier
    identical target assignment on every common record
    identical n_celllines on every common record
    max |D_rebuilt - D_deposited| < 1e-6 on C
    max |P_rebuilt - P_deposited| < 1e-6 on C
    max |E_rebuilt - E_deposited| < 1e-6 on C

The deposited artifact holds 795 records, 795 unique drug identifiers and no
drug carrying more than one target annotation, so the resampling unit is the
drug and `n_records == n_unique_drug_ids` is asserted before anything runs. Were
that to fail on rebuild, the unit stays the drug and the run stops.

**The sensitivity analyses are not interpreted unless every condition above
holds.** The 700-record floor keeps at least 88% of the deposited cohort and
limits how far reconstruction attrition can change which drugs are analyzed; it
is not derived from a power calculation. Deposited and rebuilt correlations are
compared on `C` alone; reproducing the full-cohort 0.3756 after dropping records
is not expected. A rebuild failure voids the analyses, and the failure and its
diagnostics are still reported.

Pinned here before execution. The run verifies these rather than recording new
ones; a mismatch stops the run.

| input | pin |
|---|---|
| deposited H3 artifact | sha256 `65e5d10e272037987384f89e6208de478fa17f6b8fb946add893e5a24c2d80c4` |
| `GSE92742_Broad_LINCS_sig_info.txt.gz` | sha256 `19da29c0ee12ddf27f9698cd0da40beaff58657dcde9d382aae068737e831299` |
| `lincs_shrna_siginfo.csv.gz` | sha256 `bd396fa0e1a2f00c1b5f2c8d2b35f9a056f5e5353382475655869038037ec014` |
| `frozen_drug_labels.json` | sha256 `e6d73384c5bbc501ee15f621aaa258f02f38ee0d7a9603b6c820bec3e0489d52` |
| loader | `drug-perturbation-geometry@1dc20a2`, `data/lincs_loader.py` |
| statistic and consensus construction | `direction-instability-drug-validity@ccb6223`, `geometry/direction_instability.py` |
| landmark genes | 978 |
| GCTX | `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz` from the GEO GSE92742 supplement |

The GCTX is not held locally and its checksum cannot be computed before
download. Its sha256 is recorded on first retrieval, written to the output
alongside the feature-order hash, and verified on every subsequent run. That one
input is therefore pinned by accession and filename at freeze and by checksum
from first download, not by checksum at freeze. The run also records NumPy and
SciPy versions and a hash of the final common-cohort identifier list.

## Sample size

795 drugs with an annotated target, a target with at least 3 shRNA hairpins, and
at least 5 cell lines. This is fixed by the deposited artifact and is not a
design choice.

## Maximum claim under this registration

If all three hold, the paper may say that the association between target-axis
projected cross-context dispersion and signless drug-target axis alignment
persists after
rank-based adjustment for pairwise-difference magnitude and cell-line coverage,
and exceeds the association generated when drug-target assignments are permuted
while both quantities continue to share the same target-direction axis, and
remains separated from a raw-score association practically equivalent to zero.

Regardless of outcome, the statistic remains a signless projection of pairwise
differences onto the target axis. The paper may **not** say that phenotype
projection isolates on-target consistency, that the consistency points toward
the target, that the statistic decomposes `D` into on-target and off-target
components, or that a high value means stronger mechanism conservation. Several
of those claims are in the present manuscript and are withdrawn irrespective of
these results.

If `C` holds fewer than all 795 deposited records, every estimate and claim
here is restricted to the rebuilt common cohort `C`, and the identifiers missing
from the rebuild are written to a provenance file rather than replaced.

If S1 fails: the adjusted association fell below the prespecified minimum
effect, so robustness to response magnitude was not established. If S2 fails:
correct drug-target matching could not be distinguished from coupling induced by
constructing both quantities from the same axis, and the abstract and scorecard
claim that phenotype projection separates on-target from off-target consistency
is withdrawn. If S3 fails: the raw-score contrast was not practically equivalent
to zero after adjustment, weakening the claimed separation between the raw and
projected quantities.

## A signed variant is not registered here

Dropping the absolute value would not recover movement toward versus away from a
target. For unordered context pairs the sign of `(s_i - s_j) . u` flips when the
two contexts are exchanged, so its mean depends on the arbitrary order of
contexts. A signed alignment statistic would need an ordered reference and
target-action annotations distinguishing inhibitors from agonists, and the
current enrichment is a squared cosine, which also discards direction. That is a
different estimand and is not appended to this registration.

## Design and procedure

Rebuild `lincs_subset.npz` and `lincs_shrna.npz` from GEO GSE92742 Level 5
(MODZ, 978 landmark genes) with the existing loader, on one Modal CPU worker
with a volume, checkpointing each shard and committing the volume. Recompute
`P`, `E` and `D` with the unmodified functions in
`geometry/direction_instability.py`, validate as above, then run S1-S3 once.

**Existing data:** yes. The GCTX input is pinned at freeze by GEO accession and
exact filename; its checksum is recorded at first retrieval and required to match
on every subsequent execution. All locally held metadata and artifacts are pinned
by checksum before execution.
**Data collection:** N/A — no new data are collected.
**Blinding:** no blinding is possible; the original H3 result and the
artifact-level associations above are known. No magnitude-adjusted or
target-permuted quantity has been computed or inspected.
**Missing data:** drugs that fail to rebuild are counted and reported; none are replaced.
**Outliers:** none removed.
**Exploratory analyses:** any analysis not named here is exploratory and is labeled as such.
