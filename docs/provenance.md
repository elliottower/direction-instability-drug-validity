# Provenance: the gene-axis defect, and the audit of what it reaches

Deviation 11 and Deviation 12 in `DEVIATION_LOG.md` are the registration record.
This file carries the cross-repository consequence and the paper audit, which belong
to neither registration.

## Where the defect comes from

Both LINCS extractions this project pins, `lincs_subset.npz` and `lincs_shrna.npz`,
were written by `drug-perturbation-geometry`. Three sites there parse the pinned
GCTX with a requested row order and label the values with the request rather than
with the file, and nothing in that repository calls `reindex`. cmapPy ignores the
`rid` order it is given and returns the file's own, measured in
`results/03c_h3_sensitivity/reader_agreement.json`. Both files therefore carry the
same permutation of the gene axis: correct values, wrong labels.

`experiments/modal_h3_rebuild.py:218` is the one line that avoids it, which is why
the rebuild passes the gate and the extractions do not.
`drug-perturbation-geometry/docs/provenance.md` carries that side of the story.

## The audit of this project's papers

Every gene symbol in every paper version was classified by where its identity comes
from, because the two sources are affected differently:

- **a drug-target annotation**, from `frozen_drug_labels.json` — unaffected, the
  gene axis plays no part
- **an index into a signature column** — affected, the label is wrong

| paper | gene identities | verdict |
|---|---|---|
| `direction_instability_confound_audit_v1b.tex` (current) | PSMA3, PSMB1, PSMB6, PSMC4, TOP2A, SF3B1, SNRPG, DDX23 | **clean**: the rank-shift analysis ranks the 1,676 Perturb-seq knockdown genes (`:163`, `:429`), whose identities come from the h5ad's own labels, and the remainder are drug-target annotations |
| `direction_instability_zenodo.tex` (deposited) | the same set | **clean**, same reasoning |
| `direction_instability_bmcmrm_v1-v8`, `bioinfadv_v1`, `bmcbioinfo_v1`, `psb_v3`, `psb_v4` | target annotations only | **clean** |
| `direction_instability_v1-v8`, `psb.tex`, `psb_v2.tex` (superseded drafts) | SUV39H1, MYC, CDK6 named as "chromatin remodeling genes" of a conserved signature (`v8:170-173`) | **affected**, and the same claim as the drug-transport paper. Superseded, so no correction is owed; recorded so the sentence is not revived |

A scan by symbol alone is not sufficient and was not relied on. BIRC5, CDK6 and MYC
are both landmark genes and annotated drug targets, so membership in either list
says nothing about how a sentence uses them; each occurrence was read in context.

## Open

The rank-shift analysis names Perturb-seq genes, whose labels are sound. Whether the
*quantity* being ranked is computed in the LINCS landmark space, and so inherits the
permutation even though the names do not, has not been established. It does not
affect the gene identities and it may affect the shift magnitudes.

The CRISPRi arms are the live question and are tracked in Deviation 12:
`03d:604` builds the symbol list from the extraction's declared ids and uses it to
place five external references, and `03d:950` selects R7f's removed coordinate by
gene name.

## 2026-10-06 — gate v2, after review

**What the v1 gate actually reported.** Gate v1 failed on both registered cohort
invocations. Its production compound comparison passed on both cohorts at max
3.12e-06. Its shRNA comparison passed on the 795 cohort (0.0 signatures, 8.14e-07
consensus, 7.51e-08 direction, 1.11e-16 stored) and was **not evaluated** on the 812
cohort, because the gate treated the compound-complete target universe as
shRNA-complete. The retained legacy extraction failed the v1 equality predicate on
both cohorts at max 20.0. No v1 failure is erased or relabeled; both reports and the
v1 code hash (`9dc69bea8cf9e563…`) are preserved.

**Correction to how the gap is described.** The seventeen deposited records that
carry no shRNA quantities **omit the two keys** rather than holding null. There are
zero explicit nulls in the file. An earlier description of the same records said
`proj_shrna = null`, which is wrong about the file while being right about the
boundary: the same seventeen records, the same eight targets. The hardened
measurement surfaced it because `dict.get` had collapsed absence and null into one
case. `shrna_coverage_identity_v2.json` records both counts per field and the three
record shapes (681 shRNA only, 114 both, 17 CRISPRi only).

**The eligibility boundary, asserted rather than observed.** 258 targets and 795
drug-target records are shRNA-eligible; the complement in the 812-drug compound
cohort is exactly 17 records over EEF2, EIF2S1, FNTA, HCRTR1, KCNA10, RPL3, RPS2 and
TUBB. The rebuild's shRNA target set and the deposited analysis's valued target set
are the same set, and both hash to `3f6b38fd6e282620…`. The unvalued record set and
the pinned excluded manifest both hash to `573e0486a4e33bf8…`. Comparisons are on
sets and canonical hashes, never counts, because a set with one target dropped and
another substituted has the same count.

**Cohorts are now minimal manifests.** `registry/cohorts/` holds
`cohort_795_shrna_paired.json`, `cohort_812_compound.json`,
`shrna_eligible_targets.json` and `shrna_excluded_records.json`, each with canonical
pair, drug and target hashes, counts, uniqueness constraints and a schema version.
The 795 cohort had been `phenotype_projection_results.json`, which carries projected
instability and on-target enrichment, so cohort eligibility was coupled to derived
values.

**Gate v2.** Five named predicates — cohort identity, compound source reconstruction,
shRNA reconstruction on the pinned eligible set, the eligibility partition, and
legacy-artifact custody — each evaluated in its own `try` and recorded under its own
name, then one conjunction. v1 wrapped every comparison in a single `try`, so a
structural failure in one half left the half it had already computed unreported. The
predicate set is frozen in `RELEASE_PREDICATES` and a run that evaluates a different
set is refused, so the gate's scope cannot be narrowed by an invocation. Production
authorization and legacy custody are reported separately.

**The legacy predicate is positive in both directions.** It holds when the file's
bytes match the pinned sha256, the coordinate map is a total bijection over 978
landmarks, the map composed with its inverse is the identity, the map's canonical hash
is the pinned hash, applying the map reproduces the source exactly, and the declared
labels fail the production tolerance. Bare disagreement would not satisfy it, because
almost any corruption disagrees.

Reviewed in Perplexity round 17, which asked for exact-set rather than cardinality
semantics, record-level identity in the coverage measurement, a positive legacy
identity, independent per-predicate evaluation, and pinned cohort manifests. 86 tests
pass across the gate and discordance suites.

## 2026-10-06 — a registration edited after its freeze, and why it is reverted

`prereg check` reported `experiments/PREREGISTRATION_EXTENDED.md` as CHANGED against
its pin `f4bc40b` (2026-07-08): 64 lines added, none removed.

**What was added.** Commit `ccb6223` (2026-09-11), whose message describes a
terminology scrub, also inserted a section headed "Experiment R3: HDAC-Removal
Sensitivity and Drug Pair Search for H3" into that pinned file. The section carries
`**Commit SHA:** f4bc40b` and `**Pre-registered:** 2026-07-08`, and `f4bc40b`'s tree
does not contain it.

**Where R3 is really registered.** `PREREG_ROBUSTNESS_CHECKS.md`, section "R3:
HDAC-removal sensitivity for H3", whose content was fixed at `5bb6080` on 2026-07-08
and has not changed since. `experiments/r3_hdac_sensitivity_and_drug_pair.py` was
committed the same day in `7c0e8bf`, before any R3 output existed. So the plan and
its script do predate the results; the September insertion is a restatement in the
wrong file, not a registration created after the fact.

**How the restatement differs from the registration.** The July file sets the
criterion at rho > 0.20 for HDAC removal and for both top-20 removals. The inserted
text sets 0.3, and names analyses A1 (full-set reproduction), A4 (1000-resample
bootstrap interval) and A5 (a drug-pair search) that the July file does not state in
that form. The tightened threshold is a stricter test, and the measured values clear
both: rho_proj = 0.3756 on the full 795 drugs and 0.3759 after removing the 20 HDAC
drugs (`results/r3_hdac_sensitivity/r3_results.json`).

**What depends on it.** Nothing. No manuscript reports R3. In
`paper/direction_instability_confound_audit_v1c.tex` the HDAC-removal paragraph is
commented out, and the commented text cites the registered threshold as rho > 0.20.
If R3 is ever restored to a manuscript, 0.20 is the registered criterion.

**The correction.** `experiments/PREREGISTRATION_EXTENDED.md` is restored to its
content at `f4bc40b`, which is a pure deletion of the inserted 64 lines and removes
no terminology work: the insertion was additive, and its own text still uses the
pre-scrub `projected_bracket` wording, so it reads as an older draft committed late.
The file now matches its pin. No formal amendment is filed, because no reported
quantity rests on the 0.3 criterion or on A1, A4 and A5, and the substance of all
three is recorded here.

**Why the guard did not catch it.** The `repro` mod refuses an edit to a frozen
registration, but it intercepts tool calls, so it governs edits made through an agent
and not a hand edit or another editor. A pre-commit hook refusing a commit that
modifies any file carrying a `Commit SHA:` line, unless the same commit logs the
change, would close that path.

Whether A1, A4 and A5 are reported as confirmatory under the July freeze is Elliot's
call and is not decided here.

## 2026-10-08 — the rank-residualization bias was a missing intercept

**What was measured.** Four partial-correlation estimators on the generating model
of `direction-instability-atlas`, `tests/test_partial_spearman.py`,
`test_rank_then_residualize_is_biased`: Z standard normal, two outcomes each a
linear function of Z plus independent noise, so the true partial correlation is
zero. 200 replicates at n = 5000, seed 20261008,
`results/03k_partial_spearman_estimators/estimator_comparison.json`, run
`partial_spearman_estimators_2026_10_08`.

| estimator | mean | sd |
|---|---|---|
| rank residualization with an intercept, registered at f288507 | -0.0029 | 0.0147 |
| rank residualization without an intercept, the atlas comparator | +0.2646 | 0.1108 |
| Pearson residualization without an intercept, `power_analysis.partial_spearman` | -0.0022 | 0.0141 |
| Pearson residualization with an intercept | -0.0022 | 0.0141 |

**The mechanism.** The atlas comparator regresses rank vectors on ranked
covariates with no intercept column. Ranks are strictly positive with mean
(n+1)/2, so the fit is forced through the origin and cannot represent the level:
the residual vectors it returns have means of 357 and 350 where the ranks average
2500.5, and that surviving common offset is what correlates. The atlas reference
arm omits the intercept as well, and is unaffected only because the simulated raw
variables are already centered. The test therefore compares a no-intercept fit on
uncentered data against a no-intercept fit on centered data, which is not a
comparison between the two procedures.

**What it settles here.** Two things. The estimator registered for the H3
sensitivity analyses is unbiased under conditional independence, so the floor of
0.3671 under the unique-target permutation and 0.3910 under the drug-record
permutation is not an artifact of the estimator. It is a property of the
construction: `P` and `E` are built from the same axis, and that coupling
survives reassignment. And the cross-paper contradiction recorded earlier the
same day does not exist.

**What it leaves open.** Whether the estimator is centered on this cohort's own
covariate structure at n = 795 rather than on the test's model at n = 5000. That
is H2 of `experiments/PREREG_H3_S2T_POWER_ENVELOPE_v1a_DRAFT.md` and it is a
narrower question than the one H2 was first written for.

**What depends on it elsewhere.** `direction-instability-atlas` states the bias
claim in `paper/atlas_paper_v7_plosone.tex` Methods at lines 411-416, again at
line 1430, and in `PREREGISTRATION.md` at lines 125-130. The defect affects the
regression test and the categorical methodological claim, not the estimator used
by the six reported partial-correlation analyses. The categorical claim that
rank-then-residualize is biased, and that Pearson-then-rank is the correct
method, is what needs correcting. The correction itself is Elliot's call and is
not decided here.

**How it was found.** A code review of the draft diagnostic registration read
the atlas test's comparator and observed the missing intercept. The numbers above
were measured afterwards rather than taken from the review.

### Addendum, same day: the atlas audit, done rather than inferred

The entry above asserted that no atlas result moves. That was an inference from a
simulation on centered variables and did not establish it, so the code paths were
audited. Two implementations of `partial_spearman` exist in that repository.
`experiments/utils.py` line 97 forms `Z_int = column_stack([ones, Z])` and is the
one imported by all six experiments that report a partial correlation, 05, 07,
09, 10, 11 and 14, so every reported value uses the intercept and the
manuscript's methods description is accurate for its results.
`power_analysis.py` line 18 omits the intercept and is reached only by its own
`simulate_power_ci`, whose generator draws standard normal covariates and
mean-zero outcomes, and by the regression test. The manuscript's single power
sentence, 80% power for true effects at or above 0.35 at n = 76, comes from that
simulation, where the measurement above gives -0.0022 for Pearson residualization
with and without the intercept, so it is unaffected.

So the narrow claim holds, on an audit of every path rather than on analogy: the
six reported partial-correlation analyses are unaffected. The unrestricted phrase
"no atlas number moves" is not supported and is not used. The manuscript's power
sentence came from the no-intercept copy; under its centered generator the
estimator comparison found no detectable mean shift, but that power result was
not regenerated, so it is reported with that caveat or regenerated, which is
Elliot's call. What remains wrong is the categorical methods claim, and that the
function certified by the test is not the function that produced the results.
