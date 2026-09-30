# Pre-registration amendment: Does the corrected H3 association survive magnitude, coverage and target-level shared-axis controls?

**Date:** 2026-09-22
**Status:** DESIGN FROZEN. No statistic below has been computed, and the
analysis code was not written when this was frozen.
**Commit SHA:** filled in the commit that follows this freeze
**Analysis code SHA:** pending. It is recorded in an implementation manifest,
after code review and before any run on real data. This document is not edited
to insert it.
**Amends:** `experiments/PREREG_H3_MAGNITUDE_AND_SHARED_AXIS.md` (frozen `f288507`).
**Linked:** `experiments/PREREG_H3_REFERENCE_DISCORDANCE.md`, frozen in the same
commit. That registration covers the disagreement between the shRNA and CRISPRi
references; this one covers only the threats intrinsic to H3.

`f288507` registered S1–S3 and a rebuild gate requiring the rebuilt values to
reproduce the deposited artifact. The gate failed because the deposited artifact
carried an indexing defect (Deviation 9). S1–S3 were therefore void before any of
them was computed.

This amendment makes three changes:

1. **Validation target.** It replaces the validation target with the corrected
   artifact.
2. **Validation route.** It replaces the validation route with an
   identifier-aligned reconstruction plus an integrity check.
3. **Target-level versions.** It adds target-level versions of S1–S3.

The third change is a post-error design amendment made after the target
concentration of the cohort was known. It is not part of the original registration.

The registered estimands, cutoffs and seeds of S1–S3 are unchanged, and S1–S3 are
run and reported as registered.

## Foreknowledge of data or evidence

The foreknowledge recorded in `f288507` is retained as history. Its values were
computed under the defect: projected against alignment 0.3756, raw against
alignment −0.0433, partial correlation given coverage 0.3397. The threshold
comparisons made in that registration refer to those values. The corrected
artifact is the observed baseline from here on.

Since `f288507`, the following have been computed and seen:

| quantity | value |
|---|---|
| Spearman ρ(P, E), corrected | +0.3172 (p = 4.8e-20), n = 795 |
| Spearman ρ(D, E), corrected | −0.0580 (p = 0.10), n = 795 |
| drugs over targets | 795 over 258; median 1 drug per target, maximum 37 |
| single-target and multi-target annotations | 363 and 432; the first listed target is used |

`E` is the signless target-axis alignment `cos^2(mean signature, u)`. It is
historically labeled on-target enrichment, and the deposited field keeps that
name. `P` is the target-axis projected dispersion and `D` raw direction
instability.

The linked registration lists the CRISPRi results, which are also known.

None of the following has been computed on the corrected baseline: `M_delta`, any
adjusted partial correlation, any bootstrap interval, and any permutation null.

## Research questions or hypotheses

S1, S2 and S3 are as worded in `f288507` and are not restated.

**S1-T.** The adjusted association of S1 holds with a bootstrap interval that
resamples targets instead of drugs.

**S2-T.** The adjusted association exceeds the association produced when target
directions are reassigned among unique targets, with every drug of a target
receiving the same reassigned direction.

**S3-T.** The adjusted raw association of S3 remains practically equivalent to
zero with a target-resampling interval.

S2-T carries the design. S2 as registered permutes direction vectors across drug
records. That splits drugs that share a target across different permuted targets,
which breaks the dependence the observed mapping has. S2-T preserves target group
sizes and within-target dependence. Without it, S1-T and S3-T are consistent with
an association that sharing the axis `u` produces for any assignment of targets.

S1-T and S3-T can return results in H3's favor. Each of the three can also weaken
H3.

## Inference criteria

The partial Spearman procedure, `M_delta`, `K`, the complete-case cohort and the
covariate matrix are those of `f288507`.

| hypothesis | holds when |
|---|---|
| S1 | as `f288507`: partial ρ(P, E given M_delta, K) ≥ 0.20 and the 95% drug-bootstrap interval lies above 0 |
| S1-T | partial ρ(P, E given M_delta, K) ≥ 0.20 and the 95% target-bootstrap interval lies above 0 |
| S2 | as `f288507`: p_perm < 0.01 over 10,000 record-level permutations |
| S2-T | p_perm < 0.01 over 10,000 unique-target permutations |
| S3 | as `f288507`: the 90% drug-bootstrap interval for partial ρ(D, E given M_delta, K) lies inside (−0.15, 0.15) |
| S3-T | the 90% target-bootstrap interval for the same quantity lies inside (−0.15, 0.15) |

**H3's interpretation survives only if S1-T, S2-T and S3-T all hold.** S1, S2 and S3
are reported with their verdicts for continuity. Where a drug-level verdict and
its target-level counterpart disagree, the target-level verdict governs the
manuscript. No component compensates for another, and no multiplicity correction
is applied across an intersection-union gate.

**If S2-T fails, the claim that H3 reflects annotated drug–target matching is
withdrawn, whatever S1-T shows.**

**Target bootstrap.** Each replicate draws 258 targets with replacement from the
cohort's targets and carries every drug of each drawn target, with duplicates
retained. Within each replicate:
- every variable is reranked with average ranks;
- both residual regressions are refitted with `numpy.linalg.lstsq`;
- the residual correlation is recomputed.

There are 10,000 replicates, with seed 20260926, distinct from the registered drug
bootstrap's 20260913 so that the two sets of draws are independent. Intervals are
percentile intervals. A nonfinite replicate invalidates the run and is
investigated, never discarded or redrawn.

**Unique-target permutation.** Let `π` be a uniformly random permutation of the
258 unique targets. Every drug `c` receives direction `u_{π(t(c))}`. The procedure
then:
- recomputes `P` and `E` from the reassigned direction;
- holds `D`, `M_delta` and `K` fixed;
- reranks, residualizes and recomputes the partial correlation.

Fixed points are allowed. There are 10,000 permutations, with seed 20260927,
distinct from the registered permutation's 20260914. p_perm is
`(1 + #{ρ_b ≥ ρ_obs}) / 10001`. The output reports the 2.5th, 50th and 97.5th
percentiles of the null, the full null distribution, and `ρ_obs − median(null)`.

Holding `M_delta` and `K` fixed under permutation is correct because neither
depends on the target direction: `M_delta` is the mean pairwise
signature-difference norm and `K` the number of cell lines. The run asserts that
both are unchanged when the direction is reassigned, so a future covariate that
did depend on `u` could not be held fixed by accident.

The secondary quantities of `f288507` are reported without a criterion: mean
signature norm as a covariate, and `P/M_delta`. Distinct-hairpin count per target
is reported as a record, and no stratified permutation is registered.

## Validation

Two gates are run in order. Both are fixed now.

**Reconstruction gate.** The GCTX is reprocessed by `experiments/modal_h3_rebuild.py`
(stages fetch, extract, shrna), with each per-drug matrix carrying its cell-line
identifiers. The result is then joined to the pinned extraction by identifier, not
by position:

    signature matrices joined on (drug, cell_id); genes joined on Entrez id
    both identifier sets unique, and equal as sets
    rows and columns reordered by identifier into extraction order before comparison
    the ordered identifier sequences hashed separately from the numeric arrays
    every per-drug matrix and every target consensus equal to the extraction
      within rtol = 1e-5 and atol = 1e-5

The GCTX holds float32 values whose magnitude reaches about 10, so an absolute
tolerance of 1e-6 could fail on rounding alone while the reconstruction is
correct. The tolerance above is the float32-appropriate one, and a failure here
voids the analyses.

**Integrity gate.** The cohort bundle is built by `experiments/build_h3_bundle.py`
from the pinned extraction. It must reproduce the corrected artifact on every
record:

    795 records, 795 unique drug identifiers, the same identifiers as the corrected artifact
    identical target and n_celllines on every record
    max |D|, |P|, |E| differences from the corrected artifact < 1e-6

The outcome of the two gates determines how the run proceeds:

- **Both gates hold:** S1–S3 and S1-T to S3-T run.
- **The reconstruction gate cannot be executed** because the GCTX or the volume is
  unavailable: the analyses run on the integrity gate alone. The manuscript then
  states that they operate on the pinned corrected extraction.
- **The reconstruction gate runs and fails, or the integrity gate fails:** every
  analysis here is void. The failure and its diagnostics are reported.

## Sample size

The cohort is 795 drugs over 258 targets, fixed by the corrected artifact. It is
not a design choice.

## Maximum claim under this registration

If S1-T, S2-T and S3-T hold, the paper may say:

> Within the LINCS shRNA reference, projected cross-context dispersion covaries
> with signless alignment to the annotated target axis after adjustment for
> response magnitude and coverage, and exceeds the association produced when
> target directions are reassigned among targets.

The reference is a genetic loss-of-function axis measured by shRNA on L1000. The
claim is about similarity to that axis, not about biochemical target engagement.

Regardless of outcome, the prohibitions of `f288507` stand. The paper may not say
any of the following:
- that projection isolates on-target consistency;
- that the consistency points toward the target;
- that the statistic decomposes `D` into on-target and off-target components;
- that a high value means stronger mechanism conservation.

This registration also does not license "raw instability carries no
target-related information". The linked registration reports a reference under
which raw instability is associated with alignment.

A failed component carries a set consequence:

- **S1-T fails:** robustness to response magnitude was not established.
- **S2-T fails:** correct drug–target matching could not be distinguished from
  coupling through the shared axis, and the separation claim is withdrawn.
- **S3-T fails:** the raw contrast was not practically null after adjustment.

## Inputs

| input | pin |
|---|---|
| corrected H3 artifact | sha256 `fd69e26fc9a3917323065b631688baeab8b283f735c8bf5b16210ba67bd21425` |
| superseded H3 artifact, retained for provenance | sha256 `65e5d10e272037987384f89e6208de478fa17f6b8fb946add893e5a24c2d80c4` |
| `lincs_subset.npz`, Modal volume `drug-perturbation-vol` | sha256 `2ad0f5d30ab826f9ec0cfe37f6b53b2829bfef1d73b7920c3441de6407adb4ec` |
| `lincs_shrna.npz`, same volume | sha256 `4a990e5072a43f59fdda60a2ff040f355f06be947c14bcfa87056c66332b58e7` |
| `GSE92742_Broad_LINCS_sig_info.txt.gz` | sha256 `19da29c0ee12ddf27f9698cd0da40beaff58657dcde9d382aae068737e831299` |
| `GSE92742_Broad_LINCS_gene_info.txt.gz` | sha256 `741216ccc53320119b47ab006de3bcad48963c57087c9e07f50f0d6cd088711a` |
| `lincs_shrna_siginfo.csv.gz` | sha256 `bd396fa0e1a2f00c1b5f2c8d2b35f9a056f5e5353382475655869038037ec014` |
| `frozen_drug_labels.json` | sha256 `e6d73384c5bbc501ee15f621aaa258f02f38ee0d7a9603b6c820bec3e0489d52` |
| GCTX | `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz`; sha256 as stamped on volume `di-h3` at first retrieval |
| analysis code | the freeze commit's `experiments/03c_h3_sensitivity.py`, `experiments/build_h3_bundle.py` and `experiments/modal_h3_rebuild.py` |

## Design and procedure

The code changes are limited to four:
1. The pins in `03c_h3_sensitivity.py` and the manifest keys it checks.
2. S1-T, S2-T and S3-T, added to that script.
3. Cell-line identifiers carried through `modal_h3_rebuild.py`.
4. The identifier-aligned comparison of the reconstruction gate.

The code is tested on synthetic data and reviewed before freeze, and its commit is
pinned here.

**Repairs of defect-era values: deterministic corrections of quantities affected by
Deviation 9, not new hypothesis tests.** The manuscript reports four quantities that were
computed on the superseded artifact and have no corrected replacement: the
HDAC-removal correlation, the correlations after removing the 20 lowest-`D` and the
20 highest-`P` drugs, and the drug-level bootstrap interval on the full set. Each
is recomputed on the corrected artifact with its original procedure, reported as a
repair, and carries no new criterion. The corrected value replaces the superseded
one in the manuscript whatever it shows, and the manuscript reports none of the
four until its repair has run.

Each run writes the full result structure, seeds, input hashes, replicate draws
and null distributions to `results/03c_h3_sensitivity/`. The analysis runs once.

- **Existing data:** yes, and every input is pinned above.
- **Data collection:** N/A — no new data are collected.
- **Blinding:** none is possible; the corrected unadjusted correlations are known.
- **Missing data:** records that fail a gate are reported, never replaced.
- **Outliers:** none removed.
- **Exploratory analyses:** any analysis not named here is exploratory.
