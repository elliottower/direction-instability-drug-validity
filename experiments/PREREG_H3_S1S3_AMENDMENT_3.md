# Amendment 3: the compound source, and the arms that read it

**Date:** 2026-10-05
**Status:** FROZEN. No R0-R7 statistic has been computed.
**Commit SHA:** 03175d0
**Amends:** `experiments/PREREG_H3_S1S3_AMENDMENT_2.md` (frozen `fd1ae8d`), whose
A2 this replaces, and the inputs `experiments/PREREG_H3_REFERENCE_DISCORDANCE.md`
pins.
**Discovery record:** Deviation 12 in `DEVIATION_LOG.md`.

## What A2 asserted, and what the gate found

A2 reads: "The compound source is unchanged. `lincs_subset.npz` follows its
declared axis, which the compound half of the gate confirmed before the shRNA half
ran."

It does not. Both extractions declare one gene axis and hold another, and it is the
same other: one permutation across all 978 columns, two columns in place, verified
against the pinned source on all 41,643 cohort signatures at a maximum absolute
difference of 0.0. The compound half appearing to pass came from a gate version
with an intersection defect that wrote no record.

## The repair

**B1. The authoritative compound source becomes the rebuild.** Drug signatures are
the per-drug, per-cell-line matrices written by `experiments/modal_h3_rebuild.py`
from the pinned GCTX in the frozen landmark order, extended on 2026-10-05 by the
seventeen drugs the CRISPRi arm needs. `lincs_subset.npz` is read by no production
loader. A1 remains unchanged for the shRNA side; B1 establishes the analogous
source-derived rule for compound signatures.

**B2. The gene axis every reference is placed on is the frozen landmark order.**
`experiments/03d_h3_reference_discordance.py:604` builds its symbol list from the
extraction's declared `gene_ids`; it builds it from the frozen order instead. That
list places C0, C1-K562, C1-RPE1, C1-GW and C1-GW-phenotype-positive, and selects
the coordinate R7f removes, so a single wrong list reaches every external
reference.

**B3. What this changes, and what it does not.** No registered hypothesis,
statistic, criterion, interval, permutation, seed or reading changes. The cohort,
the comparison rule, the practical null and the equivalence region stand as frozen.
Only the artifact the analyses read changes, and it changes to the one the gate
verified against the source.

**B4. The quantities that must be recomputed.** Every CRISPRi arm in R0-R7, the
harmonized comparators, and R7f's target-gene sensitivity. The shRNA arm of R0-R7
is built from the canonical signatures under A1.

**B5. Frozen preflight and pass condition.** Before any R0-R7 computation the
versioned gate runs on the two cohort manifests pinned in B8. Release is authorized
only if every predicate the frozen gate specification names holds: cohort identity;
compound reconstruction against the pinned source on each cohort; shRNA signature,
consensus and stored-direction reconstruction on the exact pinned shRNA-eligible
target set; the 812-to-795 coverage partition; and legacy-artifact identity under the
forensic coordinate check B6 specifies. Every predicate is evaluated and recorded
independently, so a structural failure in one does not leave the others unreported.
Any failed predicate, any predicate absent from a run, and any structural failure
voids the analyses. The predicate set is fixed in code as `RELEASE_PREDICATES` and
cannot be narrowed by an invocation. The gate v1 reports are preserved as failed
reports and are not superseded as records.

**B5a. The shRNA eligibility boundary.** The shRNA-eligible universe is the exact
258-target identifier set pinned in B8, not every target represented in the 812-drug
compound cohort. The complement in that cohort is exactly the seventeen pinned
drug-target records whose targets are EEF2, EIF2S1, FNTA, HCRTR1, KCNA10, RPL3, RPS2
and TUBB. The gate verifies the full eligible set, the full complement, that no
deposited record carries one shRNA quantity and not the other, and that the rebuild's
shRNA targets and the deposited analysis's valued targets are the same set.
Cardinality does not satisfy any of these: a set with one target dropped and another
substituted has the same count, so every comparison is on a set or on a canonical
hash. In the deposited records the two shRNA quantities are absent as keys rather
than present and null, and both forms count as missing.

**B6. The retired artifacts, and what their check means.** `lincs_subset.npz` and
`lincs_shrna.npz` keep their bytes and hashes so the deposited analysis remains
reproducible. Neither is a production operand of anything, and the recovered
permutation is never a production transformation.

The legacy check in the release condition is a custody check on the deposited
analysis, not a statement that the retired artifact is valid for any purpose. It
holds when every one of the following does: the file's bytes match the sha256 frozen
in `registry/frozen/expected_identities.json`; the coordinate map is a total
bijection over the 978 landmarks; composing the map with its inverse is the identity,
so its direction is the one its field name states; the map's canonical permutation
hash and its file hash both match their frozen values, which are held outside the map
because an artifact carrying its own expected hash attests to itself; the map
reproduces the pinned source **exactly, at the raw-signature level**, which is where
the map is recovered and where the map artifact records it; the relabeled matrices
reproduce the source **within the production tolerance at the aggregated level**, of
per-drug per-cell-line means; and the declared labels do not reproduce the source
within that tolerance, with the count of drugs outside it recorded rather than a
boolean. The two levels are named apart because exactness is true of raw signatures
and false of aggregated means, where even the production rebuild sits at 3.12e-06
against the source. Disagreement alone would not satisfy any of this, because almost
any corruption disagrees.
A legacy failure voids the run because custody of the deposited numbers has failed,
and the report states production authorization and legacy custody separately.

**B7. What is superseded, exactly.** Amendment 2's A2, in full: "The compound
source is unchanged. `lincs_subset.npz` follows its declared axis, which the
compound half of the gate confirmed before the shRNA half ran. Drug signatures are
aligned to the rebuild's axis by gene identifier, never by position." Its first two
sentences are false and are replaced by B1. Its third stands and is subsumed: the
alignment is by identifier. The input pins of
`PREREG_H3_REFERENCE_DISCORDANCE.md` naming `lincs_subset.npz` and
`lincs_shrna.npz` as production inputs are replaced by B8. A2 is not edited or
erased; every other clause of Amendment 2 remains operative.

**B8. The production identity, pinned.** The run consumes and records, each by
sha256: the rebuild's `landmark_gene_ids.json`, `shrna_consensus.npz`,
`shrna_signatures.npz`, `rebuild_manifest.json` and every compound shard; the four
cohort manifests under `registry/cohorts/` — `cohort_795_shrna_paired.json`,
`cohort_812_compound.json`, `shrna_eligible_targets.json` and
`shrna_excluded_records.json` — each with its own canonical pair, drug and target
hashes, its record and identifier counts, its uniqueness constraints and its schema
version; the pinned legacy coordinate map and its canonical hash; the coverage
measurement `shrna_coverage_identity_v2.json`; a passing gate report on both cohort
scopes, with its predicates, cohorts and tolerances; this amendment's freeze commit
and file hash; and the analysis code commit, including the gate version and the hash
of the preserved v1 gate. The cohorts are pinned as minimal `(drug, target)`
manifests rather than as results files, because
`phenotype_projection_results.json` carries projected instability and on-target
enrichment, and pinning it coupled cohort eligibility to derived values. The retired
extractions are recorded as legacy artifacts, separately from production inputs.

The values themselves are in `registry/frozen/expected_identities.json`, sha256
`3de491d6234c9dc9f56072227b50e7a513fdc25720a9bc113ef56c0939e3b9b2`,
and the record of what ran is
`results/03d_h3_reference_discordance/implementation_manifest_amendment_3.json`,
generated from the files by `experiments/03h_freeze_identities.py` and never patched
by hand. The gate reads its expectations from the frozen identities and not from the
artifacts they describe.

The two passing gate reports, by sha256, are:

| scope | drugs | report sha256 |
|---|---|---|
| compound | 812 | `add1f337fce02ca5d8e449de9f917f1f82ac96511347b9280d56d10f599628ed` |
| shRNA-paired | 795 | `67d0489b8742a39847d9e22f3225c4bb902e471cf47b7fb8c682e7ff598a84c1` |

Both hold all five predicates `RELEASE_PREDICATES` names, under gate code sha256
`4eda8019f9298867c5404018fe6439995b792c3f8ddc65c34660ee6e9a7ee9d0`
at rtol and atol 1e-05. The superseded gate v1 is
kept at `experiments/superseded/03e_reconstruction_gate_v1.py`, sha256
`9dc69bea8cf9e5635acad16084056a08c8d7fe2739f1ce329626d296dd8e24d1`, and its failed
reports stay on record as failed.

**B9. What the deposited records are a reproduction target for.** Raw direction
instability and the shRNA quantities are invariant under the permutation both
extractions carried, so they must still reproduce the deposited per-drug values to
the registered tolerance, and a failure voids the run. The C0 projected dispersion
and alignment are not a reproduction target: the deposited values placed the
reference on the declared axis while the signatures followed another, so the
corrected values are deliberately different. Their difference from the deposited
values is measured and reported as a correction. The four-route audit of Deviation
12 is the separate aggregate record and is not a gate.

**B10. Freeze discipline.** After this amendment is frozen, no code, cohort,
statistic, criterion or interpretive rule changes except through another recorded
amendment or deviation.

## Foreknowledge at the time of this amendment

Deviation 12 and everything in it, including the corrected pooled CRISPRi arm:
projected rho = 0.5410 and raw rho = -0.1920 on the 131-drug cohort of the C0
construction, measured before this amendment was drafted. That construction is the
one `03b` reports. **No R0-R7 statistic has been computed**, and the C1-K562,
C1-RPE1 and C1-GW arms this amendment governs have never been run under either
axis.

At this freeze the corrected C0 result is known. No C1-K562, C1-RPE1 or C1-GW
association has been computed, and their previously frozen analyses, cohorts,
hierarchy and readings are left unchanged. Corrected C0 is a known-outcome
correction interpreted under criteria frozen before it; the C1 arms are
prospectively evaluated under the earlier frozen design with this foreknowledge
disclosed; anything added or altered after seeing the C0 result is exploratory
unless justified independently of it.

## Maximum claim under this registration

With the required gates passing, R0-R7 may report the registered associations
between cross-cell-line drug-response geometry and target-reference alignment under
the five prespecified reference constructions, with each operand derived from its
own pinned source and both placed on a verified common coordinate axis. These
analyses cannot support the claim that the previously reported C0 disagreement was
biological, because that disagreement did not survive coordinate correction. The
C1, RPE1 and genome-wide analyses are first measurements under their respective
constructions, not replications of the C0 reversal.

## What the paper may say afterwards

The corrected arms are the registered quantities measured for the first time in one
coordinate system. The manuscript reports them and the deviation log carries the
history. Where a registered criterion is met on one condition and not another, the
paper says so in those terms rather than reporting a pass.
