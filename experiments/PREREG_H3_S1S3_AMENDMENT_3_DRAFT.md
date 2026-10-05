# Amendment 3: the compound source, and the arms that read it

**Date:** 2026-10-05
**Status:** DRAFT. Not frozen. No R0-R7 statistic has been computed.
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

**B5. The gate runs first, on both cohorts.** The registered 795-drug cohort, and
the 812-drug CRISPRi cohort the arms use. A failure on either voids the analyses
registered against it.

**B6. The retired artifacts are preserved.** `lincs_subset.npz` and
`lincs_shrna.npz` keep their bytes and hashes so the deposited analysis remains
reproducible. The recovered permutation stays a forensic record and is never a
production transformation.

**B7. What is superseded, exactly.** Amendment 2's A2, in full: "The compound
source is unchanged. `lincs_subset.npz` follows its declared axis, which the
compound half of the gate confirmed before the shRNA half ran. Drug signatures are
aligned to the rebuild's axis by gene identifier, never by position." Its first two
sentences are false and are replaced by B1. Its third stands and is subsumed: the
alignment is by identifier. The input pins of
`PREREG_H3_REFERENCE_DISCORDANCE.md` naming `lincs_subset.npz` and
`lincs_shrna.npz` as production inputs are replaced by B8. A2 is not edited or
erased; every other clause of Amendment 2 remains operative.

**B8. The production identity, pinned.** The run consumes and records: the rebuild's
`landmark_gene_ids.json`, `shrna_consensus.npz`, `shrna_signatures.npz`,
`rebuild_manifest.json` and every compound shard, each by sha256; the rebuild
manifest's cohort-identifier hash; the exact 812-drug identifier set, which is the
registered 795 plus the seventeen named in Deviation 12 and in
`results/03c_h3_sensitivity/rebuild_extension.json`; a passing gate report on both
cohorts, by sha256, with its cohorts, comparisons and tolerances; this amendment's
freeze commit and file hash; and the analysis code commit. The retired extractions
are recorded as legacy artifacts, separately from production inputs.

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
