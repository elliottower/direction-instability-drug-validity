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
loader. A1's treatment of the shRNA side is unchanged and now applies to both.

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

## Foreknowledge at the time of this amendment

Deviation 12 and everything in it, including the corrected pooled CRISPRi arm:
projected rho = 0.5410 and raw rho = -0.1920 on the 131-drug cohort of the C0
construction, measured before this amendment was drafted. That construction is the
one `03b` reports. **No R0-R7 statistic has been computed**, and the C1-K562,
C1-RPE1 and C1-GW arms this amendment governs have never been run under either
axis.

That the C0 arm moved from -0.1247 to +0.5410 is known. The direction of the C1
arms is not, and the registered readings apply to them unchanged.

## Maximum claim under this registration

With the gate passing on both cohorts, the R0-R7 arms can say what the frozen
registration lets them say about the agreement between a drug's consistency and
its target's genetic direction, under five reference constructions, with both
operands in one coordinate system verified against the pinned source. They cannot
say that the earlier disagreement between references was biological: Deviation 12
records that the C0 disagreement was a coordinate artifact, and the C1 arms have
no earlier result to be compared against.

## What the paper may say afterwards

The corrected arms are the registered quantities measured for the first time in one
coordinate system. The manuscript reports them and the deviation log carries the
history. Where a registered criterion is met on one condition and not another, the
paper says so in those terms rather than reporting a pass.
