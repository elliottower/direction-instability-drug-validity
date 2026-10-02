# Amendment 2: repairing the shRNA gene axis, and what H3 and H4 become

**Date:** 2026-10-02
**Status:** DRAFT. Not frozen. No repaired statistic has been computed.
**Amends:** `experiments/PREREG_H3_S1S3_CORRECTED_BASELINE.md` (design frozen
`7f57136`), and the inputs the linked
`experiments/PREREG_H3_REFERENCE_DISCORDANCE.md` pins.
**Discovery record:** Deviation 11 in `DEVIATION_LOG.md`.

## What the gate found

The registered reconstruction gate failed on its shRNA half. All 14,656 shRNA
signatures match the GEO rebuild exactly, maximum absolute difference 0.0, under
one bijective permutation of the gene axis. `lincs_shrna.npz` therefore holds
correct values under a wrong axis label, while `lincs_subset.npz` follows the axis
both files declare. Every quantity pairing a drug signature with a target
direction was computed between different coordinate systems.

The gate did what it was registered to do. Nothing here loosens it.

## The repair

**A1. The authoritative shRNA source becomes the GEO rebuild.** Target directions
are built from `shrna_signatures.npz`, the signatures parsed from the pinned GCTX
in the frozen landmark order, rather than from `lincs_shrna.npz`. The two carry
identical values; only the rebuild carries an axis that can be verified against its
own labels.

**A2. The compound source is unchanged.** `lincs_subset.npz` follows its declared
axis, which the compound half of the gate confirmed before the shRNA half ran. Drug
signatures are aligned to the rebuild's axis by gene identifier, never by position.

**A3. The consensus rule is unchanged.** Per target gene, the mean over its
signatures, at least three distinct signature ids, as the corrected artifact
computes it. Only the coordinates change.

**A4. What is recomputed.** H3's projected dispersion, signless target-axis
alignment and their correlations; H4's localization; and the `proj_shrna` and
`enrich_shrna` fields of the CRISPRi records. Raw direction instability, H1, H2,
H5, the 66-fold prediction and every CRISPRi quantity are untouched, because none
of them pairs the two matrices.

**A5. The artifacts are superseded, not overwritten.** Each defect-era artifact
moves to `superseded/mislabeled_shrna_axis/` with a README naming this amendment
and carrying both hashes, as Deviation 9's artifacts were.

**A6. The pins move with the artifacts.** The reconstruction gate, the integrity
gate, and the reference-discordance registration's reproduction gate are re-pinned
to the repaired artifacts. Their tolerances, procedures and readings are unchanged.

**A7. The gate must then pass.** Both halves, on the repaired inputs, before any
S1-S3 or R0-R7 statistic is computed. A failure still voids the analyses
registered against it.

## What this does not change

- Every hypothesis, statistic, criterion, interval, permutation and reading in the
  two frozen registrations and in Amendment 1.
- The R5 mapping, frozen at sha256 `152361cb3174a5fb…`, which uses neither matrix.
- The registered seeds.

## Foreknowledge at the time of this amendment

The forensics in Deviation 11, and the defect-era values: H3's projected
correlation 0.3172 and raw correlation −0.0580, H4's AUROC gap +0.0005, and the
CRISPRi comparisons. **No repaired value has been computed.** The repair is
specified before any of them exists.

## What the paper may say afterwards

The recomputed H3 is the registered quantity for the first time: a correlation
between projected dispersion and alignment, both measured in one coordinate
system. Whatever it shows is the result, and the defect-era values do not appear
in the manuscript.

If the recomputed association is substantially weaker, that is not a separate
finding to be explained away. It bears directly on the shared-axis concern the
original registration raised: `P` and `E` are built from the same direction, and a
direction that pointed nowhere in particular could still produce an association
between them. S2 and S2-T were registered to test exactly that, and they now run
against a reference that is what it claims to be.

The paper states the correction once, in the deviation log and the methods, and
reports the repaired numbers. It does not narrate the sequence of defects.
