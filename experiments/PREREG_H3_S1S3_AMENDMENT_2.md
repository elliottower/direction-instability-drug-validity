# Amendment 2: repairing the shRNA gene axis, and what H3 and H4 become

**Date:** 2026-10-02
**Status:** FROZEN. No repaired statistic has been computed, and no repaired
artifact has been produced.
**Commit SHA:** fd1ae8d
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

**A1. The authoritative shRNA source becomes the GEO rebuild, and the gate keeps a
second source.** Three objects are distinguished, because a comparison between an
artifact and the arrays it was written from tests nothing:

| role | object |
|---|---|
| canonical | `shrna_signatures.npz`, written by `experiments/modal_h3_rebuild.py` from the pinned GCTX: every retained shRNA signature, its signature id, and its explicit ordered gene axis. |
| derived | `shrna_consensus.npz`, the target consensuses and unit directions. A cache with its own hash, reproducible exactly from the canonical artifact, the pinned target-membership table, the frozen construction code and the frozen gene axis. |
| independent reconstruction | a fresh parse of the same pinned GCTX, performed by the gate in its own code at gate time and joined to the canonical artifact by signature id and gene id. |
| retired | `lincs_shrna.npz`, read by no production loader. Its recovered permutation stays a forensic record and never becomes a production transformation. |

**The authoritative object is signature-level.** A consensus file alone cannot be
the canonical source: averaging can cancel an error in an individual signature, so
a gate that sees only consensuses cannot verify the records underneath them. Every
target direction the analyses use is derived from the canonical signatures and the
pinned membership, and the persisted consensus is checked against that derivation
rather than trusted.

`lincs_shrna.npz` and the rebuild carry identical values, so the choice is not
between two numbers; only the rebuild carries an axis that can be verified against
its own labels, and only a parse of the source can verify it.

**A2. The compound source is unchanged.** `lincs_subset.npz` follows its declared
axis, which the compound half of the gate confirmed before the shRNA half ran. Drug
signatures are aligned to the rebuild's axis by gene identifier, never by position.

**A3. The consensus rule is unchanged.** Per target gene, the mean over its
signatures, at least three distinct signature ids, as the corrected artifact
computes it. Only the coordinates change.

**A4. What is recomputed.** H3's projected dispersion, signless target-axis
alignment and their correlations, from `experiments/03_phenotype_projection.py`;
H4's localization, from `experiments/04_localization.py`; and the `proj_shrna` and
`enrich_shrna` fields of the CRISPRi records, from
`experiments/03b_h3_crispri_ground_truth.py`. Raw direction instability, H1, H2,
H5, the 66-fold prediction and every quantity built from the CRISPRi reference
alone are untouched, because none of them pairs the two matrices: the CRISPRi
directions are placed in the landmark order the compound extraction declares,
which is the order it follows. The `proj_shrna` and `enrich_shrna` fields stored
in the CRISPRi comparison records are recomputed, because they use the shRNA
reference.

The shRNA comparator of `experiments/03d_h3_reference_discordance.py`
(`build_shrna_reference`) draws on the same source and switches with it. R0-R7
have never been computed, so nothing there is recomputed; they run for the first
time against the repaired reference. The R5 mapping is frozen and uses neither
matrix.

**A5. The artifacts are superseded, not overwritten.** Each defect-era artifact
moves to `superseded/mislabeled_shrna_axis/` with a README built like the one
Deviation 9 left at `results/03_phenotype_projection/superseded/`, carrying the
superseded and replacement sha256, the ordered gene-axis hash and the signature-id
hash for each, the target-membership hash, the generator's commit, the commit this
amendment is frozen at, and the reason. The superseded artifact stays readable for
forensic reproduction and is rejected by every production loader. Its values are
recorded as invalid measurements of their stated quantities, not as earlier
estimates of them.

**A6. The pins move with the artifacts.** The reconstruction gate, the integrity
gate, and the reference-discordance registration's reproduction gate are re-pinned
to the repaired artifacts. Their tolerances, procedures and readings are unchanged.

**A7. The gate must then pass, on two sources.** Both halves, on the repaired
inputs, before any S1-S3 or R0-R7 statistic is computed. The shRNA half compares
the canonical artifact against the gate's own parse of the pinned GCTX, and
separately compares the persisted consensuses against directions recomputed from
the canonical signatures and the pinned membership: the first reaches the
coordinate semantics, the second is what caught Deviation 9 and is what makes the
consensus file a checked cache rather than a second source of truth. A failure still voids the analyses registered
against it.

**The gate's tolerance is unchanged and is now applied as registered.** The
registered tolerance is the elementwise float32 rule, `|a-b| <= atol + rtol|b|`
with `rtol = atol = 1e-5`. The implementation reduced it to a single threshold of
`2e-5` on the maximum over unnormalized consensus entries, whose magnitude reaches
about ten: looser than the registered rule wherever an entry is below one and
stricter wherever it is above. The elementwise rule now decides every comparison.
This corrects the implementation to the registered tolerance rather than changing
the tolerance.

**A8. A declared axis is checked at load; the coordinates are checked against the
source.** Two layers, because the first alone would not have caught this defect:

1. **Schema and provenance, at load, fatal.** Every production loader asserts that
   the artifact's declared gene axis equals the frozen landmark order elementwise,
   and that the artifact's sha256 matches the execution manifest.
2. **Coordinate semantics, before analysis, fatal.** The reconstruction gate
   verifies the numeric columns against the gate's own parse of the pinned GCTX,
   joined by identifier.

The first layer is a statement about labels and cannot see labels that lie.
`experiments/build_h3_bundle.py:48` already asserts that the two extractions share
a gene axis, and that assertion passes on the defective file, because both files
declare the same identifiers. Only the second layer reaches the values.

**A9. `experiments/build_h3_bundle.py` is retired with the extraction.** The
manifest pins it as implementation of record, no registered command runs it, and
its governing premise is refuted: its docstring chooses `lincs_shrna.npz` over the
GCTX rebuild on the grounds that the rebuild "reproduced their values exactly but
not their coordinate order", which is now known to have the ordering backwards. Its
eligibility floor also counts metadata rows rather than distinct signature ids,
diverging from the three scripts that do run. It moves to
`experiments/superseded/` and the manifest is regenerated under A6.

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
