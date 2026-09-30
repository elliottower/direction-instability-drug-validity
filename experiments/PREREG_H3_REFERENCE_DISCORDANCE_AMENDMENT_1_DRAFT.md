# Amendment 1 to the reference-discordance registration: R5's exposure and mapping under the release as it is

**Date:** 2026-09-30
**Status:** DRAFT. Prospective: no drug-level PRISM value, and no association
involving one, has been computed or inspected.
**Amends:** `experiments/PREREG_H3_REFERENCE_DISCORDANCE.md`, design frozen at
`7f57136`. Only R5 changes. Every other module, gate, interval and reading stands
as frozen.
**Discovery record:** Deviation 10 in `DEVIATION_LOG.md`.

## What the frozen text assumed, and what the release is

R5 registered its toxicity measure as "minus the median, across cell lines, of the
PRISM Repurposing 19Q4 primary-screen log-fold-change (2.5 µM,
replicate-collapsed)", with exact-name matching against a compound list's
`Drug.Name` and `Synonyms` fields. Those details were written from the 24Q2 public
release, the file on disk when the registration was drafted, which covers 5 of the
131 CRISPRi-arm drugs.

The 19Q4 primary screen, retrieved on 2026-09-30 from figshare 9393293, differs:

- It is a wide matrix, 579 cell lines by 4,686 treatments, one treatment per
  compound.
- Its treatment table carries `name` and no synonym field.
- Its doses are not uniformly 2.5 µM: 4,147 HTS treatments span 0.03 to 5.6 µM and
  539 MTS004 treatments span 2.29 to 3.12 µM.

Taking each compound at whatever dose the release assigned it would change the
estimand from a response at a common exposure to a response at a
compound-specific exposure, which mixes the drug's effect with the dose it
happened to be screened at. That is not the registered quantity and could not
carry the registered reading.

## The amendment

**A1. Exposure window.** R5 uses only treatments whose dose lies within 2.0 to 3.0
µM inclusive: the registered 2.5 µM with a tolerance of 0.5. The window is fixed
here, from the treatment metadata alone, before any response value is summarized.
It is not the window that maximizes coverage, and no other window is tried after
an association is seen.

**A2. Name matching.** Exact-name matching, after case folding and whitespace
collapse, uses the release's `name` field. The release has no synonym field, and
no source outside the pinned files is consulted. A drug whose name is absent is
reported unmapped, never rescued by hand.

**A3. Treatment selection.** Each compound holds one treatment in the release, so
the registered tie-break applies only when one drug name matches several
compounds. It is unchanged: the treatment measured in the most cell lines, ties
broken by the lexical order of the column name.

**A4. The mapping is frozen before any response value is summarized.** R5 runs in
two stages:

1. **Mapping.** From the pinned treatment table, the missingness of the response
   matrix, and nothing else, write every accepted mapping with its identifier,
   route, dose, screen, candidate count and cell-line count, and every rejected
   drug with its reason. Hash that artifact and commit it.
2. **Response.** Accept only that committed artifact, and compute the toxicity
   measure and every R5 quantity from it. The stage refuses to run against a
   mapping it did not receive as a frozen file.

No mapping is revised after any response value is summarized.

**A5. Coverage, counted under A1 to A3 before the freeze.** 120 of the 131
CRISPRi-arm drugs map, over 37 targets, against the registered minimum of 60 drugs
over 20 targets. The gate is met, so R5 remains confirmatory. In the shRNA arm 718
of 795 drugs map. The doses of the mapped treatments have a median of 2.5 µM and a
range of 2.01 to 2.89 µM.

**A6. What R5 may claim.** The measure is the viability response at a near-common
exposure of 2.0 to 3.0 µM, not a dose-response summary and not a potency estimate.
Where the registered reading is returned, the paper may say that toxicity so
measured statistically accounts for the raw CRISPRi association. It may not say
that R5 isolates intrinsic toxicity, nor that any compound's response establishes
its potency.

## What this does not change

- R5's gate of 60 drugs over 20 targets, its quantities, its intervals, and its
  three readings.
- The secondary summary, the fraction of cell lines below −1, which stays
  criterion-free.
- Every other module of the registration frozen at `7f57136`.

## Foreknowledge at the time of this amendment

The counts in A5, the dose distribution, and the schema. No log-fold-change value
has been summarized, no toxicity variable has been constructed, and no association
between toxicity and any other quantity has been computed.

**Fallback, fixed here.** If the response stage finds that fewer than 60 mapped
drugs or fewer than 20 targets survive the complete-case cohort, R5 is reported
descriptively and labeled exploratory, and it neither supports nor counts against
the toxicity explanation. The cohort is not widened to rescue it.
