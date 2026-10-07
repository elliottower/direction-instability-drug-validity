# Does the pre-fix shRNA consensus indexing explain the residual difference from the deposited values?

**Date:** 2026-10-07
**Status:** FROZEN before the diagnostic ran
**Commit SHA:** e684d16
**Kind:** forensic diagnostic. No registered statistic is computed and no
confirmatory claim rests on it. It asks how the deposited numbers were produced.
**Relates to:** Deviation 9 in `DEVIATION_LOG.md`, which found the defect; Deviation
14, which records the residual this asks about;
`experiments/PREREG_H3_S1S3_AMENDMENT_5.md`, whose gate rule this cannot change.

## The question

**Does the exact historical indexing defect explain the residual across the
complete cohort and both reference-dependent quantities?**

Categorical, and answered over all 795 drugs rather than any subset of them.

## Foreknowledge

The residual is known in full before this is written. The authoritative route
differs from the deposited `proj_shrna` by at most 1.139e-06, which is 2.21 float32
ULPs at that magnitude and one drug of 795 over the frozen flat tolerance. Two
explanations were proposed and measured, and both were refuted: the float32 mean
accumulator, where the float32 and float64 accumulations of the same signatures
differ from each other by about 1e-9 while both sit about 1e-6 from the deposited
values; and a float32 projection path, where no dtype route reaches below a 1e-7
relative floor. Context membership is identical between routes for every drug
examined, and among the nine focus drugs the difference is spread over several
context pairs rather than localized to one.

Deviation 9 found that the deposited shRNA consensus was indexed by metadata row
rather than by signature id, and the current loader carries the fix. A reference
direction built on the wrong rows would move `P` and `E` without moving membership,
which is the shape of what is left unexplained. That is the reason for this
diagnostic and it was formed after seeing the residual.

## Hypotheses

**H1.** The pre-fix indexing reproduces the deposited `P_shrna` and `E_shrna`
better than the corrected indexing does, across the cohort and both quantities.

**H2.** It does not, and the residual remains unexplained.

H1 carries the diagnostic. Neither outcome licenses the defective indexing for
production.

## Implementation, fixed before it runs

The pre-fix path is reconstructed exactly as Deviation 9 describes it: the
consensus matrix indexed by metadata row position rather than by signature id,
with nothing else altered. No variant of the indexing is tried, no tolerance is
changed, and no parameter is searched. One reconstruction runs once.

`P_shrna` and `E_shrna` are computed for all 795 drugs of the registered cohort
against that reference. `D` is excluded, with the reason: `D` is one minus the mean
pairwise cosine of the drug's own signatures and does not depend on the reference
direction at all, so it cannot move under any change of reference and its agreement
would be uninformative.

Each quantity is compared against three things separately: the deposited records,
the corrected-reference legacy route, and the authoritative production route.
Every per-drug error is written to the artifact, not only the maxima.

## Inference criteria

| outcome | holds when |
|---|---|
| H1 holds | for both `P_shrna` and `E_shrna`, the maximum absolute error against the deposited records falls below 1e-08, and the median absolute error falls by at least a factor of ten against the corrected reference |
| H1 holds in part | one quantity meets that and the other does not, or the maxima fall by at least a factor of ten without reaching 1e-08 |
| H2 holds | neither quantity's maximum absolute error falls by as much as a factor of two |

**The residual is described as explained only under H1.** A single drug improving,
or a maximum that merely looks better, decides nothing: the criterion is over the
whole cohort and both quantities.

## What this cannot establish

A reference that reproduces the deposited values does not establish that the
deposited reference was built that way; it establishes that one historical defect
accounts for the difference to the stated precision. Nor can it establish anything
about signature amplitudes, for the reason Amendment 5's D5 gives.

## Disposition

The result is preserved whether it supports or rejects the candidate, and is
recorded in Deviation 14 and in Amendment 5's account of the residual. Amendment
5's gate rule does not change on the outcome: the pre-fix indexing is quarantined
forensic work and is not a production route under any result.

## Maximum claim under this diagnostic

Confirmatorily, nothing. The paper may report, as a statement about provenance
rather than about H3, whether one named historical indexing defect accounts for the
difference between the deposited values and the reconstructable routes, at the
precision stated, over all 795 drugs and both reference-dependent quantities.
