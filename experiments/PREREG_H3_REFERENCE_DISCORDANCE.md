# Registered follow-up: Why does H3's pattern reverse between the shRNA and CRISPRi references?

**Date:** 2026-09-22
**Status:** DESIGN FROZEN. No analysis below has been computed, and the
analysis code was not written when this was frozen.
**Commit SHA:** 7f57136
**Analysis code SHA:** pending. It is recorded in an implementation manifest,
after code review and before any run on real data. This document is not edited
to insert it.
**Kind:** a registered post hoc follow-up. It was designed after the reversal was
seen, and it specifies analyses that have not been run. It is not a confirmatory
test of H3 and cannot raise H3 above its original registration (`249abaf`).
**Linked:** `experiments/PREREG_H3_S1S3_CORRECTED_BASELINE.md`, frozen in the
same commit. That registration covers the shared-axis threat to H3 within the
shRNA reference (S1–S3, S1-T to S3-T). This registration covers the disagreement
between references.

## Terms

- **Reference.** A target-perturbation reference, or criterion measure, is the
  expression signature of reducing a target with a given technology in a given
  context. Neither reference is a ground truth for target engagement.
- **Alignment.** `E = cos^2(m, u)` is signless target-axis alignment: it cannot
  separate alignment from anti-alignment.
- **Reversal.** "The qualitative pattern reverses" means that the projected
  association is positive only under shRNA, and the raw association is positive
  only under CRISPRi.

## Foreknowledge of data or evidence

Values below come from the corrected per-drug records,
`results/03b_h3_crispri/h3_crispri_results.json` (sha256 `f00f8428…`):

| reference | cohort | n drugs | targets | ρ(P, E) | ρ(D, E) |
|---|---|---|---|---|---|
| shRNA | all | 795 | 258 | +0.3172 (p = 4.8e-20) | −0.0580 (p = 0.10) |
| CRISPRi K562, construction C0 | all | 131 | 41 | −0.1247 (p = 0.16) | +0.2723 (p = 0.0017) |
| shRNA | matched | 114 | 33 | +0.3398 (p = 2.2e-4) | not computed |
| CRISPRi, C0 | matched | 114 | 33 | −0.1440 (p = 0.13) | not computed |

Deviation-log Entry 5 reports a drug-level 95% interval of [−0.29, 0.05] for the
CRISPRi ρ(P, E) under C0.

**Also known, from metadata only.**

- *Target concentration.* The CRISPRi arm's targets are led by TOP2A (12 drugs),
  CDK1 (11), AURKA (11), MTOR (11), TUBB (8) and ATP1A1 (8).
- *Target annotations.* The first listed target is used:

  | cohort | single-target | multi-target | MOA contains "inhibitor" |
  |---|---|---|---|
  | shRNA arm | 363 | 432 | 462 |
  | CRISPRi arm | 61 | 70 | 121 |
  | matched | 50 | 64 | 106 |

- *K562 and hematopoietic profiles.* GSE92742 holds no K562 compound signature
  among its 71 compound lines. Nine compound lines have a hematopoietic primary
  site: HL60, JURKAT, NOMO1, PL21, SKM1, THP1, U266, U937 and WSUDLCL2. 29 of the
  131 CRISPRi-arm drugs (22 targets) have at least one hematopoietic signature.
- *PRISM coverage.* PRISM Repurposing 24Q2 public covers 5 of the 131 CRISPRi-arm
  drugs by name with log-fold-change in at least 100 lines. No log-fold-change
  value was summarized.
- *Official Replogle files.* The official pseudobulk files were downloaded, and
  only their structure and perturbation lists were read.
  - The K562-essential file has 2,285 rows over 8,563 genes, and the RPE1-essential
    file 2,679 rows over 8,749 genes.
  - The K562 genome-wide file has 11,258 rows over 8,248 genes.
  - All 41 CRISPRi-arm targets are perturbed in each of the three files.
  - The genome-wide file perturbs 169 of the 266 targets among the existing 812
    records.
- *H1 cytotoxicity labels.* H1's label marks 4 of the 131 CRISPRi-arm drugs and 7
  of the 795 shRNA-arm drugs as cytotoxic.
- One per-drug record (10-hydroxycamptothecin) was displayed while field names were
  checked.

**Not computed:** any quantity under a construction other than C0, any
target-cluster interval, any permutation, any own-target percentile, any
raw-instability correlation on the matched cohort, any analysis of RPE1, of the
genome-wide screen, or of hematopoietic profiles, any toxicity association, and
any reliability estimate.

## Explanations under consideration

| | explanation | primary analyses |
|---|---|---|
| A | The CRISPRi directions for these targets share a component across targets, from essential-gene selection or cytotoxicity, rather than carrying target identity | R2, R3, R4, R5 |
| B | The shRNA association arises from the shared L1000 platform and axis | linked S2-T |
| C | The CRISPRi directions are specific to the K562 context | R6 |
| D | Drug-level inference overstates independence, heavily represented targets dominate the estimates, or the arms differ in drug composition | R1 |
| — | The CRISPRi direction as constructed is not a valid criterion measure | R0 |

The explanations are not mutually exclusive, and each analysis is reported on its
own.

## Common definitions

**Drug quantities.** For drug `c` with first-listed target `t(c)`, let `S_c` be its
`(K_c, 978)` matrix of per-cell-line signatures from the pinned extraction, and
`m_c` its mean. For a reference `X` with unit target directions `u_t^X`:

- `D_c` is raw direction instability, always on the 978 landmarks.
- `P_c^X` is projected dispersion and `E_c^X = cos^2(m_c, u)` is alignment.
- `s_c^X = cos(m_c, u)` is signed alignment.

All four come from the unmodified functions in
`geometry/direction_instability.py`.

**Harmonized comparators.** A comparison between references must not also change
the gene space, so two further quantities are defined:

- `shRNA-h`: the shRNA reference's `P` and `E`, computed on C1's landmarks, with
  drug signatures restricted the same way.
- `D_h`: raw direction instability on C1's landmarks.

Every comparison between references uses `shRNA-h` as the paired comparator, and
every ρ(D, E) under a CRISPRi construction is reported beside the corresponding
ρ(D_h, E).

**Constructions of the CRISPRi direction.**

| name | construction |
|---|---|
| C0 | `03b_h3_crispri_ground_truth.py::build_crispri_signatures` on the scPerturb file: mean over perturbed cells minus the mean over all non-targeting cells, placed in the 978 landmarks with zeros where a landmark is absent. This produced the known values. |
| C0h | C0 restricted to the landmarks present in the scPerturb file, with drug signatures restricted to the same landmarks |
| C1 | the official Replogle K562-essential normalized pseudobulk, which is z-normalized against controls within each GEM group. The rows for the target gene are averaged, and the number of rows is recorded. The direction is restricted to landmarks present in the file, matched by symbol through the pinned gene-info file, and drug signatures are restricted to the same landmarks. **C1 is the primary CRISPRi construction for R1–R7.** |
| C1-RPE1, C1-GW | C1's procedure applied to the official RPE1-essential and K562 genome-wide pseudobulk files |

Four further rules apply to the gene space and to construction:
- A comparison between two files uses the landmarks present in both.
- Landmark symbols duplicated in a file's gene list are excluded and recorded.
- Where a gene has several rows in a file, its rows are averaged first, and the
  average is unit-normalized afterwards.
- Every direction is unit-normalized after restriction, and the shRNA reference
  uses all 978 landmarks unless it is the harmonized comparator.

**C1 analyses are void unless C1's landmark set holds at least 500 genes.** The
count is recorded by R0.2 before any statistic is computed. 500 is roughly half
the landmark space and is not derived from a power calculation. It is a usability
floor, not evidence that C1 covers the landmark space: the actual overlap is
reported wherever C1 results are, and if it is close to 500 the results are
described as applying to a substantially reduced landmark space.

**Cohorts.** All K562-essential and RPE1 comparisons use the same 131 drugs:

| cohort | definition |
|---|---|
| shRNA arm | the 795 drugs of the corrected H3 artifact |
| CRISPRi arm | the 131 drugs of the corrected CRISPRi records |
| matched | the 114 drugs in both |
| genome-wide arm | every drug in the pinned extraction with a target annotation, at least 5 cell lines, and a first-listed target perturbed in the genome-wide file |
| hematopoietic subset | CRISPRi-arm drugs with at least one hematopoietic and one non-hematopoietic signature |

**Own-target percentile.** Over the fixed set `T_X` of distinct targets in the arm:

    q_c^X = (|T_X| - 1)^-1 * sum over t in T_X, t != t(c), of
            [ 1{E_c^X(t) < E_c^X(t(c))} + 0.5 * 1{E_c^X(t) = E_c^X(t(c))} ]

`q̄^X` is the target-balanced mean: `q` is averaged within each target, and then
across targets with equal weight. If alignment carries no target identity, `q̄` has
expectation 0.5.

**Target-cluster bootstrap (TCB).** Each replicate draws the cohort's targets with
replacement and carries every drug of each drawn target, with duplicates retained.
- Statistics are recomputed, and ranks are recomputed, within each replicate.
- Paired comparisons compute both quantities within the same replicate.
- There are 10,000 replicates, with seed 20260922, and intervals are percentile
  intervals.
- A nonfinite replicate invalidates the run and is investigated, never discarded or
  redrawn.

**Unique-target permutation (UTP).** A uniformly random permutation `π` of the
cohort's unique targets assigns `u_{π(t(c))}` to every drug `c`.
- `P`, `E`, `s` and `q` are recomputed from the reassigned direction; `D` is fixed.
- There are 10,000 permutations, with seed 20260923.
- Every result reports the full null distribution, its median `ρ̃`, and
  `ρ_obs − ρ̃`.

**Target-balanced correlation.** `ρ_T` is the Spearman correlation across targets
of the within-target medians of the two variables.

**Comparison rule.** This compares an alternative estimate `ρ_alt` with a
reference estimate `ρ_ref` of the same association on the same drugs, using a
paired TCB. It applies only when `|ρ_ref| ≥ 0.15`.

| reading | condition |
|---|---|
| attenuated | `|ρ_alt| ≤ 0.5 |ρ_ref|` with the sign unchanged, **and** the paired 95% interval for `ρ_alt − ρ_ref` excludes 0 |
| reversed | the sign changes, **and** the paired 95% interval for `ρ_alt − ρ_ref` excludes 0 |
| retained | the sign is the same, `|ρ_alt| ≥ 0.5 |ρ_ref|`, **and** the 95% interval for `ρ_alt` excludes 0 |
| inconclusive | otherwise |

Every reading requires the paired interval, so a large apparent drop with an
interval that spans 0 is inconclusive rather than attenuated. Where a reading
below asks for attenuation, "reversed" satisfies it and is reported by name.

When `|ρ_ref| < 0.15`, only the paired difference and its interval are reported.

**Practical null.** An association is practically null when its 90% TCB interval
lies inside (−0.15, 0.15). In the shRNA arm, with 258 target clusters, that
condition is reachable. In the CRISPRi arm, with 41 unequal clusters, the interval
for a correlation is expected to be wider than 0.30, so the condition will rarely
be met there and a failure to meet it carries no information.

**An interval that includes 0 is never read as support for anything.** It is
reported as inconclusive unless the practical-null condition, or a stated
equivalence condition, is met.

## R0. Audit of the criterion measure

R0 carries no criterion that can upgrade a manuscript claim. It fixes which
CRISPRi construction the other analyses use, and it records what that construction
is. R0.1, R0.3 and R0.7 are audit records rather than analyses: they are written to
the output and carry no reading.

- **R0.1 Representation.** For the scPerturb file and every official file, record:
  - the dtype, minimum, maximum, fraction nonzero and fraction integer-valued of
    `X`, with single-cell files measured on 10,000 cells sampled with seed 20260925;
  - the layers, and whether `raw` is present;
  - the observation fields used.
- **R0.2 Gene space.** For each file, record:
  - landmarks matched by symbol;
  - landmarks left unmatched;
  - duplicated gene symbols;
  - a hash of the ordered landmark list used.
- **R0.3 Design,** from the single-cell files for K562-essential and RPE1: the
  control label and count, the GEM groups, the cells per target, and the number of
  GEM groups containing each target.
- **R0.4 Verification.** Compute the mean of the normalized single-cell file over
  each target's cells, and report the full per-target distribution of its cosine
  with C1, and its maximum absolute difference from C1. A failure voids nothing,
  but it changes what the reliability estimates describe: **if the median
  per-target cosine is below 0.99, R0.6 and R7c are labeled as describing the
  single-cell construction rather than C1.** Crossing 0.99 licenses no claim that
  C1 is validated; it fixes only which object the reliability estimates describe.
- **R0.5 Agreement.** Per target, compute the cosine between each pair of
  constructions on shared landmarks: C0 and C1; C1 and C1-RPE1; C1 and C1-GW.
  Report the distributions.
- **R0.6 Reliability.** Estimate split-half reliability per target:
  - For CRISPRi, randomly assign GEM groups to halves 50 times (seed 20260924).
    Take the median cosine between the half-means of the normalized single-cell
    data, for K562-essential and RPE1.
  - For shRNA, assign distinct hairpins (`pert_id`) to halves where a target has at
    least two, and otherwise split signatures at random, 50 times. Take the median
    cosine between half consensus directions.
  - A target with cells in fewer than two GEM groups has missing reliability.
  - Genome-wide reliability is not estimated; its single-cell file is 66 GB.
- **R0.7 Knockdown.** Report the target gene's own value in its C1 pseudobulk where
  the gene is measured, together with the released per-perturbation fields
  `fold_expr`, `energy_test_p_value` and `num_cells_filtered`. This is descriptive
  only.
- **R0.8 Concentration.** Compute the share of the leading eigenvalue in the
  cosine-similarity matrix among target directions, on the 33 matched targets,
  under C1, C1-RPE1 and shRNA. Estimate the difference C1 − shRNA with a TCB over
  the 33 targets. **CRISPRi directions are more concentrated than shRNA directions**
  if the 95% interval lies above 0.
- **R0.9 Headline correlations under each construction.** On the 131 drugs, compute
  ρ(P, E), ρ(D, E) and `q̄` under C0, C0h and C1, with TCB intervals. Compare C1
  with C0 by the comparison rule. The **raw reversal is construction-dependent** if
  ρ(D, E) under C1 reads attenuated relative to C0. It is
  **construction-robust** if it reads retained. Otherwise the reading is
  inconclusive. ρ(P, E) is compared the same way.

## R1. Dependence, weighting, composition (explanation D)

- **R1-dependence.** TCB 95% intervals for the drug-level ρ(P, E) and ρ(D, E), in
  the shRNA arm and in the CRISPRi arm under C1 and C0. These are estimation only.
  A widened interval is not evidence for D.
- **R1-weighting.** Compare `ρ_T` with the drug-level ρ by the comparison rule, for
  each of the correlations above.
  - Attenuated means the drug-level estimate is carried by heavily represented
    targets.
  - Retained counts against that.
- **R1-composition.** On the matched 114 drugs, report the full 2×2:
  - ρ(P, E) and ρ(D, E) under shRNA, under shRNA-h and under C1;
  - the paired differences C1 − shRNA-h for both, with TCB intervals over the 33
    targets.

  The comparator is shRNA-h, so the difference isolates the reference and not the
  gene space. The difference from the unharmonized shRNA quantities is reported
  beside it, and the gap between the two is the gene-space contribution.

  **The references differ on identical drugs** if a difference's 95% interval
  excludes 0; otherwise that comparison is inconclusive. The shRNA estimates on the
  matched drugs and on all 795 are reported side by side without a criterion,
  because the drug sets differ.
- **R1-influence** is reported under R7d.

## R2. Target identity of the alignment

For shRNA, C1, C1-RPE1 and C0, compute `q̄` with the UTP p-value, one-sided on
`q̄_b ≥ q̄_obs`, and with TCB 95% and 90% intervals.

| reading | condition |
|---|---|
| carries target identity | p < 0.01 **and** the 95% interval lies above 0.5 |
| lacks target identity | the 90% interval lies inside (0.40, 0.60) |
| inconclusive | otherwise |

The equivalence margin is 0.10 on each side rather than 0.05. If per-target `q`
were uniform, its standard deviation would be near 0.29, so with 41 targets the
standard error of `q̄` is near 0.045 and a 90% interval is about 0.15 wide. A
margin of 0.05 on each side could not be met in the CRISPRi arm at any outcome.
Even at 0.10 the reading is marginal there, and it is reachable in the shRNA arm.
The margin is prespecified and not derived from a power calculation.

Two further quantities are reported without a criterion:
- the Spearman correlation across targets between R0.6 reliability and target-mean
  `q`;
- `q̄` recomputed after `E` is z-scored across drugs within each target, which
  removes the advantage a target direction gains by carrying a component common to
  many drugs.

## R3. Permutation resistance (explanation A)

Run the UTP null for ρ(D, E) and ρ(P, E) in the CRISPRi arm, under C1 as primary
and C0 as secondary, and in the shRNA arm (unadjusted; the linked S2-T is the
adjusted version).

| reading | condition |
|---|---|
| assignment-dependent | ρ_obs lies beyond the null's 99th percentile in the direction of its sign, **and** `|ρ̃| ≤ 0.5 |ρ_obs|` |
| permutation-resistant | ρ_obs lies inside the central 95% of the null, **and** `ρ̃` has the sign of ρ_obs with `|ρ̃| ≥ 0.5 |ρ_obs|` |
| indeterminate | otherwise |

Permutation resistance alone is not attributed to any program.

One sensitivity is reported without a criterion:
- Common-direction removal. Let `v` be the leading right singular vector of the
  unit-normalized C1 directions of every perturbation in the K562-essential file,
  restricted to C1's landmarks. Project `v` out of the drug signatures and target
  directions, then recompute ρ(P, E) and ρ(D, E). This deliberately removes shared
  biology and is diagnostic only.

## R4. Essential-gene composition, from the genome-wide screen (explanation A)

Split the genome-wide arm into essential-screen targets (perturbed in the
K562-essential file) and other targets. For each subset, compute ρ(P, E) and
ρ(D, E) under C1-GW with TCB intervals. Estimate the difference essential − other
with a TCB that resamples targets within each subset independently.

**Essentiality and phenotype strength are confounded in this split.** Essential
knockdowns produce strong transcriptional phenotypes; many genome-wide
perturbations produce almost none, and their direction is then mostly noise, which
dilutes any correlation toward zero. The primary version of R4 therefore restricts
the other-target subset to targets with a transcriptional phenotype. The
restriction is fixed here, before any outcome is seen, and uses a field R0.7
already records.

A gene can hold several rows in the released file, so eligibility and construction
are both specified:

- A target enters the phenotype-positive subset when **at least one** of its
  released rows has `energy_test_p_value < 0.05`.
- Its direction is built from **only the qualifying rows**, averaged before unit
  normalization.
- The total and qualifying row counts are recorded per target.

Building the direction from qualifying rows alone keeps one phenotype definition
throughout: a target does not become phenotype-positive on the strength of one row
and then have that row averaged with unresponsive ones. Two further versions are
reported beside the primary: the direction built from all rows of the eligible
targets, and the unrestricted subset.

| reading | condition |
|---|---|
| raw association depends on essential composition | the 95% interval for the difference in ρ(D, E) lies above 0 |
| raw association practically null outside the essential set | the other-target subset's ρ(D, E) is practically null |
| H3 pattern outside the essential set | in the other-target subset, the ρ(P, E) 95% interval lies above 0 **and** ρ(D, E) is practically null |

**R4's readings apply only if the other-target subset holds at least 60 drugs over at
least 30 targets, after the energy-test restriction.** If the gate is not met, R4
is still run and reported, with the sentence: the registered minimum cohort was not
reached, so the estimates are descriptive and neither support nor count against an
essential-gene-composition explanation. R4 is never dropped after its count is
seen, and a failed sample-size gate is not a negative result.

Only the positive reading, H3's pattern outside the essential set, has consequences
for H3; a weak result in that subset is compatible with weak directions.
Genome-wide reliability is not estimated, and that limitation stands even with the
energy-test restriction.

On the 41 CRISPRi-arm targets, C1 and C1-GW are compared by the comparison rule
without a criterion. That comparison changes screen design while holding the cell
type fixed.

## R5. Toxicity (explanation A; diagnostic, not causal)

**Toxicity measure.** `T_c` is minus the median, across cell lines, of the PRISM
Repurposing 19Q4 primary-screen log-fold-change (2.5 µM, replicate-collapsed).
The fraction of lines with log-fold-change below −1 is a secondary summary without
a criterion.

**Mapping.** The only permitted sources of identifiers and names are the pinned
PRISM release files themselves. No external synonym dictionary, ontology or manual
alias is used, so the set of mapped drugs is a function of the pinned inputs alone.
- A drug is mapped first by Broad identifier stem.
- Failing that, it is mapped by exact name match after normalization, against the
  `Drug.Name` and `Synonyms` fields of the pinned compound list. Normalization is
  case folding, trimming, and collapsing internal whitespace, and nothing else.
- When several entries map, the entry with the most lines is used, with ties broken
  by lexical identifier order.
- Only entries with log-fold-change in at least 100 lines count.
- **A mapping table listing every accepted and rejected candidate, with the reason
  for each, is written and hashed before any log-fold-change value is summarized.**
  No match is revised afterwards.
- Coverage is reported by target, and by single- versus multi-target annotation.

**R5 is not run unless at least 60 CRISPRi-arm drugs over at least 20 targets map.**

On the complete-case cohort, under C1, with TCB intervals, R5 reports:
- ρ(E, T);
- the unadjusted ρ(D, E);
- the partial ρ(D, E given T), by the rank residualization of S1 with the single
  covariate `T`;
- `Δ = ρ(D, E) − partial`;
- the ratio of the partial to the unadjusted correlation, only if
  `ρ(D, E) ≥ 0.15`.

| reading | condition |
|---|---|
| toxicity accounts for the raw association | the ρ(E, T) 95% interval lies above 0, **and** either the partial correlation is practically null or it reads attenuated relative to the unadjusted correlation |
| counts against | ρ(E, T) is practically null, **or** the partial correlation reads retained |
| inconclusive | otherwise |

For the shRNA arm, ρ(E_shRNA, T) is reported without a criterion.

## R6. Cell context (explanation C)

**R6a, K562 against RPE1 (primary).** On the 131 drugs, using landmarks present in
both files, compute ρ(P, E), ρ(D, E) and `q̄` under C1 and C1-RPE1. Estimate the
paired TCB differences RPE1 − K562.

| reading | condition |
|---|---|
| raw reversal replicates in RPE1 | RPE1's ρ(D, E) reads retained relative to K562's |
| raw reversal is context-dependent | RPE1's ρ(D, E) reads attenuated relative to K562's |
| projected association is context-dependent | the 95% interval for RPE1 − K562 in ρ(P, E) lies above 0 |
| inconclusive | otherwise |

**R6b, the hematopoietic proxy: exploratory.** In the hematopoietic subset, compute
the target-balanced mean of `q_heme − q_nonheme` under C1, where `q_heme` uses the
mean of a drug's hematopoietic signatures and `q_nonheme` the mean of its other
signatures. A sensitivity analysis restricts the proxy to the six myeloid lines:
HL60, THP1, U937, NOMO1, PL21 and SKM1.

R6b carries no reading in either direction and is labeled exploratory wherever it
is reported. With 29 drugs over 22 targets it cannot change what the paper says,
and R6a is the registered test of cell context.

No K562 compound profile exists, so matched drug profiles cannot test cell type
directly.

## R7. Sensitivities, without a criterion except where gated

- **R7a Single-target annotations.** Compute ρ(P, E) and ρ(D, E) in each arm on
  single-target drugs, with TCB intervals. Estimate the single − multi difference
  with a TCB that resamples targets within each subset. No target is substituted.
- **R7b Signed alignment.** On drugs whose frozen MOA contains "inhibitor", report
  the target-balanced mean of `s` with a TCB interval, together with ρ(D, s) and
  ρ(P, s).
- **R7c Reliability.** Rerun the primary correlations and `q̄` on targets whose R0.6
  reliability is at least 0.5. This runs only if at least 20 CRISPRi and 100 shRNA
  targets remain.
- **R7d Leave one target out.** Re-estimate every primary correlation with each
  target omitted. Report the range, the target producing the largest change, and
  whether the sign is stable.
- **R7e Exclusions.** Report the eligible drugs lacking each reference, counted by
  target, together with the reliability and knockdown of included and excluded
  targets where available. No drug is filtered on alignment.
- **R7f The target's own gene.** Both references lower the target's own transcript,
  and in single-cell CRISPRi that coordinate is often the largest in `u`. A drug
  that inhibits a protein need not lower its transcript at all, so `E` partly
  measures a coordinate the drug cannot move, to a degree that may differ between
  references. Recompute `P`, `E`, `q̄` and the primary correlations with the
  target's own gene removed from both `u` and the drug signatures, under every
  reference, and report them beside the primary values.

  For each reference and each harmonized comparison, the target coordinate is
  removed only where it is present in that analysis gene space. Paired comparisons
  fix the shared gene space before removal. Where the target is absent, the
  sensitivity equals the primary analysis and is recorded as such. The numbers of
  drugs and targets whose vectors change are reported.

## Interpretation grid

| outcome | reading |
|---|---|
| R0.9: raw reversal construction-dependent | the raw half of the reversal is attributable to how the CRISPRi direction was built; the remaining analyses describe the validated C1 reference |
| R2 under C1 carries target identity, and the 95% interval for ρ(P, E) under C1 does not lie above 0 | H3 is not replicated against a CRISPRi reference that carries target identity |
| R2 under C1 carries target identity, and ρ(P, E) under C1 is practically null | H3's proposed generalization is contradicted by that reference, without any reversed association being established |
| R2 under C1 carries target identity, and the 95% interval for ρ(P, E) under C1 lies below 0 | that reference establishes a directional contradiction: the association runs the other way |
| R3 under C1: raw permutation-resistant; R2 under C1 does not carry identity; R0.8 concentrated | the raw CRISPRi association reflects a component shared across the arm's target directions, not target identity |
| the above, together with R5 accounting for the association | that shared component is statistically accounted for by broad cytotoxicity |
| R4: raw association depends on essential composition | the raw association is specific to essential-screen targets |
| R4: H3 pattern outside the essential set | H3's pattern replicates against a CRISPRi reference outside the essential-gene set |
| R6a: raw reversal replicates in RPE1 | the reversal is not specific to K562 |
| R1-weighting attenuated | the drug-level estimates are carried by heavily represented targets, and the target-balanced estimates govern claims about targets in general |
| linked S2-T fails | H3 cannot distinguish target matching from shared-axis coupling, and the separation claim is withdrawn regardless of anything here |

## Sample size

| cohort | drugs | targets |
|---|---|---|
| shRNA arm | 795 | 258 |
| CRISPRi arm, K562 and RPE1 | 131 | 41 |
| matched | 114 | 33 |
| hematopoietic subset | 29 | 22 |
| genome-wide arm | counted at build; at least 60 drugs over 30 targets in the other-target subset for R4's readings | — |
| R5 complete-case cohort | counted at mapping; at least 60 drugs over 20 targets | — |

These sizes are fixed by the data and are not a design choice.
- Intervals in the CRISPRi arm rest on 41 unequal clusters. A correlation's 95%
  interval there is expected to be wider than 0.30, and the 90% interval for `q̄`
  about 0.15 wide.
- With 258 targets, the shRNA arm's intervals are narrow enough for the
  practical-null and equivalence conditions to be met.
- The hematopoietic subset, at 29 drugs, is exploratory.
- No analysis supports a claim of no difference except through a stated
  equivalence condition, and an inconclusive reading in the CRISPRi arm is
  expected often enough that it must not be reported as a failure.

## Maximum claim under this registration

For each analysis, the paper may state the reading it returned, with its interval,
in the terms of the interpretation grid. It may describe H3 as:

> Registered shRNA criterion met after correction; convergent validity against the
> CRISPRi reference unresolved

Where the grid supports it, the paper may replace that description with the grid's
reading. The broadest construct the paper may name is similarity to a genetic
loss-of-function axis measured under a stated reference.

The paper may not:
- name a single cause of the reversal;
- call either reference a ground truth;
- say that either statistic measures on-target consistency or biochemical target
  engagement;
- say that raw instability carries no target-related information;
- attribute the CRISPRi association to cytotoxicity unless R5 ran and returned
  "toxicity accounts for the raw association";
- say that the K562 cell type is confirmed, or draw any conclusion from R6b;
- say H3 is contradicted by a reference where the grid's reading is "not
  replicated", which requires only that the interval fail to lie above 0.

Everything else is exploratory and labeled as such.

## Inputs

| input | pin |
|---|---|
| corrected CRISPRi records | sha256 `f00f8428071178bd2317139c5a1f534e458931684fb6da65c750f9f1e5456eb3` |
| corrected H3 artifact | sha256 `fd69e26fc9a3917323065b631688baeab8b283f735c8bf5b16210ba67bd21425` |
| `lincs_subset.npz` | sha256 `2ad0f5d30ab826f9ec0cfe37f6b53b2829bfef1d73b7920c3441de6407adb4ec` |
| `lincs_shrna.npz` | sha256 `4a990e5072a43f59fdda60a2ff040f355f06be947c14bcfa87056c66332b58e7` |
| scPerturb `ReplogleWeissman2022_K562_essential.h5ad`, Zenodo 10044268 | sha256 `412fd0df8c4ccea9f4db91cd88033c49200838b29d40945e48574be588b48789` |
| `K562_essential_normalized_bulk_01.h5ad`, figshare 20029387 | sha256 `c1ca6456c9c9f1aa2b02c496eb64d1dc3e6a852edbd744d682b8d2c95fd36829` |
| `rpe1_normalized_bulk_01.h5ad`, same record | sha256 `a3c5bfd0f15d63938bc80c9b8874b9cd761e3a23caf5ffe7966bae4e887ec89d` |
| `K562_gwps_normalized_bulk_01.h5ad`, same record | sha256 `37e48c474d8b5dead4151f96ea8f5fe7bbe6beb10eeea48685b740c3f74490a2` |
| `K562_essential_normalized_singlecell_01.h5ad`, same record | md5 `f1e221fbf6eac774c21c4242ed440c3f` as published; sha256 recorded at first retrieval |
| `rpe1_normalized_singlecell_01.h5ad`, same record | md5 `2c36a053960f3fae157adacdbccd4485` as published; sha256 recorded at first retrieval |
| `GSE92742_Broad_LINCS_sig_info.txt.gz` | sha256 `19da29c0ee12ddf27f9698cd0da40beaff58657dcde9d382aae068737e831299` |
| `GSE92742_Broad_LINCS_gene_info.txt.gz` | sha256 `741216ccc53320119b47ab006de3bcad48963c57087c9e07f50f0d6cd088711a` |
| `GSE92742_Broad_LINCS_cell_info.txt.gz` | sha256 `88c6ca5ff8bb3cf3af22c2eb8cbe769f9ae8e305a760959b18711877d730d61c` |
| `lincs_shrna_siginfo.csv.gz` | sha256 `bd396fa0e1a2f00c1b5f2c8d2b35f9a056f5e5353382475655869038037ec014` |
| `frozen_drug_labels.json` | sha256 `e6d73384c5bbc501ee15f621aaa258f02f38ee0d7a9603b6c820bec3e0489d52` |
| `geometry/direction_instability.py` | sha256 `e827ea0a739b207f9a736063862689b23161def07b929bbab91f973154184064` |
| PRISM Repurposing 19Q4 primary screen, log-fold-change and treatment files | by release and filename at freeze; sha256 recorded at first retrieval |

Hematopoietic lines are those whose `primary_site` in the pinned cell-info file is
"haematopoietic and lymphoid tissue" or "blood".

## Design and procedure

The work is split across two scripts:

- **`experiments/03d_h3_reference_discordance.py`** builds every cohort and
  construction from the pinned inputs, including the harmonized comparators
  `shRNA-h` and `D_h` and the target-gene-removed variant of R7f. Before any
  analysis, it checks that its C0 and shRNA quantities reproduce the corrected
  CRISPRi records to within 1e-6 for `D`, `P` and `E`, and exactly for target and
  `n_celllines`; both are computed in float64 from the same extraction, so the
  tolerance is not the float32 one used by the linked reconstruction gate. It then
  runs R0.1, R0.2, R0.5 and R0.7 to R7.
- **`experiments/modal_03d_single_cell.py`** runs R0.3, R0.4 and R0.6 on one Modal
  CPU worker. It checkpoints within each file by target batch and commits the
  volume after each batch.

**R0–R7 are void unless that reproduction holds.**

Both scripts are tested on synthetic data and reviewed before freeze, and their
commit is pinned here. Each run writes the full result structure, seeds, input
hashes, mapping tables, replicate draws and null distributions to
`results/03d_h3_reference_discordance/`. The analyses run once.

- **Existing data:** yes, and every input is pinned above or at first retrieval.
- **Data collection:** N/A — no new data are collected.
- **Blinding:** none is possible, because the reversal is known.
- **Missing data:** drugs or targets that fail to map, reproduce or meet a gate are
  counted and reported, never replaced.
- **Outliers:** none removed.
- **Exploratory analyses:** any analysis not named here is exploratory.
