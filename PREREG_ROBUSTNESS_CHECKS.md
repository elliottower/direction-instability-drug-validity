# Pre-registration: Robustness Checks for BMC MRM Revision

**Date:** 2026-07-08
**Author:** Elliot Tower
**Context:** Addressing reviewer concerns on bracket_norm_validity_bmcmrm_v3.tex.
These are robustness/sensitivity analyses on already-computed results.
All outcomes will be reported regardless of direction.

---

## R2: Lambda sensitivity analysis for H5 (transport-stable bracket)

**Motivation:** The penalty weight λ=1.0 in TS = D − λ·Var_Fréchet was
inherited from the companion paper. A reviewer asks whether the H5 result
(66/66 folds, Wilcoxon p=8.2e-13) is robust to λ choice.

**Analysis:** Re-run the full 66-fold leave-one-cell-line-out cross-validation
at λ ∈ {0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0}.

**Sanity check (λ=0):** At λ=0, TS ≡ raw DI, so TS should win ~50% of folds
(tie-breaking noise). If λ=0 shows TS winning ≥80% of folds, there is a bug
in the computation (e.g., ties broken in TS's favor). Assert this before
interpreting other λ values.

**Report:** For each λ: number of folds TS wins, Wilcoxon p-value, mean Δρ.
Present as a table or supplementary figure.

**Success criterion (pre-registered):** TS wins ≥80% of folds AND Wilcoxon
p < 0.05 for at least 5 of the 8 non-zero λ values tested. If fewer than 5
pass, we report this as evidence that H5 is λ-sensitive and the specific
λ=1.0 result should be interpreted cautiously. (λ=0 is a sanity check only
and is excluded from this count.)

**No-go criterion:** If only λ=1.0 passes, we add a limitation paragraph
stating the result is not robust to penalty weight choice.

---

## R3: HDAC-removal sensitivity for H3 (phenotype projection)

**Motivation:** Pan-HDAC inhibitors are the calibration standard for direction
instability. A reviewer asks whether the H3 correlation (ρ=0.376) is driven
by these drugs.

**Analysis:**
- (a) Remove all drugs with "HDAC" in the target field from the 795-drug H3
  dataset and recompute Spearman ρ between projected_bracket and
  on_target_enrichment.
- (b) Remove the top 20 drugs with lowest raw_bracket (most directionally
  consistent drugs regardless of class) and recompute.
- (b2) Remove the top 20 drugs with highest projected_bracket (most
  phenotype-aligned drugs) and recompute. This directly addresses Reviewer 3's
  concern: whether the projected correlation is carried by a handful of
  ultra-clean drugs rather than a broad signal.
- (c) Bootstrap 95% CI (1000 resamples) on the full 795-drug set.

**HDAC filter validation:** The HDAC filter uses substring match on the
free-text target field ("HDAC" in target.upper()). Before computing:
1. Assert the number of removed drugs is in [10, 40]. If outside this range,
   the filter is silently under- or over-selecting and results are invalid.
2. Log the name and target annotation of every removed drug so a reviewer
   can verify the filter caught the intended class.
3. Check for known synonyms ("histone deacetylase", "class I HDAC",
   "pan-HDAC") and log any drugs matching these patterns that were NOT
   caught by the primary filter.

**Report:** All four sensitivity analyses (a, b, b2, c) reported in full.

**Success criterion for (a):** ρ > 0.20 after HDAC removal. If ρ drops
below 0.20, we report that H3 is partially driven by the HDAC calibration
standard and add a limitation.

**Success criterion for (b):** ρ_proj > 0.20 after removing top-20 by
raw_bracket. If ρ drops below 0.20, we report that H3 is driven by
outlier-consistent drugs.

**Success criterion for (b2):** ρ_proj > 0.20 after removing top-20 by
projected_bracket. If ρ drops below 0.20, we report that H3 is driven by
a small number of phenotype-aligned drugs rather than a broad signal.

---

## R6: JUMP-CP stratification by source count

**Motivation:** The 78% attrition retains well-replicated compounds. A reviewer
asks whether ρ=0.987 holds for marginal compounds (5-6 sources) vs
well-replicated ones (9-10 sources).

**Analysis:** If Experiment 8 raw data (per-compound DI with source counts) is
available, stratify into bins: 5-6 sources, 7-8 sources, 9-10 sources.
Compute rank correlation between raw and corrected DI within each bin.

**Report:** Per-bin ρ values and N. If data is unavailable, report this as a
limitation that cannot be addressed without re-running the JUMP-CP pipeline.

**Note:** Experiment 8 code does not exist in this repo. If the underlying
per-compound data is not available, we will add this as an explicit limitation
rather than running new experiments.

---

## R8: Concrete worked drug pair example

**Motivation:** The "warrant upgrade" discussion uses hypothetical Drug X/Drug Y.
A reviewer asks for a real non-HDAC example.

**Analysis:** From the 795-drug H3 dataset, identify two drugs with:
- Similar raw_bracket (within 0.05)
- One has high projected_bracket AND high on_target_enrichment (passes Rung 5)
- One has low projected_bracket AND low/zero on_target_enrichment (fails Rung 5)
- Neither is HDAC-class

**Report:** The best candidate pair with full metadata. If no clean pair exists
(all similar-DI drugs have similar projection), report this as evidence that
the distinction is distributed rather than item-level.

**Selection criterion:** Top pair ranked by |Δprojected_bracket| with
|Δraw_bracket| < 0.05. No cherry-picking on drug identity.

---

## Writing-only changes (no pre-registration needed)

The following are purely expository changes that do not involve new data analysis:

- R1: Messick facet mapping table (theoretical, no data)
- R4: Clarify K=2 limitation in Perturb-seq section (rewrite, no new analysis)
- R5: Convert DEVIATION_LOG.md to supplementary table + add reporting guideline
  statement (reformatting existing content)
- R7: Formalize causal estimand and discuss DAG limitations (theoretical)
- R9: Refine scorecard table columns (reformatting existing results)
- R10: Expand Related Work section (literature, no data)

---

## Commitment

All quantitative results from R2, R3, R6, R8 will be reported in the paper
regardless of whether they support or weaken the original claims. Negative
results will be reported as limitations.
