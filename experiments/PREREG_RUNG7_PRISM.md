# Pre-registration: Rung 7 — Holdout Prediction via PRISM

**Date:** 2026-07-09
**Status:** PRE-REGISTERED (frozen before any PRISM data is downloaded or analyzed)
**Intended use:** Future study, NOT included in the current BMC MRM submission.
**Commit SHA:** 468cfdd

## Rationale

Rungs 1-6 of the validity ladder are tested within or across LINCS L1000,
Perturb-seq, and JUMP-CP. Rung 7 requires predicting an external outcome
in a genuinely independent dataset. PRISM (Profiling Relative Inhibition
Simultaneously in Mixtures; Corsello et al. 2020, Nature Cancer) provides
drug sensitivity measurements across ~900 cancer cell lines using pooled
barcoded screening — a completely independent assay modality from LINCS
transcriptomics (gene expression), with no shared parameters or data.

## Estimand bridge

Direction instability measures whether a drug's *transcriptomic* response
points in the same direction across cell lines. PRISM measures whether a
drug *inhibits growth* across cell lines. These are different modalities.
The bridge hypothesis is: drugs with conserved transcriptomic mechanism
(low DI) act through a pathway that is consistently engaged across cell
types, and consistently engaged pathways produce more predictable
sensitivity profiles. This is a hypothesis, not a given — the link
between transcriptomic direction conservation and viability consistency
is indirect and may be weak. A null result would indicate that
transcriptomic geometry does not predict viability consistency, which is
a genuine scope limitation of direction instability, not a failure of
the validity framework.

## Dataset

**PRISM Repurposing Secondary Screen** (DepMap portal, public)
- ~1,500 drugs screened across ~900 cell lines
- Measurement: log2 fold-change viability (lower = more sensitive)
- Drug identifiers: Broad IDs, matchable to LINCS via Drug Repurposing Hub

**LINCS validity scores** (already computed, frozen in results/)
- Direction instability (D), transport-stable bracket (TS), phenotype-projected
  bracket, for 8,949 drugs

## Drug matching

Match LINCS drugs to PRISM drugs via Drug Repurposing Hub broad_id
(the same source used for H3 target annotations). The matching procedure:

1. Load LINCS drug list (pert_iname from sig_info)
2. Load PRISM drug list (compound names / broad_ids from DepMap)
3. Join on Drug Repurposing Hub broad_id (exact match only, no fuzzy matching)
4. **Report matched-N before computing any correlation.** If matched-N < 200,
   flag as underpowered and report but do not interpret.
5. Retain only drugs with PRISM data in >= 50 cell lines

No manual curation of drug name synonyms. If the automatic broad_id join
misses drugs, they are missed — no post-hoc recovery.

## Predictions

### H7a: Direction instability predicts PRISM sensitivity consistency

**Prediction:** Drugs with low direction instability (conserved transcriptomic
mechanism across cell lines) show more consistent sensitivity profiles
across PRISM cell lines.

**Operationalization:**
- PRISM sensitivity consistency = 1 - coefficient of variation of
  log2 fold-change across cell lines (higher = more consistent)
- Compute Spearman correlation between LINCS direction instability and
  PRISM sensitivity consistency

**Criterion:** Spearman |rho| > 0.10 and p < 0.05/3 (Bonferroni-corrected).

**Direction:** Negative correlation expected (low DI → high consistency).

### H7b: Transport-stable bracket outpredicts raw bracket for PRISM consistency

**Prediction:** Transport-stable bracket (TS) correlates more strongly with
PRISM sensitivity consistency than raw direction instability (D).

**Operationalization:**
- Compute Spearman rho(TS, PRISM_consistency) and rho(D, PRISM_consistency)
- Compare via Steiger's test for dependent correlations

**Criterion:** |rho_TS| > |rho_D| with Steiger p < 0.05.

### H7c: Phenotype-projected bracket predicts PRISM target-specific sensitivity

**Prediction:** Among drugs with annotated targets (Drug Repurposing Hub),
phenotype-projected bracket correlates with the fraction of cell lines
where the drug's sensitivity is explained by expression of its target gene.

**Operationalization:**
- For each drug with an annotated target, compute Spearman correlation
  between target gene expression (DepMap Expression) and PRISM sensitivity
  across cell lines
- This yields a per-drug "target-expression-sensitivity" score
- Compute Spearman correlation between phenotype-projected bracket and
  this target-expression-sensitivity score

**Criterion:** Spearman |rho| > 0.10 and p < 0.05/3 (Bonferroni-corrected).

## Multiplicity correction

H7a, H7b, and H7c form one family. Bonferroni correction at alpha = 0.05
within this family (effective threshold p < 0.0167 for H7a and H7c;
H7b uses Steiger's test at nominal 0.05).

## Combination rule

- **Rung 7 confirmed:** H7a passes AND at least one of H7b/H7c passes
  (the primary prediction must hold; secondary predictions strengthen it).
- **Partial support:** H7a passes alone, or H7b/H7c pass without H7a.
  Report as suggestive.
- **Rung 7 null:** H7a fails. Report as informative null regardless of
  H7b/H7c. Interpret as a scope limitation of transcriptomic geometry
  for viability prediction.

## Null interpretation commitment

If H7 fails, the result is reported as "Rung 7: tested, not satisfied"
in any future paper that includes this analysis. The experiment is NOT
dropped or reverted to "not tested." Once run, it is reported — pass
or fail. This commitment is the purpose of pre-registration.

## What this does NOT test

- Clinical efficacy or in-vivo translatability
- Whether the validity ladder improves drug discovery outcomes
- Whether PRISM viability is the optimal external outcome
  (it is one measurable external outcome; others may be more appropriate)
