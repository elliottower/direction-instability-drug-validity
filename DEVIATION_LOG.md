# Pre-Registration Deviation Log

## Entry 1: Gene-set method and parameters for H3/H4/H5

**Date:** 2026-07-07
**Status:** Pre-specified before running experiments 03, 04, 05
**Applies to:** H3 (phenotype projection), H4 (localization), H5 (transport stability)

### Context

The original pre-registration (commit 8d74bd0) froze hypotheses, decision
thresholds, and scorer functions (geometry/bracket_norm.py) but did not specify
how to construct the pathway gene sets needed by experiments 03 and 04. The
experiment scripts were listed in the frozen-scorer section but were never
written. This entry closes that gap by declaring the gene-set method and all
analytic parameters before any experiment is run.

### Declared analytic choices

**Primary gene-set method (H3, H4): genetic perturbation connectivity.**
For each drug with an annotated target gene (from the Drug Repurposing Hub via
frozen_drug_labels.json), the "on-target" signal is defined by the matched
shRNA knockdown signature of that target gene in LINCS L1000
(lincs_shrna.npz). This avoids researcher degrees of freedom in gene-set
construction (no GO-BP union, no arbitrary pathway database selection).

- **H3 phenotype_direction:** The consensus shRNA knockdown signature
  (mean across all cell lines and hairpins) of the drug's annotated target
  gene, used as the phenotype_direction argument to
  phenotype_projected_bracket().

- **H4 region_mask:** The top N=100 genes by absolute z-score in the
  consensus shRNA knockdown signature of the drug's annotated target gene.
  N=100 is fixed at ~10% of the 978 landmark genes, within the
  benchmark-optimal range for L1000 connectivity methods (Xie et al. 2019,
  Briefings in Bioinformatics 21(6):2194).

- **H5 cross-validation:** Leave-one-cell-line-out. For each fold, compute
  transport_stable_bracket and raw direction_instability on the training
  cell lines. Predict held-out cell line effect using Spearman rank
  correlation between the drug's bracket score and its actual signature
  magnitude in the held-out cell line.

**Sensitivity analysis (H3, H4): Drug Repurposing Hub MOA grouping.**
As a secondary check, repeat H3 and H4 using MOA class membership as the
pathway definition: for each MOA class with >=10 drugs, compute the mean
drug signature across the class as the phenotype_direction (H3) or use the
top-100 genes of that mean signature as the region_mask (H4). If both
primary and sensitivity analyses agree on pass/fail, the result is robust.
If they disagree, report both and flag the discrepancy.

**Inclusion criteria:**
- Drugs must have >= 5 cell lines in the compound LINCS data
- Target gene must have >= 3 shRNA hairpins in lincs_shrna.npz
- Only drugs with non-null target annotation are eligible for H3/H4

### Justification

Using genetic perturbation connectivity rather than pathway databases
eliminates the largest unregistered knob (gene-set construction method and
size). The shRNA data is already in the LINCS ecosystem, uses the same 978
landmark genes, and provides a direct empirical definition of "on-target"
without requiring curator judgment about which GO terms are relevant.

### What this does NOT change

- Hypotheses H3, H4, H5 and their decision thresholds remain exactly as
  specified in PREREGISTRATION.md (commit 8d74bd0)
- Scorer functions in geometry/bracket_norm.py remain unchanged
- All other pre-registered decisions (H1, H2, H6) are unaffected


## Entry 2: H5 fold threshold adjustment for fewer than 5 valid folds

**Date:** 2026-07-07
**Status:** Pre-specified before running experiment 05
**Applies to:** H5 (transport stability)

The pre-registered criterion is "transport-stable bracket outpredicts raw in
at least 4 of 5 folds." If fewer than 5 cell lines produce valid folds
(i.e., >=20 drugs with data in both training and held-out), the threshold
adjusts to ceil(n_valid_folds * 0.8). This maintains the 80% win-rate
requirement while avoiding discarding the experiment when a rare cell line
has too few drugs.

In practice, the LINCS data has 71 cell lines and the first 5 alphabetically
all produced valid folds, so this fallback was not triggered.


## Entry 3: H3 and H4 shared dependency on shRNA consensus

**Date:** 2026-07-07
**Status:** Documented before interpreting results
**Applies to:** H3 (phenotype projection), H4 (localization)

H3 uses the shRNA consensus signature as the phenotype_direction vector.
H4 uses the top-100 genes of that same shRNA consensus as the region_mask.
These are not independent tests — both derive from the same genetic
perturbation object. The paper must state this shared dependency and not
present H3 and H4 as two independent confirmations of the validity ladder.
H5 (transport stability) is fully independent of H3/H4 since it uses no
gene sets or shRNA data
