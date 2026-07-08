# Pre-registration: Extended Experiments for PSB Submission

Pre-registered before any real data is analyzed for these experiments.
Commit SHA of this file serves as the timestamp.

## Experiment H5-full: Full Leave-One-Out Cross-Validation

### Motivation

The original H5 tested 5 alphabetically-selected cell lines (A375, A549,
A673, AGS, ASC). Reviewers requested full leave-one-out over all ~50 LINCS
cell lines to test whether the rare-cell-line sign reversal generalizes
beyond the alphabetical selection.

### Hypothesis

**H5-full:** Transport-stable bracket outpredicts raw bracket (higher
Spearman correlation with held-out cosine consistency) in at least 80% of
valid cell-line folds.

The 80% threshold matches the original pre-registered criterion (4/5 = 80%).

### Definitions (frozen before running)

- **Eligible drug (per fold):** A drug that (a) has a signature in the
  held-out cell line AND (b) has signatures in >= 4 other cell lines
  after removing the held-out cell line (preserving the K >= 5 inclusion
  criterion from the main study).
- **Valid fold:** A held-out cell line with >= 20 eligible drugs.
  Folds with < 20 drugs are excluded (too few for a stable Spearman
  correlation) and reported as excluded.
- **Rare vs. common cell line:** Continuous variable (number of eligible
  drugs per fold). The exploratory analysis uses Spearman correlation
  between fold drug-count and TS advantage (Delta-rho), not a binary
  split. No "rare" threshold is imposed.

### Decision criteria

- **CO-PRIMARY 1:** Fraction of valid folds where TS outpredicts raw >= 0.80.
- **CO-PRIMARY 2:** Wilcoxon signed-rank test on per-fold Delta-rho
  (rho_TS - rho_raw) across all valid folds, one-sided (H_a: median
  Delta-rho > 0), p < 0.05. This does not assume fold independence
  (folds share the drug population; the signed-rank test is on the
  paired differences, not a binomial).
- SECONDARY: Mean Delta-rho across folds, reported with 95% bootstrap CI.
- EXPLORATORY: Spearman correlation between fold drug-count (continuous)
  and fold Delta-rho. We predict a negative correlation (rarer cell lines
  show larger TS advantage), consistent with the original 5-fold finding.

### Method

Identical to Experiment 05 (05_transport_stability.py) except:
- All cell lines with >= 20 eligible drugs are used as held-out folds
  (not just the first 5 alphabetically)
- Per-fold: compute raw DI and TS on training set, correlate with
  held-out cosine consistency
- No permutation null per fold (too expensive for ~50 folds; the
  original 5-fold permutation null already established that the TS
  advantage is not algebraic, p_perm < 10^-3)

### Data

Same LINCS L1000 Level 5 signatures (GSE92742) used in all drug experiments.
Already on disk at drug-perturbation-geometry/data/lincs_subset.npz.

---

## Experiment H3-CRISPRi: Convergent Validity Check with CRISPRi Ground Truth

**This is a convergent-validity check, not a blind confirmatory test.**
The comparison value (rho_shRNA = 0.376) is known from the original H3
analysis. This experiment tests whether an independent perturbation
modality (CRISPRi) converges on the same finding, using a different
ground-truth source with different noise characteristics.

### Motivation

The original H3 used shRNA knockdown signatures as the "on-target direction"
for phenotype-projected bracket. shRNA has well-documented off-target effects
that add noise. CRISPRi (CRISPR interference) provides cleaner loss-of-function
signatures. This experiment tests whether CRISPRi ground truth produces a
consistent or stronger H3 result.

### Hypothesis

**H3-CRISPRi:** Spearman correlation between phenotype-projected bracket
(using CRISPRi signatures as the on-target direction) and on-target
enrichment is positive and > 0.

### Decision criteria (three-outcome logic)

**(a) rho_CRISPRi >= 0.376 (original shRNA value):** Convergent validity
confirmed. CRISPRi ground truth produces the same or stronger phenotype
projection signal. shRNA noise likely attenuated the original result.

**(b) 0 < rho_CRISPRi < 0.376:** Both ground truths produce positive
correlations — robustness confirmed. The attenuation most likely reflects
the CRISPRi K562-only confound: CRISPRi signatures derive from a single
cell type (K562), while shRNA consensus averages across ~20 cell lines.
A cell-type-specific on-target direction is a noisier match for drugs
tested across many cell lines, predicting a lower but still positive
correlation. This outcome supports H3's validity — two independent
perturbation modalities converge on the same qualitative finding.

**(c) rho_CRISPRi <= 0:** Interpret via the K562-only confound first:
CRISPRi signatures from a single hematopoietic cell line may not capture
the on-target direction for drugs acting in epithelial or mesenchymal
contexts. This interpretation is tested by the secondary criterion below.
Only if the secondary criterion also fails (raw bracket also uncorrelated)
should shRNA inflation be considered.

- **PRIMARY:** rho_CRISPRi > 0 (one-sided test, p < 0.05).
  Outcome (a), (b), or (c) is then determined by the point estimate.
- **SECONDARY:** Raw bracket correlation with CRISPRi enrichment should
  remain near zero (|rho| < 0.15), same as the shRNA case. If raw bracket
  correlates positively with CRISPRi enrichment but projected bracket does
  not, this would suggest CRISPRi ground truth is confounded rather than
  the shRNA result being inflated.
- **MATCHED COMPARISON:** Report both shRNA and CRISPRi results side by
  side for the intersection of drugs with both ground truth types.
  Differences are interpreted in the context of the K562-only confound.

### Method

1. Extract per-gene CRISPRi expression signatures from Replogle et al.
   2022 Perturb-seq data (K562 essential screen). For each gene, compute
   the mean expression change across perturbed cells vs control cells
   over the LINCS landmark genes.

2. Gene mapping: use the LINCS gene_info file
   (GSE92742_Broad_LINCS_gene_info.txt.gz) to map Entrez IDs to gene
   symbols. Only genes present in BOTH the Perturb-seq variable names
   AND the LINCS landmark panel (via symbol matching) are used.
   Report the number of shared genes and the number of eligible drugs
   after mapping.

3. On-target enrichment uses cosine^2 (sign-invariant), so a sign flip
   between CRISPRi repression and drug inhibition does not affect the
   enrichment measure. Phenotype-projected bracket uses the signed
   cosine projection, which is appropriate because both CRISPRi and
   drug inhibition are loss-of-function perturbations (directionally
   aligned).

4. For each drug with an annotated target gene (Drug Repurposing Hub),
   if the target gene has a CRISPRi signature, compute:
   - On-target enrichment: cosine^2 between drug mean signature and
     CRISPRi signature of target gene
   - Phenotype-projected bracket: projection of instability vector
     onto CRISPRi direction

5. Compute Spearman correlation between phenotype-projected bracket
   and on-target enrichment, separately for CRISPRi and shRNA ground
   truths.

### Data

- Drug signatures: LINCS L1000 (same as all drug experiments)
- CRISPRi signatures: ReplogleWeissman2022_K562_essential.h5ad
  (already on disk at genetic-perturbation-holonomy/data/)
- Drug target annotations: Drug Repurposing Hub (same as original H3)
- Gene name mapping: LINCS GSE92742_Broad_LINCS_gene_info.txt.gz
  (pr_gene_id -> pr_gene_symbol for landmark genes)

### Known limitations

- CRISPRi signatures are from K562 only (one cell type), while shRNA
  consensus signatures average across ~20 cell lines. The CRISPRi
  "on-target direction" is cell-type-specific, which likely reduces
  the correlation for drugs acting in non-hematopoietic contexts.
  This is the most probable explanation for outcome (b) or (c) and
  should be interpreted before considering shRNA inflation.
- Gene name matching between LINCS landmark genes and Perturb-seq
  gene names may reduce the eligible drug set. The number of shared
  genes and eligible drugs is reported before any correlation analysis.
- CRISPRi (transcriptional repression) and drug inhibition are both
  loss-of-function but may produce quantitatively different downstream
  effects for the same target gene. Cosine^2 enrichment is
  sign-invariant, mitigating directional mismatches.

---

## NOT pre-registered for PSB (future work)

### sci-Plex (dose-as-context)

sci-Plex (Srivatsan et al. 2020, Science) profiles ~188 compounds
across 3 cell lines at multiple doses via single-cell transcriptomics.
With only K=3 cell lines, this falls below the K>=5 threshold used in
LINCS and JUMP-CP, making cross-cell-line DI poorly estimated.

The better application: use the DOSE AXIS as the context dimension.
For a drug tested at doses d1, d2, ..., dk, compute direction
instability across doses within a single cell line. A drug with a
consistent mechanism should rotate its signature direction minimally
across doses (scaling magnitude but preserving direction), while a
drug that activates different pathways at different doses should show
high direction instability across doses.

This reframes DI from "consistency across cell lines" to "consistency
across doses" — a novel application that sci-Plex is uniquely suited
for and that addresses a different biological question (dose-dependent
mechanism switching vs. cross-context mechanism conservation).

Planned for the post-PSB journal expansion.
