# Pre-Registration: Combined Paper Experiments

**Status:** DRAFT — awaiting review before commit freeze.
**Date:** 2026-07-05
**Context:** Three experiments for the combined direction-instability + magnitude-confound paper. These build on two existing preprints (Zenodo, both July 5 2026):
  - "Direction Instability Predicts Cross-Cell-Line Drug Mechanism Transport in LINCS L1000"
  - "Variance- and Subspace-Based Perturbation Divergence Metrics Are Magnitude-Confounded"

**Integrity protocol:** Pre-registration document, analysis code, and decision criteria will be committed together in a single commit before any of the three experiments are run on real data. The commit SHA will be recorded here after freeze. Scripts will be complete and runnable before commit; no placeholder functions, no `TODO` blocks, no data-peeking in commented code. This inherits the standard established in the confound paper's frozen commits (249abaf, 98897ed, 92e3276).

**Commit SHA:** 8d74bd0

---

## Experiment 1: Cosine direction instability on Perturb-seq same-gene transport

### Background

The confound paper demonstrated that geodesic distance (d = -0.14) and earth mover's distance (d = -0.172) showed apparent same-gene transport signal between K562 and RPE1 cell types, but this signal inverted across cell-count bins (geodesic: +1.06 at low n to -1.24 at high n), revealing estimation-quality confounding. The confound paper claims cosine-based metrics are "immune by construction" but never tested whether direction instability finds or fails to find a same-gene transport signal on the same data.

### Methodological note

The geodesic and EMD results were computed on single-cell distributions (subsampled 200 cells, full PCA subspaces). To make a valid comparison, direction instability must be computed on a compatible representation, not pseudobulk means (which would collapse per-gene distributions to single vectors, changing the data object and confounding the metric comparison with a representation comparison).

**Approach:** Compute per-cell-type pseudobulk signatures (mean expression per gene per cell type), then compute direction instability as D = 1 - |cos(mean_K562, mean_RPE1)| per gene. This IS a different representation than what geodesic/EMD used, so the comparison is not apples-to-apples on the representation. Instead, the key diagnostic is whether cosine DI shows the same bin-inversion pattern as geodesic/EMD. If DI is truly magnitude-immune, it should show a FLAT profile across cell-count bins, unlike geodesic (+1.06 → -1.24) and EMD (+0.44 → -0.71).

### Hypothesis

Direction instability will show NO same-gene transport signal and NO bin inversion across cell-count bins.

### Decision criteria

1. **Aggregate effect:** |d| < 0.1 for same-gene vs random-gene DI comparison.
2. **Bin-flatness test (primary result):** Across the same cell-count bins used in Table 2 of the confound paper ([20,40), [40,60), [60,100), [100,200), [200,700)), all bin-level Cohen's d values will fall within [-0.15, +0.15]. This is the diagnostic that distinguishes immunity from confounding.
3. **Correlation with cell count (secondary diagnostic):** |Spearman rho| < 0.20 between direction instability and min(n_K562, n_RPE1). If |rho| > 0.20, investigate whether bin-inversion appears; the bin-flatness test (criterion 2) remains the primary immunity check regardless. Some residual correlation is expected even for immune metrics: noisier unit vectors from low-cell-count pseudobulk inflate DI as a finite-sample effect, analogous to the rho = -0.41 second-order correlation documented in the LINCS analysis.

**Interpretation guide:**
- |d| < 0.1 AND flat bins → Confirms immunity. The confounded metrics' signal was entirely estimation-quality artifact.
- |d| < 0.1 BUT bin inversion → Cosine is confounded in the same direction but with smaller magnitude. Partial immunity only.
- |d| > 0.2 AND flat bins → Real same-gene transport signal that confounded metrics missed. Would complicate the "cell-type dominates" narrative.
- |d| > 0.2 AND bin inversion → Cosine is not immune. Would undermine the combined paper's thesis.

### Data

Replogle et al. 2022 Perturb-seq data (K562 essential + RPE1), downloaded from scPerturb Zenodo (record 10044268). Pseudobulk signatures computed as mean expression per gene per cell type among cells with the same guide RNA, restricted to genes with >= 20 cells in both cell types. HVG intersection of K562 and RPE1 used as feature space (matching the confound paper's subspace extraction).

---

## Experiment 2: Multi-metric benchmark on LINCS LOO prediction

### Background

The drug paper showed that direction instability predicts held-out cell-line consistency with AUROC = 0.986 and Spearman rho = -0.98. However, this was the only metric tested. A combined paper must show that DI outperforms alternatives on the same task, and the comparison must be fair.

### Methodological note (circularity)

The drug paper's LOO outcome is cosine-based: held-out cosine between the omitted signature and the N-1 consensus. Direction instability is also cosine-based (mean pairwise cosine distance). This shared construction gives DI a structural advantage over non-cosine metrics. A fair benchmark must include a construction-neutral outcome.

**Approach:** Run the same LOO analysis with two outcomes:
- **Cosine outcome:** held-out cosine with N-1 consensus (as in the drug paper)
- **Jaccard outcome:** held-out top-gene Jaccard overlap with N-1 consensus (construction-neutral; already computed per drug)

Four predictor metrics, each tested on both outcomes:
1. Direction instability (D, cosine-based)
2. Magnitude CV (coefficient of variation of L2 norms)
3. Frechet variance (geodesic variance on the unit sphere)
4. Mean top-gene Jaccard overlap (gene-level consistency)

Plus magnitude-corrected versions: each metric residualized against mean L2 norm.

### Hypothesis

Direction instability will achieve the highest AUROC on both outcomes. On the construction-neutral Jaccard outcome, DI's advantage over alternatives will be smaller than on the cosine outcome (quantifying the construction advantage) but will remain positive.

### Decision criteria

1. **DI wins on cosine outcome:** DI AUROC > all alternatives by at least 0.02.
2. **DI wins on Jaccard outcome (primary):** DI AUROC > all alternatives on held-out Jaccard prediction. This is the construction-neutral test.
3. **Construction advantage quantified:** (DI AUROC_cosine - next-best AUROC_cosine) > (DI AUROC_jaccard - next-best AUROC_jaccard). The gap shrinks on the neutral outcome.
4. **Magnitude correction convergence (directional):** After residualizing against mean L2 norm, magnitude-corrected alternatives will converge toward DI's performance (the AUROC gap between DI and each alternative will shrink compared to uncorrected versions). No specific shrinkage threshold is pre-registered; the prediction is directional only, since the noise floor of AUROC estimates at n = 2,694 makes a precise percentage threshold uninterpretable.

**Interpretation guide:**
- DI wins on both outcomes → DI is genuinely superior, not just construction-matched.
- DI wins on cosine but not Jaccard → DI's advantage is partly/entirely a construction artifact. Would weaken the combined paper's thesis.
- Corrected alternatives match DI → Magnitude is the primary driver; cosine normalization is the fix. Strengthens the confound paper's thesis.

### Data

Same LOO analysis as the drug paper: 2,694 drugs with >= 10 cell lines in LINCS L1000. All four metrics are already computed per drug in the cross-cellline results file.

---

## Experiment 3: MOA class stratification in JUMP-CP morphological profiles

### Background

The drug paper showed that MOA classes targeting fundamental cellular machinery (HDAC inhibitors, topoisomerase inhibitors, CDK inhibitors) have lower direction instability than receptor-mediated classes in LINCS L1000 transcriptomic profiles. Replicating this stratification in a different measurement modality (Cell Painting morphological features) would demonstrate that the finding reflects biological mechanism conservation, not a transcriptomics-specific artifact.

### Preflight metadata check (performed before pre-registration, NOT peeking at DI values)

258 LINCS drugs with MOA annotations are present in JUMP-CP with >= 5 source contexts. These include:
- 19 drugs targeting fundamental cellular machinery (HDAC, topoisomerase, CDK, tubulin, proteasome, mTOR, DNA/RNA, protein synthesis, ATPase inhibitors)
- 107 drugs targeting receptors (dopamine, serotonin, histamine, adrenergic, acetylcholine, glutamate, opioid receptor agonists/antagonists, calcium/sodium channel blockers)
- 17 kinase inhibitors (EGFR, MEK, MAPK, SRC, PI3K)

The original HDAC gradient cannot be replicated (only 3 HDAC inhibitors in JUMP-CP, 1 pan-isoform, 0 selective). The experiment is therefore broadened to MOA-class-level stratification.

### Hypothesis

Drugs targeting fundamental cellular machinery will have lower direction instability in JUMP-CP morphological profiles than receptor-mediated drugs, replicating the LINCS L1000 pattern in a different measurement modality.

### Decision criteria

1. **Group separation:** Cohen's d > 0.3 between machinery-targeting drugs (n = 19) and receptor-mediated drugs (n = 107), with machinery having lower mean DI. Calibration: the LINCS effect for HDAC inhibitors (mean D = 0.62) vs the population mean (D = 0.92) corresponds to approximately d = 0.5-0.8 for machinery vs receptor-mediated classes. We set d > 0.3 as a conservative minimum, expecting attenuation in morphological space where the feature-to-mechanism mapping is less direct than in transcriptomics.
2. **Effect direction:** The rank ordering of MOA class means (machinery < kinase < receptor) in JUMP-CP matches the ordering observed in LINCS L1000.
3. **Statistical significance:** Mann-Whitney U p < 0.05 (two-sided) for machinery vs receptor comparison.

**Interpretation guide:**
- d > 0.3 and correct ordering → Cross-modal replication confirmed. Mechanism conservation is a property of the drug-target interaction, not the measurement technology.
- d > 0.3 but reversed ordering → Morphological and transcriptomic instability measure different things. Would limit generalizability claims.
- d < 0.3 → No cross-modal replication. The LINCS finding may be transcriptomics-specific. Would still be publishable as a null but weakens the combined paper.

### Power note

With n = 19 vs n = 107 and assuming d = 0.5 (half the LINCS effect), power is approximately 0.65 (one-sided alpha = 0.05). With d = 0.3 (the minimum criterion), power drops to approximately 0.35. The experiment is underpowered for small effects. A non-significant result should be interpreted cautiously.

### Data

Direction instability for 25,254 JUMP-CP compounds (already computed, stored in results/08_jump_cp/jump_cp_instabilities.npz). MOA labels matched via InChIKey cross-reference between LINCS Repurposing Hub and JUMP-CP compound metadata. DI values for the 258 matched drugs will be extracted AFTER this pre-registration is committed.

---

## Analysis code

All analysis scripts will be committed alongside this document before execution.

- `experiments/12_perturbseq_cosine_transport.py` — Experiment 1
- `experiments/13_lincs_metric_benchmark.py` — Experiment 2
- `experiments/14_jump_cp_moa_stratification.py` — Experiment 3

## What would falsify the combined paper's thesis

The combined paper argues that (i) cosine-based direction instability is immune to the magnitude confounds that plague variance- and subspace-based metrics, and (ii) this immunity makes DI a genuinely informative measure of mechanism conservation across biological contexts. This thesis is falsified if any of the following occur:

1. **DI fails to achieve flat bins in Experiment 1**, indicating that cosine normalization does not confer immunity to estimation-quality confounding. The bin-flatness test is the decisive check: if DI shows the same inversion pattern as geodesic (+1.06 to -1.24) or EMD (+0.44 to -0.71), the "immune by construction" claim is wrong.

2. **DI loses to alternatives on the Jaccard outcome in Experiment 2**, indicating that its AUROC = 0.986 advantage is entirely a construction artifact (cosine predictor matched to cosine outcome). If a non-cosine metric achieves higher AUROC on the construction-neutral Jaccard benchmark, DI's predictive value is overstated.

3. **Experiment 3 shows reversed ordering** (receptor-mediated drugs more stable than machinery in JUMP-CP), indicating that Cell Painting and transcriptomics measure different stability axes. A reversed ordering would mean the LINCS finding is modality-specific and cannot anchor a general claim about mechanism conservation.

Any single falsification would require substantial revision of the combined paper's narrative. All three together would pivot the paper from confirmatory to null-result: a pre-registered demonstration that DI's claimed advantages do not generalize.

## Exploratory analyses (not pre-registered)

Any additional analyses beyond those specified above will be clearly labeled as EXPLORATORY in both code and manuscript, consistent with the integrity protocol used throughout the confound paper.
