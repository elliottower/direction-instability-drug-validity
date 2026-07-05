# Pre-Registration Addendum: Perturb-seq Direction Instability

**Date frozen**: 2026-07-05
**Status**: Correction test frozen before running. Descriptive distances explored (see integrity boundary below).
**Depends on**: genetic-perturbation-holonomy Phase 1 results (d = -0.14, ambiguous)

## Motivation

The drug paper (Tower, 2026) showed that toxicity correction has negligible
effect on direction instability rankings in LINCS L1000: cytotoxic drugs are
already high-instability, not false-positive transporters, so removing stress
gene axes moves instability values by <0.01 (raw vs corrected rho = 0.83).

The validity ladder thesis predicts that different domains face different
failure modes, and therefore different validity conditions bite. Perturb-seq
genetic knockdowns provide a test: the analogous "violence" category is
universal-essential gene knockdowns (DepMap Chronos < -0.5 in both K562 and
RPE1), which produce lethal disruption in all cell types.

## Integrity boundary: what has been examined

The following are known before this pre-registration and cannot be un-seen:

- **Phase 1 distances** (from genetic-perturbation-holonomy): 1,676 genes with
  Grassmannian geodesic distances between K562 and RPE1 perturbation response
  subspaces. Same-gene distance effect: Cohen's d = -0.14 (ambiguous).
- **Phase 1 differential essentiality null**: Spearman r = 0.005 between
  |DepMap_K562 - DepMap_RPE1| and geodesic distance. Flat null.
- **Universal-essential direction** (explored 2026-07-05 before freezing this
  document): 855 genes with Chronos < -0.5 in both K562 and RPE1 have HIGHER
  geodesic distances (mean 2.80 vs 2.68, Cohen's d = +0.55). Universal-essential
  knockdowns DIVERGE across cell types --- the predicted "violence mimics
  transport" failure mode is INVERTED in genetic perturbations.

The following have NOT been examined:

- Subspace-level correction (projecting out essential-response subspace)
- Rank disruption between raw and corrected distances
- Any analysis of corrected distances vs biological function

## Hypotheses

### HP1 (primary): Essential-response subspace correction changes rankings

Compute the Frechet mean subspace of universal-essential gene knockdowns
(DepMap Chronos < -0.5 in both K562 and RPE1) in each cell type. Project this
out of each gene's perturbation response subspace. Recompute Grassmannian
geodesic distances on the residual subspaces.

**Pre-registered prediction**: The correction changes the distance ranking,
measured as Spearman rho between raw and corrected geodesic distances.

**Decision criterion**:
- "Bites": rho < 0.70 (substantially more rank disruption than the drug paper's 0.83)
- "Marginal": 0.70 <= rho < 0.83
- "Doesn't bite": rho >= 0.83 (comparable to the drug case)

All three outcomes are informative:
- "Bites": different validity condition matters in genetics vs pharmacology (thesis confirmed directly)
- "Marginal": intermediate evidence
- "Doesn't bite": the validity ladder's failure-mode catalog is wrong for this domain --- universal-essential knockdowns don't create a shared response axis that needs correction, because their responses DIVERGE rather than converge. The inverted failure mode (violence produces divergence, not false transport) is itself a domain-specific insight.

### HP2 (secondary): Corrected distances better predict pathway function

After essential-response correction, geodesic distance should better predict
Gene Ontology pathway coherence between K562 and RPE1 knockdown effects.
GO pathway annotations from MSigDB (C5 collection) frozen before analysis.

**Decision criterion**: Spearman |rho| between corrected distance and GO
pathway Jaccard overlap > |rho| for raw distance, by at least 0.05.

### HP3 (descriptive): The failure-mode inversion is robust

The exploratory finding that universal-essential knockdowns have HIGHER
distance (d = +0.55) will be confirmed by:
- Bootstrap resampling (200 iterations) of the Cohen's d estimate
- Stratification by essentiality severity (Chronos < -0.5, -0.7, -1.0)

This hypothesis has no pass/fail threshold; it characterizes the domain-specific
failure mode for the paper's discussion.

## Connection to Phase 1 differential-essentiality null

Phase 1 tested **differential** essentiality: |DepMap_K562 - DepMap_RPE1|
predicting geodesic distance. This was a flat null (r = 0.005).

This addendum tests **universal** essentiality: genes essential in BOTH cell
types creating a shared response subspace that inflates or deflates apparent
transport. These are orthogonal quantities:
- Differential essentiality asks whether cell-type-SPECIFIC dependency
  predicts transport failure.
- Universal essentiality asks whether shared dependency creates a confounding
  axis in the response subspace.

The Phase 1 null does not bear on the hypotheses here.

## K = 2 limitation

With two cell types, direction instability collapses to cosine distance
(1 - cos) or its subspace generalization (Grassmannian geodesic). There is no
Frechet variance, no holonomy, no cross-context transport stability. This
addendum tests only the **correction-changes-ranking** claim. Transport
stability (H5 from the main pre-registration) requires >= 3 contexts and
belongs in Phase 2 of the holonomy repo, not here.

## Analysis plan

1. Load existing subspaces from genetic-perturbation-holonomy data/extracted/
2. Identify universal-essential genes (Chronos < -0.5 in both cell types)
3. Compute Frechet mean subspace of universal-essential subspaces per cell type
4. Project out essential-response subspace from all genes' subspaces
5. Recompute Grassmannian geodesic distances on residual subspaces
6. Test HP1 (rank disruption: rho between raw and corrected)
7. Test HP2 (GO pathway prediction improvement)
8. Report HP3 (failure-mode inversion characterization)
