# Pre-registration: DepMap Confirmatory Replication of Essential-Gene Divergence

**Timestamp:** 2026-07-05T19:42:00Z
**Status:** FROZEN — do not edit after running DepMap analysis
**Prior commit (prereg chain):** 249abaf

## Exploratory prior (disclosed)

In Perturb-seq (Replogle et al. 2022, K562/RPE1), we observed that
essential gene knockouts (DepMap CERES < -0.5) have HIGHER Grassmannian
geodesic distance (more direction-unstable) than non-essential knockouts:

- Essential mean distance: 2.798
- Non-essential mean distance: 2.684
- Cohen's d = 0.548, 95% bootstrap CI [0.452, 0.650]
- p = 8.7e-31 (Mann-Whitney)

This is opposite to the "violence mimics transport" prediction (that
cell death produces uniform signatures, hence LOW instability for
essentials). Instead, essential knockouts appear to dysregulate
context-dependent pathways, producing DIVERGENT transcriptomic
directions across cell types.

**The present hypothesis is a confirmatory replication of this
exploratory finding in an independent dataset (DepMap dependency
profiles). The direction of prediction is informed by the Perturb-seq
result. This is not a blind prediction.**

## Hypothesis (H_DepMap)

Genes classified as broadly essential (common essentials per DepMap, or
CERES < -0.5 in >50% of cell lines) will show MORE context-dependent
variability in their dependency profiles than non-essential genes
(CERES > -0.3 in >80% of lines).

**Directional prediction:** Essential > Non-essential (one-tailed).

**Operationalization:** Each gene's dependency profile is its vector of
CERES/Chronos scores across all screened cell lines (~1000 dimensions).
Context-dependence is measured as the coefficient of variation (CV) or
standard deviation of this profile across cell lines, after removing
the first principal component (which captures general cell fitness /
growth rate, analogous to removing the "essential-gene subspace" in
Perturb-seq).

**Alternative operationalization (if primary is underpowered):**
Compute pairwise correlations between essential genes' profiles vs.
between non-essential genes' profiles. If essentials show LOWER
within-group correlation (more heterogeneous effects), that confirms
context-dependence.

## Decision criteria

- **Confirms (replicates):** Cohen's d > 0.35 in predicted direction
  (essential genes more variable), p < 0.01 (one-tailed Mann-Whitney)
- **Attenuated replication:** 0.20 < d < 0.35 in predicted direction,
  p < 0.05. Report honestly as partial replication.
- **Fails to replicate:** d < 0.20, or wrong sign, or p > 0.05.
  Report as non-replication.

## Effect-size floor justification

The floor of d > 0.35 is set below the exploratory CI lower bound
(0.452) to allow for cross-dataset attenuation (different assay,
different dimensionality, different definition of "instability") while
remaining above the trivially-small range. We chose NOT to set d > 0.20
because clearing that bar given a known prior of 0.55 would be
uninformative — it's too easy.

## Falsification

The hypothesis is falsified if:
1. Essential genes show LESS variability than non-essentials (sign
   reversal), OR
2. d < 0.20 (effect too small to be meaningful), OR
3. The first-PC removal changes the sign or magnitude by >50%
   (suggesting the effect is driven entirely by the fitness confound
   we're trying to control for)

## Data source

DepMap Public 24Q4 (or latest available release):
- CRISPR gene effect (Chronos): DepMap portal
- Common essential / non-essential gene lists: DepMap portal
- Cell line metadata for PC removal: sample_info.csv

## What this does and does not prove

**Does prove (if confirmed):** The essential-gene divergence pattern
replicates in an independent high-throughput screen with different
readout (viability vs. transcriptomics). Context-dependent
dysregulation is a general property of essential gene perturbation,
not an artifact of Perturb-seq technology or K562/RPE1 cell types.

**Does not prove:** Mechanism. We cannot distinguish between "essential
genes participate in different pathways in different contexts" vs.
"essential gene loss triggers different compensatory programs in
different backgrounds" from this data alone.

## Analysis plan

1. Download DepMap Chronos gene effect matrix
2. Classify genes as essential (common essential list OR Chronos < -0.5
   in >50% of lines) vs. non-essential (Chronos > -0.3 in >80% of lines)
3. For each gene, compute its profile vector across all cell lines
4. Remove first PC (fitness component) from all profiles
5. Compute SD or CV of each gene's corrected profile
6. Compare essential vs. non-essential distributions (Mann-Whitney +
   Cohen's d)
7. Report effect size, CI, and verdict per decision criteria above
