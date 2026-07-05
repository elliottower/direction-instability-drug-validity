# Pre-Registration: Bracket Norm Validity Conditions

**Registered:** 2026-07-05T12:00:00-04:00
**Status:** Analysis code frozen before running experiments on real data.

## Paper thesis

Bracket norm (direction instability) is not a universal mechanism detector.
It becomes scientifically meaningful only after conditioning on the domain's
mechanism view, phenotype target, and admissible intervention path. We
demonstrate five failure modes and propose validity-conditioned variants that
fix each one.

## Hypotheses

### H1: Toxic perturbations have high raw bracket but fail validity filtering

Known cytotoxic compounds (stress inducers, apoptosis triggers) will rank
in the top 20% by raw direction instability, but drop below median after
toxicity correction (removing stress/apoptosis gene axes).

**Decision criterion:** >=60% of known cytotoxic drugs move from top-20%
raw bracket to below-median toxicity-corrected bracket.

### H2: Successful therapeutics can have LOW raw bracket

Known clinically successful drugs with broad mechanisms (e.g., statins,
metformin, aspirin) will have raw direction instability below the 50th
percentile. Low bracket does not imply irrelevance — it implies the
mechanism acts smoothly without strong noncommutativity.

**Decision criterion:** >=3 of 5 well-known broad-mechanism drugs have
direction instability below the population median.

### H3: Phenotype-projected bracket separates mechanism from noise

For drugs with known target pathways, phenotype-projected bracket (instability
projected onto the target pathway direction) will correlate with on-target
activity, while raw bracket will not discriminate on-target from off-target
instability.

**Decision criterion:** Spearman |rho| > 0.3 between phenotype-projected
bracket and on-target gene enrichment, AND |rho| < 0.15 for raw bracket
vs same enrichment.

### H4: Localized bracket predicts mechanism better than global bracket

For drugs targeting specific pathways, the instability should concentrate
in pathway-relevant genes. The localization score (pathway bracket / global
bracket) should predict MOA annotation better than raw global bracket.

**Decision criterion:** AUROC for MOA prediction from localization score >
AUROC from raw bracket by at least 0.05.

### H5: Transport-stable bracket outpredicts raw bracket for cross-context claims

Bracket norm that replicates across contexts (low Frechet variance) should
predict held-out context effects better than raw bracket. This tests the
"don't make global claims from local evidence" principle.

**Decision criterion:** Spearman correlation between transport-stable bracket
and held-out-context prediction > correlation for raw bracket, in at least
4/5 cross-validation folds.

### H6: Neural bracket norm matches optogenetic ablation ONLY after size correction

In the Steinmetz neural data, raw bracket norm correlates with optogenetic
silencing importance — but this correlation should survive only after
neuron-count correction (BN/sqrt(n)). Without correction, the signal is
confounded by recording yield. This replicates Tower (2026) and demonstrates
the validity ladder in a non-pharmacological domain.

**Decision criterion:** raw BN vs silencing: rho > 0.5 (expected from prior
work). BN without size correction on matched-n subsets: rho < 0.3.
BN/sqrt(n): rho > 0.6.

## Datasets

### LINCS L1000 (already extracted)
- 8,949 drugs across 5+ cell lines (from drug-perturbation-geometry)
- 978 landmark genes
- MOA annotations from Broad Repurposing Hub
- Toxicity labels: curated from known cytotoxic compounds

### Steinmetz Neuropixels (already cached)
- 73 brain regions, 39 sessions, 10 mice
- Bracket norm computed in prior work (Tower 2026)
- Optogenetic silencing ground truth from Zatka-Haas et al. 2021

### scPerturb (if needed for stretch experiments)
- Cross-context CRISPR perturbations
- Only used if H1-H6 results warrant deeper investigation

## Frozen scorer

SHA-256 of analysis code will be computed at commit time:
- `geometry/bracket_norm.py`
- `experiments/01_toxicity_failure.py`
- `experiments/02_smooth_therapeutics.py`
- `experiments/03_phenotype_projection.py`
- `experiments/04_localization.py`
- `experiments/05_transport_stability.py`
- `experiments/06_neural_replication.py`

## What counts as success

The paper succeeds if:
1. At least 3/5 drug-domain hypotheses (H1-H5) are confirmed
2. The neural replication (H6) confirms prior work
3. At least one validity-conditioned variant outperforms raw bracket

## What counts as failure

The paper fails if:
1. Raw bracket already discriminates mechanism quality (validity filtering adds nothing)
2. All corrected variants perform equivalently to raw bracket
3. The toxicity confound is not real (cytotoxic drugs don't have elevated raw bracket)

## What we report either way

All results including nulls. Pre-registered misses are reported explicitly
as "we predicted X, observed Y." No cherry-picking.
