# Geometric Perturbation Metrics Measure Magnitude, Not Mechanism

## Thesis (one sentence)

Variance-based and subspace-based perturbation-divergence metrics are
dominated by effect magnitude; we demonstrate this across three data
modalities with two mechanistically distinct artifacts, provide the
causal DAG explaining why, and give the normalization that fixes it.

## Timeline structure (honest, as-it-happened)

### Part 1: Pre-registered, blind

"We hypothesized that biological confounds (toxicity genes, essential-gene
subspaces, cell-health features) would substantially reorganize
direction-instability rankings (rho < 0.70). We pre-registered decision
criteria and ran blind."

**Result:** Corrections don't bite. Rankings are robust across 3 domains
(rho = 0.83–0.99). Pre-registered prediction falsified (H1–H6, HP1–HP3,
H_tissue all null or marginal).

### Part 2: Pre-registered divergence prediction

"We observed in Perturb-seq (exploratory) that essential gene knockdowns
show higher geodesic distance (d=0.55). We pre-registered a directional
replication in DepMap (d > 0.35, SHA 98897ed) before running."

**Result:** DepMap SD-based metric 'confirms' at d=2.65. Pre-registered
criterion technically satisfied.

### Part 3: Exploratory — artifact discovered (NOT pre-registered)

"The d=2.65 was implausibly large. Post-hoc normalization checks revealed
two distinct magnitude artifacts:"

- DepMap: CV-normalized d = -1.33 (sign flip). Essential genes have bigger
  effects → mechanically more variance. Tautology.
- Perturb-seq: magnitude-matched d = -0.06 (p=0.62). After matching on
  effect size, divergence vanishes completely.

"We override our own pre-registered 'confirmation' on honesty grounds:
the primary metric was fooled by a confound our decision criteria did not
anticipate."

### Part 4: The DAG and the fix (exploratory, actionable)

Present causal DAG per domain. Show normalization removes confounded
association. Deliver rule for the field.

---

## Causal DAG (per domain)

### DepMap

```
Effect magnitude (CERES score)
    ├──→ Essentiality classification (confounder: essentiality IS magnitude)
    └──→ Score variance across cell lines (more extreme scores → more variance)
              └──→ Apparent "divergence" (SD metric)
```

**Why confounder, not mediator:** Essentiality is *operationally defined*
by dependency-score magnitude. The association between essentiality and
variance is near-tautological. Adjusting is correct.

**Test:** CV normalization (SD / |mean|): d goes 2.65 → -1.33. Fixed.

### Perturb-seq

```
Effect magnitude (transcriptomic response size)
    ├──→ Essentiality status (essentials produce larger expression changes)
    └──→ PCA subspace estimation quality (larger signals → cleaner subspaces)
              └──→ Geodesic distance (better-estimated subspaces shift distances)
```

**Why confounder, not mediator:** If magnitude were a *mediator* (essential →
big response → genuinely different directions), then matching on magnitude
would *reduce* the effect but not eliminate it. Instead, matching on
magnitude completely eliminates it (d=-0.06). The signal was entirely
artifact of estimation quality, not biology.

**Test:** Magnitude-matched comparison: d goes 0.60 → -0.06. Fixed.

### LINCS (direction instability)

```
Effect magnitude
    └──→ Direction instability? NO — metric is cosine-based (normalized)
```

**Why no confound:** Direction instability uses unit-normalized signatures
(cosine of pairwise angles). Magnitude is divided out by construction.
This explains why the LINCS rankings are robust (rho > 0.99) — the metric
isn't susceptible to the magnitude confound affecting the other two.

**Lesson:** Cosine-based metrics are inherently safe. Variance-based and
subspace-based metrics require explicit magnitude normalization.

---

## Key results table

| Domain | Metric | Raw d | After normalization | Mechanism |
|--------|--------|-------|--------------------| ----------|
| DepMap | SD of dependency profile | 2.65 | -1.33 (CV) | Magnitude → variance |
| Perturb-seq | Grassmannian distance | 0.60 | -0.06 (matched) | Magnitude → subspace quality |
| LINCS | Direction instability | N/A | N/A (inherently normalized) | Cosine = safe |

---

## Actionable rule (the deliverable)

> If you compute divergence/instability/variability metrics on perturbation
> data and then compare groups defined by effect strength (essential vs
> non-essential, toxic vs non-toxic, high-dose vs low-dose), you MUST
> normalize for effect magnitude. Otherwise you will mistake "bigger effects
> are noisier" for "this group is biologically more divergent."
>
> Three valid normalizations:
> 1. Coefficient of variation (SD / |mean effect|)
> 2. Magnitude-matching (compare groups with equal response strength)
> 3. Residualization (regress out magnitude, test residuals)
>
> Cosine-based metrics (direction instability) are inherently safe because
> they normalize by magnitude before comparison.

---

## Venue / format

- arXiv + Zenodo DOI (primary)
- Possible: NeurIPS ML4H workshop, ICML CompBio workshop (4-page format)
- Possible: PLOS ONE / F1000Research (publish-if-sound, fits honest negatives)

## What makes this citable

Researchers using DepMap, Perturb-seq, JUMP-CP, or any perturbation atlas
who compute variability metrics will need to cite this or risk the same
artifact. The "cite defensively" paper.

---

## Pre-registration integrity chain

| Commit | Content | Status |
|--------|---------|--------|
| 249abaf | Hypotheses H1-H6, HP1-HP3, H_tissue + analysis code | Blind (author had not read results) |
| 98897ed | DepMap replication H_DepMap (d > 0.35) | Blind (DepMap not yet run) |
| 92e3276 | DepMap analysis script (frozen) | Before execution |
| 4a7f3fe | All results + exploratory magnitude analysis | Post-hoc, labeled exploratory |
