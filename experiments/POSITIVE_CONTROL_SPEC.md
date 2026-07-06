# Positive Control Validation Spec (Experiment 15)

**Status: FROZEN — commit before execution. Do not modify after commit.**

**Purpose:** Establish that direction instability (DI) detects genuine
context-dependence when it exists. The confound analysis (Part I of the paper)
shows DI does not report spurious variation (specificity). This experiment
closes the complementary gap: a variability metric is only interpretable if it
can both reject spurious variation and detect genuine variation. Absent a
positive control, any null DI value is ambiguous — it could mean "stable
direction" or "instrument too blunt to see the instability." This is the
standard sensitivity/limit-of-detection requirement from assay validation in
pharmacology, independent of any external framework.

## Design: two-pronged

1. **Biological positive control** (Section A): Known context-dependent drugs
   vs known context-stable drugs, defined by external pharmacology.
2. **Synthetic spike-in** (Section B): Constructed pseudo-perturbations with
   known ground-truth instability, providing a quantitative limit-of-detection
   curve.

The synthetic spike-in is the primary evidence (unfalsifiable because ground
truth is controlled). The biological contrast is corroboration (validates that
the instrument works on real data, not just synthetic constructions).

---

## Section A: Biological Positive Control

### Positive control set: nuclear hormone receptor agonists (known-unstable)

Nuclear hormone receptors (NHR) are intracellular ligand-activated transcription
factors whose transcriptional output is famously cell-type-specific. The same
receptor (e.g., glucocorticoid receptor) drives opposite transcriptional
programs in different tissues — anti-inflammatory in immune cells,
gluconeogenic in liver, catabolic in muscle. This context-dependence is
textbook endocrinology established decades before any transcriptomic assay.
NHR agonists are therefore an externally-defined positive control: drugs whose
perturbation signatures should vary across cell lines.

**Definition:** All drugs matching the `is_nuclear_receptor()` function from
Experiment 13f v2 (the identical keyword list: glucocorticoid receptor,
estrogen receptor, androgen receptor, progesterone receptor, PPAR receptor,
retinoid receptor, vitamin D receptor, mineralocorticoid receptor).

**Count:** 125 drugs. All have n_cell_lines >= 5; 103 have n_cell_lines >= 10.

**Note on relationship to Experiment 13f:** NHR drugs were the
*excluded/sensitivity-analysis* class in the MOA classification test. Here they
serve as a *positive control* — a different experimental role. The
classification rule and the sensitivity validation ask different questions
(MOA classification asks "does DI separate machinery from receptor?"
while this asks "does DI register elevated values on known context-dependent
drugs?"). The spec makes this role distinction explicit, and the synthetic
spike-in (Section B) provides evidence with no biological-reuse question.

### Negative control set: constitutive machinery drugs (known-stable)

**Definition:** All drugs whose MOA annotation (Drug Repurposing Hub) matches
any of: "proteasome inhibitor", "tubulin inhibitor", "tubulin polymerization
inhibitor", "HDAC inhibitor", "topoisomerase inhibitor". Excludes any drug also
matching `is_nuclear_receptor()` (none expected, but the check is mechanical).

These four target classes hit universal, constitutively expressed cellular
machinery:
- **Proteasome inhibitors** target the 26S proteasome (ubiquitous protein
  degradation).
- **Tubulin inhibitors** target alpha/beta-tubulin dimers (mitotic spindle,
  universal in dividing cells).
- **HDAC inhibitors** target histone deacetylases (chromatin remodeling,
  constitutively active across cell types).
- **Topoisomerase inhibitors** target DNA topoisomerases (required for
  replication and transcription in all dividing cells).

All four target constitutive machinery present in essentially all dividing
human cells. Inhibition triggers conserved downstream cascades (unfolded
protein response, mitotic arrest, DNA damage response, transcriptional
derepression) regardless of cell identity. Using four mechanistically distinct
target classes rather than one or two guards against the possibility that a
single mechanism's idiosyncrasy (e.g., proteostasis stress specifically)
explains the stability.

**Count:** 59 drugs. All have n_cell_lines >= 5; 51 have n_cell_lines >= 10.

### Hypothesis (stated before running)

DI is higher in NHR drugs than in constitutive machinery drugs (NHR drugs are
more context-dependent). The a priori direction is: DI(NHR) > DI(machinery).

### Decision criteria

**Primary test (pre-registered pass/fail gate):** Bootstrap CI on the mean
difference mean(DI_NHR) - mean(DI_machinery), 10,000 resamples, percentile
method. The 95% CI must exclude zero. Alpha = 0.05 for this single primary
test. This is the sole gate for claiming "DI detects known context-dependence."

**Secondary (descriptive, no pass/fail gate):**
- Mann-Whitney U test (two-sided) on DI values. Report the U statistic,
  p-value, and rank-biserial correlation as an effect size.
- AUROC separating NHR from machinery using DI as the score. Report with 95%
  bootstrap CI (same harness as 13f). The positive class is NHR (expected
  higher DI).

These secondary metrics characterize the separation but do not determine
whether the experiment passes or fails.

### Robustness checks

1. **Cell-line-matched:** Restrict to drugs with n_cell_lines >= 10 (103 NHR,
   51 machinery) to ensure DI estimates are not inflated by small-n noise.
2. **Narrower negative control:** Restrict the negative control to proteasome
   + tubulin only (n=22, the most mechanistically homogeneous stable set) to
   confirm the result is not driven by heterogeneity in the broader set.

### Honest failure case

If the bootstrap CI for mean(DI_NHR) - mean(DI_machinery) includes zero, that
is a genuine sensitivity failure: DI cannot distinguish known context-dependent
from known context-stable drugs at this sample size. Report the result either
way. If the AUROC CI includes 0.5, DI has no discriminative power on this
contrast.

---

## Section B: Synthetic Spike-in (Limit-of-Detection Curve)

### Construction

The spike-in generates pseudo-perturbations with controlled ground-truth
direction instability. Each pseudo-perturbation consists of K synthetic
"cell-line signatures" in R^978 (the 978 LINCS L1000 landmark genes).

**Procedure for one pseudo-perturbation at ground-truth spread sigma:**

1. Sample a base direction vector v uniformly from the unit sphere in R^978
   (draw from N(0, I) and normalize).
2. For each of K synthetic "cell lines":
   a. Draw a perturbation direction n_i from N(0, I), orthogonalize against v
      via Gram-Schmidt, normalize to get a unit vector orthogonal to v.
   b. Draw an angle theta_i from |N(0, sigma)| (half-normal, so theta >= 0).
   c. Set the noiseless signature: s_i^0 = cos(theta_i) * v + sin(theta_i) * n_i.
   d. Add measurement noise: s_i = s_i^0 + epsilon_i, where
      epsilon_i ~ N(0, tau^2 * I_{978}).
   e. (Do NOT re-normalize s_i — real signatures are not unit-normalized.)
3. Compute DI on {s_1, ..., s_K} using the standard formula:
   DI = 1 - mean_{i<j} cos(s_i, s_j).

### Noise calibration (tau)

A noiseless synthetic generator (tau=0) produces DI = 0 exactly at sigma = 0,
making the null threshold identically zero and the LOD a grid artifact. Real
DI values have a nonzero floor from finite sampling, estimation error, and
biological noise. To make the synthetic LOD interpretable as a proxy for the
real instrument:

**Calibration procedure:** Before running the LOD sweep, set sigma = 0 and
sweep tau over a grid (tau in {0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5}).
For each tau, generate 500 pseudo-perturbations at K = 13 (the annotated-drug
median) and compute the DI distribution. Select the tau whose sigma=0 DI
distribution (mean and standard deviation) best matches the empirical DI
distribution of the negative-control machinery drugs (the 59 constitutive
machinery drugs from Section A). This anchors the synthetic null to the real
instrument's noise floor.

The selected tau is reported alongside the LOD curve. If no tau produces a
reasonable match (defined as: synthetic sigma=0 mean DI within 0.05 of
machinery mean DI, AND synthetic sigma=0 SD within a factor of 2 of machinery
SD), report the LOD curve with the closest tau and note the calibration
mismatch as a limitation.

### K sweep (second axis)

DI's sampling variance depends on K — fewer cell lines means noisier, possibly
upward-biased DI estimates. Characterizing the LOD only at one K would hide
this interaction. The spike-in sweeps K as a second axis:

**K values:** {5, 10, 13, 20, 40}.
- K=5: lower bound of cell-line coverage in the dataset (minimum for DI
  computation in the drug paper).
- K=10: common threshold for "well-measured" drugs.
- K=13: median n_cell_lines among annotated LINCS drugs.
- K=20: 75th percentile of annotated drugs.
- K=40: well-measured drugs with high coverage.

For each (sigma, K) pair, generate N_rep = 500 pseudo-perturbations. The LOD
is reported as a function of both sigma and K: LOD(K) = the smallest sigma at
which >= 95% of pseudo-perturbations exceed the null threshold at that K.

### Sigma sweep

**Sweep:** sigma in {0.0, 0.05, 0.10, 0.15, ..., 1.50} (31 values).

Note on axis labeling: sigma is the scale parameter of the half-normal
distribution from which theta_i is drawn. It is NOT the mean angle — the
realized mean angle is sigma * sqrt(2/pi) ≈ 0.798 * sigma. The x-axis of
the LOD curve is labeled "half-normal scale sigma (radians)" and the realized
mean and 90th-percentile angles are reported separately at each sigma value.

### Output: limit-of-detection curve

For each (sigma, K) pair, report:
- Mean DI across 500 pseudo-perturbations
- 95% percentile interval of DI values (2.5th to 97.5th percentile)
- Realized mean angle (degrees) and 90th-percentile angle (degrees) of the
  theta_i draws
- The fraction of pseudo-perturbations where DI > DI_null_threshold(K)

The DI_null_threshold(K) is the 95th percentile of DI at sigma = 0 for that
specific K value (estimated from the 500 sigma=0 pseudo-perturbations at that
K, with calibrated tau). This represents the false-positive threshold: the DI
value exceeded by only 5% of noise-floor-matched pseudo-perturbations.

**Limit of detection** at each K = the smallest sigma at which >= 95% of
pseudo-perturbations exceed DI_null_threshold(K). This is the instrument's
sensitivity floor in half-normal scale units. Report the corresponding mean
angle in degrees for interpretability.

### Calibration against real data

After computing the LOD surface, overlay the empirical DI distributions of
real NHR drugs and real machinery drugs on the K=13 slice of the LOD curve
(matching the annotated-drug median K). This answers: "where do real drugs fall
on the synthetic dose-response curve?" If NHR drugs cluster at high
sigma-equivalent and machinery drugs cluster at low sigma-equivalent, the
synthetic and biological evidence are concordant.

### Parameters (frozen)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| K (cell lines per pseudo-drug) | {5, 10, 13, 20, 40} | Spans the real drug K distribution |
| N_rep (pseudo-drugs per sigma per K) | 500 | Sufficient for tight percentile intervals |
| Sigma range | [0.0, 1.50] step 0.05 | 31 values; mean angle range 0° to ~69° |
| Dimensionality | 978 | Matches LINCS L1000 landmark gene count |
| Noise calibration (tau) | Matched to machinery DI distribution | Anchors synthetic null to real instrument |
| Detection threshold | 95th pctile at sigma=0, per K | Standard 5% false-positive rate |
| Detection criterion | >= 95% power | Standard assay LOD definition |

---

## Pre-registration commitments

1. **Class definitions are external.** NHR = textbook endocrinology.
   Proteasome/tubulin/HDAC/topoisomerase = constitutive machinery. These are
   not defined by looking at DI values.
2. **Hypothesis direction is a priori.** NHR should have HIGHER DI (more
   context-dependent). This follows from the pharmacology, not from data.
3. **Report either way.** If DI does not separate NHR from machinery, that
   is a real sensitivity limitation and weakens all null DI interpretations
   in the paper.
4. **No data-dependent modifications.** The keyword lists, decision criteria,
   and spike-in parameters above are frozen. If the result is null, we do not
   add drugs, change thresholds, or redefine classes.
5. **Drug counts verified for feasibility only.** We confirmed 125 NHR and
   59 constitutive machinery drugs exist in the dataset with adequate
   cell-line coverage. No DI values were examined before freezing this spec.
6. **Single primary test.** The pass/fail gate is the bootstrap CI on
   mean(DI_NHR) - mean(DI_machinery) at alpha = 0.05. Mann-Whitney and AUROC
   are secondary/descriptive.

## Implementation notes

- Reuse `is_nuclear_receptor()` from 13f_v2 (identical keyword list).
- Reuse bootstrap harness from 13f_v2 for AUROC CIs.
- Output to `results/15_positive_control/`.
- Script: `experiments/15_positive_control_validation.py`.
- Save all results to JSON (never rely on stdout).
- Tau calibration runs first and its result (selected tau, quality of match)
  is saved before the main LOD sweep.
