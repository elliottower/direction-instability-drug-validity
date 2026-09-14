# direction-instability-drug-validity

Confound audit of a cross-context consistency score across drug, genetic and morphological perturbations.

## What this is

Direction instability measures how much a perturbation's response direction changes across contexts: `D = 1 - mean pairwise cosine` over a perturbation's context signatures. Scores like it are used to rank drugs and genes by how well a mechanism carries from one context to another. Whether those rankings survive correction for the obvious confounds — toxicity programs, off-target consistency, essentiality, generic morphology — had not been tested.

This repository holds the audit. Hypotheses and decision criteria were frozen and published before the outcomes were inspected, and every reported number is written to a file under `results/` by the script that computes it.

## Key results

| | Prediction | Result | Status |
|---|---|---|:--|
| H1 | Toxicity correction reclassifies cytotoxic drugs | Δ < 0.01 | Refuted — the criterion passes as an artifact |
| H2 | Broad-mechanism drugs fall below the median | 7 of 10 | Suggestive; does not survive multiplicity correction |
| H3 | Target-axis projection separates on- from off-target | ρ = 0.376 vs −0.043 raw | **Confirmed**; survives Holm–Bonferroni |
| H4 | Localization predicts mechanism of action | AUROC gap 0.018 | Informative null; domain-dependent |
| H5 | The variance-regularized score outpredicts raw | 10 of 66 folds | Not confirmed |
| HP1 | Essential-gene correction reorganizes the ranking | ρ = 0.91 | Correction does not bite |
| HP2 | Corrected distances predict pathway function | Δρ = −0.02 | Null |
| HP3 | Universal-essential knockdowns diverge | d = 0.60 | Directional support, inverting the expected confound |
| Exp. 8 | — | ρ = 0.998, n = 24,992 | Unregistered sensitivity analysis |

Corrections leave the population ranking largely intact while moving specific items: proteasome and spliceosome components shift 660–900 positions of 1,676 under a correction whose global effect is ρ = 0.91. Raw direction instability predicts held-out cell-line agreement in all 66 leave-one-out folds (mean ρ = −0.602); no proposed variant improves on it.

Data: 8,949 LINCS L1000 compounds, 1,676 Perturb-seq CRISPRi knockdowns, 24,992 JUMP-CP Cell Painting compounds.

## Repo structure

```
geometry/direction_instability.py   The score and its variants
experiments/                        One script per experiment, plus registrations
results/                            Every reported number, one directory per experiment
paper/                              LaTeX source, figures, and each version
DEVIATION_LOG.md                    Departures from the registrations, with what each cost
```

## Pre-registration

| Registration | Frozen | Covers |
|---|---|---|
| `PREREGISTRATION.md`, `PREREGISTRATION_PERTURBSEQ.md` | `249abaf` | H1–H6, HP1–HP3 |
| `PREREG_ROBUSTNESS_CHECKS.md` | `868695e` | Lambda and HDAC sensitivity, source-count stratification |
| `experiments/PREREG_H3_MAGNITUDE_AND_SHARED_AXIS.md` | `f288507` | Magnitude, coverage and shared-axis controls on H3 |

Deviations are recorded in `DEVIATION_LOG.md` against the registration they depart from, each with the artifact that reproduces it.

## Data

Not in git. LINCS L1000 signatures come from GEO GSE92742 Level 5 and are built by
`drug-perturbation-geometry/data/lincs_loader.py`; the JUMP-CP profiles stream from the Cell
Painting Gallery; Perturb-seq is Replogle et al. 2022.

## Setup

```bash
uv sync
```

## Reproducing

```bash
# Synthetic validation, no data required
uv run python experiments/01_toxicity_failure.py --synthetic

# Real analysis, requires the LINCS build
uv run python experiments/03_phenotype_projection.py --real --data ../drug-perturbation-geometry/data/

# JUMP-CP, streamed and run on one CPU worker
modal run --detach experiments/modal_08_jump_cp.py --stage all
```

## Status

Analysis complete; the manuscript is in revision. One registered hypothesis is confirmed, and three sensitivity analyses of it are registered and not yet run.

## License

MIT
