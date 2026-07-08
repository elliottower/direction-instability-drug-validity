"""Experiment 13: Multi-metric LOO benchmark on LINCS.

Pre-registered in PREREGISTRATION_COMBINED_PAPER.md.

Compares direction instability against three alternative metrics for predicting
held-out cell-line consistency. Uses two outcomes:
  - Cosine outcome: held-out cosine with N-1 consensus (as in drug paper)
  - Jaccard outcome: held-out top-gene Jaccard (construction-neutral)

Tests whether DI's AUROC advantage persists on the construction-neutral outcome,
and whether magnitude-corrected alternatives converge toward DI.

MIN_CELL_LINES = 10 here (vs 5 in experiments 03-05) because this benchmark
focuses on the most data-rich drugs where metric differences are meaningful,
while the hypothesis tests use 5 to maximize statistical power on smaller samples.

Decision criteria compare DI against each alternative individually (not
max-of-others, which would inflate the alternative's apparent performance).

Data: drug_instability_8949.csv from drug-perturbation-geometry/zenodo_v1/
      holdout_prediction.json from drug-perturbation-geometry/results/03_core_defenses/

Usage:
    uv run python experiments/13_lincs_metric_benchmark.py --data PATH_TO_DRUG_REPO
    uv run python experiments/13_lincs_metric_benchmark.py  # uses default sibling repo path
"""
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score


MIN_CELL_LINES = 10
CONSISTENCY_THRESHOLD = 0.5


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def auroc_for_predictor(predictor: np.ndarray, outcome: np.ndarray, threshold: float, higher_is_worse: bool = True) -> float:
    """AUROC for binary classification: outcome >= threshold -> consistent."""
    labels = (outcome >= threshold).astype(int)
    if labels.sum() == 0 or labels.sum() == len(labels):
        return np.nan
    scores = -predictor if higher_is_worse else predictor
    return roc_auc_score(labels, scores)


def residualize(x: np.ndarray, covariate: np.ndarray) -> np.ndarray:
    """Residualize x against covariate via OLS."""
    cov = np.column_stack([covariate, np.ones(len(covariate))])
    beta, _, _, _ = np.linalg.lstsq(cov, x, rcond=None)
    return x - cov @ beta


def bootstrap_auroc_gap(predictor_a, predictor_b, outcome, threshold,
                        higher_is_worse_a, higher_is_worse_b,
                        n_boot=2000, seed=42):
    """Bootstrap 95% CI on (AUROC_a - AUROC_b)."""
    rng = np.random.default_rng(seed)
    n = len(predictor_a)
    gaps = []
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        auc_a = auroc_for_predictor(predictor_a[idx], outcome[idx], threshold, higher_is_worse_a)
        auc_b = auroc_for_predictor(predictor_b[idx], outcome[idx], threshold, higher_is_worse_b)
        if np.isnan(auc_a) or np.isnan(auc_b):
            continue
        gaps.append(auc_a - auc_b)
    if len(gaps) < 100:
        return np.nan, np.nan
    gaps = np.array(gaps)
    return float(np.percentile(gaps, 2.5)), float(np.percentile(gaps, 97.5))


def main(drug_repo: Path):
    log("=== EXPERIMENT 13: MULTI-METRIC LOO BENCHMARK ===")
    log("Pre-registered in PREREGISTRATION_COMBINED_PAPER.md")

    instability_csv = drug_repo / "zenodo_v1" / "drug_instability_8949.csv"
    holdout_json = drug_repo / "results" / "03_core_defenses" / "holdout_prediction.json"
    output_dir = Path("results/13_lincs_benchmark")

    if not instability_csv.exists():
        log(f"ERROR: Per-drug CSV not found at {instability_csv}")
        sys.exit(1)
    if not holdout_json.exists():
        log(f"ERROR: Holdout results not found at {holdout_json}")
        sys.exit(1)

    df = pd.read_csv(instability_csv)
    log(f"  Loaded {len(df)} drugs from instability CSV")
    log(f"  Columns: {list(df.columns)}")

    with open(holdout_json) as f:
        holdout = json.load(f)
    log(f"  Loaded {len(holdout)} drugs from holdout prediction")

    holdout_df = pd.DataFrame(holdout)
    log(f"  Holdout columns: {list(holdout_df.columns)}")

    df_10 = df[df["n_cell_lines"] >= MIN_CELL_LINES].copy()
    log(f"  Drugs with >= {MIN_CELL_LINES} cell lines: {len(df_10)}")

    holdout_df = holdout_df.rename(columns={"drug": "drug_name"})
    merged = df_10.merge(holdout_df[["drug_name", "mean_heldout_cosine"]], on="drug_name", how="inner")
    log(f"  Merged (instability + holdout): {len(merged)}")

    predictors = {
        "direction_instability": merged["direction_instability"].values,
        "magnitude_cv": merged["magnitude_cv"].values,
        "frechet_variance": merged["frechet_variance"].values,
        "mean_top_gene_jaccard": merged["mean_top_gene_jaccard"].values,
    }
    HIGHER_IS_WORSE = {
        "direction_instability": True,
        "magnitude_cv": True,
        "frechet_variance": True,
        "mean_top_gene_jaccard": False,
    }

    cosine_outcome = merged["mean_heldout_cosine"].values
    jaccard_outcome = merged["mean_top_gene_jaccard"].values
    mean_norm = merged["mean_norm"].values

    log("\n--- Predictor correlations ---")
    for name, pred in predictors.items():
        rho_cos, _ = stats.spearmanr(pred, cosine_outcome)
        rho_jac, _ = stats.spearmanr(pred, jaccard_outcome)
        log(f"  {name:30s} rho_cosine={rho_cos:+.4f}  rho_jaccard={rho_jac:+.4f}")

    log("\n--- AUROC on cosine outcome ---")
    cosine_aurocs = {}
    for name, pred in predictors.items():
        auc = auroc_for_predictor(pred, cosine_outcome, CONSISTENCY_THRESHOLD, HIGHER_IS_WORSE[name])
        cosine_aurocs[name] = auc
        log(f"  {name:30s} AUROC={auc:.4f}")

    log("\n--- AUROC on Jaccard outcome (PRIMARY — construction-neutral) ---")
    jaccard_aurocs = {}
    jaccard_threshold = np.median(jaccard_outcome)
    log(f"  Jaccard threshold (median): {jaccard_threshold:.4f}")
    for name, pred in predictors.items():
        auc = auroc_for_predictor(pred, jaccard_outcome, jaccard_threshold, HIGHER_IS_WORSE[name])
        jaccard_aurocs[name] = auc
        circular = " (CIRCULAR — same column)" if name == "mean_top_gene_jaccard" else ""
        log(f"  {name:30s} AUROC={auc:.4f}{circular}")

    di_pred = predictors["direction_instability"]

    log("\n--- Decision criteria (DI vs each alternative individually) ---")
    cosine_gaps = {}
    for name, pred in predictors.items():
        if name == "direction_instability":
            continue
        gap = cosine_aurocs["direction_instability"] - cosine_aurocs[name]
        ci_lo, ci_hi = bootstrap_auroc_gap(
            di_pred, pred, cosine_outcome, CONSISTENCY_THRESHOLD,
            HIGHER_IS_WORSE["direction_instability"], HIGHER_IS_WORSE[name])
        cosine_gaps[name] = {"gap": float(gap), "ci": [ci_lo, ci_hi]}
        passes = gap >= 0.02
        log(f"  DI vs {name:25s} cosine gap={gap:+.4f} 95% CI [{ci_lo:+.4f}, {ci_hi:+.4f}] -> {'PASS' if passes else 'FAIL'} (>=0.02)")

    jaccard_gaps = {}
    for name, pred in predictors.items():
        if name in ("direction_instability", "mean_top_gene_jaccard"):
            continue
        gap = jaccard_aurocs["direction_instability"] - jaccard_aurocs[name]
        ci_lo, ci_hi = bootstrap_auroc_gap(
            di_pred, pred, jaccard_outcome, jaccard_threshold,
            HIGHER_IS_WORSE["direction_instability"], HIGHER_IS_WORSE[name])
        jaccard_gaps[name] = {"gap": float(gap), "ci": [ci_lo, ci_hi]}
        passes = gap > 0
        log(f"  DI vs {name:25s} Jaccard gap={gap:+.4f} 95% CI [{ci_lo:+.4f}, {ci_hi:+.4f}] -> {'PASS' if passes else 'FAIL'} (>0)")

    c1_all = all(g["gap"] >= 0.02 for g in cosine_gaps.values())
    c2_all = all(g["gap"] > 0 for g in jaccard_gaps.values())

    log(f"\n  Criterion 1 (DI wins on cosine vs ALL alternatives by >= 0.02): {'PASS' if c1_all else 'FAIL'}")
    log(f"  Criterion 2 (DI wins on Jaccard vs ALL non-circular alternatives): {'PASS' if c2_all else 'FAIL'}")

    log("\n--- Magnitude-corrected alternatives ---")
    corrected_jaccard_aurocs = {}
    for name, pred in predictors.items():
        corrected = residualize(pred, mean_norm)
        auc = auroc_for_predictor(corrected, jaccard_outcome, jaccard_threshold, HIGHER_IS_WORSE[name])
        corrected_jaccard_aurocs[name] = auc
        circular = " (CIRCULAR)" if name == "mean_top_gene_jaccard" else ""
        log(f"  {name:30s} corrected_AUROC={auc:.4f} (raw={jaccard_aurocs[name]:.4f}){circular}")

    log("\n--- Spearman correlations with mean_norm ---")
    for name, pred in predictors.items():
        rho, p = stats.spearmanr(pred, mean_norm)
        log(f"  {name:30s} rho(mean_norm)={rho:+.4f} (p={p:.2e})")

    log("\n--- Summary ---")
    if c1_all and c2_all:
        log("  RESULT: DI wins on both outcomes vs all alternatives.")
        log("  Advantage is genuine, not construction-matched.")
    elif c1_all and not c2_all:
        log("  RESULT: DI wins on cosine but not Jaccard for all alternatives.")
        log("  Advantage is partly/entirely a construction artifact.")
    elif not c1_all and c2_all:
        log("  RESULT: DI does not dominate on cosine but wins on Jaccard.")
    else:
        log("  RESULT: DI does not dominate on either outcome.")

    output_dir.mkdir(parents=True, exist_ok=True)
    results = {
        "experiment": "13_lincs_metric_benchmark",
        "preregistered": "PREREGISTRATION_COMBINED_PAPER.md",
        "n_drugs": len(merged),
        "min_cell_lines": MIN_CELL_LINES,
        "cosine_aurocs": {k: float(v) for k, v in cosine_aurocs.items()},
        "jaccard_aurocs": {k: float(v) for k, v in jaccard_aurocs.items()},
        "corrected_jaccard_aurocs": {k: float(v) for k, v in corrected_jaccard_aurocs.items()},
        "pairwise_cosine_gaps": cosine_gaps,
        "pairwise_jaccard_gaps": jaccard_gaps,
        "criteria": {
            "c1_all_cosine_pass": bool(c1_all),
            "c2_all_jaccard_pass": bool(c2_all),
        },
        "jaccard_threshold": float(jaccard_threshold),
        "cosine_threshold": float(CONSISTENCY_THRESHOLD),
    }

    with open(output_dir / "lincs_benchmark_results.json", "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {output_dir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).resolve().parent.parent.parent / "drug-perturbation-geometry",
                        help="Path to drug-perturbation-geometry repo (default: sibling directory)")
    args = parser.parse_args()
    main(args.data)
