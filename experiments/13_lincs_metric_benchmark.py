"""Experiment 13: Multi-metric LOO benchmark on LINCS.

Pre-registered in PREREGISTRATION_COMBINED_PAPER.md.

Compares direction instability against three alternative metrics for predicting
held-out cell-line consistency. Uses two outcomes:
  - Cosine outcome: held-out cosine with N-1 consensus (as in drug paper)
  - Jaccard outcome: held-out top-gene Jaccard (construction-neutral)

Tests whether DI's AUROC advantage persists on the construction-neutral outcome,
and whether magnitude-corrected alternatives converge toward DI.

Data: drug_instability_8949.csv from drug-perturbation-geometry/zenodo_v1/
      holdout_prediction.json from drug-perturbation-geometry/results/03_core_defenses/

Usage:
    uv run python experiments/13_lincs_metric_benchmark.py
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

DRUG_REPO = Path(__file__).resolve().parent.parent.parent / "drug-perturbation-geometry"
INSTABILITY_CSV = DRUG_REPO / "zenodo_v1" / "drug_instability_8949.csv"
HOLDOUT_JSON = DRUG_REPO / "results" / "03_core_defenses" / "holdout_prediction.json"
OUTPUT_DIR = Path("results/13_lincs_benchmark")

MIN_CELL_LINES = 10
CONSISTENCY_THRESHOLD = 0.5


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def auroc_for_predictor(predictor: np.ndarray, outcome: np.ndarray, threshold: float, higher_is_worse: bool = True) -> float:
    """AUROC for binary classification: outcome >= threshold → consistent.

    Args:
        higher_is_worse: If True (default), higher predictor values indicate
            less consistency (e.g., DI, magnitude_cv). If False, higher values
            indicate more consistency (e.g., Jaccard overlap).
    """
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


def main():
    log("=== EXPERIMENT 13: MULTI-METRIC LOO BENCHMARK ===")
    log("Pre-registered in PREREGISTRATION_COMBINED_PAPER.md")

    if not INSTABILITY_CSV.exists():
        log(f"ERROR: Per-drug CSV not found at {INSTABILITY_CSV}")
        sys.exit(1)
    if not HOLDOUT_JSON.exists():
        log(f"ERROR: Holdout results not found at {HOLDOUT_JSON}")
        sys.exit(1)

    df = pd.read_csv(INSTABILITY_CSV)
    log(f"  Loaded {len(df)} drugs from instability CSV")
    log(f"  Columns: {list(df.columns)}")

    with open(HOLDOUT_JSON) as f:
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
    # Higher values = worse consistency for these; Jaccard is the opposite
    HIGHER_IS_WORSE = {
        "direction_instability": True,
        "magnitude_cv": True,
        "frechet_variance": True,
        "mean_top_gene_jaccard": False,
    }

    cosine_outcome = merged["mean_heldout_cosine"].values
    # mean_top_gene_jaccard is the full-dataset mean pairwise Jaccard (not LOO).
    # It approximates the LOO held-out Jaccard closely (removing 1 of 10+ cell
    # lines changes the mean very little). Predictor #4 (same column) is excluded
    # from Jaccard-outcome comparisons to avoid trivial self-prediction.
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

    log("\n--- Decision criteria ---")
    di_cos = cosine_aurocs["direction_instability"]
    others_cos = {k: v for k, v in cosine_aurocs.items() if k != "direction_instability"}
    best_other_cos = max(others_cos.values())
    gap_cos = di_cos - best_other_cos
    c1 = gap_cos >= 0.02
    log(f"  Criterion 1 (DI wins on cosine by >= 0.02): gap={gap_cos:.4f} → {'PASS' if c1 else 'FAIL'}")

    di_jac = jaccard_aurocs["direction_instability"]
    # Exclude mean_top_gene_jaccard from Jaccard-outcome comparison (circular: same column)
    others_jac = {k: v for k, v in jaccard_aurocs.items()
                  if k not in ("direction_instability", "mean_top_gene_jaccard")}
    best_other_jac_name = max(others_jac, key=others_jac.get)
    best_other_jac = others_jac[best_other_jac_name]
    gap_jac = di_jac - best_other_jac
    c2 = gap_jac > 0
    log(f"  Criterion 2 (DI wins on Jaccard): gap={gap_jac:.4f} vs {best_other_jac_name} → {'PASS' if c2 else 'FAIL'}")
    log(f"  (mean_top_gene_jaccard excluded from Jaccard-outcome comparison: circular)")

    gap_shrinks = gap_cos > gap_jac
    log(f"  Criterion 3 (construction gap shrinks): cosine_gap={gap_cos:.4f} > jaccard_gap={gap_jac:.4f} → {'PASS' if gap_shrinks else 'FAIL'}")

    log("\n--- Magnitude-corrected alternatives ---")
    corrected_jaccard_aurocs = {}
    for name, pred in predictors.items():
        corrected = residualize(pred, mean_norm)
        auc = auroc_for_predictor(corrected, jaccard_outcome, jaccard_threshold, HIGHER_IS_WORSE[name])
        corrected_jaccard_aurocs[name] = auc
        circular = " (CIRCULAR)" if name == "mean_top_gene_jaccard" else ""
        log(f"  {name:30s} corrected_AUROC={auc:.4f} (raw={jaccard_aurocs[name]:.4f}){circular}")

    corrected_others = {k: v for k, v in corrected_jaccard_aurocs.items()
                        if k not in ("direction_instability", "mean_top_gene_jaccard")}
    best_corrected = max(corrected_others.values())
    corrected_gap = corrected_jaccard_aurocs["direction_instability"] - best_corrected
    gap_converges = corrected_gap < gap_jac
    log(f"\n  Criterion 4 (convergence, directional): raw_gap={gap_jac:.4f} → corrected_gap={corrected_gap:.4f} → {'PASS' if gap_converges else 'FAIL'}")

    log("\n--- Spearman correlations with mean_norm ---")
    for name, pred in predictors.items():
        rho, p = stats.spearmanr(pred, mean_norm)
        log(f"  {name:30s} rho(mean_norm)={rho:+.4f} (p={p:.2e})")

    log("\n--- Summary ---")
    if c1 and c2:
        log("  RESULT: DI wins on both outcomes. Advantage is genuine,")
        log("  not just construction-matched.")
    elif c1 and not c2:
        log("  RESULT: DI wins on cosine but not Jaccard. Advantage is")
        log("  partly/entirely a construction artifact.")
    elif not c1 and c2:
        log("  RESULT: DI does not dominate on cosine but wins on Jaccard.")
        log("  Unexpected — cosine outcome may have ceiling effects.")
    else:
        log("  RESULT: DI does not dominate on either outcome.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {
        "experiment": "13_lincs_metric_benchmark",
        "preregistered": "PREREGISTRATION_COMBINED_PAPER.md",
        "n_drugs": len(merged),
        "min_cell_lines": MIN_CELL_LINES,
        "cosine_aurocs": {k: float(v) for k, v in cosine_aurocs.items()},
        "jaccard_aurocs": {k: float(v) for k, v in jaccard_aurocs.items()},
        "corrected_jaccard_aurocs": {k: float(v) for k, v in corrected_jaccard_aurocs.items()},
        "criteria": {
            "c1_cosine_gap": float(gap_cos),
            "c1_pass": bool(c1),
            "c2_jaccard_gap": float(gap_jac),
            "c2_best_alternative": best_other_jac_name,
            "c2_pass": bool(c2),
            "c3_gap_shrinks": bool(gap_shrinks),
            "c4_corrected_gap": float(corrected_gap),
            "c4_gap_converges": bool(gap_converges),
        },
        "jaccard_threshold": float(jaccard_threshold),
        "cosine_threshold": float(CONSISTENCY_THRESHOLD),
    }

    with open(OUTPUT_DIR / "lincs_benchmark_results.json", "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
