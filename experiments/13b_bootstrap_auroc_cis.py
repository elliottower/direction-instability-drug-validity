"""Experiment 13b: Bootstrap confidence intervals for all AUROCs in Experiment 13.

Adds 95% bootstrap CIs to every AUROC reported in the multi-metric benchmark.
Also computes CIs for AUROC *differences* (DI minus alternatives) to test
whether gaps are statistically meaningful.

Usage:
    uv run python experiments/13b_bootstrap_auroc_cis.py
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
N_BOOTSTRAP = 10_000


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def auroc_for_predictor(predictor: np.ndarray, outcome: np.ndarray, threshold: float, higher_is_worse: bool = True) -> float:
    labels = (outcome >= threshold).astype(int)
    if labels.sum() == 0 or labels.sum() == len(labels):
        return np.nan
    scores = -predictor if higher_is_worse else predictor
    return roc_auc_score(labels, scores)


def residualize(x: np.ndarray, covariate: np.ndarray) -> np.ndarray:
    cov = np.column_stack([covariate, np.ones(len(covariate))])
    beta, _, _, _ = np.linalg.lstsq(cov, x, rcond=None)
    return x - cov @ beta


def bootstrap_auroc(predictor: np.ndarray, outcome: np.ndarray, threshold: float,
                    higher_is_worse: bool, n_bootstrap: int = N_BOOTSTRAP) -> dict:
    """Bootstrap CI for a single AUROC."""
    labels = (outcome >= threshold).astype(int)
    scores = -predictor if higher_is_worse else predictor
    point = roc_auc_score(labels, scores)

    n = len(labels)
    boot_aurocs = np.zeros(n_bootstrap)
    rng = np.random.default_rng(42)
    valid = 0
    attempts = 0
    while valid < n_bootstrap and attempts < n_bootstrap * 3:
        idx = rng.integers(0, n, size=n)
        b_labels = labels[idx]
        b_scores = scores[idx]
        if b_labels.sum() == 0 or b_labels.sum() == len(b_labels):
            attempts += 1
            continue
        boot_aurocs[valid] = roc_auc_score(b_labels, b_scores)
        valid += 1
        attempts += 1

    boot_aurocs = boot_aurocs[:valid]
    ci_lo = np.percentile(boot_aurocs, 2.5)
    ci_hi = np.percentile(boot_aurocs, 97.5)
    return {"point": float(point), "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
            "n_boot": valid, "se": float(np.std(boot_aurocs))}


def bootstrap_auroc_difference(pred_a: np.ndarray, pred_b: np.ndarray,
                                outcome: np.ndarray, threshold: float,
                                hiw_a: bool, hiw_b: bool,
                                n_bootstrap: int = N_BOOTSTRAP) -> dict:
    """Bootstrap CI for AUROC(A) - AUROC(B), paired on same resamples."""
    labels = (outcome >= threshold).astype(int)
    scores_a = -pred_a if hiw_a else pred_a
    scores_b = -pred_b if hiw_b else pred_b

    point_a = roc_auc_score(labels, scores_a)
    point_b = roc_auc_score(labels, scores_b)
    point_diff = point_a - point_b

    n = len(labels)
    boot_diffs = np.zeros(n_bootstrap)
    rng = np.random.default_rng(42)
    valid = 0
    attempts = 0
    while valid < n_bootstrap and attempts < n_bootstrap * 3:
        idx = rng.integers(0, n, size=n)
        b_labels = labels[idx]
        if b_labels.sum() == 0 or b_labels.sum() == len(b_labels):
            attempts += 1
            continue
        auc_a = roc_auc_score(b_labels, scores_a[idx])
        auc_b = roc_auc_score(b_labels, scores_b[idx])
        boot_diffs[valid] = auc_a - auc_b
        valid += 1
        attempts += 1

    boot_diffs = boot_diffs[:valid]
    ci_lo = np.percentile(boot_diffs, 2.5)
    ci_hi = np.percentile(boot_diffs, 97.5)
    p_value = float(np.mean(boot_diffs <= 0)) if point_diff > 0 else float(np.mean(boot_diffs >= 0))

    return {"point_diff": float(point_diff), "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
            "p_value": p_value, "n_boot": valid,
            "auroc_a": float(point_a), "auroc_b": float(point_b)}


def main():
    log("=== EXPERIMENT 13b: BOOTSTRAP CIs FOR ALL AUROCs ===")

    df = pd.read_csv(INSTABILITY_CSV)
    with open(HOLDOUT_JSON) as f:
        holdout = json.load(f)
    holdout_df = pd.DataFrame(holdout).rename(columns={"drug": "drug_name"})
    df_10 = df[df["n_cell_lines"] >= MIN_CELL_LINES].copy()
    merged = df_10.merge(holdout_df[["drug_name", "mean_heldout_cosine"]], on="drug_name", how="inner")
    log(f"  {len(merged)} drugs")

    predictors = {
        "direction_instability": merged["direction_instability"].values,
        "magnitude_cv": merged["magnitude_cv"].values,
        "frechet_variance": merged["frechet_variance"].values,
        "mean_top_gene_jaccard": merged["mean_top_gene_jaccard"].values,
    }
    HIGHER_IS_WORSE = {
        "direction_instability": True, "magnitude_cv": True,
        "frechet_variance": True, "mean_top_gene_jaccard": False,
    }

    cosine_outcome = merged["mean_heldout_cosine"].values
    jaccard_outcome = merged["mean_top_gene_jaccard"].values
    jaccard_threshold = float(np.median(jaccard_outcome))
    mean_norm = merged["mean_norm"].values

    results = {"n_drugs": len(merged), "n_bootstrap": N_BOOTSTRAP, "cosine_threshold": CONSISTENCY_THRESHOLD,
               "jaccard_threshold": jaccard_threshold}

    # --- Individual AUROCs with CIs ---
    log("\n--- Cosine outcome AUROCs with 95% CI ---")
    results["cosine_aurocs"] = {}
    for name, pred in predictors.items():
        r = bootstrap_auroc(pred, cosine_outcome, CONSISTENCY_THRESHOLD, HIGHER_IS_WORSE[name])
        results["cosine_aurocs"][name] = r
        log(f"  {name:30s} {r['point']:.4f} [{r['ci_lo']:.4f}, {r['ci_hi']:.4f}]")

    log("\n--- Jaccard outcome AUROCs with 95% CI ---")
    results["jaccard_aurocs"] = {}
    for name, pred in predictors.items():
        if name == "mean_top_gene_jaccard":
            log(f"  {name:30s} SKIPPED (circular)")
            continue
        r = bootstrap_auroc(pred, jaccard_outcome, jaccard_threshold, HIGHER_IS_WORSE[name])
        results["jaccard_aurocs"][name] = r
        log(f"  {name:30s} {r['point']:.4f} [{r['ci_lo']:.4f}, {r['ci_hi']:.4f}]")

    log("\n--- Corrected Jaccard outcome AUROCs with 95% CI ---")
    results["corrected_jaccard_aurocs"] = {}
    for name, pred in predictors.items():
        if name == "mean_top_gene_jaccard":
            continue
        corrected = residualize(pred, mean_norm)
        r = bootstrap_auroc(corrected, jaccard_outcome, jaccard_threshold, HIGHER_IS_WORSE[name])
        results["corrected_jaccard_aurocs"][name] = r
        log(f"  {name:30s} {r['point']:.4f} [{r['ci_lo']:.4f}, {r['ci_hi']:.4f}]")

    # --- Paired AUROC differences (DI vs each alternative) ---
    log("\n--- DI vs alternatives: paired AUROC differences ---")
    di_pred = predictors["direction_instability"]
    di_hiw = HIGHER_IS_WORSE["direction_instability"]

    results["paired_differences"] = {}
    for outcome_name, outcome, threshold in [
        ("cosine", cosine_outcome, CONSISTENCY_THRESHOLD),
        ("jaccard", jaccard_outcome, jaccard_threshold),
    ]:
        results["paired_differences"][outcome_name] = {}
        for name, pred in predictors.items():
            if name == "direction_instability":
                continue
            if name == "mean_top_gene_jaccard" and outcome_name == "jaccard":
                continue
            r = bootstrap_auroc_difference(di_pred, pred, outcome, threshold,
                                            di_hiw, HIGHER_IS_WORSE[name])
            results["paired_differences"][outcome_name][name] = r
            sig = "*" if (r["ci_lo"] > 0 or r["ci_hi"] < 0) else "ns"
            log(f"  {outcome_name:8s} DI - {name:25s} = {r['point_diff']:+.4f} "
                f"[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] p={r['p_value']:.4f} {sig}")

    # --- Also bootstrap the original 0.986 AUROC from the drug paper ---
    log("\n--- Original drug paper AUROC (frac_consistent > 0.5 outcome) ---")
    holdout_full = pd.DataFrame(holdout).rename(columns={"drug": "drug_name"})
    merged_full = df_10.merge(holdout_full[["drug_name", "mean_heldout_cosine", "frac_consistent"]],
                               on="drug_name", how="inner")
    frac_outcome = merged_full["frac_consistent"].values
    di_full = merged_full["direction_instability"].values
    r986 = bootstrap_auroc(di_full, frac_outcome, 0.5, higher_is_worse=True)
    results["original_paper_auroc"] = r986
    log(f"  DI → frac_consistent>0.5: {r986['point']:.4f} [{r986['ci_lo']:.4f}, {r986['ci_hi']:.4f}]")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outpath = OUTPUT_DIR / "bootstrap_auroc_cis.json"
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {outpath}")


if __name__ == "__main__":
    main()
