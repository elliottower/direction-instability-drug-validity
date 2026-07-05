"""Experiment 13e: Bootstrap AUPRC for all predictor-outcome combinations.

Complements 13b (AUROC) with Average Precision (AUPRC), which is more
informative when class balance is skewed. Uses same data, same predictors,
same outcomes, same paired bootstrap approach.

Usage:
    uv run python experiments/13e_bootstrap_auprc.py
"""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

DRUG_REPO = Path(__file__).resolve().parent.parent.parent / "drug-perturbation-geometry"
INSTABILITY_CSV = DRUG_REPO / "zenodo_v1" / "drug_instability_8949.csv"
HOLDOUT_JSON = DRUG_REPO / "results" / "03_core_defenses" / "holdout_prediction.json"
OUTPUT_DIR = Path("results/13_lincs_benchmark")

MIN_CELL_LINES = 10
CONSISTENCY_THRESHOLD = 0.5
N_BOOTSTRAP = 10_000


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def bootstrap_auprc(predictor: np.ndarray, outcome: np.ndarray, threshold: float,
                    higher_is_worse: bool, n_bootstrap: int = N_BOOTSTRAP) -> dict:
    labels = (outcome >= threshold).astype(int)
    scores = -predictor if higher_is_worse else predictor
    point = average_precision_score(labels, scores)

    n = len(labels)
    boot_aps = np.zeros(n_bootstrap)
    rng = np.random.default_rng(42)
    valid = 0
    attempts = 0
    while valid < n_bootstrap and attempts < n_bootstrap * 3:
        idx = rng.integers(0, n, size=n)
        b_labels = labels[idx]
        if b_labels.sum() == 0 or b_labels.sum() == len(b_labels):
            attempts += 1
            continue
        boot_aps[valid] = average_precision_score(b_labels, scores[idx])
        valid += 1
        attempts += 1

    boot_aps = boot_aps[:valid]
    return {
        "point": float(point),
        "ci_lo": float(np.percentile(boot_aps, 2.5)),
        "ci_hi": float(np.percentile(boot_aps, 97.5)),
        "n_boot": valid,
        "se": float(np.std(boot_aps)),
    }


def bootstrap_auprc_difference(pred_a: np.ndarray, pred_b: np.ndarray,
                               outcome: np.ndarray, threshold: float,
                               hiw_a: bool, hiw_b: bool,
                               n_bootstrap: int = N_BOOTSTRAP) -> dict:
    labels = (outcome >= threshold).astype(int)
    scores_a = -pred_a if hiw_a else pred_a
    scores_b = -pred_b if hiw_b else pred_b

    point_a = average_precision_score(labels, scores_a)
    point_b = average_precision_score(labels, scores_b)
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
        ap_a = average_precision_score(b_labels, scores_a[idx])
        ap_b = average_precision_score(b_labels, scores_b[idx])
        boot_diffs[valid] = ap_a - ap_b
        valid += 1
        attempts += 1

    boot_diffs = boot_diffs[:valid]
    ci_lo = float(np.percentile(boot_diffs, 2.5))
    ci_hi = float(np.percentile(boot_diffs, 97.5))
    p_value = float(np.mean(boot_diffs <= 0)) if point_diff > 0 else float(np.mean(boot_diffs >= 0))

    return {
        "point_diff": float(point_diff),
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "p_value": p_value,
        "n_boot": valid,
        "auprc_a": float(point_a),
        "auprc_b": float(point_b),
    }


def main():
    log("=== EXPERIMENT 13e: BOOTSTRAP AUPRC FOR ALL PREDICTOR-OUTCOME PAIRS ===")

    df = pd.read_csv(INSTABILITY_CSV)
    with open(HOLDOUT_JSON) as f:
        holdout = json.load(f)
    holdout_df = pd.DataFrame(holdout).rename(columns={"drug": "drug_name"})
    df_10 = df[df["n_cell_lines"] >= MIN_CELL_LINES].copy()
    merged = df_10.merge(holdout_df[["drug_name", "mean_heldout_cosine", "frac_consistent"]],
                         on="drug_name", how="inner")
    log(f"  {len(merged)} drugs with >= {MIN_CELL_LINES} cell lines")

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
    frac_outcome = merged["frac_consistent"].values

    # Report class balance
    n_cos_pos = (cosine_outcome >= CONSISTENCY_THRESHOLD).sum()
    n_jac_pos = (jaccard_outcome >= jaccard_threshold).sum()
    n_frac_pos = (frac_outcome >= 0.5).sum()
    n = len(merged)
    log(f"  Class balance — cosine: {n_cos_pos}/{n} ({n_cos_pos/n:.1%}), "
        f"jaccard: {n_jac_pos}/{n} ({n_jac_pos/n:.1%}), "
        f"frac_consistent: {n_frac_pos}/{n} ({n_frac_pos/n:.1%})")

    results = {
        "n_drugs": n,
        "n_bootstrap": N_BOOTSTRAP,
        "cosine_threshold": CONSISTENCY_THRESHOLD,
        "jaccard_threshold": jaccard_threshold,
        "class_balance": {
            "cosine_positive_frac": float(n_cos_pos / n),
            "jaccard_positive_frac": float(n_jac_pos / n),
            "frac_consistent_positive_frac": float(n_frac_pos / n),
        },
    }

    # --- Cosine outcome ---
    log("\n--- Cosine outcome AUPRC with 95% CI ---")
    results["cosine_auprc"] = {}
    for name, pred in predictors.items():
        r = bootstrap_auprc(pred, cosine_outcome, CONSISTENCY_THRESHOLD, HIGHER_IS_WORSE[name])
        results["cosine_auprc"][name] = r
        log(f"  {name:30s} {r['point']:.4f} [{r['ci_lo']:.4f}, {r['ci_hi']:.4f}]")

    # --- Jaccard outcome ---
    log("\n--- Jaccard outcome AUPRC with 95% CI ---")
    results["jaccard_auprc"] = {}
    for name, pred in predictors.items():
        if name == "mean_top_gene_jaccard":
            log(f"  {name:30s} SKIPPED (circular)")
            continue
        r = bootstrap_auprc(pred, jaccard_outcome, jaccard_threshold, HIGHER_IS_WORSE[name])
        results["jaccard_auprc"][name] = r
        log(f"  {name:30s} {r['point']:.4f} [{r['ci_lo']:.4f}, {r['ci_hi']:.4f}]")

    # --- Original paper outcome (frac_consistent > 0.5) ---
    log("\n--- frac_consistent > 0.5 AUPRC with 95% CI ---")
    results["frac_consistent_auprc"] = {}
    for name, pred in predictors.items():
        if name == "mean_top_gene_jaccard":
            continue
        r = bootstrap_auprc(pred, frac_outcome, 0.5, HIGHER_IS_WORSE[name])
        results["frac_consistent_auprc"][name] = r
        log(f"  {name:30s} {r['point']:.4f} [{r['ci_lo']:.4f}, {r['ci_hi']:.4f}]")

    # --- Paired AUPRC differences (DI vs each alternative) ---
    log("\n--- DI vs alternatives: paired AUPRC differences ---")
    di_pred = predictors["direction_instability"]
    di_hiw = HIGHER_IS_WORSE["direction_instability"]

    results["paired_differences"] = {}
    for outcome_name, outcome, threshold in [
        ("cosine", cosine_outcome, CONSISTENCY_THRESHOLD),
        ("jaccard", jaccard_outcome, jaccard_threshold),
        ("frac_consistent", frac_outcome, 0.5),
    ]:
        results["paired_differences"][outcome_name] = {}
        for name, pred in predictors.items():
            if name == "direction_instability":
                continue
            if name == "mean_top_gene_jaccard" and outcome_name == "jaccard":
                continue
            r = bootstrap_auprc_difference(di_pred, pred, outcome, threshold,
                                           di_hiw, HIGHER_IS_WORSE[name])
            sig = "*" if (r["ci_lo"] > 0 or r["ci_hi"] < 0) else "ns"
            results["paired_differences"][outcome_name][name] = r
            log(f"  {outcome_name:16s} DI - {name:25s} = {r['point_diff']:+.4f} "
                f"[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] p={r['p_value']:.4f} {sig}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outpath = OUTPUT_DIR / "bootstrap_auprc.json"
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {outpath}")


if __name__ == "__main__":
    main()
