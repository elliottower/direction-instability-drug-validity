"""Experiment 13c: Bootstrap CI for the original drug paper AUROC (0.986).

The drug transport paper reported LOO AUROC = 0.986 using mean_loo_di as the
predictor and frac_consistent > 0.5 as the binary outcome. This script
bootstraps that exact setup to produce a 95% CI.

Usage:
    uv run python experiments/13c_bootstrap_original_auroc.py
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

DRUG_REPO = Path(__file__).resolve().parent.parent.parent / "drug-perturbation-geometry"
INSTABILITY_CSV = DRUG_REPO / "zenodo_v1" / "drug_instability_8949.csv"
HOLDOUT_JSON = DRUG_REPO / "results" / "03_core_defenses" / "holdout_prediction.json"
OUTPUT_DIR = Path("results/13_lincs_benchmark")

MIN_CELL_LINES = 10
N_BOOTSTRAP = 10_000


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    log("=== EXPERIMENT 13c: BOOTSTRAP CI FOR ORIGINAL 0.986 AUROC ===")

    df = pd.read_csv(INSTABILITY_CSV)
    with open(HOLDOUT_JSON) as f:
        holdout = json.load(f)
    holdout_df = pd.DataFrame(holdout).rename(columns={"drug": "drug_name"})

    df_10 = df[df["n_cell_lines"] >= MIN_CELL_LINES].copy()
    merged = df_10.merge(
        holdout_df[["drug_name", "mean_loo_di", "frac_consistent"]],
        on="drug_name", how="inner"
    )
    log(f"  {len(merged)} drugs")

    labels = (merged["frac_consistent"].values > 0.5).astype(int)
    scores = -merged["mean_loo_di"].values
    n_pos = labels.sum()
    n_neg = len(labels) - n_pos
    log(f"  Positive (frac_consistent > 0.5): {n_pos}")
    log(f"  Negative: {n_neg}")

    point = roc_auc_score(labels, scores)
    log(f"  Point estimate: {point:.4f}")

    rng = np.random.default_rng(42)
    n = len(labels)
    boot_aurocs = []
    attempts = 0
    while len(boot_aurocs) < N_BOOTSTRAP and attempts < N_BOOTSTRAP * 3:
        idx = rng.integers(0, n, size=n)
        b_labels = labels[idx]
        if b_labels.sum() == 0 or b_labels.sum() == len(b_labels):
            attempts += 1
            continue
        boot_aurocs.append(roc_auc_score(b_labels, scores[idx]))
        attempts += 1

    boot_aurocs = np.array(boot_aurocs)
    ci_lo = float(np.percentile(boot_aurocs, 2.5))
    ci_hi = float(np.percentile(boot_aurocs, 97.5))
    se = float(np.std(boot_aurocs))
    log(f"  Bootstrap ({len(boot_aurocs)} resamples): {np.mean(boot_aurocs):.4f} [{ci_lo:.4f}, {ci_hi:.4f}]")
    log(f"  SE: {se:.4f}")

    result = {
        "experiment": "13c_bootstrap_original_auroc",
        "description": "Bootstrap CI for original drug paper AUROC (LOO DI -> frac_consistent > 0.5)",
        "n_drugs": len(merged),
        "n_positive": int(n_pos),
        "n_negative": int(n_neg),
        "predictor": "mean_loo_di",
        "outcome": "frac_consistent > 0.5",
        "point_auroc": float(point),
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "se": se,
        "n_bootstrap": len(boot_aurocs),
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outpath = OUTPUT_DIR / "bootstrap_original_986.json"
    with open(outpath, "w") as f:
        json.dump(result, f, indent=2)
    log(f"  Saved to {outpath}")


if __name__ == "__main__":
    main()
