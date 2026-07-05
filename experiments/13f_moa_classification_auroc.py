"""Experiment 13f: MOA classification AUROC — the third construction-neutral outcome.

Can each metric predict a drug's mechanism-of-action class (machinery vs receptor)?
This outcome is completely independent of cosine consistency and Jaccard overlap:
the label comes from pharmacology, not geometry.

This addresses the devil's advocate critique that DI and Jaccard share
magnitude-invariance, so DI's Jaccard advantage could be structural rather than
biological. MOA labels have no geometric construction — if DI predicts them,
it captures real biology.

Computes:
- AUROC and AUPRC for each metric predicting machinery vs receptor
- Paired bootstrap differences (DI vs alternatives)
- Three-class ordinal prediction (machinery < kinase < receptor) via Spearman

Usage:
    uv run python experiments/13f_moa_classification_auroc.py
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

DRUG_REPO = Path(__file__).resolve().parent.parent.parent / "drug-perturbation-geometry"
INSTABILITY_CSV = DRUG_REPO / "zenodo_v1" / "drug_instability_8949.csv"
OUTPUT_DIR = Path("results/13_lincs_benchmark")

N_BOOTSTRAP = 10_000

sys.path.insert(0, str(Path(__file__).resolve().parent))
mod14 = __import__("importlib").import_module("14_jump_cp_moa_stratification")
classify_moa = mod14.classify_moa


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def bootstrap_binary_clf(scores: np.ndarray, labels: np.ndarray,
                         n_bootstrap: int = N_BOOTSTRAP, seed: int = 42) -> dict:
    point_auroc = roc_auc_score(labels, scores)
    point_auprc = average_precision_score(labels, scores)

    rng = np.random.default_rng(seed)
    n = len(labels)
    boot_aurocs = np.zeros(n_bootstrap)
    boot_auprcs = np.zeros(n_bootstrap)
    valid = 0
    attempts = 0
    while valid < n_bootstrap and attempts < n_bootstrap * 3:
        idx = rng.integers(0, n, size=n)
        b_labels = labels[idx]
        if b_labels.sum() == 0 or b_labels.sum() == len(b_labels):
            attempts += 1
            continue
        boot_aurocs[valid] = roc_auc_score(b_labels, scores[idx])
        boot_auprcs[valid] = average_precision_score(b_labels, scores[idx])
        valid += 1
        attempts += 1

    boot_aurocs = boot_aurocs[:valid]
    boot_auprcs = boot_auprcs[:valid]
    return {
        "auroc": float(point_auroc),
        "auroc_ci_lo": float(np.percentile(boot_aurocs, 2.5)),
        "auroc_ci_hi": float(np.percentile(boot_aurocs, 97.5)),
        "auroc_se": float(np.std(boot_aurocs)),
        "auprc": float(point_auprc),
        "auprc_ci_lo": float(np.percentile(boot_auprcs, 2.5)),
        "auprc_ci_hi": float(np.percentile(boot_auprcs, 97.5)),
        "auprc_se": float(np.std(boot_auprcs)),
        "n_boot": valid,
    }


def bootstrap_clf_difference(scores_a: np.ndarray, scores_b: np.ndarray,
                             labels: np.ndarray,
                             n_bootstrap: int = N_BOOTSTRAP, seed: int = 42) -> dict:
    point_a = roc_auc_score(labels, scores_a)
    point_b = roc_auc_score(labels, scores_b)
    point_diff = point_a - point_b

    rng = np.random.default_rng(seed)
    n = len(labels)
    boot_diffs = np.zeros(n_bootstrap)
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
    ci_lo = float(np.percentile(boot_diffs, 2.5))
    ci_hi = float(np.percentile(boot_diffs, 97.5))
    p_value = float(np.mean(boot_diffs <= 0)) if point_diff > 0 else float(np.mean(boot_diffs >= 0))

    return {
        "point_diff": float(point_diff),
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "p_value": p_value,
        "n_boot": valid,
        "auroc_a": float(point_a),
        "auroc_b": float(point_b),
    }


def bootstrap_spearman(x: np.ndarray, y: np.ndarray,
                       n_bootstrap: int = N_BOOTSTRAP, seed: int = 42) -> dict:
    point_r, point_p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    n = len(x)
    boot_rs = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        r, _ = stats.spearmanr(x[idx], y[idx])
        boot_rs[i] = r
    return {
        "rho": float(point_r),
        "p_value": float(point_p),
        "ci_lo": float(np.percentile(boot_rs, 2.5)),
        "ci_hi": float(np.percentile(boot_rs, 97.5)),
        "se": float(np.std(boot_rs)),
    }


def main():
    log("=== EXPERIMENT 13f: MOA CLASSIFICATION AUROC (THIRD OUTCOME) ===")
    log("  This outcome is pharmacological, not geometric — independent of")
    log("  cosine consistency and Jaccard overlap.")

    df = pd.read_csv(INSTABILITY_CSV)
    df["moa_class"] = df["moa"].apply(lambda x: classify_moa(x) if pd.notna(x) else "other")

    machinery = df[df["moa_class"] == "machinery"]
    receptor = df[df["moa_class"] == "receptor"]
    kinase = df[df["moa_class"] == "kinase"]
    log(f"\n  Drug counts — machinery: {len(machinery)}, kinase: {len(kinase)}, "
        f"receptor: {len(receptor)}, other/unlabeled: {len(df) - len(machinery) - len(receptor) - len(kinase)}")

    # --- Binary: machinery (1) vs receptor (0) ---
    binary_df = pd.concat([machinery, receptor])
    labels = (binary_df["moa_class"] == "machinery").astype(int).values
    log(f"\n  Binary classification: {labels.sum()} machinery vs {len(labels) - labels.sum()} receptor")
    log(f"  Base rate (machinery): {labels.mean():.3f}")

    predictors = {
        "direction_instability": -binary_df["direction_instability"].values,
        "frechet_variance": -binary_df["frechet_variance"].values,
        "magnitude_cv": -binary_df["magnitude_cv"].values,
    }

    results = {
        "experiment": "13f_moa_classification_auroc",
        "description": "Can each metric predict machinery vs receptor MOA class?",
        "n_machinery": int(labels.sum()),
        "n_receptor": int(len(labels) - labels.sum()),
        "base_rate_machinery": float(labels.mean()),
    }

    log("\n--- Binary AUROC + AUPRC (machinery vs receptor) ---")
    results["binary_classification"] = {}
    for name, scores in predictors.items():
        r = bootstrap_binary_clf(scores, labels, seed=42)
        results["binary_classification"][name] = r
        log(f"  {name:25s} AUROC={r['auroc']:.4f} [{r['auroc_ci_lo']:.4f}, {r['auroc_ci_hi']:.4f}]  "
            f"AUPRC={r['auprc']:.4f} [{r['auprc_ci_lo']:.4f}, {r['auprc_ci_hi']:.4f}]")

    # --- Paired differences ---
    log("\n--- DI vs alternatives: paired AUROC differences ---")
    di_scores = predictors["direction_instability"]
    results["paired_differences"] = {}
    for name, scores in predictors.items():
        if name == "direction_instability":
            continue
        r = bootstrap_clf_difference(di_scores, scores, labels, seed=42)
        sig = "*" if (r["ci_lo"] > 0 or r["ci_hi"] < 0) else "ns"
        results["paired_differences"][name] = r
        log(f"  DI - {name:25s} = {r['point_diff']:+.4f} [{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] "
            f"p={r['p_value']:.4f} {sig}")

    # --- Three-class ordinal: machinery < kinase < receptor ---
    log("\n--- Three-class ordinal (Spearman: metric rank vs MOA ordinal) ---")
    ordinal_df = pd.concat([machinery, kinase, receptor])
    ordinal_map = {"machinery": 0, "kinase": 1, "receptor": 2}
    ordinal_labels = ordinal_df["moa_class"].map(ordinal_map).values

    ordinal_predictors = {
        "direction_instability": ordinal_df["direction_instability"].values,
        "frechet_variance": ordinal_df["frechet_variance"].values,
        "magnitude_cv": ordinal_df["magnitude_cv"].values,
    }

    results["ordinal_spearman"] = {}
    results["ordinal_n"] = {
        "machinery": int(len(machinery)),
        "kinase": int(len(kinase)),
        "receptor": int(len(receptor)),
    }
    for name, pred in ordinal_predictors.items():
        r = bootstrap_spearman(pred, ordinal_labels, seed=42)
        results["ordinal_spearman"][name] = r
        log(f"  {name:25s} rho={r['rho']:.4f} [{r['ci_lo']:.4f}, {r['ci_hi']:.4f}] p={r['p_value']:.2e}")

    # --- Save ---
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outpath = OUTPUT_DIR / "moa_classification_auroc.json"
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {outpath}")


if __name__ == "__main__":
    main()
