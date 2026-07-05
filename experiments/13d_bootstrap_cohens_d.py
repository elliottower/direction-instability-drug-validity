"""Experiment 13d: Bootstrap CIs for all Cohen's d values reported in the paper.

Bootstraps every group-comparison d value using actual per-unit data:
- LINCS MOA stratification (machinery vs receptor, from 8949-drug CSV)
- JUMP-CP MOA stratification (machinery vs receptor, from matched_drugs.csv)
- Perturb-seq cosine transport (same-gene vs random-gene, from npz)

DepMap values (d=2.65, d=-1.33) already have bootstrap CIs from the
magnitude confound paper and are not re-bootstrapped here.

Usage:
    uv run python experiments/13d_bootstrap_cohens_d.py
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

OUTPUT_DIR = Path("results/13_lincs_benchmark")
N_BOOTSTRAP = 10_000

DRUG_REPO = Path(__file__).resolve().parent.parent.parent / "drug-perturbation-geometry"
INSTABILITY_CSV = DRUG_REPO / "zenodo_v1" / "drug_instability_8949.csv"
JUMP_MATCHED = Path("results/14_jump_cp_moa/matched_drugs.csv")
PERTURBSEQ_NPZ = Path("results/12_perturbseq_cosine/perturbseq_cosine_di.npz")
PERTURBSEQ_JSON = Path("results/12_perturbseq_cosine/perturbseq_cosine_transport.json")


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def cohens_d(g1: np.ndarray, g2: np.ndarray) -> float:
    n1, n2 = len(g1), len(g2)
    var1, var2 = np.var(g1, ddof=1), np.var(g2, ddof=1)
    pooled_sd = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_sd < 1e-10:
        return 0.0
    return (np.mean(g1) - np.mean(g2)) / pooled_sd


def bootstrap_d(g1: np.ndarray, g2: np.ndarray, n_bootstrap: int = N_BOOTSTRAP,
                seed: int = 42) -> dict:
    """Bootstrap CI for Cohen's d between two groups using actual data."""
    point = cohens_d(g1, g2)
    rng = np.random.default_rng(seed)
    boot_ds = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        b1 = g1[rng.integers(0, len(g1), size=len(g1))]
        b2 = g2[rng.integers(0, len(g2), size=len(g2))]
        boot_ds[i] = cohens_d(b1, b2)
    return {
        "point": float(point),
        "ci_lo": float(np.percentile(boot_ds, 2.5)),
        "ci_hi": float(np.percentile(boot_ds, 97.5)),
        "se": float(np.std(boot_ds)),
        "n1": len(g1),
        "n2": len(g2),
    }


def main():
    log("=== EXPERIMENT 13d: BOOTSTRAP CIs FOR ALL COHEN'S d VALUES ===")
    results = {}

    # --- LINCS MOA stratification (actual per-drug DI values) ---
    log("\n--- LINCS MOA stratification (raw data from 8949-drug CSV) ---")
    if INSTABILITY_CSV.exists():
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        mod14 = __import__("importlib").import_module("14_jump_cp_moa_stratification")
        classify_moa = mod14.classify_moa

        df = pd.read_csv(INSTABILITY_CSV)
        df["moa_class"] = df["moa"].apply(lambda x: classify_moa(x) if pd.notna(x) else "other")

        machinery = df[df["moa_class"] == "machinery"]["direction_instability"].values
        receptor = df[df["moa_class"] == "receptor"]["direction_instability"].values
        kinase = df[df["moa_class"] == "kinase"]["direction_instability"].values
        log(f"  Machinery n={len(machinery)}, Receptor n={len(receptor)}, Kinase n={len(kinase)}")

        r = bootstrap_d(machinery, receptor, seed=42)
        r["description"] = "Machinery vs receptor DI in LINCS (lower = more consistent)"
        results["lincs_machinery_vs_receptor"] = r
        log(f"  Machinery vs receptor: d = {r['point']:.3f} [{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]")

        r2 = bootstrap_d(machinery, kinase, seed=43)
        r2["description"] = "Machinery vs kinase DI in LINCS"
        results["lincs_machinery_vs_kinase"] = r2
        log(f"  Machinery vs kinase:   d = {r2['point']:.3f} [{r2['ci_lo']:.3f}, {r2['ci_hi']:.3f}]")

        r3 = bootstrap_d(kinase, receptor, seed=44)
        r3["description"] = "Kinase vs receptor DI in LINCS"
        results["lincs_kinase_vs_receptor"] = r3
        log(f"  Kinase vs receptor:    d = {r3['point']:.3f} [{r3['ci_lo']:.3f}, {r3['ci_hi']:.3f}]")
    else:
        log(f"  SKIPPED: {INSTABILITY_CSV} not found")

    # --- JUMP-CP MOA stratification (actual per-drug DI values) ---
    log("\n--- JUMP-CP MOA stratification (raw data from matched_drugs.csv) ---")
    if JUMP_MATCHED.exists():
        jdf = pd.read_csv(JUMP_MATCHED)
        j_machinery = jdf[jdf["moa_class"] == "machinery"]["jump_di"].values
        j_receptor = jdf[jdf["moa_class"] == "receptor"]["jump_di"].values
        j_kinase = jdf[jdf["moa_class"] == "kinase"]["jump_di"].values
        log(f"  Machinery n={len(j_machinery)}, Receptor n={len(j_receptor)}, Kinase n={len(j_kinase)}")

        r = bootstrap_d(j_machinery, j_receptor, seed=45)
        r["description"] = "Machinery vs receptor DI in JUMP-CP"
        results["jump_cp_machinery_vs_receptor"] = r
        log(f"  Machinery vs receptor: d = {r['point']:.3f} [{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]")

        r2 = bootstrap_d(j_machinery, j_kinase, seed=46)
        r2["description"] = "Machinery vs kinase DI in JUMP-CP"
        results["jump_cp_machinery_vs_kinase"] = r2
        log(f"  Machinery vs kinase:   d = {r2['point']:.3f} [{r2['ci_lo']:.3f}, {r2['ci_hi']:.3f}]")

        r3 = bootstrap_d(j_kinase, j_receptor, seed=47)
        r3["description"] = "Kinase vs receptor DI in JUMP-CP"
        results["jump_cp_kinase_vs_receptor"] = r3
        log(f"  Kinase vs receptor:    d = {r3['point']:.3f} [{r3['ci_lo']:.3f}, {r3['ci_hi']:.3f}]")
    else:
        log(f"  SKIPPED: {JUMP_MATCHED} not found")

    # --- Perturb-seq cosine transport (actual per-gene DI values) ---
    log("\n--- Perturb-seq cosine transport (raw data from npz) ---")
    if PERTURBSEQ_NPZ.exists() and PERTURBSEQ_JSON.exists():
        data = np.load(PERTURBSEQ_NPZ, allow_pickle=True)
        same_di = data["same_gene_di"]
        n_genes = len(same_di)

        with open(PERTURBSEQ_JSON) as f:
            exp12 = json.load(f)
        random_mean = exp12["aggregate"]["random_gene_di_mean"]
        random_std = exp12["aggregate"]["random_gene_di_std"]

        rng = np.random.default_rng(48)
        boot_ds = np.zeros(N_BOOTSTRAP)
        for i in range(N_BOOTSTRAP):
            b_same = same_di[rng.integers(0, n_genes, size=n_genes)]
            b_random = rng.normal(random_mean, random_std, size=n_genes)
            boot_ds[i] = cohens_d(b_same, b_random)

        results["perturbseq_cosine_transport"] = {
            "description": "Same-gene vs random-gene cosine DI (random from permutation distribution)",
            "point": float(exp12["aggregate"]["cohens_d"]),
            "ci_lo": float(np.percentile(boot_ds, 2.5)),
            "ci_hi": float(np.percentile(boot_ds, 97.5)),
            "se": float(np.std(boot_ds)),
            "n_same": n_genes,
            "note": "Random group bootstrapped from N(mean, std) of 1000-permutation null",
        }
        r = results["perturbseq_cosine_transport"]
        log(f"  Same vs random: d = {r['point']:.3f} [{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]")
    else:
        log(f"  SKIPPED: npz or json not found")

    # --- Save ---
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outpath = OUTPUT_DIR / "bootstrap_cohens_d.json"
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {outpath}")


if __name__ == "__main__":
    main()
