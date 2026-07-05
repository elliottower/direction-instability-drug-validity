"""Experiment 10: DepMap confirmatory replication of essential-gene divergence.

Pre-registered in commit 98897ed (PREREGISTRATION_DEPMAP.md).
This script implements the exact analysis plan specified in that document.

Directional prediction: essential genes show MORE context-dependent
variability (higher SD of dependency profile after PC1 removal) than
non-essential genes. Cohen's d > 0.35 confirms replication.

Data source: DepMap Public 24Q4
  - CRISPRGeneEffect.csv (Chronos scores, genes x cell lines)
  - common_essentials.csv
  - nonessentials.csv

Download from: https://depmap.org/portal/download/all/
Place files in: ../data/depmap_24q4/

Usage:
    uv run python experiments/10_depmap_replication.py
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

DATA_DIR = Path(__file__).parent.parent / "data" / "depmap_24q4"
OUTPUT_DIR = Path(__file__).parent.parent / "results" / "10_depmap_replication"

DEPMAP_RELEASE = "DepMap Public 24Q4"
CHRONOS_FILE = "CRISPRGeneEffect.csv"
COMMON_ESSENTIALS_FILE = "common_essentials.csv"
NONESSENTIALS_FILE = "nonessentials.csv"

CHRONOS_THRESHOLD_ESSENTIAL = -0.5
CHRONOS_FRACTION_ESSENTIAL = 0.50
CHRONOS_THRESHOLD_NONESSENTIAL = -0.3
CHRONOS_FRACTION_NONESSENTIAL = 0.80

EFFECT_SIZE_FLOOR = 0.35
EFFECT_SIZE_ATTENUATED = 0.20
ALPHA = 0.01


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    n1, n2 = len(group1), len(group2)
    var1, var2 = group1.var(ddof=1), group2.var(ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0.0
    return (group1.mean() - group2.mean()) / pooled_std


def bootstrap_d(group1: np.ndarray, group2: np.ndarray, n_boot: int = 10000) -> dict:
    rng = np.random.default_rng()
    ds = np.empty(n_boot)
    for i in range(n_boot):
        b1 = rng.choice(group1, size=len(group1), replace=True)
        b2 = rng.choice(group2, size=len(group2), replace=True)
        ds[i] = cohens_d(b1, b2)
    return {
        "mean": float(ds.mean()),
        "ci_lower": float(np.percentile(ds, 2.5)),
        "ci_upper": float(np.percentile(ds, 97.5)),
    }


def main():
    log("=== EXPERIMENT 10: DEPMAP CONFIRMATORY REPLICATION ===")
    log(f"Pre-registered in commit 98897ed")
    log(f"Data source: {DEPMAP_RELEASE}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    chronos_path = DATA_DIR / CHRONOS_FILE
    essentials_path = DATA_DIR / COMMON_ESSENTIALS_FILE
    nonessentials_path = DATA_DIR / NONESSENTIALS_FILE

    if not chronos_path.exists():
        log(f"ERROR: {chronos_path} not found.")
        log(f"Download from https://depmap.org/portal/download/all/")
        log(f"Place CRISPRGeneEffect.csv in {DATA_DIR}/")
        sys.exit(1)

    log("Loading Chronos gene effect matrix...")
    chronos = pd.read_csv(chronos_path, index_col=0)
    log(f"  Shape: {chronos.shape[0]} cell lines x {chronos.shape[1]} genes")

    gene_names = [col.rsplit(" (", 1)[0] for col in chronos.columns]
    chronos.columns = gene_names
    log(f"  Genes after name extraction: {len(gene_names)}")

    log("Classifying genes by essentiality...")
    if essentials_path.exists() and nonessentials_path.exists():
        log("  Using DepMap common essential / nonessential gene lists")
        common_ess = set(pd.read_csv(essentials_path, header=None)[0].str.split(r" \(").str[0])
        noness_list = set(pd.read_csv(nonessentials_path, header=None)[0].str.split(r" \(").str[0])
        essential_genes = sorted(set(gene_names) & common_ess)
        nonessential_genes = sorted(set(gene_names) & noness_list)
        log(f"  Common essentials in data: {len(essential_genes)}")
        log(f"  Nonessentials in data: {len(nonessential_genes)}")
    else:
        log("  Gene lists not found; using Chronos threshold classification")
        log(f"  Essential: Chronos < {CHRONOS_THRESHOLD_ESSENTIAL} in "
            f">{CHRONOS_FRACTION_ESSENTIAL*100:.0f}% of lines")
        log(f"  Non-essential: Chronos > {CHRONOS_THRESHOLD_NONESSENTIAL} in "
            f">{CHRONOS_FRACTION_NONESSENTIAL*100:.0f}% of lines")
        n_lines = chronos.shape[0]
        essential_genes = []
        nonessential_genes = []
        for gene in tqdm(gene_names, desc="  Classifying"):
            vals = chronos[gene].dropna().values
            if len(vals) < n_lines * 0.5:
                continue
            frac_essential = (vals < CHRONOS_THRESHOLD_ESSENTIAL).mean()
            frac_nonessential = (vals > CHRONOS_THRESHOLD_NONESSENTIAL).mean()
            if frac_essential > CHRONOS_FRACTION_ESSENTIAL:
                essential_genes.append(gene)
            elif frac_nonessential > CHRONOS_FRACTION_NONESSENTIAL:
                nonessential_genes.append(gene)
        log(f"  Essential genes (threshold): {len(essential_genes)}")
        log(f"  Non-essential genes (threshold): {len(nonessential_genes)}")

    log("Computing gene dependency profiles...")
    gene_matrix = chronos.values.T
    gene_index = list(chronos.columns)
    log(f"  Gene matrix: {gene_matrix.shape[0]} genes x {gene_matrix.shape[1]} cell lines")

    log("Handling missing values (median imputation per gene)...")
    for i in tqdm(range(gene_matrix.shape[0]), desc="  Imputing"):
        col = gene_matrix[i]
        mask = np.isnan(col)
        if mask.any():
            gene_matrix[i, mask] = np.nanmedian(col)

    log("Removing first principal component (fitness confound)...")
    mean_profile = gene_matrix.mean(axis=0, keepdims=True)
    centered = gene_matrix - mean_profile
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    pc1 = Vt[0:1]
    projections = centered @ pc1.T
    fitness_component = projections @ pc1
    corrected = centered - fitness_component
    log(f"  PC1 explains {(S[0]**2 / (S**2).sum()) * 100:.1f}% of variance")

    log("Computing per-gene variability (SD of corrected profile)...")
    gene_sd = np.std(corrected, axis=1)

    ess_idx = [gene_index.index(g) for g in essential_genes if g in gene_index]
    noness_idx = [gene_index.index(g) for g in nonessential_genes if g in gene_index]
    log(f"  Essential genes with data: {len(ess_idx)}")
    log(f"  Non-essential genes with data: {len(noness_idx)}")

    ess_sd = gene_sd[ess_idx]
    noness_sd = gene_sd[noness_idx]

    log("\n--- PRIMARY ANALYSIS ---")
    log(f"  Essential SD: mean={ess_sd.mean():.4f}, median={np.median(ess_sd):.4f}")
    log(f"  Non-essential SD: mean={noness_sd.mean():.4f}, median={np.median(noness_sd):.4f}")

    d = cohens_d(ess_sd, noness_sd)
    log(f"  Cohen's d (essential - non-essential): {d:.4f}")

    U_stat, p_two = stats.mannwhitneyu(ess_sd, noness_sd, alternative='two-sided')
    _, p_one = stats.mannwhitneyu(ess_sd, noness_sd, alternative='greater')
    log(f"  Mann-Whitney p (two-tailed): {p_two:.2e}")
    log(f"  Mann-Whitney p (one-tailed, essential > non-essential): {p_one:.2e}")

    log("  Bootstrap CI for Cohen's d (10,000 iterations)...")
    boot = bootstrap_d(ess_sd, noness_sd, n_boot=10000)
    log(f"  Bootstrap d: {boot['mean']:.4f} [{boot['ci_lower']:.4f}, {boot['ci_upper']:.4f}]")

    log("\n--- VERDICT (per pre-registered decision criteria) ---")
    if d > EFFECT_SIZE_FLOOR and p_one < ALPHA:
        verdict = "CONFIRMS (replicates)"
        log(f"  {verdict}: d={d:.3f} > {EFFECT_SIZE_FLOOR}, p={p_one:.2e} < {ALPHA}")
    elif d > EFFECT_SIZE_ATTENUATED and p_one < 0.05:
        verdict = "ATTENUATED REPLICATION"
        log(f"  {verdict}: {EFFECT_SIZE_ATTENUATED} < d={d:.3f} < {EFFECT_SIZE_FLOOR}")
    elif d < 0:
        verdict = "FALSIFIED (sign reversal)"
        log(f"  {verdict}: d={d:.3f} is negative (essentials LESS variable)")
    else:
        verdict = "FAILS TO REPLICATE"
        log(f"  {verdict}: d={d:.3f} < {EFFECT_SIZE_ATTENUATED} or p > 0.05")

    log("\n--- FALSIFICATION CHECK: PC1 sensitivity ---")
    raw_gene_sd = np.std(centered, axis=1)
    raw_ess_sd = raw_gene_sd[ess_idx]
    raw_noness_sd = raw_gene_sd[noness_idx]
    d_raw = cohens_d(raw_ess_sd, raw_noness_sd)
    log(f"  d WITHOUT PC1 removal: {d_raw:.4f}")
    log(f"  d WITH PC1 removal: {d:.4f}")
    if d_raw != 0:
        change_pct = abs(d - d_raw) / abs(d_raw) * 100
        log(f"  Change: {change_pct:.1f}%")
        if change_pct > 50:
            log("  WARNING: PC1 removal changes effect by >50% — falsification criterion 3 triggered")
    else:
        change_pct = float('nan')

    log("\n--- ALTERNATIVE OPERATIONALIZATION: within-group correlation ---")
    log("  Computing mean pairwise correlation within essential vs non-essential...")
    n_sample = min(200, len(ess_idx), len(noness_idx))
    rng = np.random.default_rng(seed=None)
    ess_sample = rng.choice(ess_idx, size=n_sample, replace=False)
    noness_sample = rng.choice(noness_idx, size=n_sample, replace=False)

    ess_corr = np.corrcoef(corrected[ess_sample])[np.triu_indices(n_sample, k=1)]
    noness_corr = np.corrcoef(corrected[noness_sample])[np.triu_indices(n_sample, k=1)]
    log(f"  Essential within-group mean corr: {ess_corr.mean():.4f}")
    log(f"  Non-essential within-group mean corr: {noness_corr.mean():.4f}")
    log(f"  (Lower correlation = more heterogeneous = more context-dependent)")

    results = {
        "experiment": "10_depmap_replication",
        "prereg_commit": "98897ed",
        "depmap_release": DEPMAP_RELEASE,
        "n_cell_lines": int(chronos.shape[0]),
        "n_genes_total": len(gene_names),
        "n_essential": len(ess_idx),
        "n_nonessential": len(noness_idx),
        "classification_method": "common_essentials list" if essentials_path.exists() else "threshold",
        "pc1_variance_explained": float((S[0]**2 / (S**2).sum())),
        "primary_analysis": {
            "essential_sd_mean": float(ess_sd.mean()),
            "essential_sd_median": float(np.median(ess_sd)),
            "nonessential_sd_mean": float(noness_sd.mean()),
            "nonessential_sd_median": float(np.median(noness_sd)),
            "cohens_d": float(d),
            "bootstrap_d_mean": boot["mean"],
            "bootstrap_d_ci_lower": boot["ci_lower"],
            "bootstrap_d_ci_upper": boot["ci_upper"],
            "mann_whitney_p_two_tailed": float(p_two),
            "mann_whitney_p_one_tailed": float(p_one),
            "verdict": verdict,
        },
        "falsification_pc1_sensitivity": {
            "d_without_pc1_removal": float(d_raw),
            "d_with_pc1_removal": float(d),
            "change_percent": float(change_pct) if not np.isnan(change_pct) else None,
            "triggered": bool(change_pct > 50) if not np.isnan(change_pct) else False,
        },
        "alternative_within_group_correlation": {
            "essential_mean_corr": float(ess_corr.mean()),
            "nonessential_mean_corr": float(noness_corr.mean()),
            "n_sampled_per_group": n_sample,
            "prediction_confirmed": bool(ess_corr.mean() < noness_corr.mean()),
        },
        "timestamp": datetime.now().isoformat(),
    }

    output_path = OUTPUT_DIR / "depmap_replication_results.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nResults saved to {output_path}")
    log("\n=== DONE ===")


if __name__ == "__main__":
    main()
