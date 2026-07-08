"""HP3 permutation null: shuffle essentiality labels to test background divergence.

Tests whether the essential > non-essential geodesic distance difference
(Cohen's d = +0.55) exceeds what you'd get from random label assignment.
This separates essential-gene-specific divergence from the baseline
K562/RPE1 phenotypic divergence.

Uses pre-computed results from genetic-perturbation-holonomy repo.
Output: results/07_perturb_seq/hp3_permutation_null.json
"""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats

HOLONOMY_REPO = Path(__file__).parent.parent.parent / "genetic-perturbation-holonomy"
RESULTS_PATH = HOLONOMY_REPO / "results" / "01_subspace_transport" / "results.json"
DEPMAP_PATH = HOLONOMY_REPO / "data" / "extracted" / "frozen_depmap_labels.json"
OUTPUT_DIR = Path(__file__).parent.parent / "results" / "07_perturb_seq"
CHRONOS_THRESHOLD = -0.5
N_PERM = 10000


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    log("HP3 permutation null test")

    transport_results = json.load(open(RESULTS_PATH))
    per_gene = transport_results["per_gene"]
    depmap = json.load(open(DEPMAP_PATH))

    genes = []
    geodesic_distances = []
    is_essential = []

    for gene, data in per_gene.items():
        if gene not in depmap:
            continue
        k562_eff = depmap[gene]["k562_effect"]
        rpe1_eff = depmap[gene]["rpe1_effect"]
        essential = (k562_eff < CHRONOS_THRESHOLD) and (rpe1_eff < CHRONOS_THRESHOLD)

        genes.append(gene)
        geodesic_distances.append(data["geodesic_distance"])
        is_essential.append(essential)

    geodesic_distances = np.array(geodesic_distances)
    is_essential = np.array(is_essential)
    n_total = len(genes)
    n_essential = int(is_essential.sum())
    n_nonessential = n_total - n_essential

    log(f"  {n_total} genes with both geodesic distance and DepMap labels")
    log(f"  {n_essential} universal-essential, {n_nonessential} non-essential")

    ess_dists = geodesic_distances[is_essential]
    noness_dists = geodesic_distances[~is_essential]

    observed_diff = float(ess_dists.mean() - noness_dists.mean())
    observed_d = observed_diff / np.sqrt(
        ((len(ess_dists) - 1) * ess_dists.var() + (len(noness_dists) - 1) * noness_dists.var())
        / (len(ess_dists) + len(noness_dists) - 2)
    )
    observed_mw = stats.mannwhitneyu(ess_dists, noness_dists, alternative="two-sided")

    log(f"  Observed: essential mean={ess_dists.mean():.4f}, non-essential mean={noness_dists.mean():.4f}")
    log(f"  Observed diff={observed_diff:.4f}, Cohen's d={observed_d:.4f}")
    log(f"  Mann-Whitney p={observed_mw.pvalue:.2e}")

    rng = np.random.default_rng(seed=2026070707)
    null_diffs = np.zeros(N_PERM)

    log(f"  Running {N_PERM} permutations...")
    for i in range(N_PERM):
        shuffled = rng.permutation(is_essential)
        null_ess = geodesic_distances[shuffled]
        null_noness = geodesic_distances[~shuffled]
        null_diffs[i] = null_ess.mean() - null_noness.mean()

    p_value = float(np.mean(np.abs(null_diffs) >= abs(observed_diff)))
    null_mean = float(null_diffs.mean())
    null_std = float(null_diffs.std())

    result = {
        "n_genes": n_total,
        "n_essential": n_essential,
        "n_nonessential": n_nonessential,
        "chronos_threshold": CHRONOS_THRESHOLD,
        "observed_diff": round(observed_diff, 4),
        "observed_cohens_d": round(float(observed_d), 4),
        "observed_mann_whitney_p": float(observed_mw.pvalue),
        "essential_mean": round(float(ess_dists.mean()), 4),
        "nonessential_mean": round(float(noness_dists.mean()), 4),
        "n_permutations": N_PERM,
        "null_mean_diff": round(null_mean, 6),
        "null_std_diff": round(null_std, 4),
        "permutation_p_value": p_value,
        "null_quantiles": {
            "q025": round(float(np.quantile(null_diffs, 0.025)), 4),
            "q975": round(float(np.quantile(null_diffs, 0.975)), 4),
        },
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_DIR / "hp3_permutation_null.json", "w") as f:
        json.dump(result, f, indent=2)

    log(f"\n=== HP3 PERMUTATION NULL ===")
    log(f"Observed essential-nonessential diff: {observed_diff:.4f} (d={observed_d:.4f})")
    log(f"Null distribution: mean={null_mean:.6f}, std={null_std:.4f}")
    log(f"Null 95% interval: [{np.quantile(null_diffs, 0.025):.4f}, {np.quantile(null_diffs, 0.975):.4f}]")
    log(f"Permutation p-value (two-sided): {p_value:.4f}")

    if p_value < 0.05:
        log("PASS: essential-gene divergence exceeds background K562/RPE1 divergence")
    else:
        log("FAIL: essential-gene divergence is within background K562/RPE1 divergence")

    log(f"Saved to {OUTPUT_DIR / 'hp3_permutation_null.json'}")


if __name__ == "__main__":
    main()
