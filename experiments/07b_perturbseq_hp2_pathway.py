"""Experiment 07b: HP2 — Corrected distances predict GO pathway function.

Tests whether essential-gene-corrected Grassmannian geodesic distances
better predict Gene Ontology pathway coherence than raw distances.

Prediction: Spearman |rho| between corrected distance and GO pathway
Jaccard overlap > |rho| for raw distance, by at least 0.05.

Uses MSigDB C5 (GO Biological Process) gene sets. Downloads the GMT
file from MSigDB public FTP on first run, caches locally.

Usage:
    uv run python experiments/07b_perturbseq_hp2_pathway.py \
        --data ../genetic-perturbation-holonomy/data/extracted/ \
        --results results/07_perturbseq_correction/
"""
import argparse
import json
import sys
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats
from tqdm import tqdm

MSIGDB_C5_BP_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/"
    "c5.go.bp.v2023.2.Hs.symbols.gmt"
)


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def download_gmt(cache_path: Path) -> Path:
    """Download GO BP GMT file if not cached."""
    if cache_path.exists():
        log(f"  Using cached GMT: {cache_path}")
        return cache_path
    log(f"  Downloading GO BP gene sets from MSigDB...")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(MSIGDB_C5_BP_URL, cache_path)
    log(f"  Saved to {cache_path}")
    return cache_path


def parse_gmt(gmt_path: Path) -> dict[str, set[str]]:
    """Parse GMT file into {pathway_name: set of gene symbols}."""
    pathways = {}
    with open(gmt_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = set(parts[2:])
            pathways[name] = genes
    return pathways


def build_gene_pathway_sets(
    gene_names: list[str], pathways: dict[str, set[str]]
) -> dict[str, set[str]]:
    """For each gene in our dataset, collect which pathways it belongs to."""
    gene_to_pathways = {}
    for gene in gene_names:
        memberships = set()
        for pw_name, pw_genes in pathways.items():
            if gene in pw_genes:
                memberships.add(pw_name)
        if memberships:
            gene_to_pathways[gene] = memberships
    return gene_to_pathways


def pairwise_pathway_jaccard(
    gene_names: list[str],
    gene_to_pathways: dict[str, set[str]],
    sample_n: int = 50000,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample gene pairs and compute pathway Jaccard overlap.

    Returns (pair_indices, jaccard_values) where pair_indices is (N, 2)
    and jaccard_values is (N,).
    """
    if rng is None:
        rng = np.random.default_rng()

    annotated_indices = [i for i, g in enumerate(gene_names) if g in gene_to_pathways]
    n_annotated = len(annotated_indices)
    log(f"  {n_annotated}/{len(gene_names)} genes have GO BP annotations")

    if n_annotated < 100:
        log("  ERROR: Too few annotated genes for meaningful analysis")
        sys.exit(1)

    n_possible_pairs = n_annotated * (n_annotated - 1) // 2
    actual_n = min(sample_n, n_possible_pairs)
    log(f"  Sampling {actual_n} pairs from {n_possible_pairs} possible")

    pair_idx_a = rng.choice(n_annotated, size=actual_n * 2, replace=True)
    pair_idx_b = rng.choice(n_annotated, size=actual_n * 2, replace=True)
    mask = pair_idx_a != pair_idx_b
    pair_idx_a = pair_idx_a[mask][:actual_n]
    pair_idx_b = pair_idx_b[mask][:actual_n]

    pairs = np.column_stack([
        np.array(annotated_indices)[pair_idx_a],
        np.array(annotated_indices)[pair_idx_b],
    ])

    jaccards = np.zeros(len(pairs))
    for k in tqdm(range(len(pairs)), desc="Jaccard"):
        g1 = gene_names[pairs[k, 0]]
        g2 = gene_names[pairs[k, 1]]
        s1 = gene_to_pathways[g1]
        s2 = gene_to_pathways[g2]
        intersection = len(s1 & s2)
        union = len(s1 | s2)
        jaccards[k] = intersection / union if union > 0 else 0.0

    return pairs, jaccards


def run_hp2(data_dir: Path, results_dir: Path, output_dir: Path):
    """Test HP2: corrected distances better predict pathway function."""
    log("=== HP2: GO Pathway Function Prediction ===")

    cache_dir = output_dir.parent / "data" / "go_cache"
    gmt_path = download_gmt(cache_dir / "c5.go.bp.v2023.2.Hs.symbols.gmt")
    pathways = parse_gmt(gmt_path)
    log(f"  Loaded {len(pathways)} GO BP pathways")

    with open(data_dir / "metadata.json") as f:
        meta = json.load(f)
    gene_names = meta["gene_names"]

    corr_data = np.load(results_dir / "corrected_distances.npz")
    raw_dists = corr_data["raw_dists"]
    corr_dists = corr_data["corr_dists"]
    log(f"  Loaded distances for {len(raw_dists)} genes")

    log("Building gene-pathway membership sets...")
    gene_to_pathways = build_gene_pathway_sets(gene_names, pathways)

    log("Computing pairwise pathway Jaccard (sampled)...")
    rng = np.random.default_rng(42)
    pairs, jaccards = pairwise_pathway_jaccard(
        gene_names, gene_to_pathways, sample_n=50000, rng=rng
    )

    pair_raw_dists = np.abs(raw_dists[pairs[:, 0]] - raw_dists[pairs[:, 1]])
    pair_corr_dists = np.abs(corr_dists[pairs[:, 0]] - corr_dists[pairs[:, 1]])

    pair_raw_mean = (raw_dists[pairs[:, 0]] + raw_dists[pairs[:, 1]]) / 2
    pair_corr_mean = (corr_dists[pairs[:, 0]] + corr_dists[pairs[:, 1]]) / 2

    log("\nCorrelation: pathway Jaccard vs distance difference (|d_i - d_j|)")
    rho_raw_diff, p_raw_diff = stats.spearmanr(pair_raw_dists, jaccards)
    rho_corr_diff, p_corr_diff = stats.spearmanr(pair_corr_dists, jaccards)
    log(f"  Raw distance diff vs Jaccard:       rho={rho_raw_diff:.4f}, p={p_raw_diff:.2e}")
    log(f"  Corrected distance diff vs Jaccard: rho={rho_corr_diff:.4f}, p={p_corr_diff:.2e}")

    log("\nCorrelation: pathway Jaccard vs mean distance ((d_i + d_j)/2)")
    rho_raw_mean, p_raw_mean = stats.spearmanr(pair_raw_mean, jaccards)
    rho_corr_mean, p_corr_mean = stats.spearmanr(pair_corr_mean, jaccards)
    log(f"  Raw mean distance vs Jaccard:       rho={rho_raw_mean:.4f}, p={p_raw_mean:.2e}")
    log(f"  Corrected mean distance vs Jaccard: rho={rho_corr_mean:.4f}, p={p_corr_mean:.2e}")

    improvement_diff = abs(rho_corr_diff) - abs(rho_raw_diff)
    improvement_mean = abs(rho_corr_mean) - abs(rho_raw_mean)
    log(f"\n  Improvement (diff): {improvement_diff:+.4f}")
    log(f"  Improvement (mean): {improvement_mean:+.4f}")

    bites = improvement_diff >= 0.05 or improvement_mean >= 0.05
    log(f"\n  HP2 criterion (improvement >= 0.05): {'CONFIRMED' if bites else 'DOES NOT BITE'}")

    results = {
        "n_pathways": len(pathways),
        "n_annotated_genes": len(gene_to_pathways),
        "n_pairs_sampled": len(pairs),
        "distance_difference_correlation": {
            "raw_rho": float(rho_raw_diff),
            "raw_p": float(p_raw_diff),
            "corrected_rho": float(rho_corr_diff),
            "corrected_p": float(p_corr_diff),
            "improvement": float(improvement_diff),
        },
        "mean_distance_correlation": {
            "raw_rho": float(rho_raw_mean),
            "raw_p": float(p_raw_mean),
            "corrected_rho": float(rho_corr_mean),
            "corrected_p": float(p_corr_mean),
            "improvement": float(improvement_mean),
        },
        "hp2_bites": bool(bites),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "hp2_pathway_results.json", "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {output_dir / 'hp2_pathway_results.json'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--data", type=Path, required=True,
        help="Path to genetic-perturbation-holonomy/data/extracted/",
    )
    parser.add_argument(
        "--results", type=Path,
        default=Path("results/07_perturbseq_correction"),
        help="Path to HP1 results (corrected_distances.npz)",
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/07b_hp2_pathway"),
    )
    args = parser.parse_args()
    run_hp2(args.data, args.results, args.output)
