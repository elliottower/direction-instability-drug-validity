"""Experiment 12: Cosine direction instability on Perturb-seq same-gene transport.

Pre-registered in PREREGISTRATION_COMBINED_PAPER.md.

Tests whether direction instability (cosine-based) shows the same bin-inversion
confound pattern as geodesic distance and EMD on Perturb-seq same-gene transport.

Computes pseudobulk signatures per gene per cell type (K562, RPE1), then
direction instability as D = 1 - |cos(mean_K562, mean_RPE1)|. Evaluates:
  1. Aggregate Cohen's d for same-gene vs random-gene DI
  2. Bin-flatness: per-bin d across cell-count bins
  3. Correlation between DI and min(n_K562, n_RPE1)

Data: Replogle et al. 2022 Perturb-seq, downloaded from scPerturb Zenodo
  (record 10044268). Requires .h5ad files in
  ../genetic-perturbation-holonomy/data/raw/.

Usage:
    uv run python experiments/12_perturbseq_cosine_transport.py
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import stats
from tqdm import tqdm

HOLONOMY_DIR = Path(__file__).resolve().parent.parent.parent / "genetic-perturbation-holonomy"
SUBSPACE_PATH = HOLONOMY_DIR / "data" / "extracted" / "subspaces.npz"
METADATA_PATH = HOLONOMY_DIR / "data" / "extracted" / "metadata.json"
RAW_DIR = HOLONOMY_DIR / "data"
OUTPUT_DIR = Path("results/12_perturbseq_cosine")

K562_FILENAME = "ReplogleWeissman2022_K562_essential.h5ad"
RPE1_FILENAME = "ReplogleWeissman2022_rpe1.h5ad"

BINS = [(20, 40), (40, 60), (60, 100), (100, 200), (200, 700)]


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def load_pseudobulk_signatures(k562_path: Path, rpe1_path: Path, gene_names: list[str], hvg_features: list[str] | None = None):
    """Compute pseudobulk mean expression per gene per cell type.

    Returns:
        sigs_k562: (n_genes, n_features) pseudobulk for K562
        sigs_rpe1: (n_genes, n_features) pseudobulk for RPE1
        valid_genes: list of gene names with >= 20 cells in both
        counts_k562: per-gene cell counts in K562
        counts_rpe1: per-gene cell counts in RPE1
    """
    log("Loading K562 h5ad...")
    adata_k562 = ad.read_h5ad(k562_path)
    log(f"  K562 shape: {adata_k562.shape}")

    log("Loading RPE1 h5ad...")
    adata_rpe1 = ad.read_h5ad(rpe1_path)
    log(f"  RPE1 shape: {adata_rpe1.shape}")

    guide_col = None
    for candidate in ["gene", "perturbation", "guide_id", "gene_name"]:
        if candidate in adata_k562.obs.columns:
            guide_col = candidate
            break
    if guide_col is None:
        log(f"  Available obs columns (K562): {list(adata_k562.obs.columns)}")
        log("ERROR: No gene/perturbation column found in K562 obs")
        sys.exit(1)
    log(f"  Guide column: {guide_col}")

    if hvg_features is not None:
        k562_var = list(adata_k562.var_names)
        rpe1_var = list(adata_rpe1.var_names)
        shared_hvg = [g for g in hvg_features if g in k562_var and g in rpe1_var]
        log(f"  Using {len(shared_hvg)} shared HVG features")
        adata_k562 = adata_k562[:, shared_hvg]
        adata_rpe1 = adata_rpe1[:, shared_hvg]
    else:
        shared_var = sorted(set(adata_k562.var_names) & set(adata_rpe1.var_names))
        log(f"  Shared var features: {len(shared_var)}")
        adata_k562 = adata_k562[:, shared_var]
        adata_rpe1 = adata_rpe1[:, shared_var]

    sigs_k562 = {}
    counts_k562 = {}
    log("Computing K562 pseudobulk...")
    k562_genes = adata_k562.obs[guide_col]
    for gene in tqdm(gene_names, desc="K562 pseudobulk"):
        mask = k562_genes == gene
        n = mask.sum()
        if n >= 20:
            X = adata_k562[mask].X
            if hasattr(X, 'toarray'):
                X = X.toarray()
            sigs_k562[gene] = np.mean(X.astype(np.float64), axis=0)
            counts_k562[gene] = int(n)

    sigs_rpe1 = {}
    counts_rpe1 = {}
    log("Computing RPE1 pseudobulk...")
    rpe1_genes = adata_rpe1.obs[guide_col]
    for gene in tqdm(gene_names, desc="RPE1 pseudobulk"):
        mask = rpe1_genes == gene
        n = mask.sum()
        if n >= 20:
            X = adata_rpe1[mask].X
            if hasattr(X, 'toarray'):
                X = X.toarray()
            sigs_rpe1[gene] = np.mean(X.astype(np.float64), axis=0)
            counts_rpe1[gene] = int(n)

    valid_genes = sorted(set(sigs_k562) & set(sigs_rpe1))
    log(f"  Genes with >= 20 cells in both: {len(valid_genes)}")

    n_features = adata_k562.shape[1]
    arr_k562 = np.zeros((len(valid_genes), n_features), dtype=np.float64)
    arr_rpe1 = np.zeros((len(valid_genes), n_features), dtype=np.float64)
    arr_counts_k562 = np.zeros(len(valid_genes), dtype=int)
    arr_counts_rpe1 = np.zeros(len(valid_genes), dtype=int)

    for i, gene in enumerate(valid_genes):
        arr_k562[i] = sigs_k562[gene]
        arr_rpe1[i] = sigs_rpe1[gene]
        arr_counts_k562[i] = counts_k562[gene]
        arr_counts_rpe1[i] = counts_rpe1[gene]

    return arr_k562, arr_rpe1, valid_genes, arr_counts_k562, arr_counts_rpe1


def cosine_direction_instability(v1: np.ndarray, v2: np.ndarray) -> float:
    """D = 1 - |cos(v1, v2)| for two vectors."""
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < 1e-10 or n2 < 1e-10:
        return np.nan
    return 1.0 - abs(np.dot(v1, v2) / (n1 * n2))


def compute_same_gene_di(sigs_k562: np.ndarray, sigs_rpe1: np.ndarray) -> np.ndarray:
    """Compute DI for same-gene pairs."""
    n = sigs_k562.shape[0]
    di = np.zeros(n)
    for i in range(n):
        di[i] = cosine_direction_instability(sigs_k562[i], sigs_rpe1[i])
    return di


def compute_random_gene_di(sigs_k562: np.ndarray, sigs_rpe1: np.ndarray, n_permutations: int = 1000) -> np.ndarray:
    """Compute DI for random-gene pairs (permuted RPE1 assignment)."""
    n = sigs_k562.shape[0]
    all_di = np.zeros(n * n_permutations)
    for p in range(n_permutations):
        perm = np.random.permutation(n)
        for i in range(n):
            all_di[p * n + i] = cosine_direction_instability(sigs_k562[i], sigs_rpe1[perm[i]])
    return all_di


def cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """Cohen's d (pooled SD)."""
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_sd = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_sd < 1e-10:
        return 0.0
    return (np.mean(group1) - np.mean(group2)) / pooled_sd


def main():
    log("=== EXPERIMENT 12: COSINE DI ON PERTURB-SEQ TRANSPORT ===")
    log("Pre-registered in PREREGISTRATION_COMBINED_PAPER.md")

    if not METADATA_PATH.exists():
        log(f"ERROR: Metadata not found at {METADATA_PATH}")
        sys.exit(1)

    with open(METADATA_PATH) as f:
        metadata = json.load(f)

    gene_names = metadata["gene_names"]
    cell_counts_k562 = np.array(metadata["cell_counts_k562"])
    cell_counts_rpe1 = np.array(metadata["cell_counts_rpe1"])
    log(f"  {len(gene_names)} genes from metadata")

    k562_path = RAW_DIR / K562_FILENAME
    rpe1_path = RAW_DIR / RPE1_FILENAME

    if not k562_path.exists() or not rpe1_path.exists():
        log(f"ERROR: Raw h5ad files not found.")
        log(f"  Expected K562: {k562_path}")
        log(f"  Expected RPE1: {rpe1_path}")
        log("  Run the download script first:")
        log("    cd ../genetic-perturbation-holonomy && uv run python data/download.py")
        sys.exit(1)

    sigs_k562, sigs_rpe1, valid_genes, counts_k562, counts_rpe1 = (
        load_pseudobulk_signatures(k562_path, rpe1_path, gene_names)
    )

    n_genes = len(valid_genes)
    min_counts = np.minimum(counts_k562, counts_rpe1)

    log("\n--- Same-gene direction instability ---")
    same_di = compute_same_gene_di(sigs_k562, sigs_rpe1)
    log(f"  Same-gene DI: mean={same_di.mean():.4f}, std={same_di.std():.4f}")

    log("\n--- Random-gene direction instability (1000 permutations) ---")
    random_di = compute_random_gene_di(sigs_k562, sigs_rpe1, n_permutations=1000)
    log(f"  Random-gene DI: mean={random_di.mean():.4f}, std={random_di.std():.4f}")

    aggregate_d = cohens_d(same_di, random_di)
    _, aggregate_p = stats.mannwhitneyu(same_di, random_di, alternative='two-sided')
    log(f"\n  Aggregate Cohen's d: {aggregate_d:.4f}")
    log(f"  Mann-Whitney p: {aggregate_p:.2e}")
    log(f"  Decision criterion 1: |d| < 0.1 → {'PASS' if abs(aggregate_d) < 0.1 else 'FAIL'}")

    log("\n--- Bin-flatness test (PRIMARY RESULT) ---")
    bin_results = []
    for lo, hi in BINS:
        mask = (min_counts >= lo) & (min_counts < hi)
        n_in_bin = mask.sum()
        if n_in_bin < 5:
            log(f"  [{lo},{hi}): n={n_in_bin} — skipping (too few)")
            bin_results.append({
                "range": f"[{lo},{hi})", "n_genes": int(n_in_bin),
                "d_vs_random": None, "skipped": True
            })
            continue

        bin_same_di = same_di[mask]
        bin_random_di = compute_random_gene_di(
            sigs_k562[mask], sigs_rpe1[mask], n_permutations=500
        )
        bin_d = cohens_d(bin_same_di, bin_random_di)
        log(f"  [{lo},{hi}): n={n_in_bin}, d={bin_d:+.4f}")
        bin_results.append({
            "range": f"[{lo},{hi})", "n_genes": int(n_in_bin),
            "d_vs_random": float(bin_d), "skipped": False
        })

    valid_bin_ds = [b["d_vs_random"] for b in bin_results if not b.get("skipped")]
    all_in_band = all(abs(d) <= 0.15 for d in valid_bin_ds)
    log(f"\n  All bins within [-0.15, +0.15]: {'YES' if all_in_band else 'NO'}")
    log(f"  Decision criterion 2 (bin-flatness): {'PASS' if all_in_band else 'FAIL'}")

    if not all_in_band:
        max_d = max(abs(d) for d in valid_bin_ds)
        log(f"  Largest |d| across bins: {max_d:.4f}")
        ds_sorted = sorted(valid_bin_ds)
        log(f"  Bin d range: {ds_sorted[0]:+.4f} to {ds_sorted[-1]:+.4f}")
        if ds_sorted[0] < -0.15 and ds_sorted[-1] > 0.15:
            log("  WARNING: Bin inversion detected — same pattern as geodesic/EMD")

    log("\n--- Cell-count correlation (SECONDARY) ---")
    rho_count, p_count = stats.spearmanr(same_di, min_counts)
    log(f"  Spearman rho(DI, min_cell_count): {rho_count:.4f} (p={p_count:.2e})")
    log(f"  Decision criterion 3: |rho| < 0.20 → {'PASS' if abs(rho_count) < 0.20 else 'INVESTIGATE'}")
    if abs(rho_count) >= 0.20:
        log("  Investigating: does bin-inversion appear?")
        log("  (Bin-flatness test above is the primary immunity check)")

    log("\n--- Summary ---")
    c1 = abs(aggregate_d) < 0.1
    c2 = all_in_band
    c3 = abs(rho_count) < 0.20

    if c1 and c2:
        log("  RESULT: Confirms immunity. Cosine DI shows no transport signal")
        log("  and no bin-inversion. Confounded metrics' signal was entirely")
        log("  estimation-quality artifact.")
    elif c1 and not c2:
        log("  RESULT: Partial immunity. Aggregate d is small but bins are")
        log("  not flat — cosine is confounded in the same direction as")
        log("  geodesic/EMD but with smaller magnitude.")
    elif not c1 and c2:
        log("  RESULT: Real transport signal with flat bins. This would")
        log("  indicate genuine same-gene transport that confounded metrics missed.")
    else:
        log("  RESULT: Cosine is not immune. Shows both transport signal and")
        log("  bin-inversion. Undermines the combined paper's thesis.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {
        "experiment": "12_perturbseq_cosine_transport",
        "preregistered": "PREREGISTRATION_COMBINED_PAPER.md",
        "n_genes": n_genes,
        "aggregate": {
            "cohens_d": float(aggregate_d),
            "mann_whitney_p": float(aggregate_p),
            "same_gene_di_mean": float(same_di.mean()),
            "same_gene_di_std": float(same_di.std()),
            "random_gene_di_mean": float(random_di.mean()),
            "random_gene_di_std": float(random_di.std()),
            "criterion_1_pass": bool(c1),
        },
        "bin_flatness": {
            "bins": bin_results,
            "all_in_band": bool(all_in_band),
            "criterion_2_pass": bool(c2),
        },
        "cell_count_correlation": {
            "spearman_rho": float(rho_count),
            "p_value": float(p_count),
            "criterion_3_pass": bool(c3),
        },
    }

    with open(OUTPUT_DIR / "perturbseq_cosine_transport.json", "w") as f:
        json.dump(results, f, indent=2)

    np.savez(
        OUTPUT_DIR / "perturbseq_cosine_di.npz",
        genes=np.array(valid_genes, dtype=object),
        same_gene_di=same_di,
        min_cell_counts=min_counts,
        counts_k562=counts_k562,
        counts_rpe1=counts_rpe1,
    )
    log(f"\nSaved to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
