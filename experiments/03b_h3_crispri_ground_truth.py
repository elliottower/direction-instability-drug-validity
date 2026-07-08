"""Experiment 03b: Phenotype projection with CRISPRi ground truth (H3-CRISPRi).

Replaces shRNA consensus signatures with CRISPRi (Replogle et al. 2022)
expression signatures as the "on-target direction" for phenotype-projected
bracket. CRISPRi has fewer off-target effects than shRNA, so if H3
strengthens, this confirms shRNA noise attenuated the original result.

Pre-registered in PREREGISTRATION_EXTENDED.md before running on real data.

Decision criterion: rho_CRISPRi >= rho_shRNA (0.376).

Reports BOTH shRNA and CRISPRi results side by side for the same drug set
(intersection of drugs with both ground truth types).

Usage:
    uv run python experiments/03b_h3_crispri_ground_truth.py --synthetic
    uv run python experiments/03b_h3_crispri_ground_truth.py --real \
        --data PATH_TO_LINCS --perturbseq PATH_TO_H5AD
"""
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from geometry.bracket_norm import direction_instability, phenotype_projected_bracket


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def cosine_squared(mean_sig: np.ndarray, direction: np.ndarray) -> float:
    norm_sig = np.linalg.norm(mean_sig)
    norm_dir = np.linalg.norm(direction)
    if norm_sig < 1e-10 or norm_dir < 1e-10:
        return 0.0
    cos = float(mean_sig @ direction / (norm_sig * norm_dir))
    return cos ** 2


def build_crispri_signatures(h5ad_path: Path, landmark_symbols: list[str]) -> dict[str, np.ndarray]:
    """Extract per-gene CRISPRi expression signatures projected onto LINCS landmarks.

    For each perturbed gene, computes mean expression change (perturbed - control)
    across all perturbed cells, then subsets to the LINCS 978 landmark genes.
    """
    log(f"Loading Perturb-seq data from {h5ad_path}...")
    adata = ad.read_h5ad(h5ad_path)
    log(f"  {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    perturbseq_genes = list(adata.var_names)
    landmark_set = set(landmark_symbols)
    shared_genes = [g for g in perturbseq_genes if g in landmark_set]
    log(f"  {len(shared_genes)} genes shared between Perturb-seq and LINCS landmarks")

    if len(shared_genes) < 100:
        log(f"  WARNING: only {len(shared_genes)} shared genes, results may be noisy")

    shared_idx_perturbseq = [perturbseq_genes.index(g) for g in shared_genes]
    landmark_order = {g: i for i, g in enumerate(landmark_symbols)}
    shared_landmark_idx = [landmark_order[g] for g in shared_genes]

    control_mask = adata.obs["gene"] == "non-targeting"
    if control_mask.sum() == 0:
        control_genes = ["non-targeting", "control", ""]
        for cg in control_genes:
            control_mask = adata.obs["gene"] == cg
            if control_mask.sum() > 0:
                break
    log(f"  {control_mask.sum():,} control cells")

    if control_mask.sum() == 0:
        log("  ERROR: no control cells found")
        sys.exit(1)

    X = adata.X
    if hasattr(X, 'toarray'):
        log("  Converting sparse matrix to dense (control subset)...")
        control_expr = X[control_mask][:, shared_idx_perturbseq].toarray()
    else:
        control_expr = X[control_mask][:, shared_idx_perturbseq]
    control_mean = control_expr.mean(axis=0).A1 if hasattr(control_expr, 'A1') else np.asarray(control_expr.mean(axis=0)).flatten()

    unique_genes = [g for g in adata.obs["gene"].unique() if g not in ("non-targeting", "control", "")]
    log(f"  Computing CRISPRi signatures for {len(unique_genes)} genes...")

    signatures = {}
    for gene in tqdm(unique_genes, desc="CRISPRi signatures"):
        gene_mask = adata.obs["gene"] == gene
        n_cells = gene_mask.sum()
        if n_cells < 10:
            continue

        if hasattr(X, 'toarray'):
            gene_expr = X[gene_mask][:, shared_idx_perturbseq].toarray()
        else:
            gene_expr = X[gene_mask][:, shared_idx_perturbseq]
        gene_mean = gene_expr.mean(axis=0).A1 if hasattr(gene_expr, 'A1') else np.asarray(gene_expr.mean(axis=0)).flatten()

        diff = gene_mean - control_mean

        full_sig = np.zeros(len(landmark_symbols))
        for i, lidx in enumerate(shared_landmark_idx):
            full_sig[lidx] = diff[i]

        norm = np.linalg.norm(full_sig)
        if norm > 1e-10:
            signatures[gene] = full_sig

    log(f"  {len(signatures)} CRISPRi signatures with nonzero norm")
    return signatures


def build_shrna_consensus(shrna_sigs, shrna_siginfo):
    """Build consensus shRNA signature per target gene."""
    sig_ids = list(shrna_siginfo["sig_id"])
    sig_id_to_idx = {sid: i for i, sid in enumerate(sig_ids)}
    assert len(set(sig_ids)) == len(sig_ids), "Duplicate sig_ids in shRNA data"

    filtered = shrna_siginfo[shrna_siginfo["sig_id"].isin(sig_id_to_idx)].copy()
    filtered["_idx"] = filtered["sig_id"].map(sig_id_to_idx)

    consensus = {}
    for gene, group in filtered.groupby("pert_iname"):
        indices = group["_idx"].values
        if len(indices) < 3:
            continue
        consensus[gene] = shrna_sigs[indices].mean(axis=0)
    return consensus


def run_synthetic():
    """Validate: CRISPRi ground truth should give same or stronger correlation."""
    log("=== SYNTHETIC MODE ===")
    rng = np.random.default_rng(seed=2026070803)
    n_genes = 200
    n_contexts = 8
    n_drugs = 100

    results = []
    for i in tqdm(range(n_drugs), desc="Drugs"):
        target_dir = rng.standard_normal(n_genes)
        target_dir /= np.linalg.norm(target_dir)

        on_target_strength = rng.uniform(0, 1)
        sigs = np.array([
            target_dir * on_target_strength * 2
            + rng.standard_normal(n_genes) * 0.5
            for _ in range(n_contexts)
        ])
        crispri_dir = target_dir + rng.standard_normal(n_genes) * 0.1
        shrna_dir = target_dir + rng.standard_normal(n_genes) * 0.4

        raw = direction_instability(sigs)
        proj_crispri = phenotype_projected_bracket(sigs, crispri_dir)
        proj_shrna = phenotype_projected_bracket(sigs, shrna_dir)
        enrich_crispri = cosine_squared(sigs.mean(axis=0), crispri_dir)
        enrich_shrna = cosine_squared(sigs.mean(axis=0), shrna_dir)

        results.append({
            "raw": raw,
            "proj_crispri": proj_crispri,
            "proj_shrna": proj_shrna,
            "enrich_crispri": enrich_crispri,
            "enrich_shrna": enrich_shrna,
        })

    rho_crispri = stats.spearmanr(
        [r["proj_crispri"] for r in results],
        [r["enrich_crispri"] for r in results],
    ).statistic
    rho_shrna = stats.spearmanr(
        [r["proj_shrna"] for r in results],
        [r["enrich_shrna"] for r in results],
    ).statistic

    log(f"\nCRISPRi projected rho: {rho_crispri:.4f}")
    log(f"shRNA projected rho:   {rho_shrna:.4f}")
    log(f"CRISPRi >= shRNA: {'PASS' if rho_crispri >= rho_shrna else 'FAIL'}")


def run_real(data_dir: Path, perturbseq_path: Path, output_dir: Path):
    """Run H3-CRISPRi on real data, reporting both ground truths side by side."""
    log("=== REAL MODE: H3-CRISPRi ===")

    sigs_path = data_dir / "lincs_subset.npz"
    siginfo_path = data_dir / "GSE92742_Broad_LINCS_sig_info.txt.gz"
    gene_info_path = data_dir / "GSE92742_Broad_LINCS_gene_info.txt.gz"
    shrna_path = data_dir / "lincs_shrna.npz"
    shrna_siginfo_path = data_dir / "lincs_shrna_siginfo.csv.gz"
    labels_path = data_dir / "frozen_drug_labels.json"

    for p in [sigs_path, gene_info_path, labels_path]:
        if not p.exists():
            log(f"ERROR: {p} not found.")
            sys.exit(1)

    log("Loading LINCS gene info for symbol mapping...")
    gene_info = pd.read_csv(gene_info_path, sep="\t")
    landmark = gene_info[gene_info["pr_is_lm"] == 1].copy()
    entrez_to_symbol = dict(zip(landmark["pr_gene_id"].astype(str), landmark["pr_gene_symbol"]))

    log("Loading compound signatures...")
    data = np.load(sigs_path, allow_pickle=True)
    compound_sigs = data["signatures"]
    compound_sig_ids = list(data["sig_ids"])
    compound_gene_ids = list(data["gene_ids"])
    log(f"  {compound_sigs.shape[0]:,} signatures x {compound_sigs.shape[1]} genes")

    landmark_symbols = [entrez_to_symbol.get(gid, gid) for gid in compound_gene_ids]
    n_mapped = sum(1 for s in landmark_symbols if not s.isdigit())
    log(f"  {n_mapped}/{len(landmark_symbols)} gene IDs mapped to symbols")

    log("Building CRISPRi signatures...")
    crispri_consensus = build_crispri_signatures(perturbseq_path, landmark_symbols)

    has_shrna = shrna_path.exists() and shrna_siginfo_path.exists()
    shrna_consensus = {}
    if has_shrna:
        log("Loading shRNA signatures for comparison...")
        shrna_data = np.load(shrna_path, allow_pickle=True)
        shrna_sigs = shrna_data["signatures"]
        shrna_siginfo = pd.read_csv(shrna_siginfo_path)
        shrna_consensus = build_shrna_consensus(shrna_sigs, shrna_siginfo)
        log(f"  {len(shrna_consensus)} shRNA consensus signatures")
    else:
        log("  shRNA data not found, reporting CRISPRi only")

    log("Loading drug labels...")
    with open(labels_path) as f:
        labels = json.load(f)
    drug_targets = {}
    for d in labels["drugs"]:
        target = d.get("target")
        if target:
            drug_targets[d["pert_iname"]] = target.split("|")[0].strip()

    log("Loading compound sig info and building drug-cell matrix...")
    siginfo = pd.read_csv(siginfo_path, sep="\t", low_memory=False)
    sig_id_to_idx = {sid: i for i, sid in enumerate(compound_sig_ids)}
    filtered = siginfo[siginfo["sig_id"].isin(sig_id_to_idx) & siginfo["pert_iname"].notna()].copy()
    filtered["_idx"] = filtered["sig_id"].map(sig_id_to_idx)

    drug_cell_map = {}
    for (drug, cell), group in filtered.groupby(["pert_iname", "cell_id"]):
        if drug not in drug_cell_map:
            drug_cell_map[drug] = {}
        indices = group["_idx"].values
        drug_cell_map[drug][cell] = compound_sigs[indices].mean(axis=0)

    n_drugs_with_target = sum(1 for d in drug_cell_map if d in drug_targets and len(drug_cell_map[d]) >= 5)
    n_crispri_eligible = sum(1 for d in drug_cell_map if d in drug_targets and drug_targets[d] in crispri_consensus and len(drug_cell_map[d]) >= 5)
    n_shrna_eligible = sum(1 for d in drug_cell_map if d in drug_targets and drug_targets[d] in shrna_consensus and len(drug_cell_map[d]) >= 5) if has_shrna else 0
    log(f"\n=== POST-MAPPING ELIGIBLE DRUG COUNTS ===")
    log(f"  Drugs with target annotation and >= 5 cell lines: {n_drugs_with_target}")
    log(f"  Drugs with CRISPRi ground truth available: {n_crispri_eligible}")
    log(f"  Drugs with shRNA ground truth available:   {n_shrna_eligible}")
    assert n_crispri_eligible >= 10, (
        f"Only {n_crispri_eligible} drugs have CRISPRi ground truth — too few for correlation. "
        f"Check gene-name mapping (shared genes: {len([g for g in landmark_symbols if not g.isdigit()])})"
    )

    log("Computing phenotype-projected brackets with both ground truths...")
    results = []
    for drug, cells in tqdm(drug_cell_map.items(), desc="Drugs"):
        if drug not in drug_targets or len(cells) < 5:
            continue
        target = drug_targets[drug]

        has_crispri = target in crispri_consensus
        has_shrna_gt = target in shrna_consensus

        if not has_crispri and not has_shrna_gt:
            continue

        sigs_matrix = np.array(list(cells.values()))
        mean_sig = sigs_matrix.mean(axis=0)
        raw = direction_instability(sigs_matrix)

        entry = {
            "drug": drug,
            "target": target,
            "n_celllines": len(cells),
            "raw_bracket": float(raw),
        }

        if has_crispri:
            crispri_dir = crispri_consensus[target]
            entry["proj_crispri"] = float(phenotype_projected_bracket(sigs_matrix, crispri_dir))
            entry["enrich_crispri"] = float(cosine_squared(mean_sig, crispri_dir))

        if has_shrna_gt:
            shrna_dir = shrna_consensus[target]
            entry["proj_shrna"] = float(phenotype_projected_bracket(sigs_matrix, shrna_dir))
            entry["enrich_shrna"] = float(cosine_squared(mean_sig, shrna_dir))

        results.append(entry)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "h3_crispri_results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    crispri_results = [r for r in results if "proj_crispri" in r]
    shrna_results = [r for r in results if "proj_shrna" in r]
    both_results = [r for r in results if "proj_crispri" in r and "proj_shrna" in r]

    log(f"\n=== H3-CRISPRi RESULTS ===")
    log(f"Drugs with CRISPRi ground truth: {len(crispri_results)}")
    log(f"Drugs with shRNA ground truth:   {len(shrna_results)}")
    log(f"Drugs with BOTH:                 {len(both_results)}")

    if len(crispri_results) >= 10:
        rho_proj_c, p_proj_c = stats.spearmanr(
            [r["proj_crispri"] for r in crispri_results],
            [r["enrich_crispri"] for r in crispri_results],
        )
        rho_raw_c, p_raw_c = stats.spearmanr(
            [r["raw_bracket"] for r in crispri_results],
            [r["enrich_crispri"] for r in crispri_results],
        )
        log(f"\nCRISPRi ground truth (n={len(crispri_results)}):")
        log(f"  Projected bracket rho: {rho_proj_c:.4f} (p={p_proj_c:.2e})")
        log(f"  Raw bracket rho:       {rho_raw_c:.4f} (p={p_raw_c:.2e})")
    else:
        log("\n  Too few drugs with CRISPRi ground truth for correlation")
        rho_proj_c = None

    if len(shrna_results) >= 10:
        rho_proj_s, p_proj_s = stats.spearmanr(
            [r["proj_shrna"] for r in shrna_results],
            [r["enrich_shrna"] for r in shrna_results],
        )
        rho_raw_s, p_raw_s = stats.spearmanr(
            [r["raw_bracket"] for r in shrna_results],
            [r["enrich_shrna"] for r in shrna_results],
        )
        log(f"\nshRNA ground truth (n={len(shrna_results)}):")
        log(f"  Projected bracket rho: {rho_proj_s:.4f} (p={p_proj_s:.2e})")
        log(f"  Raw bracket rho:       {rho_raw_s:.4f} (p={p_raw_s:.2e})")
    else:
        rho_proj_s = None

    if len(both_results) >= 10:
        rho_both_c, _ = stats.spearmanr(
            [r["proj_crispri"] for r in both_results],
            [r["enrich_crispri"] for r in both_results],
        )
        rho_both_s, _ = stats.spearmanr(
            [r["proj_shrna"] for r in both_results],
            [r["enrich_shrna"] for r in both_results],
        )
        log(f"\nMatched comparison (same {len(both_results)} drugs, both ground truths):")
        log(f"  CRISPRi projected rho: {rho_both_c:.4f}")
        log(f"  shRNA projected rho:   {rho_both_s:.4f}")
        log(f"  CRISPRi >= shRNA: {'YES' if rho_both_c >= rho_both_s else 'NO'}")

    log(f"\n=== CRITERION (three-outcome convergent validity check) ===")
    log(f"  NOTE: This is a convergent-validity check, not a blind confirmatory test.")
    log(f"  The comparison value (rho_shRNA = 0.376) is known from the original H3.")
    if rho_proj_c is not None:
        original_shrna_rho = 0.376
        if rho_proj_c >= original_shrna_rho:
            log(f"  OUTCOME (a): rho_CRISPRi ({rho_proj_c:.4f}) >= rho_shRNA ({original_shrna_rho})")
            log(f"  Convergent validity confirmed. CRISPRi ground truth produces same or")
            log(f"  stronger signal. shRNA noise likely attenuated the original result.")
        elif rho_proj_c > 0:
            log(f"  OUTCOME (b): 0 < rho_CRISPRi ({rho_proj_c:.4f}) < rho_shRNA ({original_shrna_rho})")
            log(f"  Both ground truths produce positive correlations — robustness confirmed.")
            log(f"  Attenuation likely reflects K562-only confound (CRISPRi from single cell type")
            log(f"  vs shRNA averaged across ~20 cell lines).")
        else:
            log(f"  OUTCOME (c): rho_CRISPRi ({rho_proj_c:.4f}) <= 0")
            log(f"  Interpret via K562-only confound first: CRISPRi signatures from a single")
            log(f"  hematopoietic cell line may not capture on-target direction for drugs acting")
            log(f"  in non-hematopoietic contexts.")

    log(f"\nResults saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--synthetic", action="store_true")
    mode.add_argument("--real", action="store_true")
    parser.add_argument("--data", type=Path,
                        help="Path to drug-perturbation-geometry/data/")
    parser.add_argument("--perturbseq", type=Path,
                        help="Path to ReplogleWeissman2022_K562_essential.h5ad")
    parser.add_argument("--output", type=Path,
                        default=Path("results/03b_h3_crispri"))
    args = parser.parse_args()

    if args.synthetic:
        run_synthetic()
    elif args.real:
        if not args.data:
            parser.error("--real requires --data")
        if not args.perturbseq:
            parser.error("--real requires --perturbseq")
        run_real(args.data, args.perturbseq, args.output)
