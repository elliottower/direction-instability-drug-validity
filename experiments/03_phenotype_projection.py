"""Experiment 03: Phenotype-projected bracket separates mechanism from noise.

Tests H3: for drugs with known target pathways, phenotype-projected bracket
(instability projected onto the target pathway direction) correlates with
on-target activity, while raw bracket does not.

Decision criterion (pre-registered): Spearman |rho| > 0.3 between
phenotype-projected bracket and on-target gene enrichment, AND |rho| < 0.15
for raw bracket vs same enrichment.

Gene-set method (declared in DEVIATION_LOG.md entry 1):
  Primary: genetic perturbation connectivity — phenotype_direction is the
  consensus shRNA knockdown signature of the drug's annotated target gene.
  Sensitivity: Drug Repurposing Hub MOA grouping.

Usage:
    uv run python experiments/03_phenotype_projection.py --synthetic
    uv run python experiments/03_phenotype_projection.py --real --data PATH_TO_LINCS
"""
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from geometry.bracket_norm import direction_instability, phenotype_projected_bracket


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def build_shrna_consensus(shrna_sigs, shrna_siginfo):
    """Build consensus shRNA signature per target gene (mean across hairpins and cells)."""
    sig_id_to_idx = {sid: i for i, sid in enumerate(shrna_siginfo["sig_id"])}
    gene_sigs = {}
    for _, row in shrna_siginfo.iterrows():
        gene = row["pert_iname"]
        sid = row["sig_id"]
        if sid in sig_id_to_idx:
            if gene not in gene_sigs:
                gene_sigs[gene] = []
            gene_sigs[gene].append(sig_id_to_idx[sid])

    consensus = {}
    for gene, indices in gene_sigs.items():
        if len(indices) < 3:
            continue
        consensus[gene] = shrna_sigs[indices].mean(axis=0)
    return consensus


def compute_on_target_enrichment(drug_sigs_matrix, shrna_direction):
    """Fraction of drug signature variance along the on-target direction."""
    mean_sig = drug_sigs_matrix.mean(axis=0)
    total_var = float(np.var(mean_sig))
    if total_var < 1e-10:
        return 0.0
    direction = shrna_direction / (np.linalg.norm(shrna_direction) + 1e-10)
    proj = float(mean_sig @ direction)
    return proj**2 / (total_var * len(mean_sig))


def run_synthetic():
    """Validate: phenotype-projected bracket correlates with on-target enrichment."""
    log("=== SYNTHETIC MODE ===")
    rng = np.random.default_rng()
    n_genes = 200
    n_contexts = 8
    n_drugs = 100
    target_size = 30

    results = []
    for i in tqdm(range(n_drugs), desc="Drugs"):
        target_start = rng.integers(0, n_genes - target_size)
        on_target_strength = rng.uniform(0.5, 5.0)

        phenotype_dir = np.zeros(n_genes)
        phenotype_dir[target_start:target_start + target_size] = rng.standard_normal(target_size)
        phenotype_dir /= np.linalg.norm(phenotype_dir) + 1e-10

        sigs = rng.standard_normal((n_contexts, n_genes)) * 0.5
        for ctx in range(n_contexts):
            on_target_noise = rng.standard_normal(target_size) * on_target_strength
            sigs[ctx, target_start:target_start + target_size] += on_target_noise
            sigs[ctx] += rng.standard_normal(n_genes) * rng.uniform(0, 2)

        raw = direction_instability(sigs)
        projected = phenotype_projected_bracket(sigs, phenotype_dir)
        enrichment = on_target_strength

        results.append({
            "raw_bracket": raw,
            "projected_bracket": projected,
            "on_target_enrichment": enrichment,
        })

    raw_vals = [r["raw_bracket"] for r in results]
    proj_vals = [r["projected_bracket"] for r in results]
    enrich_vals = [r["on_target_enrichment"] for r in results]

    raw_rho, raw_p = stats.spearmanr(raw_vals, enrich_vals)
    proj_rho, proj_p = stats.spearmanr(proj_vals, enrich_vals)

    log(f"\nSpearman correlations with on-target enrichment:")
    log(f"  Raw bracket:       rho={raw_rho:.4f} (p={raw_p:.2e})")
    log(f"  Projected bracket: rho={proj_rho:.4f} (p={proj_p:.2e})")

    pass_proj = abs(proj_rho) > 0.3
    pass_raw = abs(raw_rho) < 0.15

    log(f"\nH3 criteria:")
    log(f"  |projected rho| > 0.3: {abs(proj_rho):.4f} -> {'PASS' if pass_proj else 'FAIL'}")
    log(f"  |raw rho| < 0.15:      {abs(raw_rho):.4f} -> {'PASS' if pass_raw else 'FAIL'}")

    if pass_proj and pass_raw:
        log("  H3 PASS: projected bracket correlates with on-target, raw does not")
    else:
        log("  H3 FAIL")


def run_real(data_dir: Path, output_dir: Path):
    """Run H3 on real LINCS data using shRNA genetic perturbation connectivity."""
    log("=== REAL MODE ===")

    sigs_path = data_dir / "lincs_subset.npz"
    siginfo_path = data_dir / "GSE92742_Broad_LINCS_sig_info.txt.gz"
    shrna_path = data_dir / "lincs_shrna.npz"
    shrna_siginfo_path = data_dir / "lincs_shrna_siginfo.csv.gz"
    labels_path = data_dir / "frozen_drug_labels.json"

    for p in [sigs_path, shrna_path, shrna_siginfo_path, labels_path]:
        if not p.exists():
            log(f"ERROR: {p} not found.")
            sys.exit(1)

    log("Loading compound signatures...")
    data = np.load(sigs_path, allow_pickle=True)
    compound_sigs = data["signatures"]
    compound_sig_ids = list(data["sig_ids"])
    gene_ids = list(data["gene_ids"])
    log(f"  {compound_sigs.shape[0]:,} signatures x {compound_sigs.shape[1]} genes")

    log("Loading shRNA signatures...")
    shrna_data = np.load(shrna_path, allow_pickle=True)
    shrna_sigs = shrna_data["signatures"]
    log(f"  {shrna_sigs.shape[0]:,} shRNA signatures")

    log("Loading shRNA sig info...")
    shrna_siginfo = pd.read_csv(shrna_siginfo_path)

    log("Building shRNA consensus signatures...")
    shrna_consensus = build_shrna_consensus(shrna_sigs, shrna_siginfo)
    log(f"  {len(shrna_consensus)} target genes with >= 3 hairpins")

    log("Loading drug labels...")
    with open(labels_path) as f:
        labels = json.load(f)
    drugs_list = labels["drugs"]
    drug_targets = {}
    for d in drugs_list:
        target = d.get("target")
        if target:
            primary_target = target.split("|")[0].strip()
            drug_targets[d["pert_iname"]] = primary_target

    log("Loading compound sig info...")
    siginfo = pd.read_csv(siginfo_path, sep="\t", low_memory=False)

    log("Building drug-cell matrix...")
    sig_id_to_idx = {sid: i for i, sid in enumerate(compound_sig_ids)}
    drug_cell_map = {}
    for _, row in siginfo.iterrows():
        sid = row.get("sig_id")
        drug = row.get("pert_iname", "")
        cell = row.get("cell_id", "")
        if sid in sig_id_to_idx and drug:
            if drug not in drug_cell_map:
                drug_cell_map[drug] = {}
            if cell not in drug_cell_map[drug]:
                drug_cell_map[drug][cell] = []
            drug_cell_map[drug][cell].append(sig_id_to_idx[sid])

    log("Computing phenotype-projected brackets...")
    results = []
    skipped_no_target = 0
    skipped_no_shrna = 0
    skipped_few_cells = 0

    for drug, cells in tqdm(drug_cell_map.items(), desc="Drugs"):
        if drug not in drug_targets:
            skipped_no_target += 1
            continue
        target = drug_targets[drug]
        if target not in shrna_consensus:
            skipped_no_shrna += 1
            continue
        if len(cells) < 5:
            skipped_few_cells += 1
            continue

        cell_means = []
        for cell, indices in cells.items():
            cell_means.append(compound_sigs[indices].mean(axis=0))
        sigs_matrix = np.array(cell_means)

        phenotype_dir = shrna_consensus[target]
        raw = direction_instability(sigs_matrix)
        projected = phenotype_projected_bracket(sigs_matrix, phenotype_dir)
        enrichment = compute_on_target_enrichment(sigs_matrix, phenotype_dir)

        results.append({
            "drug": drug,
            "target": target,
            "n_celllines": len(cells),
            "raw_bracket": float(raw),
            "projected_bracket": float(projected),
            "on_target_enrichment": float(enrichment),
        })

    log(f"\nComputed for {len(results)} drugs")
    log(f"Skipped: {skipped_no_target} no target, {skipped_no_shrna} no shRNA, {skipped_few_cells} < 5 cells")

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "phenotype_projection_results.json", "w") as f:
        json.dump(results, f, indent=2)

    if len(results) < 10:
        log("ERROR: too few drugs to compute correlations")
        sys.exit(1)

    raw_vals = [r["raw_bracket"] for r in results]
    proj_vals = [r["projected_bracket"] for r in results]
    enrich_vals = [r["on_target_enrichment"] for r in results]

    raw_rho, raw_p = stats.spearmanr(raw_vals, enrich_vals)
    proj_rho, proj_p = stats.spearmanr(proj_vals, enrich_vals)

    log(f"\n=== H3 RESULTS (PRIMARY: genetic perturbation connectivity) ===")
    log(f"N drugs: {len(results)}")
    log(f"Spearman correlations with on-target enrichment:")
    log(f"  Raw bracket:       rho={raw_rho:.4f} (p={raw_p:.2e})")
    log(f"  Projected bracket: rho={proj_rho:.4f} (p={proj_p:.2e})")

    pass_proj = abs(proj_rho) > 0.3
    pass_raw = abs(raw_rho) < 0.15

    log(f"\nH3 criteria:")
    log(f"  |projected rho| > 0.3: {abs(proj_rho):.4f} -> {'PASS' if pass_proj else 'FAIL'}")
    log(f"  |raw rho| < 0.15:      {abs(raw_rho):.4f} -> {'PASS' if pass_raw else 'FAIL'}")

    if pass_proj and pass_raw:
        log("  H3 CONFIRMED: projected bracket correlates with on-target, raw does not")
    else:
        log("  H3 NOT CONFIRMED")
        if not pass_proj:
            log(f"    Projected rho too low: {abs(proj_rho):.4f} (needed > 0.3)")
        if not pass_raw:
            log(f"    Raw rho too high: {abs(raw_rho):.4f} (needed < 0.15)")

    log(f"\nResults saved to {output_dir / 'phenotype_projection_results.json'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--synthetic", action="store_true")
    mode.add_argument("--real", action="store_true")
    parser.add_argument("--data", type=Path,
                        help="Path to drug-perturbation-geometry/data/")
    parser.add_argument("--output", type=Path,
                        default=Path("results/03_phenotype_projection"))
    args = parser.parse_args()

    if args.synthetic:
        run_synthetic()
    elif args.real:
        if not args.data:
            parser.error("--real requires --data")
        run_real(args.data, args.output)
