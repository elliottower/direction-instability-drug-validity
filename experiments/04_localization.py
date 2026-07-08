"""Experiment 04: Localized bracket predicts mechanism better than global bracket.

Tests H4: for drugs targeting specific pathways, the instability concentrates
in pathway-relevant genes. The localization score (pathway bracket / global
bracket) predicts MOA annotation better than raw global bracket.

Decision criterion (pre-registered): AUROC for MOA prediction from
localization score > AUROC from raw bracket by at least 0.05.

Gene-set method (declared in DEVIATION_LOG.md entry 1):
  Primary: region_mask = top N=100 genes by absolute z-score in the consensus
  shRNA knockdown signature of the drug's annotated target gene.
  Sensitivity: Drug Repurposing Hub MOA grouping.

Usage:
    uv run python experiments/04_localization.py --synthetic
    uv run python experiments/04_localization.py --real --data PATH_TO_LINCS
"""
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from geometry.bracket_norm import direction_instability, localization_score


REGION_MASK_SIZE = 100


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def build_shrna_consensus(shrna_sigs, shrna_siginfo):
    """Build consensus shRNA signature per target gene."""
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


def make_region_mask(shrna_sig, n_genes, top_n=REGION_MASK_SIZE):
    """Top-N genes by absolute z-score as the pathway region mask."""
    abs_z = np.abs(shrna_sig)
    top_idx = np.argsort(abs_z)[-top_n:]
    mask = np.zeros(n_genes, dtype=bool)
    mask[top_idx] = True
    return mask


def run_synthetic():
    """Validate: localization score predicts MOA better than raw bracket."""
    log("=== SYNTHETIC MODE ===")
    rng = np.random.default_rng()
    n_genes = 200
    n_contexts = 8
    n_moa_classes = 5
    n_drugs_per_moa = 30
    pathway_size = 40

    drug_results = []
    drug_labels = []

    for moa in range(n_moa_classes):
        pathway_start = moa * pathway_size
        for _ in range(n_drugs_per_moa):
            sigs = rng.standard_normal((n_contexts, n_genes)) * 0.3
            strength = rng.uniform(1.0, 4.0)
            for ctx in range(n_contexts):
                sigs[ctx, pathway_start:pathway_start + pathway_size] += (
                    rng.standard_normal(pathway_size) * strength
                )

            region_mask = np.zeros(n_genes, dtype=bool)
            region_mask[pathway_start:pathway_start + pathway_size] = True

            raw = direction_instability(sigs)
            loc = localization_score(sigs, region_mask)

            drug_results.append({"raw": raw, "localization": loc, "moa": moa})
            drug_labels.append(moa)

    log(f"Generated {len(drug_results)} drugs across {n_moa_classes} MOA classes")

    raw_aurocs = []
    loc_aurocs = []
    for moa in range(n_moa_classes):
        y_true = [1 if d["moa"] == moa else 0 for d in drug_results]
        raw_scores = [d["raw"] for d in drug_results]
        loc_scores = [d["localization"] for d in drug_results]

        raw_auc = roc_auc_score(y_true, raw_scores)
        loc_auc = roc_auc_score(y_true, loc_scores)
        raw_aurocs.append(raw_auc)
        loc_aurocs.append(loc_auc)

    mean_raw = np.mean(raw_aurocs)
    mean_loc = np.mean(loc_aurocs)
    gap = mean_loc - mean_raw

    log(f"\nMean AUROC for MOA prediction:")
    log(f"  Raw bracket:        {mean_raw:.4f}")
    log(f"  Localization score:  {mean_loc:.4f}")
    log(f"  Gap:                 {gap:+.4f}")

    if gap > 0.05:
        log("  H4 PASS: localization AUROC > raw AUROC by > 0.05")
    else:
        log(f"  H4 FAIL: gap = {gap:.4f} (needed > 0.05)")


def run_real(data_dir: Path, output_dir: Path):
    """Run H4 on real LINCS data using shRNA-derived region masks."""
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
    n_genes = compound_sigs.shape[1]
    log(f"  {compound_sigs.shape[0]:,} signatures x {n_genes} genes")

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

    drug_meta = {}
    for d in drugs_list:
        target = d.get("target")
        moa = d.get("moa")
        if target and moa:
            primary_target = target.split("|")[0].strip()
            drug_meta[d["pert_iname"]] = {"target": primary_target, "moa": moa}

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

    log("Computing localization scores...")
    results = []
    skipped = {"no_meta": 0, "no_shrna": 0, "few_cells": 0}

    for drug, cells in tqdm(drug_cell_map.items(), desc="Drugs"):
        if drug not in drug_meta:
            skipped["no_meta"] += 1
            continue
        meta = drug_meta[drug]
        if meta["target"] not in shrna_consensus:
            skipped["no_shrna"] += 1
            continue
        if len(cells) < 5:
            skipped["few_cells"] += 1
            continue

        cell_means = []
        for cell, indices in cells.items():
            cell_means.append(compound_sigs[indices].mean(axis=0))
        sigs_matrix = np.array(cell_means)

        region_mask = make_region_mask(shrna_consensus[meta["target"]], n_genes)
        raw = direction_instability(sigs_matrix)
        loc = localization_score(sigs_matrix, region_mask)

        results.append({
            "drug": drug,
            "target": meta["target"],
            "moa": meta["moa"],
            "n_celllines": len(cells),
            "raw_bracket": float(raw),
            "localization_score": float(loc),
        })

    log(f"\nComputed for {len(results)} drugs")
    log(f"Skipped: {skipped}")

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "localization_results.json", "w") as f:
        json.dump(results, f, indent=2)

    from collections import Counter
    moa_counts = Counter(r["moa"] for r in results)
    eligible_moas = [moa for moa, count in moa_counts.items() if count >= 10]
    log(f"\nMOA classes with >= 10 drugs: {len(eligible_moas)}")

    if len(eligible_moas) < 2:
        log("ERROR: too few MOA classes for AUROC computation")
        sys.exit(1)

    eligible_drugs = [r for r in results if r["moa"] in eligible_moas]

    raw_aurocs = []
    loc_aurocs = []
    for moa in eligible_moas:
        y_true = [1 if d["moa"] == moa else 0 for d in eligible_drugs]
        raw_scores = [d["raw_bracket"] for d in eligible_drugs]
        loc_scores = [d["localization_score"] for d in eligible_drugs]

        if sum(y_true) < 2 or sum(y_true) == len(y_true):
            continue

        raw_auc = roc_auc_score(y_true, raw_scores)
        loc_auc = roc_auc_score(y_true, loc_scores)
        raw_aurocs.append(raw_auc)
        loc_aurocs.append(loc_auc)
        log(f"  {moa:40s}: raw={raw_auc:.3f}, loc={loc_auc:.3f}, gap={loc_auc - raw_auc:+.3f}")

    mean_raw = np.mean(raw_aurocs)
    mean_loc = np.mean(loc_aurocs)
    gap = mean_loc - mean_raw

    log(f"\n=== H4 RESULTS (PRIMARY: shRNA-derived region masks, N={REGION_MASK_SIZE}) ===")
    log(f"N drugs: {len(eligible_drugs)}, N MOA classes: {len(raw_aurocs)}")
    log(f"Mean AUROC for MOA prediction:")
    log(f"  Raw bracket:        {mean_raw:.4f}")
    log(f"  Localization score:  {mean_loc:.4f}")
    log(f"  Gap:                 {gap:+.4f}")

    if gap > 0.05:
        log(f"  H4 CONFIRMED: localization AUROC > raw AUROC by {gap:.4f} (> 0.05)")
    else:
        log(f"  H4 NOT CONFIRMED: gap = {gap:.4f} (needed > 0.05)")

    log(f"\nResults saved to {output_dir / 'localization_results.json'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--synthetic", action="store_true")
    mode.add_argument("--real", action="store_true")
    parser.add_argument("--data", type=Path,
                        help="Path to drug-perturbation-geometry/data/")
    parser.add_argument("--output", type=Path,
                        default=Path("results/04_localization"))
    args = parser.parse_args()

    if args.synthetic:
        run_synthetic()
    elif args.real:
        if not args.data:
            parser.error("--real requires --data")
        run_real(args.data, args.output)
