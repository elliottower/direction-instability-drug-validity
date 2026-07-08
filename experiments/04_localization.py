"""Experiment 04: Localized bracket predicts mechanism better than global bracket.

Tests H4: for drugs targeting specific pathways, the instability concentrates
in pathway-relevant genes. The localization score (pathway bracket / global
bracket) predicts MOA annotation better than raw global bracket.

Decision criterion (pre-registered): AUROC for MOA prediction from
localization score > AUROC from raw bracket by at least 0.05.

AUROC is macro-averaged: one-vs-rest AUROC per MOA class, then unweighted mean.
Bootstrap 95% CI on the gap is reported.

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


def make_region_mask(shrna_sig, n_genes, top_n=REGION_MASK_SIZE):
    """Top-N genes by absolute z-score as the pathway region mask."""
    abs_z = np.abs(shrna_sig)
    top_idx = np.argsort(abs_z)[-top_n:]
    mask = np.zeros(n_genes, dtype=bool)
    mask[top_idx] = True
    return mask


def compute_macro_auroc(drugs, eligible_moas):
    """Macro-averaged one-vs-rest AUROC across MOA classes.

    Returns (mean_raw_auroc, mean_loc_auroc, per_moa_raw, per_moa_loc).
    """
    raw_aurocs = []
    loc_aurocs = []
    for moa in eligible_moas:
        y_true = [1 if d["moa"] == moa else 0 for d in drugs]
        if sum(y_true) < 2 or sum(y_true) == len(y_true):
            continue
        raw_scores = [d["raw_bracket"] for d in drugs]
        loc_scores = [d["localization_score"] for d in drugs]
        raw_aurocs.append(roc_auc_score(y_true, raw_scores))
        loc_aurocs.append(roc_auc_score(y_true, loc_scores))
    return np.mean(raw_aurocs), np.mean(loc_aurocs), raw_aurocs, loc_aurocs


def bootstrap_auroc_gap(drugs, eligible_moas, n_boot=2000, seed=42):
    """Bootstrap 95% CI on (macro_loc_auroc - macro_raw_auroc)."""
    rng = np.random.default_rng(seed)
    gaps = []
    n = len(drugs)
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        boot_drugs = [drugs[i] for i in idx]
        raw, loc, _, _ = compute_macro_auroc(boot_drugs, eligible_moas)
        gaps.append(loc - raw)
    gaps = np.array(gaps)
    return float(np.percentile(gaps, 2.5)), float(np.percentile(gaps, 97.5))


def run_synthetic():
    """Validate: localization score predicts MOA better than raw bracket."""
    log("=== SYNTHETIC MODE ===")
    rng = np.random.default_rng(seed=2026070704)
    n_genes = 200
    n_contexts = 8
    n_moa_classes = 5
    n_drugs_per_moa = 30
    pathway_size = 40

    drug_results = []

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

            drug_results.append({"raw_bracket": raw, "localization_score": loc, "moa": moa})

    log(f"Generated {len(drug_results)} drugs across {n_moa_classes} MOA classes")

    moas = list(range(n_moa_classes))
    mean_raw, mean_loc, _, _ = compute_macro_auroc(drug_results, moas)
    gap = mean_loc - mean_raw

    log(f"\nMacro-averaged AUROC for MOA prediction:")
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
    compound_gene_ids = list(data["gene_ids"])
    n_genes = compound_sigs.shape[1]
    assert len(set(compound_sig_ids)) == len(compound_sig_ids), "Duplicate sig_ids in compound data"
    log(f"  {compound_sigs.shape[0]:,} signatures x {n_genes} genes")

    log("Loading shRNA signatures...")
    shrna_data = np.load(shrna_path, allow_pickle=True)
    shrna_sigs = shrna_data["signatures"]
    shrna_gene_ids = list(shrna_data["gene_ids"])
    log(f"  {shrna_sigs.shape[0]:,} shRNA signatures")

    assert compound_gene_ids == shrna_gene_ids, (
        f"Gene-axis mismatch: compound has {len(compound_gene_ids)} genes, "
        f"shRNA has {len(shrna_gene_ids)} genes"
    )
    log("  Gene-axis alignment verified: compound and shRNA use identical gene ordering")

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
    filtered = siginfo[siginfo["sig_id"].isin(sig_id_to_idx) & siginfo["pert_iname"].notna()].copy()
    filtered["_idx"] = filtered["sig_id"].map(sig_id_to_idx)

    drug_cell_map = {}
    for (drug, cell), group in filtered.groupby(["pert_iname", "cell_id"]):
        if drug not in drug_cell_map:
            drug_cell_map[drug] = {}
        indices = group["_idx"].values
        drug_cell_map[drug][cell] = compound_sigs[indices].mean(axis=0)

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

        sigs_matrix = np.array(list(cells.values()))
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

    mean_raw, mean_loc, per_moa_raw, per_moa_loc = compute_macro_auroc(eligible_drugs, eligible_moas)
    gap = mean_loc - mean_raw

    for moa, raw_auc, loc_auc in zip(eligible_moas, per_moa_raw, per_moa_loc):
        log(f"  {moa:40s}: raw={raw_auc:.3f}, loc={loc_auc:.3f}, gap={loc_auc - raw_auc:+.3f}")

    log(f"\nBootstrapping 95% CI on AUROC gap (2000 resamples)...")
    ci_lo, ci_hi = bootstrap_auroc_gap(eligible_drugs, eligible_moas)

    log(f"\n=== H4 RESULTS (PRIMARY: shRNA-derived region masks, N={REGION_MASK_SIZE}) ===")
    log(f"N drugs: {len(eligible_drugs)}, N MOA classes: {len(per_moa_raw)}")
    log(f"Macro-averaged AUROC for MOA prediction:")
    log(f"  Raw bracket:        {mean_raw:.4f}")
    log(f"  Localization score:  {mean_loc:.4f}")
    log(f"  Gap:                 {gap:+.4f}  95% CI [{ci_lo:+.4f}, {ci_hi:+.4f}]")

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
