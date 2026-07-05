"""Experiment 09: Same-tissue direction instability in LINCS L1000.

Tests whether restricting to cell lines from the same tissue causes the
lineage-correction to bite (rho < 0.70), unlike the all-tissue toxicity
correction (rho = 0.83).

Hypothesis: within a single tissue, lineage-shared programs create false
consistency. Removing high-variance "lineage genes" should reorganize the
ranking more than removing stress genes did across all tissues.

Usage:
    uv run python experiments/09_same_tissue_lincs.py
"""
import json
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

DRUG_GEOM_DIR = Path("/Users/elliottower/Documents/GitHub/drug-perturbation-geometry")
DATA_DIR = DRUG_GEOM_DIR / "data"
OUTPUT_DIR = Path(__file__).parent.parent / "results" / "09_same_tissue"

TISSUES = {
    "large_intestine": ["CL34", "HCT116", "HT115", "HT29", "LOVO", "MDST8",
                        "NCIH508", "NCIH716", "RKO", "SNU1040", "SNUC5",
                        "SW480", "SW620", "SW948"],
    "lung": ["A549", "CORL23", "DV90", "H1299", "HCC15", "HCC515",
             "NCIH1694", "NCIH1836", "NCIH2073", "NCIH596", "SKLU1", "T3M10"],
    "heme": ["HL60", "JURKAT", "NOMO1", "PL21", "SKM1", "THP1", "U937", "WSUDLCL2"],
    "breast": ["BT20", "HS578T", "MCF10A", "MCF7", "MDAMB231", "SKBR3"],
    "ovary": ["COV644", "EFO27", "OV7", "RMGI", "RMUGS", "TYKNU"],
    "skin": ["A375", "FIBRNPC", "MCH58", "SKMEL1", "SKMEL28"],
}


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def direction_instability(signatures: np.ndarray) -> float:
    norms = np.linalg.norm(signatures, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    unit = signatures / norms
    K = unit.shape[0]
    if K < 2:
        return np.nan
    cos_matrix = unit @ unit.T
    mask = np.triu(np.ones((K, K), dtype=bool), k=1)
    return 1.0 - float(cos_matrix[mask].mean())


def main():
    log("=== EXPERIMENT 09: SAME-TISSUE DIRECTION INSTABILITY ===")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    log("Loading LINCS subset...")
    data = np.load(DATA_DIR / "lincs_subset.npz", allow_pickle=True)
    signatures = data["signatures"]
    sig_ids = data["sig_ids"]
    gene_ids = list(data["gene_ids"])
    log(f"  {signatures.shape[0]:,} signatures x {signatures.shape[1]} genes")

    log("Loading siginfo for drug/cell mapping...")
    siginfo = pd.read_csv(DATA_DIR / "GSE92742_Broad_LINCS_sig_info.txt.gz",
                          sep="\t", low_memory=False)
    sig_to_info = siginfo.set_index("sig_id")[["pert_iname", "cell_id"]].to_dict("index")

    log("Building drug x cell-line signature matrix...")
    drug_cell_sigs = defaultdict(dict)
    for i, sid in enumerate(tqdm(sig_ids, desc="Mapping")):
        info = sig_to_info.get(sid)
        if info is None:
            continue
        drug = info["pert_iname"]
        cell = info["cell_id"]
        if cell not in drug_cell_sigs[drug]:
            drug_cell_sigs[drug][cell] = []
        drug_cell_sigs[drug][cell].append(i)

    log(f"  {len(drug_cell_sigs)} drugs mapped")

    log("Aggregating per drug x cell (mean signature)...")
    drug_cell_means = {}
    for drug, cells in tqdm(drug_cell_sigs.items(), desc="Aggregating"):
        drug_cell_means[drug] = {}
        for cell, indices in cells.items():
            drug_cell_means[drug][cell] = signatures[indices].mean(axis=0)

    all_results = {}

    for tissue_name, tissue_cells in TISSUES.items():
        log(f"\n--- Tissue: {tissue_name} ({len(tissue_cells)} cell lines) ---")
        tissue_set = set(tissue_cells)

        drugs_in_tissue = []
        for drug, cells in drug_cell_means.items():
            present = [c for c in cells if c in tissue_set]
            if len(present) >= 5:
                drugs_in_tissue.append(drug)

        log(f"  Drugs tested in >= 5 {tissue_name} cell lines: {len(drugs_in_tissue)}")

        if len(drugs_in_tissue) < 50:
            log(f"  SKIP: too few drugs (need >= 50)")
            all_results[tissue_name] = {
                "n_drugs": len(drugs_in_tissue),
                "skipped": True,
                "reason": "fewer than 50 drugs with K>=5 in this tissue",
            }
            continue

        log("  Computing within-tissue direction instability...")
        raw_instabilities = {}
        for drug in tqdm(drugs_in_tissue, desc=f"  {tissue_name} raw"):
            cells = drug_cell_means[drug]
            tissue_sigs = np.array([cells[c] for c in tissue_cells if c in cells])
            if tissue_sigs.shape[0] >= 5:
                raw_instabilities[drug] = direction_instability(tissue_sigs)

        log(f"  Computed raw D for {len(raw_instabilities)} drugs")

        log("  Identifying lineage genes (top 100 by within-tissue variance)...")
        all_tissue_sigs = []
        for drug in drugs_in_tissue[:500]:
            cells = drug_cell_means[drug]
            for c in tissue_cells:
                if c in cells:
                    all_tissue_sigs.append(cells[c])
        all_tissue_sigs = np.array(all_tissue_sigs)
        gene_variances = all_tissue_sigs.var(axis=0)
        lineage_gene_idx = np.argsort(gene_variances)[-100:]
        non_lineage_idx = np.array([i for i in range(978) if i not in set(lineage_gene_idx)])
        log(f"  Lineage genes: top 100 by variance (var range: "
            f"{gene_variances[lineage_gene_idx].min():.3f} - "
            f"{gene_variances[lineage_gene_idx].max():.3f})")

        log("  Computing lineage-corrected direction instability...")
        corrected_instabilities = {}
        for drug in tqdm(drugs_in_tissue, desc=f"  {tissue_name} corrected"):
            cells = drug_cell_means[drug]
            tissue_sigs = np.array([cells[c] for c in tissue_cells if c in cells])
            if tissue_sigs.shape[0] >= 5:
                corrected_instabilities[drug] = direction_instability(
                    tissue_sigs[:, non_lineage_idx])

        shared_drugs = sorted(set(raw_instabilities) & set(corrected_instabilities))
        raw_vals = np.array([raw_instabilities[d] for d in shared_drugs])
        corr_vals = np.array([corrected_instabilities[d] for d in shared_drugs])

        rho, p = stats.spearmanr(raw_vals, corr_vals)
        log(f"  Spearman rho (raw vs lineage-corrected): {rho:.4f}, p={p:.2e}")
        log(f"  Raw D: mean={raw_vals.mean():.4f}, median={np.median(raw_vals):.4f}")
        log(f"  Corrected D: mean={corr_vals.mean():.4f}, median={np.median(corr_vals):.4f}")

        if rho < 0.70:
            verdict = "BITES"
        elif rho < 0.83:
            verdict = "MARGINAL"
        else:
            verdict = "DOESN'T BITE"
        log(f"  Verdict: {verdict}")

        raw_ranks = stats.rankdata(raw_vals)
        corr_ranks = stats.rankdata(corr_vals)
        rank_changes = np.abs(corr_ranks - raw_ranks)
        top_movers_idx = np.argsort(rank_changes)[-10:][::-1]

        log(f"  Top 10 rank movers:")
        top_movers = []
        for idx in top_movers_idx:
            drug = shared_drugs[idx]
            log(f"    {drug:30s} raw_rank={int(raw_ranks[idx]):5d} -> "
                f"corr_rank={int(corr_ranks[idx]):5d} (delta={int(rank_changes[idx]):+5d})")
            top_movers.append({
                "drug": drug,
                "raw_rank": int(raw_ranks[idx]),
                "corrected_rank": int(corr_ranks[idx]),
                "rank_change": int(rank_changes[idx]),
            })

        all_results[tissue_name] = {
            "n_cell_lines": len(tissue_cells),
            "n_drugs": len(shared_drugs),
            "raw_vs_corrected_rho": float(rho),
            "raw_vs_corrected_p": float(p),
            "verdict": verdict,
            "raw_mean": float(raw_vals.mean()),
            "raw_median": float(np.median(raw_vals)),
            "corrected_mean": float(corr_vals.mean()),
            "corrected_median": float(np.median(corr_vals)),
            "n_lineage_genes_removed": 100,
            "top_movers": top_movers,
        }

    log("\n\n=== SUMMARY ===")
    log(f"{'Tissue':<20s} {'N drugs':>8s} {'rho':>6s} {'Verdict':<15s}")
    log("-" * 55)
    for tissue, res in all_results.items():
        if res.get("skipped"):
            log(f"{tissue:<20s} {'N/A':>8s} {'N/A':>6s} SKIPPED ({res['reason'][:30]})")
        else:
            log(f"{tissue:<20s} {res['n_drugs']:>8d} {res['raw_vs_corrected_rho']:>6.3f} "
                f"{res['verdict']:<15s}")

    log(f"\nCf. all-tissue toxicity correction: rho = 0.83")
    log(f"Cf. Perturb-seq essential correction: rho = 0.91")
    log(f"Cf. JUMP-CP cell-health correction: rho = 0.987")

    with open(OUTPUT_DIR / "same_tissue_results.json", "w") as f:
        json.dump(all_results, f, indent=2)
    log(f"\nResults saved to {OUTPUT_DIR / 'same_tissue_results.json'}")

    log("\n=== DONE ===")


if __name__ == "__main__":
    main()
