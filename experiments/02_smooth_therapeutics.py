"""Experiment 02: Successful therapeutics can have low raw bracket.

Tests H2: known broad-mechanism drugs (statins, metformin, aspirin,
dexamethasone, sirolimus) have raw direction instability below the population
median. Low bracket does not imply irrelevance.

This demonstrates failure mode 2: "useful interventions can be smooth."

Usage:
    uv run python experiments/02_smooth_therapeutics.py --data PATH_TO_LINCS
    uv run python experiments/02_smooth_therapeutics.py --synthetic
"""
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from geometry.bracket_norm import direction_instability, magnitude_cv


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


BROAD_MECHANISM_DRUGS = {
    "atorvastatin": "HMG-CoA reductase inhibitor — universally expressed target",
    "simvastatin": "HMG-CoA reductase inhibitor — universally expressed target",
    "lovastatin": "HMG-CoA reductase inhibitor — universally expressed target",
    "metformin": "AMPK activator / mitochondrial complex I — fundamental metabolism",
    "dexamethasone": "Glucocorticoid receptor agonist — ubiquitous receptor",
    "sirolimus": "mTOR inhibitor — universal growth pathway",
    "methotrexate": "DHFR inhibitor — fundamental folate metabolism",
    "aspirin": "COX inhibitor — ubiquitous prostaglandin pathway",
    "vorinostat": "Pan-HDAC inhibitor — universal chromatin modifier",
    "trichostatin-a": "Pan-HDAC inhibitor — universal chromatin modifier",
}

CONTEXT_SPECIFIC_DRUGS = {
    "vemurafenib": "BRAF V600E — only mutant melanoma lines",
    "imatinib": "BCR-ABL — CML-specific fusion protein",
    "erlotinib": "EGFR — lung/epithelial cancers",
    "crizotinib": "ALK/ROS1 — rare lung cancer subtypes",
    "lapatinib": "HER2 — amplified breast cancer only",
    "gefitinib": "EGFR — mutant-dependent efficacy",
    "PCI-34051": "HDAC8-selective — tissue-restricted expression",
    "tubastatin-a": "HDAC6-selective — cell-type dependent effects",
}


def run_synthetic():
    """Validate: smooth therapeutics have lower bracket than specific drugs."""
    log("=== SYNTHETIC MODE ===")
    rng = np.random.default_rng()
    n_genes = 978
    n_contexts = 10

    broad_brackets = []
    specific_brackets = []

    for _ in range(50):
        sigs = rng.standard_normal((n_contexts, n_genes)) * 0.2
        sigs += rng.standard_normal((1, n_genes)) * 1.5
        broad_brackets.append(direction_instability(sigs))

    for _ in range(50):
        sigs = rng.standard_normal((n_contexts, n_genes)) * 0.5
        active_contexts = rng.choice(n_contexts, size=n_contexts // 2, replace=False)
        sigs[active_contexts] += rng.standard_normal((len(active_contexts), n_genes)) * 2.0
        specific_brackets.append(direction_instability(sigs))

    log(f"Broad mechanism: {np.mean(broad_brackets):.4f} +/- {np.std(broad_brackets):.4f}")
    log(f"Context-specific: {np.mean(specific_brackets):.4f} +/- {np.std(specific_brackets):.4f}")

    u, p = stats.mannwhitneyu(broad_brackets, specific_brackets, alternative="less")
    log(f"Mann-Whitney U (broad < specific): p={p:.2e}")

    median_all = np.median(broad_brackets + specific_brackets)
    below = sum(1 for b in broad_brackets if b < median_all)
    log(f"Broad below median: {below}/50 ({100*below/50:.0f}%)")

    if below >= 30:
        log("PASS: broad-mechanism drugs have lower bracket (synthetic)")
    else:
        log("FAIL: separation not achieved in synthetic data")
        sys.exit(1)


def run_real(data_dir: Path, output_dir: Path):
    """Run H2 on real LINCS data."""
    log("=== REAL MODE ===")

    sigs_path = data_dir / "lincs_subset.npz"
    siginfo_path = data_dir / "GSE92742_Broad_LINCS_sig_info.txt.gz"

    if not sigs_path.exists():
        log(f"ERROR: {sigs_path} not found")
        sys.exit(1)

    log("Loading signatures...")
    data = np.load(sigs_path, allow_pickle=True)
    signatures = data["signatures"]
    sig_ids = data["sig_ids"]
    log(f"  {signatures.shape[0]:,} x {signatures.shape[1]}")

    import pandas as pd
    log("Loading siginfo...")
    siginfo = pd.read_csv(siginfo_path, sep="\t", low_memory=False)

    sig_id_to_idx = {sid: i for i, sid in enumerate(sig_ids)}
    drug_sig_map = {}
    for _, row in siginfo.iterrows():
        sid = row.get("sig_id")
        drug = row.get("pert_iname", "")
        cell = row.get("cell_id", "")
        if sid in sig_id_to_idx and drug:
            if drug not in drug_sig_map:
                drug_sig_map[drug] = {}
            if cell not in drug_sig_map[drug]:
                drug_sig_map[drug][cell] = []
            drug_sig_map[drug][cell].append(sig_id_to_idx[sid])

    log("Computing direction instability for all drugs with 5+ cell lines...")
    all_brackets = {}
    for drug_name, cell_sigs in tqdm(drug_sig_map.items(), desc="Drugs"):
        if len(cell_sigs) < 5:
            continue
        cell_means = []
        for cell, indices in cell_sigs.items():
            cell_means.append(signatures[indices].mean(axis=0))
        sigs_matrix = np.array(cell_means)
        all_brackets[drug_name] = {
            "raw_bracket": float(direction_instability(sigs_matrix)),
            "magnitude_cv": float(magnitude_cv(sigs_matrix)),
            "n_celllines": len(cell_sigs),
        }

    all_values = [v["raw_bracket"] for v in all_brackets.values()]
    population_median = np.median(all_values)
    log(f"\nPopulation: {len(all_brackets)} drugs, median bracket = {population_median:.4f}")

    log(f"\n=== H2: BROAD-MECHANISM DRUGS ===")
    broad_found = []
    for drug, reason in BROAD_MECHANISM_DRUGS.items():
        if drug in all_brackets:
            b = all_brackets[drug]
            below = b["raw_bracket"] < population_median
            broad_found.append({"drug": drug, "reason": reason, **b, "below_median": below})
            marker = "BELOW" if below else "ABOVE"
            log(f"  {drug:20s} bracket={b['raw_bracket']:.4f} [{marker}] — {reason}")
        else:
            log(f"  {drug:20s} NOT FOUND in 5+ cell-line set")

    log(f"\n=== CONTEXT-SPECIFIC DRUGS (expected HIGH bracket) ===")
    specific_found = []
    for drug, reason in CONTEXT_SPECIFIC_DRUGS.items():
        if drug in all_brackets:
            b = all_brackets[drug]
            above = b["raw_bracket"] > population_median
            specific_found.append({"drug": drug, "reason": reason, **b, "above_median": above})
            marker = "ABOVE" if above else "BELOW"
            log(f"  {drug:20s} bracket={b['raw_bracket']:.4f} [{marker}] — {reason}")

    n_broad_below = sum(1 for d in broad_found if d["below_median"])
    n_specific_above = sum(1 for d in specific_found if d["above_median"])

    log(f"\n=== H2 SUMMARY ===")
    log(f"Broad drugs below median: {n_broad_below}/{len(broad_found)}")
    log(f"Specific drugs above median: {n_specific_above}/{len(specific_found)}")

    if len(broad_found) >= 5:
        fraction = n_broad_below / len(broad_found)
        if fraction >= 0.6:
            log(f"H2 CONFIRMED: {fraction:.0%} of broad-mechanism drugs below median")
        else:
            log(f"H2 NOT CONFIRMED: only {fraction:.0%} below median (needed >=60%)")

    output_dir.mkdir(parents=True, exist_ok=True)

    def make_serializable(obj):
        if isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        return obj

    results = {
        "population_median": float(population_median),
        "n_drugs_total": len(all_brackets),
        "broad_mechanism_drugs": [
            {k: make_serializable(v) for k, v in d.items()} for d in broad_found
        ],
        "context_specific_drugs": [
            {k: make_serializable(v) for k, v in d.items()} for d in specific_found
        ],
        "h2_broad_below_median": n_broad_below,
        "h2_broad_total": len(broad_found),
        "h2_specific_above_median": n_specific_above,
        "h2_specific_total": len(specific_found),
    }
    with open(output_dir / "smooth_therapeutics_results.json", "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--synthetic", action="store_true")
    mode.add_argument("--real", action="store_true")
    parser.add_argument("--data", type=Path)
    parser.add_argument("--output", type=Path, default=Path("results/02_smooth_therapeutics"))
    args = parser.parse_args()

    if args.synthetic:
        run_synthetic()
    elif args.real:
        if not args.data:
            parser.error("--real requires --data")
        run_real(args.data, args.output)
