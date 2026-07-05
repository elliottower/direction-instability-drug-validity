"""Experiment 01: Toxic perturbations have high raw bracket but fail validity filtering.

Tests H1: known cytotoxic compounds rank highly by raw direction instability
but drop after toxicity correction. This demonstrates failure mode 1:
"violence is not control."

Uses LINCS L1000 data already extracted for the drug transport paper.
Curates a toxicity gene set from known stress/apoptosis markers, then
shows raw vs corrected bracket for known toxic vs therapeutic drugs.

Usage:
    uv run python experiments/01_toxicity_failure.py --data PATH_TO_LINCS
    uv run python experiments/01_toxicity_failure.py --synthetic
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
from geometry.bracket_norm import (
    direction_instability,
    toxicity_corrected_bracket,
)


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


STRESS_APOPTOSIS_GENES = [
    "HSPA1A", "HSPA1B", "HSP90AA1", "HSPA5", "HSPB1", "DNAJB1",
    "BAX", "BCL2", "CASP3", "CASP8", "CASP9", "CYCS", "APAF1",
    "TP53", "MDM2", "CDKN1A", "GADD45A", "GADD45B",
    "ATF4", "DDIT3", "XBP1", "ATF6", "EIF2AK3",
    "NFKB1", "RELA", "TNF", "IL6", "IL1B",
    "FOS", "JUN", "EGR1", "MYC",
    "HMOX1", "NQO1", "GCLM", "TXNRD1",
    "CCND1", "CCNE1", "CDK4", "CDK2", "RB1",
]

KNOWN_CYTOTOXIC = [
    "staurosporine", "camptothecin", "doxorubicin", "etoposide",
    "cisplatin", "vincristine", "paclitaxel", "actinomycin-d",
    "thapsigargin", "tunicamycin", "brefeldin-a", "rotenone",
    "antimycin-a", "oligomycin-a", "menadione", "hydrogen-peroxide",
]

KNOWN_BROAD_THERAPEUTICS = [
    "vorinostat", "trichostatin-a", "panobinostat",
    "sirolimus", "dexamethasone", "methotrexate",
    "atorvastatin", "simvastatin", "lovastatin",
    "metformin",
]


def run_synthetic():
    """Validate pipeline on synthetic data with known toxicity structure.

    The planted structure: toxic drugs have a CONSISTENT stress signature
    across contexts (low raw instability = looks like transport). But when
    we remove stress genes, they're exposed — nothing consistent remains.
    Therapeutic drugs have consistent signal in NON-stress genes, so
    correction doesn't hurt them.
    """
    log("=== SYNTHETIC MODE ===")
    rng = np.random.default_rng()
    n_genes = 200
    n_contexts = 8
    n_toxic = 50
    n_therapeutic = 50
    n_stress = 30

    stress_idx = np.arange(n_stress)

    results_toxic = []
    results_therapeutic = []

    log("Generating synthetic toxic drugs (consistent stress, noise elsewhere)...")
    for _ in tqdm(range(n_toxic), desc="Toxic"):
        shared_stress = rng.standard_normal(n_stress) * 3.0
        sigs = rng.standard_normal((n_contexts, n_genes)) * 0.5
        for ctx in range(n_contexts):
            sigs[ctx, :n_stress] = shared_stress + rng.standard_normal(n_stress) * 0.3
        raw = direction_instability(sigs)
        corrected = toxicity_corrected_bracket(sigs, stress_idx)
        results_toxic.append({"raw": raw, "corrected": corrected})

    log("Generating synthetic therapeutic drugs (consistent mechanism in non-stress genes)...")
    for _ in tqdm(range(n_therapeutic), desc="Therapeutic"):
        target_start = n_stress + rng.integers(0, n_genes - n_stress - 20)
        shared_mech = rng.standard_normal(20) * 3.0
        sigs = rng.standard_normal((n_contexts, n_genes)) * 0.5
        for ctx in range(n_contexts):
            sigs[ctx, target_start:target_start + 20] = shared_mech + rng.standard_normal(20) * 0.3
        raw = direction_instability(sigs)
        corrected = toxicity_corrected_bracket(sigs, stress_idx)
        results_therapeutic.append({"raw": raw, "corrected": corrected})

    raw_toxic = [r["raw"] for r in results_toxic]
    raw_ther = [r["raw"] for r in results_therapeutic]
    corr_toxic = [r["corrected"] for r in results_toxic]
    corr_ther = [r["corrected"] for r in results_therapeutic]

    log(f"\nRaw instability:  toxic={np.mean(raw_toxic):.4f}, therapeutic={np.mean(raw_ther):.4f}")
    log(f"Corrected:        toxic={np.mean(corr_toxic):.4f}, therapeutic={np.mean(corr_ther):.4f}")
    log(f"(Low = stable/transporting. High = unstable/context-dependent.)")

    log(f"\nKey insight: BOTH have low raw bracket (both look like they transport).")
    log(f"But after removing stress genes:")
    log(f"  Toxic: instability RISES (stress was the only consistent signal)")
    log(f"  Therapeutic: instability stays LOW (mechanism is in non-stress genes)")

    corr_increase_toxic = np.mean(corr_toxic) - np.mean(raw_toxic)
    corr_increase_ther = np.mean(corr_ther) - np.mean(raw_ther)
    log(f"\nCorrection effect (instability change):")
    log(f"  Toxic: {corr_increase_toxic:+.4f} (should increase)")
    log(f"  Therapeutic: {corr_increase_ther:+.4f} (should stay same or decrease)")

    u_corr, p_corr = stats.mannwhitneyu(corr_toxic, corr_ther, alternative="greater")
    log(f"\nCorrected bracket: toxic > therapeutic p={p_corr:.2e}")

    toxic_exposed = sum(1 for t, c in zip(raw_toxic, corr_toxic) if c > t + 0.05)
    ther_preserved = sum(1 for t, c in zip(raw_ther, corr_ther) if c < t + 0.05)

    log(f"\nH1 check:")
    log(f"  Toxic drugs exposed by correction (instability rises): "
        f"{toxic_exposed}/{n_toxic} ({100*toxic_exposed/n_toxic:.0f}%)")
    log(f"  Therapeutic drugs preserved (instability stable): "
        f"{ther_preserved}/{n_therapeutic} ({100*ther_preserved/n_therapeutic:.0f}%)")

    if toxic_exposed > n_toxic * 0.5 and ther_preserved > n_therapeutic * 0.5:
        log("  PASS: toxicity correction exposes false-positive transporters")
    else:
        log("  FAIL: correction does not separate toxic from therapeutic")
        sys.exit(1)


def run_real(data_dir: Path, output_dir: Path):
    """Run H1 on real LINCS data."""
    log("=== REAL MODE ===")

    sigs_path = data_dir / "lincs_subset.npz"
    labels_path = data_dir / "frozen_drug_labels.json"
    siginfo_path = data_dir / "GSE92742_Broad_LINCS_sig_info.txt.gz"
    geneinfo_path = data_dir / "GSE92742_Broad_LINCS_gene_info.txt.gz"

    if not sigs_path.exists():
        log(f"ERROR: {sigs_path} not found. Run drug-perturbation-geometry extraction first.")
        sys.exit(1)

    log("Loading signatures...")
    data = np.load(sigs_path, allow_pickle=True)
    signatures = data["signatures"]
    sig_ids = data["sig_ids"]
    gene_ids = list(data["gene_ids"])
    log(f"  {signatures.shape[0]:,} signatures x {signatures.shape[1]} genes")

    log("Loading drug labels...")
    import pandas as pd
    siginfo = pd.read_csv(siginfo_path, sep="\t", low_memory=False)
    with open(labels_path) as f:
        drug_labels = json.load(f)

    log("Loading gene info for stress gene mapping...")
    geneinfo = pd.read_csv(geneinfo_path, sep="\t", low_memory=False)
    gene_symbols = geneinfo.set_index("pr_gene_id")["pr_gene_symbol"].to_dict()

    gene_id_to_symbol = {}
    for gid in gene_ids:
        gid_int = int(gid) if gid.isdigit() else gid
        if gid_int in gene_symbols:
            gene_id_to_symbol[gid] = gene_symbols[gid_int]

    stress_indices = []
    for i, gid in enumerate(gene_ids):
        symbol = gene_id_to_symbol.get(gid, "")
        if symbol in STRESS_APOPTOSIS_GENES:
            stress_indices.append(i)
    stress_indices = np.array(stress_indices)
    log(f"  Mapped {len(stress_indices)} stress/apoptosis genes")

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

    log("Computing bracket norms...")
    results = []
    for drug_name, cell_sigs in tqdm(drug_sig_map.items(), desc="Drugs"):
        if len(cell_sigs) < 5:
            continue
        cell_means = []
        for cell, indices in cell_sigs.items():
            cell_means.append(signatures[indices].mean(axis=0))
        sigs_matrix = np.array(cell_means)

        raw = direction_instability(sigs_matrix)
        corrected = toxicity_corrected_bracket(sigs_matrix, stress_indices)

        is_toxic = drug_name.lower() in [d.lower() for d in KNOWN_CYTOTOXIC]
        is_therapeutic = drug_name.lower() in [d.lower() for d in KNOWN_BROAD_THERAPEUTICS]

        results.append({
            "drug": drug_name,
            "n_celllines": len(cell_sigs),
            "raw_bracket": float(raw),
            "corrected_bracket": float(corrected),
            "is_toxic": is_toxic,
            "is_therapeutic": is_therapeutic,
            "moa": drug_labels.get(drug_name, {}).get("moa", "unknown"),
        })

    log(f"Computed for {len(results)} drugs")

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "toxicity_results.json", "w") as f:
        json.dump(results, f, indent=2)

    toxic_drugs = [r for r in results if r["is_toxic"]]
    therapeutic_drugs = [r for r in results if r["is_therapeutic"]]
    all_raw = [r["raw_bracket"] for r in results]
    all_corr = [r["corrected_bracket"] for r in results]

    log(f"\n=== H1 RESULTS ===")
    log(f"Known cytotoxic found: {len(toxic_drugs)}")
    log(f"Known broad therapeutic found: {len(therapeutic_drugs)}")

    if toxic_drugs:
        p80 = np.percentile(all_raw, 80)
        median_corr = np.median(all_corr)

        toxic_in_top20 = sum(1 for d in toxic_drugs if d["raw_bracket"] > p80)
        toxic_below_median = sum(1 for d in toxic_drugs if d["corrected_bracket"] < median_corr)

        log(f"\nToxic drugs in top 20% raw bracket: {toxic_in_top20}/{len(toxic_drugs)}")
        log(f"Toxic drugs below median corrected: {toxic_below_median}/{len(toxic_drugs)}")
        raw_strs = [f"{d['drug']}={d['raw_bracket']:.3f}" for d in toxic_drugs]
        corr_strs = [f"{d['drug']}={d['corrected_bracket']:.3f}" for d in toxic_drugs]
        log(f"  Raw brackets: {raw_strs}")
        log(f"  Corrected:    {corr_strs}")

        fraction_pass = toxic_below_median / len(toxic_drugs) if toxic_drugs else 0
        if fraction_pass >= 0.6:
            log(f"\n  H1 CONFIRMED: {fraction_pass:.0%} of toxic drugs drop below median after correction")
        else:
            log(f"\n  H1 NOT CONFIRMED: only {fraction_pass:.0%} drop (needed >=60%)")

    if therapeutic_drugs:
        log(f"\nTherapeutic drugs raw brackets:")
        for d in sorted(therapeutic_drugs, key=lambda x: x["raw_bracket"]):
            log(f"  {d['drug']:20s} raw={d['raw_bracket']:.4f} corrected={d['corrected_bracket']:.4f}")

    log(f"\nResults saved to {output_dir / 'toxicity_results.json'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--synthetic", action="store_true")
    mode.add_argument("--real", action="store_true")
    parser.add_argument("--data", type=Path,
                        help="Path to drug-perturbation-geometry/data/")
    parser.add_argument("--output", type=Path,
                        default=Path("results/01_toxicity_failure"))
    args = parser.parse_args()

    if args.synthetic:
        run_synthetic()
    elif args.real:
        if not args.data:
            parser.error("--real requires --data")
        run_real(args.data, args.output)
