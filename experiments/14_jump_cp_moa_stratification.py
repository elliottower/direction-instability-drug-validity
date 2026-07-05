"""Experiment 14: MOA class stratification in JUMP-CP morphological profiles.

Pre-registered in PREREGISTRATION_COMBINED_PAPER.md.

Tests whether the LINCS finding (machinery-targeting drugs have lower direction
instability than receptor-mediated drugs) replicates in JUMP-CP Cell Painting
morphological profiles.

Compares three MOA groups:
  - Fundamental machinery (HDAC, topoisomerase, CDK, tubulin, proteasome, etc.)
  - Kinase inhibitors (EGFR, MEK, MAPK, SRC, PI3K)
  - Receptor-mediated (dopamine, serotonin, histamine, adrenergic, etc.)

Data: Direction instability from results/08_jump_cp/jump_cp_instabilities.npz
      MOA labels from drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv
      InChIKey matching via data/jump_cp/compound_metadata.csv.gz

Usage:
    uv run python experiments/14_jump_cp_moa_stratification.py
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

DRUG_REPO = Path(__file__).resolve().parent.parent.parent / "drug-perturbation-geometry"
INSTABILITY_CSV = DRUG_REPO / "zenodo_v1" / "drug_instability_8949.csv"
LINCS_PERTINFO = DRUG_REPO / "data" / "GSE92742_Broad_LINCS_pert_info.txt.gz"
REPURPOSING_HUB = DRUG_REPO / "data" / "repurposing_hub_drugs.txt"

JUMP_DI_PATH = Path("results/08_jump_cp/jump_cp_instabilities.npz")
JUMP_META_PATH = Path("data/jump_cp/compound_metadata.csv.gz")

OUTPUT_DIR = Path("results/14_jump_cp_moa")

MACHINERY_MOAS = [
    "HDAC inhibitor",
    "topoisomerase inhibitor",
    "CDK inhibitor",
    "tubulin inhibitor",
    "proteasome inhibitor",
    "mTOR inhibitor",
    "DNA synthesis inhibitor",
    "RNA synthesis inhibitor",
    "protein synthesis inhibitor",
    "ATPase inhibitor",
    "DNA alkylating agent",
    "ribonucleotide reductase inhibitor",
]

RECEPTOR_MOAS = [
    "dopamine receptor agonist",
    "dopamine receptor antagonist",
    "serotonin receptor agonist",
    "serotonin receptor antagonist",
    "histamine receptor antagonist",
    "adrenergic receptor agonist",
    "adrenergic receptor antagonist",
    "acetylcholine receptor agonist",
    "acetylcholine receptor antagonist",
    "glutamate receptor antagonist",
    "glutamate receptor agonist",
    "opioid receptor agonist",
    "opioid receptor antagonist",
    "calcium channel blocker",
    "sodium channel blocker",
    "potassium channel blocker",
    "GABA receptor agonist",
    "GABA receptor antagonist",
    "cannabinoid receptor agonist",
    "cannabinoid receptor antagonist",
    "benzodiazepine receptor agonist",
]

KINASE_MOAS = [
    "EGFR inhibitor",
    "MEK inhibitor",
    "MAPK inhibitor",
    "SRC inhibitor",
    "PI3K inhibitor",
    "JAK inhibitor",
    "FLT3 inhibitor",
    "VEGFR inhibitor",
    "PDGFR inhibitor",
    "ABL inhibitor",
    "RAF inhibitor",
    "Aurora kinase inhibitor",
    "PKC inhibitor",
]


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def classify_moa(moa_str: str) -> str | None:
    """Classify a MOA string into machinery/kinase/receptor/other."""
    if pd.isna(moa_str) or moa_str == "":
        return None
    moa_lower = moa_str.lower()
    for m in MACHINERY_MOAS:
        if m.lower() in moa_lower:
            return "machinery"
    for m in KINASE_MOAS:
        if m.lower() in moa_lower:
            return "kinase"
    for m in RECEPTOR_MOAS:
        if m.lower() in moa_lower:
            return "receptor"
    return "other"


def cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """Cohen's d (pooled SD)."""
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_sd = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_sd < 1e-10:
        return 0.0
    return (np.mean(group1) - np.mean(group2)) / pooled_sd


def main():
    log("=== EXPERIMENT 14: MOA STRATIFICATION IN JUMP-CP ===")
    log("Pre-registered in PREREGISTRATION_COMBINED_PAPER.md")

    if not JUMP_DI_PATH.exists():
        log(f"ERROR: JUMP-CP DI not found at {JUMP_DI_PATH}")
        log("  Run experiments/08_jump_cp_instability.py first")
        sys.exit(1)
    if not JUMP_META_PATH.exists():
        log(f"ERROR: JUMP-CP metadata not found at {JUMP_META_PATH}")
        sys.exit(1)
    if not INSTABILITY_CSV.exists():
        log(f"ERROR: LINCS instability CSV not found at {INSTABILITY_CSV}")
        sys.exit(1)

    log("Loading JUMP-CP direction instabilities...")
    jump_data = np.load(JUMP_DI_PATH, allow_pickle=True)
    jump_compounds = jump_data["compounds"]
    jump_raw_di = jump_data["raw"]
    log(f"  {len(jump_compounds)} compounds with DI values")

    log("Loading JUMP-CP compound metadata (InChIKey mapping)...")
    jump_meta = pd.read_csv(JUMP_META_PATH)
    log(f"  {len(jump_meta)} rows in compound metadata")
    log(f"  Columns: {list(jump_meta.columns)}")

    jump_di_df = pd.DataFrame({
        "Metadata_JCP2022": jump_compounds,
        "jump_di": jump_raw_di,
    })
    jump_merged = jump_di_df.merge(
        jump_meta[["Metadata_JCP2022", "Metadata_InChIKey"]].drop_duplicates(),
        on="Metadata_JCP2022", how="inner"
    )
    log(f"  Matched to InChIKey: {len(jump_merged)}")

    jump_merged["inchikey_prefix"] = jump_merged["Metadata_InChIKey"].str[:14]

    log("Loading LINCS drug annotations...")
    lincs_df = pd.read_csv(INSTABILITY_CSV)
    log(f"  {len(lincs_df)} LINCS drugs")

    lincs_with_moa = lincs_df[lincs_df["moa"].notna()].copy()
    log(f"  With MOA annotation: {len(lincs_with_moa)}")

    log("Loading LINCS pert_info for InChIKey mapping...")
    if LINCS_PERTINFO.exists():
        pert_info = pd.read_csv(LINCS_PERTINFO, sep="\t", low_memory=False)
        log(f"  pert_info: {pert_info.shape}")
        if "inchi_key" in pert_info.columns:
            pert_inchi = pert_info[["pert_iname", "inchi_key"]].dropna().drop_duplicates()
            pert_inchi["inchikey_prefix"] = pert_inchi["inchi_key"].str[:14]
            lincs_with_moa = lincs_with_moa.merge(
                pert_inchi.rename(columns={"pert_iname": "drug_name"}),
                on="drug_name", how="left"
            )
            log(f"  LINCS drugs with InChIKey: {lincs_with_moa['inchikey_prefix'].notna().sum()}")
    elif REPURPOSING_HUB.exists():
        log("  Using Repurposing Hub for InChIKey...")
        hub = pd.read_csv(REPURPOSING_HUB, sep="\t", low_memory=False)
        if "InChIKey" in hub.columns:
            hub["inchikey_prefix"] = hub["InChIKey"].str[:14]
            lincs_with_moa = lincs_with_moa.merge(
                hub[["pert_iname", "inchikey_prefix"]].rename(
                    columns={"pert_iname": "drug_name"}
                ).drop_duplicates(),
                on="drug_name", how="left"
            )

    lincs_with_inchi = lincs_with_moa[lincs_with_moa["inchikey_prefix"].notna()].copy()
    log(f"  LINCS drugs with InChIKey for matching: {len(lincs_with_inchi)}")

    matched = jump_merged.merge(
        lincs_with_inchi[["drug_name", "moa", "direction_instability", "inchikey_prefix"]],
        on="inchikey_prefix", how="inner"
    ).drop_duplicates(subset=["drug_name"])
    log(f"  Matched LINCS→JUMP-CP drugs: {len(matched)}")

    matched["moa_class"] = matched["moa"].apply(classify_moa)
    class_counts = matched["moa_class"].value_counts()
    log(f"\n  MOA class distribution:")
    for cls, n in class_counts.items():
        log(f"    {cls}: {n}")

    machinery = matched[matched["moa_class"] == "machinery"]["jump_di"].values
    receptor = matched[matched["moa_class"] == "receptor"]["jump_di"].values
    kinase = matched[matched["moa_class"] == "kinase"]["jump_di"].values

    log(f"\n--- Group statistics ---")
    for name, arr in [("machinery", machinery), ("kinase", kinase), ("receptor", receptor)]:
        if len(arr) > 0:
            log(f"  {name:12s}: n={len(arr):3d}, mean={arr.mean():.4f}, std={arr.std():.4f}, median={np.median(arr):.4f}")
        else:
            log(f"  {name:12s}: n=0")

    log(f"\n--- Decision criteria ---")

    if len(machinery) < 3 or len(receptor) < 3:
        log("  ERROR: Insufficient sample size for comparison")
        log(f"  machinery n={len(machinery)}, receptor n={len(receptor)}")
        d_mr = np.nan
        u_stat, u_p = np.nan, np.nan
        c1 = False
        c2 = False
        c3 = False
    else:
        d_mr = cohens_d(receptor, machinery)
        u_stat, u_p = stats.mannwhitneyu(machinery, receptor, alternative='two-sided')

        c1 = d_mr > 0.3 and machinery.mean() < receptor.mean()
        log(f"  Criterion 1 (d > 0.3, machinery lower): d={d_mr:.4f} → {'PASS' if c1 else 'FAIL'}")

        if len(kinase) > 0:
            ordering_correct = machinery.mean() < kinase.mean() < receptor.mean()
            c2 = ordering_correct
            log(f"  Criterion 2 (ordering machinery < kinase < receptor):")
            log(f"    machinery={machinery.mean():.4f} < kinase={kinase.mean():.4f} < receptor={receptor.mean():.4f} → {'PASS' if c2 else 'FAIL'}")
        else:
            c2 = False
            log(f"  Criterion 2: No kinase drugs matched — cannot test ordering")

        c3 = u_p < 0.05
        log(f"  Criterion 3 (Mann-Whitney p < 0.05): U={u_stat:.1f}, p={u_p:.4f} → {'PASS' if c3 else 'FAIL'}")

    log(f"\n--- LINCS reference (for comparison) ---")
    lincs_machinery = lincs_with_moa[lincs_with_moa["moa"].apply(classify_moa) == "machinery"]["direction_instability"].values
    lincs_receptor = lincs_with_moa[lincs_with_moa["moa"].apply(classify_moa) == "receptor"]["direction_instability"].values
    lincs_kinase = lincs_with_moa[lincs_with_moa["moa"].apply(classify_moa) == "kinase"]["direction_instability"].values
    for name, arr in [("machinery", lincs_machinery), ("kinase", lincs_kinase), ("receptor", lincs_receptor)]:
        if len(arr) > 0:
            log(f"  LINCS {name:12s}: n={len(arr):3d}, mean={arr.mean():.4f}")
    if len(lincs_machinery) > 0 and len(lincs_receptor) > 0:
        lincs_d = cohens_d(lincs_receptor, lincs_machinery)
        log(f"  LINCS Cohen's d (receptor vs machinery): {lincs_d:.4f}")

    log(f"\n--- Summary ---")
    if c1 and c2 and c3:
        log("  RESULT: Cross-modal replication confirmed. Mechanism conservation")
        log("  is a property of the drug-target interaction, not the measurement.")
    elif c1 and not c2:
        log("  RESULT: Machinery vs receptor separation replicated but ordering")
        log("  with kinase inhibitors differs from LINCS.")
    elif not c1 and d_mr > 0:
        log("  RESULT: Correct direction but small effect. Morphological profiles")
        log("  show attenuated version of LINCS pattern.")
    elif not c1 and d_mr <= 0:
        log("  RESULT: No cross-modal replication or reversed ordering.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {
        "experiment": "14_jump_cp_moa_stratification",
        "preregistered": "PREREGISTRATION_COMBINED_PAPER.md",
        "n_matched": len(matched),
        "group_sizes": {
            "machinery": len(machinery),
            "kinase": len(kinase),
            "receptor": len(receptor),
            "other": int((matched["moa_class"] == "other").sum()),
            "unclassified": int(matched["moa_class"].isna().sum()),
        },
        "jump_cp": {
            "machinery_mean": float(machinery.mean()) if len(machinery) > 0 else None,
            "machinery_std": float(machinery.std()) if len(machinery) > 0 else None,
            "kinase_mean": float(kinase.mean()) if len(kinase) > 0 else None,
            "kinase_std": float(kinase.std()) if len(kinase) > 0 else None,
            "receptor_mean": float(receptor.mean()) if len(receptor) > 0 else None,
            "receptor_std": float(receptor.std()) if len(receptor) > 0 else None,
            "cohens_d_receptor_vs_machinery": float(d_mr) if not np.isnan(d_mr) else None,
            "mann_whitney_u": float(u_stat) if not np.isnan(u_stat) else None,
            "mann_whitney_p": float(u_p) if not np.isnan(u_p) else None,
        },
        "criteria": {
            "c1_d_gt_03_and_correct_direction": bool(c1),
            "c2_ordering_correct": bool(c2),
            "c3_significant": bool(c3),
        },
    }

    matched_out = matched[["drug_name", "moa", "moa_class", "jump_di",
                            "direction_instability", "Metadata_JCP2022"]].copy()
    matched_out = matched_out.rename(columns={"direction_instability": "lincs_di"})
    matched_out.to_csv(OUTPUT_DIR / "matched_drugs.csv", index=False)

    with open(OUTPUT_DIR / "moa_stratification_results.json", "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
