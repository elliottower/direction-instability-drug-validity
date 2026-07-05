"""Experiment 13f v2: MOA classification AUROC using the broad principled definition.

Uses the frozen classification spec (commit 31e5e7a) which defines:
- Machinery = core genomic/epigenetic/proteostatic/cytoskeletal targets
- Receptor = cell-surface receptors and ion channels
- Everything else excluded automatically by the rule

This is a POST-HOC EXPLORATORY analysis (not pre-registered).

Includes nuclear-receptor sensitivity analysis: primary excludes them,
robustness checks fold them into machinery and into receptor. The result
is considered robust only if DI > Frechet holds in all three variants.

Usage:
    uv run python experiments/13f_moa_classification_auroc_v2.py
"""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

DRUG_REPO = Path(__file__).resolve().parent.parent.parent / "drug-perturbation-geometry"
INSTABILITY_CSV = DRUG_REPO / "zenodo_v1" / "drug_instability_8949.csv"
OUTPUT_DIR = Path("results/13_lincs_benchmark")

N_BOOTSTRAP = 10_000
SPEC_COMMIT = "57867de"

MACHINERY_KEYWORDS = [
    "HDAC inhibitor", "topoisomerase inhibitor", "CDK inhibitor",
    "proteasome inhibitor", "mTOR inhibitor", "DNA synthesis inhibitor",
    "RNA synthesis inhibitor", "protein synthesis inhibitor",
    "ATPase inhibitor", "DNA alkylating agent", "ribonucleotide reductase inhibitor",
    "tubulin inhibitor", "tubulin polymerization inhibitor",
    "PARP inhibitor", "thymidylate synthase inhibitor",
    "dihydrofolate reductase inhibitor",
    "histone methyltransferase inhibitor", "histone demethylase inhibitor",
    "DNA methyltransferase inhibitor", "BET inhibitor", "sirtuin inhibitor",
    "cell cycle inhibitor", "MDM inhibitor",
    "autophagy inhibitor", "HSP inhibitor",
    "NFkB pathway inhibitor",
]

KINASE_KEYWORDS = [
    "EGFR inhibitor", "MEK inhibitor", "MAPK inhibitor", "p38 MAPK inhibitor",
    "SRC inhibitor", "PI3K inhibitor", "JAK inhibitor", "FLT3 inhibitor",
    "VEGFR inhibitor", "PDGFR inhibitor", "ABL inhibitor", "RAF inhibitor",
    "PKC inhibitor", "AKT inhibitor", "glycogen synthase kinase inhibitor",
    "IKK inhibitor", "rho associated kinase inhibitor", "tyrosine kinase inhibitor",
    "checkpoint kinase inhibitor", "ALK inhibitor", "BTK inhibitor",
    "KIT inhibitor", "Aurora kinase inhibitor",
]

RECEPTOR_KEYWORDS = [
    "dopamine receptor", "serotonin receptor", "histamine receptor",
    "adrenergic receptor", "acetylcholine receptor", "glutamate receptor",
    "opioid receptor", "cannabinoid receptor", "adenosine receptor",
    "angiotensin receptor", "prostanoid receptor", "tachykinin antagonist",
    "leukotriene receptor", "endothelin receptor", "benzodiazepine receptor",
    "calcium channel blocker", "T-type calcium channel blocker",
    "sodium channel blocker", "potassium channel blocker",
    "potassium channel activator", "GABA receptor",
]

NUCLEAR_RECEPTOR_KEYWORDS = [
    "glucocorticoid receptor", "estrogen receptor", "androgen receptor",
    "progesterone receptor", "PPAR receptor", "retinoid receptor",
    "vitamin D receptor", "mineralocorticoid receptor",
]


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def classify_broad(moa_str: str) -> str:
    if pd.isna(moa_str) or moa_str == "":
        return "no_annotation"
    if is_nuclear_receptor(moa_str):
        return "excluded"
    moa_lower = moa_str.lower()
    for kw in MACHINERY_KEYWORDS:
        if kw.lower() in moa_lower:
            return "machinery"
    for kw in KINASE_KEYWORDS:
        if kw.lower() in moa_lower:
            return "kinase"
    for kw in RECEPTOR_KEYWORDS:
        if kw.lower() in moa_lower:
            return "receptor"
    return "excluded"


def is_nuclear_receptor(moa_str: str) -> bool:
    if pd.isna(moa_str) or moa_str == "":
        return False
    moa_lower = moa_str.lower()
    return any(kw.lower() in moa_lower for kw in NUCLEAR_RECEPTOR_KEYWORDS)


def bootstrap_binary_clf(scores: np.ndarray, labels: np.ndarray,
                         n_bootstrap: int = N_BOOTSTRAP, seed: int = 42) -> dict:
    point_auroc = roc_auc_score(labels, scores)
    point_auprc = average_precision_score(labels, scores)

    rng = np.random.default_rng(seed)
    n = len(labels)
    boot_aurocs = np.zeros(n_bootstrap)
    boot_auprcs = np.zeros(n_bootstrap)
    valid = 0
    attempts = 0
    while valid < n_bootstrap and attempts < n_bootstrap * 3:
        idx = rng.integers(0, n, size=n)
        b_labels = labels[idx]
        if b_labels.sum() == 0 or b_labels.sum() == len(b_labels):
            attempts += 1
            continue
        boot_aurocs[valid] = roc_auc_score(b_labels, scores[idx])
        boot_auprcs[valid] = average_precision_score(b_labels, scores[idx])
        valid += 1
        attempts += 1

    boot_aurocs = boot_aurocs[:valid]
    boot_auprcs = boot_auprcs[:valid]
    return {
        "auroc": float(point_auroc),
        "auroc_ci_lo": float(np.percentile(boot_aurocs, 2.5)),
        "auroc_ci_hi": float(np.percentile(boot_aurocs, 97.5)),
        "auroc_se": float(np.std(boot_aurocs)),
        "auprc": float(point_auprc),
        "auprc_ci_lo": float(np.percentile(boot_auprcs, 2.5)),
        "auprc_ci_hi": float(np.percentile(boot_auprcs, 97.5)),
        "auprc_se": float(np.std(boot_auprcs)),
        "n_boot": valid,
    }


def bootstrap_clf_difference(scores_a: np.ndarray, scores_b: np.ndarray,
                             labels: np.ndarray,
                             n_bootstrap: int = N_BOOTSTRAP, seed: int = 42) -> dict:
    point_a = roc_auc_score(labels, scores_a)
    point_b = roc_auc_score(labels, scores_b)
    point_diff = point_a - point_b

    rng = np.random.default_rng(seed)
    n = len(labels)
    boot_diffs = np.zeros(n_bootstrap)
    valid = 0
    attempts = 0
    while valid < n_bootstrap and attempts < n_bootstrap * 3:
        idx = rng.integers(0, n, size=n)
        b_labels = labels[idx]
        if b_labels.sum() == 0 or b_labels.sum() == len(b_labels):
            attempts += 1
            continue
        auc_a = roc_auc_score(b_labels, scores_a[idx])
        auc_b = roc_auc_score(b_labels, scores_b[idx])
        boot_diffs[valid] = auc_a - auc_b
        valid += 1
        attempts += 1

    boot_diffs = boot_diffs[:valid]
    ci_lo = float(np.percentile(boot_diffs, 2.5))
    ci_hi = float(np.percentile(boot_diffs, 97.5))
    p_value = float(np.mean(boot_diffs <= 0)) if point_diff > 0 else float(np.mean(boot_diffs >= 0))

    return {
        "point_diff": float(point_diff),
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "p_value": p_value,
        "n_boot": valid,
        "auroc_a": float(point_a),
        "auroc_b": float(point_b),
    }


def bootstrap_spearman(x: np.ndarray, y: np.ndarray,
                       n_bootstrap: int = N_BOOTSTRAP, seed: int = 42) -> dict:
    point_r, point_p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    n = len(x)
    boot_rs = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        r, _ = stats.spearmanr(x[idx], y[idx])
        boot_rs[i] = r
    return {
        "rho": float(point_r),
        "p_value": float(point_p),
        "ci_lo": float(np.percentile(boot_rs, 2.5)),
        "ci_hi": float(np.percentile(boot_rs, 97.5)),
        "se": float(np.std(boot_rs)),
    }


def run_binary_analysis(df_subset: pd.DataFrame, label_col: str,
                        variant_name: str) -> dict:
    labels = (df_subset[label_col] == "machinery").astype(int).values
    n_mach = labels.sum()
    n_recv = len(labels) - n_mach

    log(f"\n  [{variant_name}] n_machinery={n_mach}, n_receptor={n_recv}, "
        f"base_rate={n_mach / len(labels):.3f}")

    predictors = {
        "direction_instability": -df_subset["direction_instability"].values,
        "frechet_variance": -df_subset["frechet_variance"].values,
        "magnitude_cv": -df_subset["magnitude_cv"].values,
    }

    result = {
        "n_machinery": int(n_mach),
        "n_receptor": int(n_recv),
        "base_rate_machinery": float(n_mach / len(labels)),
        "binary_classification": {},
        "paired_differences": {},
    }

    for name, scores in predictors.items():
        r = bootstrap_binary_clf(scores, labels, seed=42)
        result["binary_classification"][name] = r
        log(f"    {name:25s} AUROC={r['auroc']:.4f} [{r['auroc_ci_lo']:.4f}, {r['auroc_ci_hi']:.4f}]  "
            f"AUPRC={r['auprc']:.4f} [{r['auprc_ci_lo']:.4f}, {r['auprc_ci_hi']:.4f}]")

    di_scores = predictors["direction_instability"]
    for name, scores in predictors.items():
        if name == "direction_instability":
            continue
        r = bootstrap_clf_difference(di_scores, scores, labels, seed=42)
        sig = "*" if (r["ci_lo"] > 0 or r["ci_hi"] < 0) else "ns"
        result["paired_differences"][name] = r
        log(f"    DI - {name:25s} = {r['point_diff']:+.4f} "
            f"[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] p={r['p_value']:.4f} {sig}")

    return result


def main():
    log("=== EXPERIMENT 13f v2: MOA CLASSIFICATION (BROAD PRINCIPLED DEFINITION) ===")
    log(f"  Classification spec frozen at commit {SPEC_COMMIT}")
    log("  STATUS: POST-HOC EXPLORATORY (not pre-registered)")

    df = pd.read_csv(INSTABILITY_CSV)
    df["broad_class"] = df["moa"].apply(classify_broad)
    df["is_nuclear_receptor"] = df["moa"].apply(is_nuclear_receptor)

    counts = df["broad_class"].value_counts()
    log(f"\n  Classification counts:")
    for cls in ["machinery", "kinase", "receptor", "excluded", "no_annotation"]:
        log(f"    {cls}: {counts.get(cls, 0)}")

    # === ASSERTION BLOCK (from spec) ===
    annotated = df[df["broad_class"] != "no_annotation"]
    n_mach = (annotated["broad_class"] == "machinery").sum()
    n_kin = (annotated["broad_class"] == "kinase").sum()
    n_recv = (annotated["broad_class"] == "receptor").sum()
    n_excl = (annotated["broad_class"] == "excluded").sum()

    assert n_mach + n_kin + n_recv + n_excl == len(annotated), \
        f"Classes don't sum: {n_mach}+{n_kin}+{n_recv}+{n_excl} != {len(annotated)}"
    assert len(annotated) == df["moa"].notna().sum(), \
        "Annotated count mismatch"

    binary_df = annotated[annotated["broad_class"].isin(["machinery", "receptor"])].copy()
    assert len(binary_df) == n_mach + n_recv, \
        f"Binary set should be {n_mach}+{n_recv}, got {len(binary_df)}"
    assert not binary_df["broad_class"].isin(["kinase", "excluded"]).any(), \
        "Kinase or excluded drugs leaked into binary set"

    n_nhr_total = df["is_nuclear_receptor"].sum()
    n_nhr_excluded = ((annotated["broad_class"] == "excluded") & annotated["is_nuclear_receptor"]).sum()
    assert n_nhr_excluded == n_nhr_total, \
        f"NHR contamination: {n_nhr_total - n_nhr_excluded} NHR drugs not in excluded"
    nhr_in_binary = binary_df[binary_df["is_nuclear_receptor"]]
    assert len(nhr_in_binary) == 0, \
        f"NHR drugs in binary set: {nhr_in_binary['drug_name'].tolist()}"

    log(f"\n  NHR drugs: all {n_nhr_excluded} in excluded (deterministic, NHR checked first)")
    log(f"  Assertions passed: classes disjoint, zero NHR in binary set, "
        f"excluded removed, counts sum to {len(annotated)}")

    results = {
        "experiment": "13f_v2_moa_classification_broad",
        "spec_commit": SPEC_COMMIT,
        "status": "POST-HOC EXPLORATORY",
        "classification_principle": "core genomic/epigenetic/proteostatic/cytoskeletal "
                                    "vs cell-surface receptor/ion channel",
        "counts": {
            "machinery": int(n_mach),
            "kinase": int(n_kin),
            "receptor": int(n_recv),
            "excluded": int(n_excl),
            "no_annotation": int((df["broad_class"] == "no_annotation").sum()),
        },
    }

    # === PRIMARY: nuclear receptors excluded ===
    log("\n--- PRIMARY: nuclear receptors excluded ---")
    results["primary"] = run_binary_analysis(binary_df, "broad_class", "primary")

    # === SENSITIVITY 1: nuclear receptors folded into machinery ===
    log("\n--- SENSITIVITY: nuclear receptors -> machinery ---")
    nhr_excluded = annotated[
        (annotated["broad_class"] == "excluded") & annotated["is_nuclear_receptor"]
    ].copy()
    log(f"  Moving {len(nhr_excluded)} nuclear receptor drugs to machinery")

    nhr_as_mach = nhr_excluded.copy()
    nhr_as_mach["broad_class"] = "machinery"
    sens_mach_df = pd.concat([binary_df, nhr_as_mach])
    sens_mach_df = sens_mach_df[sens_mach_df["broad_class"].isin(["machinery", "receptor"])]
    results["sensitivity_nhr_machinery"] = run_binary_analysis(
        sens_mach_df, "broad_class", "NHR→machinery"
    )

    # === SENSITIVITY 2: nuclear receptors folded into receptor ===
    log("\n--- SENSITIVITY: nuclear receptors -> receptor ---")
    nhr_as_recv = nhr_excluded.copy()
    nhr_as_recv["broad_class"] = "receptor"
    sens_recv_df = pd.concat([binary_df, nhr_as_recv])
    sens_recv_df = sens_recv_df[sens_recv_df["broad_class"].isin(["machinery", "receptor"])]
    results["sensitivity_nhr_receptor"] = run_binary_analysis(
        sens_recv_df, "broad_class", "NHR→receptor"
    )

    # === ROBUSTNESS CHECK ===
    # Primary: require significance (paired CI excludes zero).
    # Sensitivity: require directional (point estimate DI > Frechet).
    primary_diff = results["primary"]["paired_differences"]["frechet_variance"]
    di_sig_primary = primary_diff["ci_lo"] > 0
    di_dir_nhr_mach = (results["sensitivity_nhr_machinery"]["paired_differences"]["frechet_variance"]["point_diff"] > 0)
    di_dir_nhr_recv = (results["sensitivity_nhr_receptor"]["paired_differences"]["frechet_variance"]["point_diff"] > 0)

    robust = di_sig_primary and di_dir_nhr_mach and di_dir_nhr_recv
    results["robustness"] = {
        "di_significant_primary": di_sig_primary,
        "primary_ci_lo": primary_diff["ci_lo"],
        "di_directional_nhr_machinery": di_dir_nhr_mach,
        "di_directional_nhr_receptor": di_dir_nhr_recv,
        "all_hold": robust,
        "rule": "primary CI excludes zero AND both sensitivity point estimates > 0",
        "conclusion": "ROBUST: DI significantly > Frechet in primary, directionally in both sensitivity variants" if robust
                       else "NOT ROBUST: DI advantage not established (primary not significant or sensitivity reversal)",
    }
    log(f"\n  Robustness gate: primary sig={di_sig_primary} (ci_lo={primary_diff['ci_lo']:.4f}), "
        f"NHR→mach dir={di_dir_nhr_mach}, NHR→recv dir={di_dir_nhr_recv}")
    log(f"  Robustness: {'PASS' if robust else 'FAIL'}")

    # === ORDINAL SPEARMAN (three-class, primary exclusion only) ===
    log("\n--- Ordinal Spearman (machinery < kinase < receptor) ---")
    ordinal_df = annotated[annotated["broad_class"].isin(["machinery", "kinase", "receptor"])]
    ordinal_map = {"machinery": 0, "kinase": 1, "receptor": 2}
    ordinal_labels = ordinal_df["broad_class"].map(ordinal_map).values

    results["ordinal_spearman"] = {}
    results["ordinal_n"] = {
        "machinery": int(n_mach), "kinase": int(n_kin), "receptor": int(n_recv)
    }
    for name, col in [("direction_instability", "direction_instability"),
                      ("frechet_variance", "frechet_variance"),
                      ("magnitude_cv", "magnitude_cv")]:
        r = bootstrap_spearman(ordinal_df[col].values, ordinal_labels, seed=42)
        results["ordinal_spearman"][name] = r
        log(f"  {name:25s} rho={r['rho']:.4f} [{r['ci_lo']:.4f}, {r['ci_hi']:.4f}] p={r['p_value']:.2e}")

    # === SAVE ===
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outpath = OUTPUT_DIR / "moa_classification_auroc_v2.json"
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {outpath}")


if __name__ == "__main__":
    main()
