"""Experiment 15: Positive control validation.

Two-pronged design:
  A) Biological: NHR drugs (known context-dependent) vs constitutive machinery
     drugs (known context-stable). Primary test: bootstrap CI on mean DI
     difference excludes zero.
  B) Synthetic spike-in: Noise-calibrated pseudo-perturbations with controlled
     ground-truth instability. Produces a limit-of-detection curve as a
     function of angular spread and number of cell lines.

Spec frozen at commit 49bb198. Do not modify spec after that commit.

Usage:
    uv run python experiments/15_positive_control_validation.py
"""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

DRUG_REPO = Path(__file__).resolve().parent.parent.parent / "drug-perturbation-geometry"
INSTABILITY_CSV = DRUG_REPO / "zenodo_v1" / "drug_instability_8949.csv"
OUTPUT_DIR = Path("results/15_positive_control")

N_BOOTSTRAP = 10_000
SPEC_COMMIT = "49bb198"

NUCLEAR_RECEPTOR_KEYWORDS = [
    "glucocorticoid receptor", "estrogen receptor", "androgen receptor",
    "progesterone receptor", "PPAR receptor", "retinoid receptor",
    "vitamin D receptor", "mineralocorticoid receptor",
]

CONSTITUTIVE_MACHINERY_KEYWORDS = [
    "proteasome inhibitor",
    "tubulin inhibitor", "tubulin polymerization inhibitor",
    "HDAC inhibitor",
    "topoisomerase inhibitor",
]

NARROW_MACHINERY_KEYWORDS = [
    "proteasome inhibitor",
    "tubulin inhibitor", "tubulin polymerization inhibitor",
]

SIGMA_GRID = np.arange(0.0, 1.55, 0.05)
K_VALUES = [5, 10, 13, 20, 40]
N_REP = 500
TAU_GRID = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
D_SIG = 978


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def is_nuclear_receptor(moa_str: str) -> bool:
    if pd.isna(moa_str) or moa_str == "":
        return False
    moa_lower = moa_str.lower()
    return any(kw.lower() in moa_lower for kw in NUCLEAR_RECEPTOR_KEYWORDS)


def is_constitutive_machinery(moa_str: str) -> bool:
    if pd.isna(moa_str) or moa_str == "":
        return False
    if is_nuclear_receptor(moa_str):
        return False
    moa_lower = moa_str.lower()
    return any(kw.lower() in moa_lower for kw in CONSTITUTIVE_MACHINERY_KEYWORDS)


def is_narrow_machinery(moa_str: str) -> bool:
    if pd.isna(moa_str) or moa_str == "":
        return False
    if is_nuclear_receptor(moa_str):
        return False
    moa_lower = moa_str.lower()
    return any(kw.lower() in moa_lower for kw in NARROW_MACHINERY_KEYWORDS)


def compute_di(signatures: np.ndarray) -> float:
    norms = np.linalg.norm(signatures, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-12)
    unit = signatures / norms
    cosines = unit @ unit.T
    K = len(signatures)
    mask = np.triu(np.ones((K, K), dtype=bool), k=1)
    return float(1.0 - cosines[mask].mean())


def bootstrap_mean_diff(group_a: np.ndarray, group_b: np.ndarray,
                        n_bootstrap: int = N_BOOTSTRAP, seed: int = 42) -> dict:
    point_diff = float(group_a.mean() - group_b.mean())
    rng = np.random.default_rng(seed)
    na, nb = len(group_a), len(group_b)
    boot_diffs = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        ba = group_a[rng.integers(0, na, size=na)]
        bb = group_b[rng.integers(0, nb, size=nb)]
        boot_diffs[i] = ba.mean() - bb.mean()
    ci_lo = float(np.percentile(boot_diffs, 2.5))
    ci_hi = float(np.percentile(boot_diffs, 97.5))
    p_value = float(np.mean(boot_diffs <= 0)) if point_diff > 0 else float(np.mean(boot_diffs >= 0))
    return {
        "point_diff": point_diff,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "p_value": p_value,
        "n_boot": n_bootstrap,
        "mean_a": float(group_a.mean()),
        "mean_b": float(group_b.mean()),
        "std_a": float(group_a.std()),
        "std_b": float(group_b.std()),
    }


def bootstrap_auroc(scores: np.ndarray, labels: np.ndarray,
                    n_bootstrap: int = N_BOOTSTRAP, seed: int = 42) -> dict:
    point_auroc = roc_auc_score(labels, scores)
    rng = np.random.default_rng(seed)
    n = len(labels)
    boot_aurocs = np.zeros(n_bootstrap)
    valid = 0
    attempts = 0
    while valid < n_bootstrap and attempts < n_bootstrap * 3:
        idx = rng.integers(0, n, size=n)
        b_labels = labels[idx]
        if b_labels.sum() == 0 or b_labels.sum() == len(b_labels):
            attempts += 1
            continue
        boot_aurocs[valid] = roc_auc_score(b_labels, scores[idx])
        valid += 1
        attempts += 1
    boot_aurocs = boot_aurocs[:valid]
    return {
        "auroc": float(point_auroc),
        "ci_lo": float(np.percentile(boot_aurocs, 2.5)),
        "ci_hi": float(np.percentile(boot_aurocs, 97.5)),
        "se": float(np.std(boot_aurocs)),
        "n_boot": valid,
    }


def generate_spike_in(sigma: float, K: int, tau: float, rng: np.random.Generator) -> float:
    v = rng.standard_normal(D_SIG)
    v /= np.linalg.norm(v)
    signatures = np.zeros((K, D_SIG))
    for i in range(K):
        n_i = rng.standard_normal(D_SIG)
        n_i -= np.dot(n_i, v) * v
        n_norm = np.linalg.norm(n_i)
        if n_norm > 1e-12:
            n_i /= n_norm
        theta_i = abs(rng.normal(0, sigma)) if sigma > 0 else 0.0
        s_i = np.cos(theta_i) * v + np.sin(theta_i) * n_i
        if tau > 0:
            s_i += rng.normal(0, tau, size=D_SIG)
        signatures[i] = s_i
    return compute_di(signatures)


def calibrate_tau(machinery_di: np.ndarray, K_cal: int = 13, n_rep: int = N_REP,
                  seed: int = 123) -> dict:
    target_mean = float(machinery_di.mean())
    target_std = float(machinery_di.std())
    log(f"  Calibration target: mean={target_mean:.4f}, std={target_std:.4f}")

    best_tau = TAU_GRID[0]
    best_score = float("inf")
    cal_results = {}

    for tau in TAU_GRID:
        rng = np.random.default_rng(seed)
        dis = np.array([generate_spike_in(0.0, K_cal, tau, rng) for _ in range(n_rep)])
        syn_mean = float(dis.mean())
        syn_std = float(dis.std())
        mean_err = abs(syn_mean - target_mean)
        std_ratio = max(syn_std / target_std, target_std / syn_std) if target_std > 0 and syn_std > 0 else float("inf")
        score = mean_err + 0.1 * abs(std_ratio - 1.0)
        cal_results[str(tau)] = {
            "syn_mean": syn_mean,
            "syn_std": syn_std,
            "mean_err": mean_err,
            "std_ratio": std_ratio,
        }
        if score < best_score:
            best_score = score
            best_tau = tau
        log(f"    tau={tau:.3f}: syn_mean={syn_mean:.4f}, syn_std={syn_std:.4f}, "
            f"mean_err={mean_err:.4f}, std_ratio={std_ratio:.2f}")

    best_info = cal_results[str(best_tau)]
    match_ok = (best_info["mean_err"] < 0.05 and best_info["std_ratio"] < 2.0)
    log(f"  Selected tau={best_tau:.3f} (match_ok={match_ok})")

    return {
        "selected_tau": best_tau,
        "match_ok": match_ok,
        "target_mean": target_mean,
        "target_std": target_std,
        "best_syn_mean": best_info["syn_mean"],
        "best_syn_std": best_info["syn_std"],
        "best_mean_err": best_info["mean_err"],
        "best_std_ratio": best_info["std_ratio"],
        "all_taus": cal_results,
    }


def run_lod_sweep(tau: float, seed: int = 456) -> dict:
    lod_results = {}
    for K in K_VALUES:
        log(f"\n  LOD sweep: K={K}")
        k_results = []
        for sigma in SIGMA_GRID:
            rng = np.random.default_rng(seed + int(K * 1000 + sigma * 100))
            dis = np.array([generate_spike_in(sigma, K, tau, rng) for _ in range(N_REP)])

            if sigma > 0:
                mean_angle_deg = float(np.degrees(sigma * np.sqrt(2.0 / np.pi)))
                p90_angle_deg = float(np.degrees(sigma * 1.2816))
            else:
                mean_angle_deg = 0.0
                p90_angle_deg = 0.0

            k_results.append({
                "sigma": float(sigma),
                "mean_angle_deg": mean_angle_deg,
                "p90_angle_deg": p90_angle_deg,
                "di_mean": float(dis.mean()),
                "di_p2_5": float(np.percentile(dis, 2.5)),
                "di_p97_5": float(np.percentile(dis, 97.5)),
                "di_values": dis.tolist(),
            })

        null_dis = np.array(k_results[0]["di_values"])
        threshold = float(np.percentile(null_dis, 95))

        for entry in k_results:
            entry_dis = np.array(entry["di_values"])
            entry["frac_above_threshold"] = float(np.mean(entry_dis > threshold))
            del entry["di_values"]

        lod_sigma = None
        for entry in k_results:
            if entry["frac_above_threshold"] >= 0.95:
                lod_sigma = entry["sigma"]
                break

        lod_results[str(K)] = {
            "null_threshold": threshold,
            "lod_sigma": lod_sigma,
            "lod_mean_angle_deg": float(np.degrees(lod_sigma * np.sqrt(2.0 / np.pi))) if lod_sigma else None,
            "sweep": k_results,
        }
        log(f"    null_threshold={threshold:.6f}, LOD sigma={lod_sigma}")

    return lod_results


def run_section_a(df: pd.DataFrame) -> dict:
    log("\n=== SECTION A: BIOLOGICAL POSITIVE CONTROL ===")

    nhr_mask = df["moa"].apply(is_nuclear_receptor)
    mach_mask = df["moa"].apply(is_constitutive_machinery)
    narrow_mach_mask = df["moa"].apply(is_narrow_machinery)

    nhr_di = df.loc[nhr_mask, "direction_instability"].values
    mach_di = df.loc[mach_mask, "direction_instability"].values
    narrow_mach_di = df.loc[narrow_mach_mask, "direction_instability"].values

    n_nhr = len(nhr_di)
    n_mach = len(mach_di)
    n_narrow = len(narrow_mach_di)
    log(f"  NHR (positive control): n={n_nhr}")
    log(f"  Constitutive machinery (negative control): n={n_mach}")
    log(f"  Narrow machinery (robustness): n={n_narrow}")

    assert n_nhr == 125, f"Expected 125 NHR, got {n_nhr}"
    assert n_mach == 59, f"Expected 59 constitutive machinery, got {n_mach}"
    overlap = nhr_mask & mach_mask
    assert overlap.sum() == 0, f"NHR/machinery overlap: {overlap.sum()}"

    # --- Primary test: bootstrap CI on mean difference ---
    log("\n--- Primary test: bootstrap CI on mean(DI_NHR) - mean(DI_machinery) ---")
    primary = bootstrap_mean_diff(nhr_di, mach_di, seed=42)
    ci_excludes_zero = primary["ci_lo"] > 0
    log(f"  mean(NHR)={primary['mean_a']:.4f} (std={primary['std_a']:.4f})")
    log(f"  mean(machinery)={primary['mean_b']:.4f} (std={primary['std_b']:.4f})")
    log(f"  diff={primary['point_diff']:.4f} [{primary['ci_lo']:.4f}, {primary['ci_hi']:.4f}] p={primary['p_value']:.6f}")
    log(f"  CI excludes zero: {ci_excludes_zero} {'PASS' if ci_excludes_zero else 'FAIL'}")

    # --- Secondary: Mann-Whitney U ---
    log("\n--- Secondary: Mann-Whitney U ---")
    u_stat, u_p = stats.mannwhitneyu(nhr_di, mach_di, alternative="two-sided")
    rank_biserial = 1.0 - (2.0 * u_stat) / (n_nhr * n_mach)
    log(f"  U={u_stat:.1f}, p={u_p:.2e}, rank-biserial r={rank_biserial:.4f}")

    # --- Secondary: AUROC ---
    log("\n--- Secondary: AUROC (NHR=positive class) ---")
    all_di = np.concatenate([nhr_di, mach_di])
    labels = np.concatenate([np.ones(n_nhr), np.zeros(n_mach)])
    auroc_result = bootstrap_auroc(all_di, labels, seed=42)
    log(f"  AUROC={auroc_result['auroc']:.4f} [{auroc_result['ci_lo']:.4f}, {auroc_result['ci_hi']:.4f}]")

    # --- Robustness 1: n_cell_lines >= 10 ---
    log("\n--- Robustness: n_cell_lines >= 10 ---")
    nhr_10 = df.loc[nhr_mask & (df["n_cell_lines"] >= 10), "direction_instability"].values
    mach_10 = df.loc[mach_mask & (df["n_cell_lines"] >= 10), "direction_instability"].values
    log(f"  NHR n>=10: {len(nhr_10)}, machinery n>=10: {len(mach_10)}")
    robust_matched = bootstrap_mean_diff(nhr_10, mach_10, seed=43)
    log(f"  diff={robust_matched['point_diff']:.4f} [{robust_matched['ci_lo']:.4f}, {robust_matched['ci_hi']:.4f}]")

    # --- Robustness 2: narrow machinery only ---
    log("\n--- Robustness: narrow machinery (proteasome + tubulin only) ---")
    log(f"  Narrow machinery: n={n_narrow}")
    robust_narrow = bootstrap_mean_diff(nhr_di, narrow_mach_di, seed=44)
    log(f"  diff={robust_narrow['point_diff']:.4f} [{robust_narrow['ci_lo']:.4f}, {robust_narrow['ci_hi']:.4f}]")

    return {
        "n_nhr": n_nhr,
        "n_constitutive_machinery": n_mach,
        "n_narrow_machinery": n_narrow,
        "primary_test": {
            "description": "Bootstrap CI on mean(DI_NHR) - mean(DI_machinery), alpha=0.05",
            "result": primary,
            "ci_excludes_zero": ci_excludes_zero,
            "pass": ci_excludes_zero,
        },
        "secondary_mannwhitney": {
            "U": float(u_stat),
            "p_value": float(u_p),
            "rank_biserial": float(rank_biserial),
        },
        "secondary_auroc": auroc_result,
        "robustness_n_ge_10": {
            "n_nhr": len(nhr_10),
            "n_machinery": len(mach_10),
            "result": robust_matched,
            "ci_excludes_zero": robust_matched["ci_lo"] > 0,
        },
        "robustness_narrow_machinery": {
            "n_narrow": n_narrow,
            "result": robust_narrow,
            "ci_excludes_zero": robust_narrow["ci_lo"] > 0,
        },
    }


def run_section_b(machinery_di: np.ndarray) -> dict:
    log("\n=== SECTION B: SYNTHETIC SPIKE-IN (LIMIT-OF-DETECTION) ===")

    log("\n--- Tau calibration ---")
    cal = calibrate_tau(machinery_di, K_cal=13, n_rep=N_REP, seed=123)
    tau = cal["selected_tau"]

    log("\n--- LOD sweep ---")
    lod = run_lod_sweep(tau, seed=456)

    return {
        "calibration": cal,
        "lod_by_K": lod,
    }


def main():
    log("=== EXPERIMENT 15: POSITIVE CONTROL VALIDATION ===")
    log(f"  Spec frozen at commit {SPEC_COMMIT}")

    df = pd.read_csv(INSTABILITY_CSV)
    log(f"  Loaded {len(df)} drugs")

    section_a = run_section_a(df)

    mach_mask = df["moa"].apply(is_constitutive_machinery)
    machinery_di = df.loc[mach_mask, "direction_instability"].values
    section_b = run_section_b(machinery_di)

    # --- Calibration overlay: where do real drugs fall? ---
    log("\n--- Calibration overlay ---")
    nhr_mask = df["moa"].apply(is_nuclear_receptor)
    nhr_di = df.loc[nhr_mask, "direction_instability"].values
    overlay = {
        "nhr_di_mean": float(nhr_di.mean()),
        "nhr_di_median": float(np.median(nhr_di)),
        "nhr_di_p25": float(np.percentile(nhr_di, 25)),
        "nhr_di_p75": float(np.percentile(nhr_di, 75)),
        "machinery_di_mean": float(machinery_di.mean()),
        "machinery_di_median": float(np.median(machinery_di)),
        "machinery_di_p25": float(np.percentile(machinery_di, 25)),
        "machinery_di_p75": float(np.percentile(machinery_di, 75)),
    }
    log(f"  NHR DI: mean={overlay['nhr_di_mean']:.4f}, median={overlay['nhr_di_median']:.4f}")
    log(f"  Machinery DI: mean={overlay['machinery_di_mean']:.4f}, median={overlay['machinery_di_median']:.4f}")

    results = {
        "experiment": "15_positive_control_validation",
        "spec_commit": SPEC_COMMIT,
        "section_a_biological": section_a,
        "section_b_synthetic": section_b,
        "calibration_overlay": overlay,
        "pass_fail": {
            "primary_gate": section_a["primary_test"]["pass"],
            "description": "Bootstrap CI on mean(DI_NHR) - mean(DI_machinery) excludes zero at alpha=0.05",
        },
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outpath = OUTPUT_DIR / "positive_control_validation.json"
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {outpath}")

    if section_a["primary_test"]["pass"]:
        log("\n  RESULT: PASS — DI detects known context-dependent perturbations")
    else:
        log("\n  RESULT: FAIL — DI does not significantly separate NHR from machinery")


if __name__ == "__main__":
    main()
