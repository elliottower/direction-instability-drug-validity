"""Experiment R3: HDAC-removal sensitivity and non-HDAC drug pair search for H3.

Tests whether the H3 result (Spearman rho = 0.376 between phenotype-projected
bracket and on-target enrichment) survives removal of HDAC inhibitors and
other perturbations. Finds concrete non-HDAC drug pairs illustrating the
discriminative power of phenotype projection.

Pre-registered in experiments/PREREGISTRATION_EXTENDED.md before execution.

Usage:
    uv run python experiments/r3_hdac_sensitivity_and_drug_pair.py
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats
from tqdm import tqdm

RESULTS_JSON = Path(__file__).resolve().parent.parent / "results" / "03_phenotype_projection" / "phenotype_projection_results.json"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "results" / "r3_hdac_sensitivity"

N_BOOTSTRAP = 1000
RAW_BRACKET_MATCH_TOL = 0.05
TOP_N_CONSISTENT_REMOVE = 20
N_PAIRS_REPORT = 5


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def spearman_with_pvalue(x, y):
    """Compute Spearman rho and two-sided p-value."""
    result = stats.spearmanr(x, y)
    return float(result.correlation), float(result.pvalue)


def compute_correlations(drugs: list[dict], label: str) -> dict:
    """Compute projected and raw Spearman correlations for a drug set."""
    proj = np.array([d["projected_bracket"] for d in drugs])
    raw = np.array([d["raw_bracket"] for d in drugs])
    enrich = np.array([d["on_target_enrichment"] for d in drugs])

    rho_proj, p_proj = spearman_with_pvalue(proj, enrich)
    rho_raw, p_raw = spearman_with_pvalue(raw, enrich)

    result = {
        "label": label,
        "n_drugs": len(drugs),
        "rho_proj": rho_proj,
        "p_proj": p_proj,
        "rho_raw": rho_raw,
        "p_raw": p_raw,
    }

    log(f"  {label}: n={len(drugs)}, rho_proj={rho_proj:.4f} (p={p_proj:.2e}), rho_raw={rho_raw:.4f} (p={p_raw:.2e})")
    return result


def bootstrap_spearman(x: np.ndarray, y: np.ndarray, n_bootstrap: int) -> dict:
    """Bootstrap 95% percentile CI for Spearman rho."""
    rng = np.random.default_rng(2026070803)
    n = len(x)
    boot_rhos = np.zeros(n_bootstrap)

    for i in tqdm(range(n_bootstrap), desc="Bootstrap"):
        idx = rng.integers(0, n, size=n)
        rho, _ = stats.spearmanr(x[idx], y[idx])
        boot_rhos[i] = rho

    ci_lo = float(np.percentile(boot_rhos, 2.5))
    ci_hi = float(np.percentile(boot_rhos, 97.5))
    point = float(stats.spearmanr(x, y).correlation)

    return {
        "point": point,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "se": float(np.std(boot_rhos)),
        "n_bootstrap": n_bootstrap,
    }


def find_drug_pairs(drugs: list[dict]) -> list[dict]:
    """Find non-HDAC drug pairs with similar raw_bracket but divergent projected_bracket and enrichment."""
    non_hdac = [d for d in drugs if "HDAC" not in d["target"].upper()]

    proj_vals = np.array([d["projected_bracket"] for d in non_hdac])
    enrich_vals = np.array([d["on_target_enrichment"] for d in non_hdac])

    proj_p75 = float(np.percentile(proj_vals, 75))
    proj_p25 = float(np.percentile(proj_vals, 25))
    enrich_p75 = float(np.percentile(enrich_vals, 75))
    enrich_p25 = float(np.percentile(enrich_vals, 25))

    log(f"  Thresholds: proj P25={proj_p25:.3f}, P75={proj_p75:.3f}, enrich P25={enrich_p25:.4f}, P75={enrich_p75:.4f}")

    # Classify drugs
    high_drugs = [d for d in non_hdac
                  if d["projected_bracket"] >= proj_p75
                  and d["on_target_enrichment"] >= enrich_p75]
    low_drugs = [d for d in non_hdac
                 if d["projected_bracket"] <= proj_p25
                 and d["on_target_enrichment"] <= enrich_p25]

    log(f"  High-high drugs (proj>=P75 & enrich>=P75): {len(high_drugs)}")
    log(f"  Low-low drugs (proj<=P25 & enrich<=P25): {len(low_drugs)}")

    # Find pairs with similar raw_bracket
    pairs = []
    for h in high_drugs:
        for lo in low_drugs:
            raw_diff = abs(h["raw_bracket"] - lo["raw_bracket"])
            if raw_diff <= RAW_BRACKET_MATCH_TOL:
                proj_diff = abs(h["projected_bracket"] - lo["projected_bracket"])
                pairs.append({
                    "high_drug": h["drug"],
                    "high_target": h["target"],
                    "high_raw_bracket": h["raw_bracket"],
                    "high_projected_bracket": h["projected_bracket"],
                    "high_on_target_enrichment": h["on_target_enrichment"],
                    "high_n_celllines": h["n_celllines"],
                    "low_drug": lo["drug"],
                    "low_target": lo["target"],
                    "low_raw_bracket": lo["raw_bracket"],
                    "low_projected_bracket": lo["projected_bracket"],
                    "low_on_target_enrichment": lo["on_target_enrichment"],
                    "low_n_celllines": lo["n_celllines"],
                    "raw_bracket_diff": raw_diff,
                    "projected_bracket_diff": proj_diff,
                })

    # Sort by projected_bracket_diff descending (most dramatic separation)
    pairs.sort(key=lambda p: p["projected_bracket_diff"], reverse=True)
    return pairs


def main():
    log("=== Experiment R3: HDAC-removal sensitivity and drug pair search ===")

    # Load data
    if not RESULTS_JSON.exists():
        log(f"ERROR: Results file not found: {RESULTS_JSON}")
        sys.exit(1)

    with open(RESULTS_JSON) as f:
        all_drugs = json.load(f)

    log(f"Loaded {len(all_drugs)} drugs from {RESULTS_JSON.name}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_results = {}

    # --- A1: Full-set reproduction ---
    log("\n--- A1: Full-set reproduction ---")
    a1 = compute_correlations(all_drugs, "full_set")
    all_results["a1_full_set"] = a1

    # --- A2: HDAC removal ---
    log("\n--- A2: HDAC removal ---")
    hdac_drugs = [d for d in all_drugs if "HDAC" in d["target"].upper()]
    non_hdac_drugs = [d for d in all_drugs if "HDAC" not in d["target"].upper()]

    # HDAC filter validation (pre-registered requirement)
    log(f"  HDAC filter: removed {len(hdac_drugs)} drugs")
    assert 10 <= len(hdac_drugs) <= 40, (
        f"HDAC filter removed {len(hdac_drugs)} drugs — outside expected [10, 40] range. "
        f"Filter may be silently under- or over-selecting. Investigate before proceeding."
    )
    log(f"  Removed HDAC drugs (name → target):")
    for d in sorted(hdac_drugs, key=lambda x: x["drug"]):
        log(f"    {d['drug']}: target={d['target']}, raw={d['raw_bracket']:.3f}, proj={d['projected_bracket']:.3f}")

    # Check for known HDAC synonyms missed by substring filter
    hdac_synonyms = ["histone deacetylase", "class i hdac", "class ii hdac", "pan-hdac", "sirtuin"]
    missed_by_primary = [
        d for d in non_hdac_drugs
        if any(syn in d["target"].lower() for syn in hdac_synonyms)
    ]
    if missed_by_primary:
        log(f"  WARNING: {len(missed_by_primary)} drugs match HDAC synonyms but were NOT caught by primary filter:")
        for d in missed_by_primary:
            log(f"    {d['drug']}: target={d['target']}")
    else:
        log(f"  Synonym check passed: no drugs matching {hdac_synonyms} missed by primary filter")

    log(f"  {len(non_hdac_drugs)} drugs remain after HDAC removal")

    a2 = compute_correlations(non_hdac_drugs, "hdac_removed")
    a2["n_hdac_removed"] = len(hdac_drugs)
    a2["hdac_drug_names"] = sorted(d["drug"] for d in hdac_drugs)
    a2["hdac_targets"] = sorted(set(d["target"] for d in hdac_drugs))
    a2["missed_synonyms"] = [{"drug": d["drug"], "target": d["target"]} for d in missed_by_primary]
    all_results["a2_hdac_removed"] = a2

    # Decision: does rho_proj remain > 0.3?
    passes_primary = a2["rho_proj"] > 0.3
    log(f"  PRIMARY CRITERION (rho_proj > 0.3 after HDAC removal): {'PASS' if passes_primary else 'FAIL'} (rho={a2['rho_proj']:.4f})")
    all_results["a2_hdac_removed"]["passes_primary"] = passes_primary

    # --- A3: Top-20 most consistent removal (by raw_bracket) ---
    log("\n--- A3: Top-20 most consistent drug removal (by raw_bracket) ---")
    sorted_by_raw = sorted(all_drugs, key=lambda d: d["raw_bracket"])
    top20_consistent = sorted_by_raw[:TOP_N_CONSISTENT_REMOVE]
    remaining = sorted_by_raw[TOP_N_CONSISTENT_REMOVE:]

    log(f"  Removed {TOP_N_CONSISTENT_REMOVE} drugs with lowest raw_bracket:")
    for d in top20_consistent:
        log(f"    {d['drug']}: raw={d['raw_bracket']:.3f}, proj={d['projected_bracket']:.3f}, target={d['target']}")

    a3 = compute_correlations(remaining, "top20_raw_removed")
    a3["removed_drugs"] = [d["drug"] for d in top20_consistent]
    a3["raw_bracket_range_removed"] = [
        float(top20_consistent[0]["raw_bracket"]),
        float(top20_consistent[-1]["raw_bracket"]),
    ]
    passes_a3 = a3["rho_proj"] > 0.20
    a3["passes_criterion"] = passes_a3
    log(f"  A3 CRITERION (rho_proj > 0.20 after top-20 raw removal): {'PASS' if passes_a3 else 'FAIL'} (rho={a3['rho_proj']:.4f})")
    all_results["a3_top20_raw_removed"] = a3

    # --- A3b: Top-20 most phenotype-aligned removal (by projected_bracket) ---
    log("\n--- A3b: Top-20 most phenotype-aligned drug removal (by projected_bracket) ---")
    sorted_by_proj = sorted(all_drugs, key=lambda d: d["projected_bracket"], reverse=True)
    top20_projected = sorted_by_proj[:TOP_N_CONSISTENT_REMOVE]
    remaining_proj = sorted_by_proj[TOP_N_CONSISTENT_REMOVE:]

    log(f"  Removed {TOP_N_CONSISTENT_REMOVE} drugs with highest projected_bracket:")
    for d in top20_projected:
        log(f"    {d['drug']}: proj={d['projected_bracket']:.3f}, raw={d['raw_bracket']:.3f}, enrich={d['on_target_enrichment']:.4f}, target={d['target']}")

    a3b = compute_correlations(remaining_proj, "top20_projected_removed")
    a3b["removed_drugs"] = [d["drug"] for d in top20_projected]
    a3b["projected_bracket_range_removed"] = [
        float(top20_projected[-1]["projected_bracket"]),
        float(top20_projected[0]["projected_bracket"]),
    ]
    passes_a3b = a3b["rho_proj"] > 0.20
    a3b["passes_criterion"] = passes_a3b
    log(f"  A3b CRITERION (rho_proj > 0.20 after top-20 projected removal): {'PASS' if passes_a3b else 'FAIL'} (rho={a3b['rho_proj']:.4f})")
    all_results["a3b_top20_projected_removed"] = a3b

    # --- A4: Bootstrap CI ---
    log("\n--- A4: Bootstrap CI (1000 resamples) ---")
    proj_arr = np.array([d["projected_bracket"] for d in all_drugs])
    raw_arr = np.array([d["raw_bracket"] for d in all_drugs])
    enrich_arr = np.array([d["on_target_enrichment"] for d in all_drugs])

    log("  Bootstrapping rho_proj...")
    boot_proj = bootstrap_spearman(proj_arr, enrich_arr, N_BOOTSTRAP)
    log(f"  rho_proj: {boot_proj['point']:.4f} [{boot_proj['ci_lo']:.4f}, {boot_proj['ci_hi']:.4f}]")

    log("  Bootstrapping rho_raw...")
    boot_raw = bootstrap_spearman(raw_arr, enrich_arr, N_BOOTSTRAP)
    log(f"  rho_raw: {boot_raw['point']:.4f} [{boot_raw['ci_lo']:.4f}, {boot_raw['ci_hi']:.4f}]")

    all_results["a4_bootstrap"] = {
        "rho_proj_bootstrap": boot_proj,
        "rho_raw_bootstrap": boot_raw,
        "ci_excludes_zero_proj": boot_proj["ci_lo"] > 0 or boot_proj["ci_hi"] < 0,
        "ci_excludes_zero_raw": boot_raw["ci_lo"] > 0 or boot_raw["ci_hi"] < 0,
    }

    passes_secondary = boot_proj["ci_lo"] > 0 or boot_proj["ci_hi"] < 0
    log(f"  SECONDARY CRITERION (95% CI excludes zero for rho_proj): {'PASS' if passes_secondary else 'FAIL'}")

    # --- A5: Drug pair search ---
    log("\n--- A5: Non-HDAC drug pair search ---")
    pairs = find_drug_pairs(all_drugs)
    log(f"  Found {len(pairs)} candidate pairs")

    top_pairs = pairs[:N_PAIRS_REPORT]
    for i, p in enumerate(top_pairs):
        log(f"\n  Pair {i + 1}:")
        log(f"    HIGH: {p['high_drug']} (target={p['high_target']}, n={p['high_n_celllines']})")
        log(f"      raw={p['high_raw_bracket']:.3f}, proj={p['high_projected_bracket']:.3f}, enrich={p['high_on_target_enrichment']:.4f}")
        log(f"    LOW:  {p['low_drug']} (target={p['low_target']}, n={p['low_n_celllines']})")
        log(f"      raw={p['low_raw_bracket']:.3f}, proj={p['low_projected_bracket']:.3f}, enrich={p['low_on_target_enrichment']:.4f}")
        log(f"    raw_bracket_diff={p['raw_bracket_diff']:.4f}, projected_bracket_diff={p['projected_bracket_diff']:.3f}")

    all_results["a5_drug_pairs"] = {
        "n_total_pairs": len(pairs),
        "top_pairs": top_pairs,
        "passes_exploratory": len(pairs) > 0,
    }

    # --- Summary ---
    log("\n=== SUMMARY ===")
    log(f"A1 full-set:         rho_proj={a1['rho_proj']:.4f}, rho_raw={a1['rho_raw']:.4f}")
    log(f"A2 HDAC-removed:     rho_proj={a2['rho_proj']:.4f}, rho_raw={a2['rho_raw']:.4f}  PRIMARY: {'PASS' if passes_primary else 'FAIL'}")
    log(f"A3 top20-raw:        rho_proj={a3['rho_proj']:.4f}, rho_raw={a3['rho_raw']:.4f}  CRITERION: {'PASS' if passes_a3 else 'FAIL'}")
    log(f"A3b top20-projected: rho_proj={a3b['rho_proj']:.4f}, rho_raw={a3b['rho_raw']:.4f}  CRITERION: {'PASS' if passes_a3b else 'FAIL'}")
    log(f"A4 bootstrap CI:     rho_proj [{boot_proj['ci_lo']:.4f}, {boot_proj['ci_hi']:.4f}]  SECONDARY: {'PASS' if passes_secondary else 'FAIL'}")
    log(f"A5 drug pairs:       {len(pairs)} found  EXPLORATORY: {'PASS' if len(pairs) > 0 else 'FAIL'}")

    all_results["summary"] = {
        "passes_primary": passes_primary,
        "passes_a3": passes_a3,
        "passes_a3b": passes_a3b,
        "passes_secondary": passes_secondary,
        "passes_exploratory": len(pairs) > 0,
    }

    # Save results
    out_path = OUTPUT_DIR / "r3_results.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    log(f"\nResults saved to {out_path}")

    # Save drug pairs separately for easy inspection
    pairs_path = OUTPUT_DIR / "r3_drug_pairs.json"
    with open(pairs_path, "w") as f:
        json.dump(top_pairs, f, indent=2)
    log(f"Drug pairs saved to {pairs_path}")


if __name__ == "__main__":
    main()
