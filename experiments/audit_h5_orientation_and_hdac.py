"""Audit of two defects found in direction_instability_confound_audit_v1.

Recomputes every number quoted in the 2026-09-07 audit so that none of them
lives only in a terminal. Writes the whole result structure to
results/audit_h5_orientation_and_hdac/audit.json.

Two findings:

1. H5 orientation. The leave-one-out win test in
   experiments/05c_h5_full_leave_one_out.py is `ts_rho > raw_rho`, comparing
   SIGNED Spearman correlations against held-out cosine. Direction instability
   measures inconsistency and held-out cosine measures agreement, so a working
   D must correlate negatively. TS = D - lambda*Var_Frechet subtracts a term
   that also grows with inconsistency and grows faster, inverting the ordering,
   so TS correlates positively. The signed comparison is therefore won by
   orientation rather than by prediction.

2. HDAC figure. fig5_hdac() in paper/generate_figures.py hardcodes six D values
   and six cell-line counts that reproduce no stored artifact. Treating
   selectivity as three ordered categories rather than six untied ranks gives
   the same Kendall tau_b for both candidate quantities.

Usage:
    uv run --no-project --with numpy python experiments/audit_h5_orientation_and_hdac.py
"""
import json
from itertools import combinations, permutations
from math import sqrt
from pathlib import Path

import numpy as np

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
COMPANION = Path("/Users/elliottower/Documents/GitHub/drug-perturbation-geometry")
OUT_DIR = REPO / "results" / "audit_h5_orientation_and_hdac"

HDAC_DRUGS = ["panobinostat", "vorinostat", "belinostat", "entinostat", "tubacin", "pci-34051"]
HDAC_CLASS = ["Pan-HDAC", "Pan-HDAC", "Pan-HDAC", "Class I", "HDAC6", "HDAC8"]
HDAC_CATEGORY_RANK = [1, 1, 1, 2, 3, 3]          # pan < class-I < isoform-selective
PLOTTED_D = [0.14, 0.16, 0.19, 0.27, 0.32, 0.45]  # hardcoded in fig5_hdac()
PLOTTED_N = [64, 63, 42, 38, 25, 18]


def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    return float(rx @ ry / sqrt((rx @ rx) * (ry @ ry)))


def exact_spearman_p(rho, n=6):
    """Exact permutation p-values by full enumeration of n! orderings."""
    vals = []
    for p in permutations(range(1, n + 1)):
        d2 = sum((a - b) ** 2 for a, b in zip(p, range(1, n + 1)))
        vals.append(1 - 6 * d2 / (n * (n * n - 1)))
    tol = 1e-9
    one = sum(1 for v in vals if v >= rho - tol) / len(vals)
    two = sum(1 for v in vals if abs(v) >= rho - tol) / len(vals)
    return round(one, 4), round(two, 4)


def kendall_tau_b(categories, values):
    """tau_b over cross-category pairs, with ties in the category variable."""
    n = len(values)
    concordant = discordant = cross = 0
    for i, j in combinations(range(n), 2):
        if categories[i] == categories[j]:
            continue
        cross += 1
        s = (categories[j] - categories[i]) * (values[j] - values[i])
        concordant += s > 0
        discordant += s < 0
    n0 = n * (n - 1) // 2
    sizes = {}
    for c in categories:
        sizes[c] = sizes.get(c, 0) + 1
    ties_x = sum(k * (k - 1) // 2 for k in sizes.values())
    tau = (concordant - discordant) / sqrt((n0 - ties_x) * n0)
    return {"cross_category_pairs": cross, "concordant": concordant,
            "discordant": discordant, "tau_b": round(tau, 4)}


def direction_instability(signatures):
    """1 - mean pairwise cosine over distinct pairs. Complement of PPC."""
    unit = signatures / np.maximum(np.linalg.norm(signatures, axis=1, keepdims=True), 1e-10)
    cos = unit @ unit.T
    iu = np.triu_indices(cos.shape[0], k=1)
    return 1.0 - float(cos[iu].mean())


def transport_stable(signatures, frechet_penalty=1.0):
    """D - lambda * mean squared angular deviation from the Frechet mean."""
    raw = direction_instability(signatures)
    unit = signatures / np.maximum(np.linalg.norm(signatures, axis=1, keepdims=True), 1e-10)
    mean_dir = unit.mean(axis=0)
    mean_dir = mean_dir / (np.linalg.norm(mean_dir) + 1e-10)
    deviations = np.arccos(np.clip(unit @ mean_dir, -1, 1))
    return raw - frechet_penalty * float(np.mean(deviations ** 2))


def audit_h5_orientation():
    folds = json.loads((REPO / "results/05c_h5_full/h5_full_leave_one_out_results.json").read_text())
    raw = [f["raw_spearman"] for f in folds]
    ts = [f["ts_spearman"] for f in folds]
    n = len(folds)
    largest = sorted(folds, key=lambda f: -f["n_drugs"])[:2]
    return {
        "n_folds": n,
        "win_test_in_manuscript": "ts_rho > raw_rho (signed)",
        "signed_wins": sum(1 for r, t in zip(raw, ts) if t > r),
        "magnitude_wins": sum(1 for r, t in zip(raw, ts) if abs(t) > abs(r)),
        "mean_raw_rho": round(float(np.mean(raw)), 4),
        "mean_ts_rho": round(float(np.mean(ts)), 4),
        "mean_abs_raw_rho": round(float(np.mean(np.abs(raw))), 4),
        "mean_abs_ts_rho": round(float(np.mean(np.abs(ts))), 4),
        "raw_negative_in_all_folds": all(r < 0 for r in raw),
        "ts_positive_in_n_folds": sum(1 for t in ts if t > 0),
        "two_largest_folds": [
            {"cell": f["holdout_cell"], "n_drugs": f["n_drugs"],
             "raw_rho": round(f["raw_spearman"], 4), "ts_rho": round(f["ts_spearman"], 4)}
            for f in largest
        ],
    }


def audit_score_inversion(n_trials=400, n_genes=200, seed=0):
    """Is TS an order-reversing transform of D? Synthetic data, no repo values."""
    rng = np.random.default_rng(seed)
    d_vals, ts_vals = [], []
    for _ in range(n_trials):
        k = int(rng.integers(5, 25))
        noise = float(rng.uniform(0.1, 4.0))
        shared = rng.standard_normal(n_genes)
        sigs = shared + rng.standard_normal((k, n_genes)) * noise
        d_vals.append(direction_instability(sigs))
        ts_vals.append(transport_stable(sigs))
    return {
        "n_synthetic_perturbations": n_trials,
        "seed": seed,
        "spearman_D_vs_TS": round(spearman(np.array(d_vals), np.array(ts_vals)), 4),
        "D_range": [round(min(d_vals), 4), round(max(d_vals), 4)],
        "TS_range": [round(min(ts_vals), 4), round(max(ts_vals), 4)],
        "note": "Strong negative rank correlation means TS reverses D's ordering.",
    }


def audit_hdac():
    companion = {str(r.get("drug_name", "")).lower(): r
                 for r in json.loads((COMPANION / "results/01_cross_cellline/real_results.json").read_text())}
    local = {str(r.get("drug", "")).lower(): r
             for r in json.loads((REPO / "results/01_toxicity_failure/toxicity_results.json").read_text())}

    per_drug, series = [], {"plotted_in_figure5": PLOTTED_D,
                            "companion_direction_instability": [],
                            "thispaper_corrected_instability": [],
                            "thispaper_raw_instability": []}
    for i, name in enumerate(HDAC_DRUGS):
        c, t = companion[name], local[name]
        series["companion_direction_instability"].append(round(c["direction_instability"], 4))
        series["thispaper_corrected_instability"].append(round(t["corrected_instability"], 4))
        series["thispaper_raw_instability"].append(round(t["raw_instability"], 4))
        per_drug.append({
            "drug": name, "selectivity_class": HDAC_CLASS[i],
            "category_rank": HDAC_CATEGORY_RANK[i], "untied_rank_as_used_in_paper": i + 1,
            "plotted_D": PLOTTED_D[i], "plotted_n_celllines": PLOTTED_N[i],
            "companion_direction_instability": round(c["direction_instability"], 4),
            "companion_n_cell_lines": c["n_cell_lines"],
            "thispaper_raw_instability": round(t["raw_instability"], 4),
            "thispaper_corrected_instability": round(t["corrected_instability"], 4),
            "thispaper_n_celllines": t["n_celllines"],
        })

    untied = list(range(1, 7))
    stats = {}
    for key, vals in series.items():
        rho = spearman(np.array(vals, dtype=float), np.array(untied, dtype=float))
        one, two = exact_spearman_p(rho)
        stats[key] = {
            "spearman_vs_untied_rank_1_to_6": round(rho, 4),
            "exact_one_tailed_p": one, "exact_two_tailed_p": two,
            "kendall_tau_b_vs_three_categories": kendall_tau_b(HDAC_CATEGORY_RANK, vals),
        }

    # Does any stored companion field reproduce the plotted values?
    fields = [k for k, v in companion[HDAC_DRUGS[0]].items() if isinstance(v, (int, float))]
    scan = sorted(
        ({"field": f,
          "mean_abs_error_vs_plotted": round(float(np.mean(
              [abs(companion[n][f] - t) for n, t in zip(HDAC_DRUGS, PLOTTED_D)])), 4),
          "values": [round(companion[n][f], 4) for n in HDAC_DRUGS]}
         for f in fields),
        key=lambda r: r["mean_abs_error_vs_plotted"])[:5]

    return {"per_drug": per_drug, "statistics": stats,
            "search_for_stored_field_matching_plotted_values": scan,
            "conclusion": "No stored numeric field reproduces the plotted values. "
                          "Under three ordered categories both candidate quantities give "
                          "the same tau_b, so the choice between them is provenance, not result."}


def main():
    audit = {
        "generated_by": "experiments/audit_h5_orientation_and_hdac.py",
        "finding_1_h5_orientation": audit_h5_orientation(),
        "finding_1_score_inversion_check": audit_score_inversion(),
        "finding_2_hdac_figure": audit_hdac(),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "audit.json").write_text(json.dumps(audit, indent=2))

    h5 = audit["finding_1_h5_orientation"]
    print(f"H5: signed wins {h5['signed_wins']}/{h5['n_folds']}, "
          f"magnitude wins {h5['magnitude_wins']}/{h5['n_folds']}")
    print(f"    mean |raw rho| = {h5['mean_abs_raw_rho']}, mean |TS rho| = {h5['mean_abs_ts_rho']}")
    print(f"    Spearman(D, TS) synthetic = {audit['finding_1_score_inversion_check']['spearman_D_vs_TS']}")
    for k, v in audit["finding_2_hdac_figure"]["statistics"].items():
        print(f"HDAC {k:34s} rho={v['spearman_vs_untied_rank_1_to_6']:+.4f} "
              f"tau_b={v['kendall_tau_b_vs_three_categories']['tau_b']:+.4f}")
    print(f"\nWritten to {OUT_DIR / 'audit.json'}")


if __name__ == "__main__":
    main()
