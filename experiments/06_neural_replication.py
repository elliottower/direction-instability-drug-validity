"""Experiment 06: Neural bracket norm replication with validity ladder.

Tests H6: Bracket norm predicts optogenetic silencing importance in neural
data, but ONLY after neuron-count correction (BN/sqrt(n)). Without the
correction, the signal is confounded by recording yield.

This replicates Tower (2026) and demonstrates the validity ladder applied
to a non-pharmacological domain: the same bracket norm that captures drug
mechanism transport also captures neural causal importance, but only when
the appropriate validity condition (size correction) is applied.

Data: Steinmetz et al. 2019 Neuropixels + Zatka-Haas et al. 2021 silencing.
All cached locally from prior work.

Usage:
    uv run python experiments/06_neural_replication.py --data PATH_TO_BRACKET_NORM_CACHE
    uv run python experiments/06_neural_replication.py --synthetic
"""
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent))


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


SILENCING_EFFECTS = {
    "PL": 0.3333, "ORB": 0.3085, "VISpm": 0.2248, "MOp": 0.1869,
    "SSp": 0.1862, "VISl": 0.1722, "SSs": 0.1600, "MOs": 0.1529,
    "ACA": 0.1451, "RSP": 0.1421, "VISp": 0.1414, "VISam": 0.0818,
}


def compute_bracket_norm_per_region(
    activity: dict,
    labels: dict,
    time_window: slice = slice(15, 35),
) -> dict:
    """Compute bracket norm per brain region.

    Bracket norm = how much the coding direction changes as stimulus
    evidence varies. Operationalized as 1 - mean(pairwise cosine) of
    choice-coding vectors across evidence quartiles.

    Args:
        activity: region -> (n_trials, n_neurons, n_timebins) array
        labels: region -> (n_trials,) choice labels

    Returns:
        region -> {"raw_bn": float, "n_neurons": int, "corrected_bn": float}
    """
    results = {}
    for region in activity:
        X = activity[region][:, :, time_window].mean(axis=2)
        y = labels[region]
        n_neurons = X.shape[1]

        unique_y = np.unique(y)
        if len(unique_y) < 2:
            continue

        quartiles = np.quantile(np.abs(y - y.mean()), [0.25, 0.5, 0.75])
        q_labels = np.digitize(np.abs(y - y.mean()), quartiles)

        coding_vectors = []
        for q in range(4):
            mask = q_labels == q
            if mask.sum() < 5:
                continue
            X_q = X[mask]
            y_q = y[mask]
            mean_left = X_q[y_q == unique_y[0]].mean(axis=0)
            mean_right = X_q[y_q == unique_y[1]].mean(axis=0)
            diff = mean_right - mean_left
            norm = np.linalg.norm(diff)
            if norm > 1e-10:
                coding_vectors.append(diff / norm)

        if len(coding_vectors) < 3:
            continue

        coding_vectors = np.array(coding_vectors)
        cosines = coding_vectors @ coding_vectors.T
        n = cosines.shape[0]
        triu_idx = np.triu_indices(n, k=1)
        raw_bn = 1.0 - float(cosines[triu_idx].mean())
        corrected_bn = raw_bn / np.sqrt(n_neurons)

        results[region] = {
            "raw_bn": float(raw_bn),
            "n_neurons": n_neurons,
            "corrected_bn": float(corrected_bn),
        }

    return results


def run_synthetic():
    """Estimator unit-test: validates that BN/sqrt(n) removes a planted size confound.

    This is NOT evidence for H6 — it verifies the correction arithmetic by
    injecting a known sqrt(n) + linear-n dependence and checking recovery.
    H6 evidence comes only from the --real Steinmetz path.
    """
    log("=== SYNTHETIC MODE (estimator unit-test, not H6 evidence) ===")
    rng = np.random.default_rng(seed=20260707_06)

    n_regions = 20
    neuron_counts = rng.integers(10, 200, size=n_regions)
    true_causal = rng.uniform(0, 0.4, size=n_regions)

    raw_bns = []
    corrected_bns = []
    for i in range(n_regions):
        n = neuron_counts[i]
        raw = true_causal[i] * np.sqrt(n) * 0.1 + 0.005 * n + rng.normal(0, 0.01)
        raw_bns.append(max(0.001, raw))
        corrected_bns.append(max(0.001, raw) / np.sqrt(n))

    raw_bns = np.array(raw_bns)
    corrected_bns = np.array(corrected_bns)

    r_raw_n, _ = stats.spearmanr(raw_bns, neuron_counts)
    r_corr_n, _ = stats.spearmanr(corrected_bns, neuron_counts)
    r_raw_causal, _ = stats.spearmanr(raw_bns, true_causal)
    r_corr_causal, _ = stats.spearmanr(corrected_bns, true_causal)

    log(f"Raw BN vs neuron count:  rho={r_raw_n:.3f} (should be high — confounded)")
    log(f"Corrected vs count:      rho={r_corr_n:.3f} (should be low — deconfounded)")
    log(f"Raw BN vs true causal:   rho={r_raw_causal:.3f}")
    log(f"Corrected vs true causal: rho={r_corr_causal:.3f} (should be highest)")

    if abs(r_corr_n) < abs(r_raw_n) and r_corr_causal > r_raw_causal:
        log("PASS: correction removes size confound and improves causal correlation")
    else:
        log("FAIL: correction does not work as expected")
        sys.exit(1)


def run_real(data_dir: Path, output_dir: Path):
    """Run H6 on cached Steinmetz bracket norm results."""
    log("=== REAL MODE ===")

    cache_files = list(data_dir.glob("*bracket*")) + list(data_dir.glob("*bundle*"))
    log(f"Looking for cached results in {data_dir}")
    log(f"  Found: {[f.name for f in cache_files]}")

    bn_results_path = data_dir / "bracket_norms_per_region.json"
    if not bn_results_path.exists():
        results_dirs = list(data_dir.glob("results*"))
        for rd in results_dirs:
            candidates = list(rd.glob("*bracket*")) + list(rd.glob("*bn_*"))
            if candidates:
                bn_results_path = candidates[0]
                break

    if not bn_results_path.exists():
        log("No cached bracket norm results found. Computing from raw Steinmetz data...")
        log("ERROR: Raw Steinmetz data loading not implemented in this script.")
        log("  Copy results from bracket-norm repo: results_bundle/batch3/")
        sys.exit(1)

    log(f"Loading cached results from {bn_results_path}")
    with open(bn_results_path) as f:
        cached = json.load(f)

    region_bns = {}
    if isinstance(cached, dict) and "per_region" in cached:
        region_bns = cached["per_region"]
    elif isinstance(cached, dict):
        region_bns = cached
    else:
        log(f"Unexpected format in {bn_results_path}")
        sys.exit(1)

    matched_regions = []
    for region, effect in SILENCING_EFFECTS.items():
        if region in region_bns:
            entry = region_bns[region]
            matched_regions.append({
                "region": region,
                "silencing_effect": effect,
                "raw_bn": entry.get("raw_bn", entry.get("bracket_norm", 0)),
                "n_neurons": entry.get("n_neurons", entry.get("neuron_count", 1)),
            })

    if len(matched_regions) < 6:
        log(f"Only {len(matched_regions)} regions matched — insufficient for test")
        sys.exit(1)

    raw_bns = np.array([r["raw_bn"] for r in matched_regions])
    n_neurons = np.array([r["n_neurons"] for r in matched_regions])
    silencing = np.array([r["silencing_effect"] for r in matched_regions])
    corrected_bns = raw_bns / np.sqrt(n_neurons)

    rho_raw, p_raw = stats.spearmanr(raw_bns, silencing)
    rho_corr, p_corr = stats.spearmanr(corrected_bns, silencing)
    rho_raw_n, _ = stats.spearmanr(raw_bns, n_neurons)
    rho_corr_n, _ = stats.spearmanr(corrected_bns, n_neurons)

    log(f"\n=== H6 RESULTS ({len(matched_regions)} matched regions) ===")
    log(f"Raw BN vs silencing:        rho={rho_raw:.3f}, p={p_raw:.3e}")
    log(f"BN/sqrt(n) vs silencing:    rho={rho_corr:.3f}, p={p_corr:.3e}")
    log(f"Raw BN vs neuron count:     rho={rho_raw_n:.3f}")
    log(f"BN/sqrt(n) vs neuron count: rho={rho_corr_n:.3f}")

    log(f"\nPer-region breakdown:")
    for r in sorted(matched_regions, key=lambda x: x["silencing_effect"], reverse=True):
        bn_corr = r["raw_bn"] / np.sqrt(r["n_neurons"])
        log(f"  {r['region']:8s} silencing={r['silencing_effect']:.3f} "
            f"raw_bn={r['raw_bn']:.4f} n={r['n_neurons']:3d} "
            f"bn/sqrt(n)={bn_corr:.5f}")

    if rho_corr > 0.6 and p_corr < 0.05:
        log(f"\nH6 CONFIRMED: corrected bracket norm predicts silencing (rho={rho_corr:.3f})")
    elif rho_raw > 0.5:
        log(f"\nH6 PARTIAL: raw bracket correlates but correction doesn't improve")
    else:
        log(f"\nH6 NOT CONFIRMED")

    output_dir.mkdir(parents=True, exist_ok=True)
    results = {
        "n_matched_regions": len(matched_regions),
        "rho_raw_vs_silencing": float(rho_raw),
        "p_raw_vs_silencing": float(p_raw),
        "rho_corrected_vs_silencing": float(rho_corr),
        "p_corrected_vs_silencing": float(p_corr),
        "rho_raw_vs_neuron_count": float(rho_raw_n),
        "rho_corrected_vs_neuron_count": float(rho_corr_n),
        "regions": matched_regions,
    }
    with open(output_dir / "neural_replication_results.json", "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--synthetic", action="store_true")
    mode.add_argument("--real", action="store_true")
    parser.add_argument("--data", type=Path,
                        help="Path to bracket-norm/data/cache/steinmetz or results_bundle/batch3")
    parser.add_argument("--output", type=Path, default=Path("results/06_neural_replication"))
    args = parser.parse_args()

    if args.synthetic:
        run_synthetic()
    elif args.real:
        if not args.data:
            parser.error("--real requires --data")
        run_real(args.data, args.output)
