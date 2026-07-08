"""H5 permutation null: variance-preserving signature shuffle.

Tests whether the transport-stable bracket's advantage over raw DI in
predicting held-out cell-line cosine is an algebraic artifact of the
coupling between Frechet variance and held-out outcome.

The null shuffles held-out SIGNATURES across drugs while holding each
drug's training data (and therefore its raw DI, Frechet variance, TS,
and consensus vector) fixed. Cosines are recomputed against each drug's
original consensus. This preserves the mechanical coupling between
consensus tightness and cosine magnitude -- if TS's advantage survives
this null, it's algebraic; if it collapses, the advantage is genuine.

Contrast with the naive label-shuffle (shuffling pre-computed cosines),
which destroys the coupling entirely and tests the wrong hypothesis.

Output: results/05_transport_stability/permutation_null.json
"""
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from geometry.bracket_norm import direction_instability, transport_stable_bracket


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def cosine_sim(a: np.ndarray, b: np.ndarray) -> float:
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-10 or nb < 1e-10:
        return 0.0
    return float(a @ b / (na * nb))


def build_drug_cell_matrix(signatures, sig_ids, siginfo):
    sig_id_set = set(sig_ids)
    sig_id_to_idx = {sid: i for i, sid in enumerate(sig_ids)}
    filtered = siginfo[siginfo["sig_id"].isin(sig_id_set) & siginfo["pert_iname"].notna()].copy()
    filtered["_idx"] = filtered["sig_id"].map(sig_id_to_idx)
    drug_cell_map = {}
    for (drug, cell), group in filtered.groupby(["pert_iname", "cell_id"]):
        if drug not in drug_cell_map:
            drug_cell_map[drug] = {}
        indices = group["_idx"].values
        drug_cell_map[drug][cell] = signatures[indices].mean(axis=0)
    return {drug: cells for drug, cells in drug_cell_map.items() if len(cells) >= 5}


def run_permutation(data_dir: Path, output_dir: Path, n_perm: int = 1000):
    log(f"H5 variance-preserving permutation null ({n_perm} permutations)")

    data = np.load(data_dir / "lincs_subset.npz", allow_pickle=True)
    signatures = data["signatures"]
    sig_ids = list(data["sig_ids"])
    log(f"  {signatures.shape[0]:,} signatures x {signatures.shape[1]} genes")

    siginfo = pd.read_csv(data_dir / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t", low_memory=False)
    drug_cell_means = build_drug_cell_matrix(signatures, sig_ids, siginfo)
    log(f"  {len(drug_cell_means)} drugs")

    all_cells = sorted(set(c for cells in drug_cell_means.values() for c in cells))
    n_folds = min(5, len(all_cells))

    fold_data = []
    for fold_idx in range(n_folds):
        holdout_cell = all_cells[fold_idx]
        fold_raw = []
        fold_ts = []
        fold_holdout_sigs = []
        fold_consensus_vecs = []

        for drug, cells in drug_cell_means.items():
            if holdout_cell not in cells:
                continue
            train_cells = {c: s for c, s in cells.items() if c != holdout_cell}
            if len(train_cells) < 4:
                continue
            train_sigs = np.array(list(train_cells.values()))
            holdout_sig = cells[holdout_cell]
            consensus = train_sigs.mean(axis=0)
            raw = direction_instability(train_sigs)
            ts = transport_stable_bracket(train_sigs)

            fold_raw.append(float(raw))
            fold_ts.append(float(ts))
            fold_holdout_sigs.append(holdout_sig)
            fold_consensus_vecs.append(consensus)

        if len(fold_raw) >= 20:
            holdout_cos = [
                cosine_sim(fold_holdout_sigs[i], fold_consensus_vecs[i])
                for i in range(len(fold_raw))
            ]
            fold_data.append({
                "holdout_cell": holdout_cell,
                "raw": fold_raw,
                "ts": fold_ts,
                "holdout_sigs": np.array(fold_holdout_sigs),
                "consensus_vecs": np.array(fold_consensus_vecs),
                "holdout_cos": holdout_cos,
            })

    log(f"  {len(fold_data)} valid folds")

    observed_ts_wins = 0
    observed_ts_deltas = []
    for fd in fold_data:
        raw_rho = stats.spearmanr(fd["raw"], fd["holdout_cos"]).statistic
        ts_rho = stats.spearmanr(fd["ts"], fd["holdout_cos"]).statistic
        observed_ts_deltas.append(ts_rho - raw_rho)
        if ts_rho > raw_rho:
            observed_ts_wins += 1

    observed_mean_delta = float(np.mean(observed_ts_deltas))
    log(f"  Observed TS wins: {observed_ts_wins}/{len(fold_data)}")
    log(f"  Observed mean rho delta (TS - raw): {observed_mean_delta:.4f}")

    rng = np.random.default_rng(seed=2026070705)
    null_ts_wins_distribution = []
    null_mean_delta_distribution = []

    for perm_i in tqdm(range(n_perm), desc="Permutations"):
        perm_ts_wins = 0
        perm_deltas = []
        for fd in fold_data:
            n = len(fd["raw"])
            perm_idx = rng.permutation(n)
            shuffled_holdout_sigs = fd["holdout_sigs"][perm_idx]
            perm_holdout_cos = [
                cosine_sim(shuffled_holdout_sigs[i], fd["consensus_vecs"][i])
                for i in range(n)
            ]
            raw_rho = stats.spearmanr(fd["raw"], perm_holdout_cos).statistic
            ts_rho = stats.spearmanr(fd["ts"], perm_holdout_cos).statistic
            perm_deltas.append(ts_rho - raw_rho)
            if ts_rho > raw_rho:
                perm_ts_wins += 1
        null_ts_wins_distribution.append(perm_ts_wins)
        null_mean_delta_distribution.append(float(np.mean(perm_deltas)))

    null_ts_wins_distribution = np.array(null_ts_wins_distribution)
    null_mean_delta_distribution = np.array(null_mean_delta_distribution)

    p_value_wins = float(np.mean(null_ts_wins_distribution >= observed_ts_wins))
    p_value_delta = float(np.mean(null_mean_delta_distribution >= observed_mean_delta))
    null_wins_mean = float(np.mean(null_ts_wins_distribution))
    null_wins_std = float(np.std(null_ts_wins_distribution))
    null_delta_mean = float(np.mean(null_mean_delta_distribution))
    null_delta_std = float(np.std(null_mean_delta_distribution))

    result = {
        "test": "variance_preserving_signature_shuffle",
        "description": (
            "Shuffles held-out signatures across drugs while holding "
            "training data (raw DI, Frechet variance, TS, consensus) fixed. "
            "Cosines recomputed against original consensus to preserve the "
            "mechanical coupling between consensus tightness and cosine magnitude."
        ),
        "observed_ts_wins": observed_ts_wins,
        "observed_mean_rho_delta": round(observed_mean_delta, 4),
        "n_folds": len(fold_data),
        "n_permutations": n_perm,
        "wins_null_mean": round(null_wins_mean, 2),
        "wins_null_std": round(null_wins_std, 2),
        "wins_p_value": p_value_wins,
        "delta_null_mean": round(null_delta_mean, 4),
        "delta_null_std": round(null_delta_std, 4),
        "delta_p_value": p_value_delta,
        "wins_null_quantiles": {
            "q05": float(np.quantile(null_ts_wins_distribution, 0.05)),
            "q50": float(np.quantile(null_ts_wins_distribution, 0.50)),
            "q95": float(np.quantile(null_ts_wins_distribution, 0.95)),
        },
        "delta_null_quantiles": {
            "q025": round(float(np.quantile(null_mean_delta_distribution, 0.025)), 4),
            "q975": round(float(np.quantile(null_mean_delta_distribution, 0.975)), 4),
        },
        "interpretation": (
            "If p_value_wins < 0.05 AND p_value_delta < 0.05: "
            "TS advantage persists under variance-preserving null → "
            "advantage is algebraic (Frechet variance mechanically predicts "
            "held-out cosine). "
            "If p_value_wins >= 0.05 OR p_value_delta >= 0.05: "
            "TS advantage collapses when drug-specific signal is broken → "
            "advantage is genuine, not an artifact of the penalization formula."
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "permutation_null.json", "w") as f:
        json.dump(result, f, indent=2)

    log(f"\n=== H5 VARIANCE-PRESERVING PERMUTATION NULL ===")
    log(f"Observed TS wins: {observed_ts_wins}/{len(fold_data)}")
    log(f"Observed mean rho delta: {observed_mean_delta:.4f}")
    log(f"Null wins distribution: mean={null_wins_mean:.2f}, std={null_wins_std:.2f}")
    log(f"Null delta distribution: mean={null_delta_mean:.4f}, std={null_delta_std:.4f}")
    log(f"p-value (wins): {p_value_wins:.4f}")
    log(f"p-value (delta): {p_value_delta:.4f}")

    if p_value_wins >= 0.05 or p_value_delta >= 0.05:
        log("PASS: TS advantage collapses under variance-preserving null")
        log("  → advantage is drug-specific, not an algebraic artifact")
    else:
        log("CAUTION: TS advantage persists under variance-preserving null")
        log("  → advantage may reflect algebraic coupling, not genuine signal")

    log(f"Saved to {output_dir / 'permutation_null.json'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).parent.parent.parent / "drug-perturbation-geometry" / "data")
    parser.add_argument("--n-perm", type=int, default=1000)
    args = parser.parse_args()
    output_dir = Path(__file__).parent.parent / "results" / "05_transport_stability"
    run_permutation(args.data, output_dir, args.n_perm)
