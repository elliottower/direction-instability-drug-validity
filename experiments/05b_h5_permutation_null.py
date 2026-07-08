"""H5 permutation null: shuffle drug labels within each fold.

Tests whether the transport-stable > raw pattern survives under the null
hypothesis that drug identity is unrelated to held-out cosine. If TS wins
are common under the null, the algebraic coupling between Frechet variance
and held-out cosine could explain the result.

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
    log(f"H5 permutation null ({n_perm} permutations)")

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
        fold_drugs = []
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
            holdout_cos = cosine_sim(holdout_sig, consensus)
            fold_drugs.append((float(raw), float(ts), holdout_cos))

        if len(fold_drugs) >= 20:
            fold_data.append({
                "holdout_cell": holdout_cell,
                "raw": [d[0] for d in fold_drugs],
                "ts": [d[1] for d in fold_drugs],
                "holdout": [d[2] for d in fold_drugs],
            })

    log(f"  {len(fold_data)} valid folds")

    observed_ts_wins = 0
    for fd in fold_data:
        raw_rho = stats.spearmanr(fd["raw"], fd["holdout"]).statistic
        ts_rho = stats.spearmanr(fd["ts"], fd["holdout"]).statistic
        if ts_rho > raw_rho:
            observed_ts_wins += 1

    log(f"  Observed TS wins: {observed_ts_wins}/{len(fold_data)}")

    rng = np.random.default_rng(seed=2026070705)
    null_ts_wins_distribution = []

    for perm_i in tqdm(range(n_perm), desc="Permutations"):
        perm_ts_wins = 0
        for fd in fold_data:
            n = len(fd["holdout"])
            shuffled_holdout = list(fd["holdout"])
            rng.shuffle(shuffled_holdout)
            raw_rho = stats.spearmanr(fd["raw"], shuffled_holdout).statistic
            ts_rho = stats.spearmanr(fd["ts"], shuffled_holdout).statistic
            if ts_rho > raw_rho:
                perm_ts_wins += 1
        null_ts_wins_distribution.append(perm_ts_wins)

    null_ts_wins_distribution = np.array(null_ts_wins_distribution)
    p_value = float(np.mean(null_ts_wins_distribution >= observed_ts_wins))
    null_mean = float(np.mean(null_ts_wins_distribution))
    null_std = float(np.std(null_ts_wins_distribution))

    result = {
        "observed_ts_wins": observed_ts_wins,
        "n_folds": len(fold_data),
        "n_permutations": n_perm,
        "null_mean_ts_wins": round(null_mean, 2),
        "null_std_ts_wins": round(null_std, 2),
        "p_value": p_value,
        "null_distribution_quantiles": {
            "q05": float(np.quantile(null_ts_wins_distribution, 0.05)),
            "q50": float(np.quantile(null_ts_wins_distribution, 0.50)),
            "q95": float(np.quantile(null_ts_wins_distribution, 0.95)),
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "permutation_null.json", "w") as f:
        json.dump(result, f, indent=2)

    log(f"\n=== H5 PERMUTATION NULL ===")
    log(f"Observed TS wins: {observed_ts_wins}/{len(fold_data)}")
    log(f"Null distribution: mean={null_mean:.2f}, std={null_std:.2f}")
    log(f"p-value: {p_value:.4f}")
    if p_value < 0.05:
        log("PASS: TS advantage is not explained by algebraic coupling")
    else:
        log("FAIL: TS advantage could be explained by algebraic coupling")

    log(f"Saved to {output_dir / 'permutation_null.json'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).parent.parent.parent / "drug-perturbation-geometry" / "data")
    parser.add_argument("--n-perm", type=int, default=1000)
    args = parser.parse_args()
    output_dir = Path(__file__).parent.parent / "results" / "05_transport_stability"
    run_permutation(args.data, output_dir, args.n_perm)
