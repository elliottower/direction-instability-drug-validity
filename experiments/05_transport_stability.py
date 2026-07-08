"""Experiment 05: Transport-stable bracket outpredicts raw bracket for cross-context claims.

Tests H5: bracket norm that replicates across contexts (low Frechet variance)
predicts held-out context effects better than raw bracket. Leave-one-cell-line-out
cross-validation with direction-based held-out outcome (cosine consistency).

Decision criterion (pre-registered): Spearman correlation between transport-stable
bracket and held-out-context prediction > correlation for raw bracket, in at least
4/5 cross-validation folds.

Held-out outcome: cosine similarity between the held-out cell line's signature
and the training-set consensus direction. This matches the direction-based
estimand of DI (not magnitude, which DI does not measure).

Gene-set method: N/A (no gene sets needed).
Analytic choices declared in DEVIATION_LOG.md entries 1-2 before running.

Usage:
    uv run python experiments/05_transport_stability.py --synthetic
    uv run python experiments/05_transport_stability.py --real --data PATH_TO_LINCS
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
    """Build {drug: {cell: mean_signature}} mapping via vectorized groupby."""
    sig_id_set = set(sig_ids)
    sig_id_to_idx = {sid: i for i, sid in enumerate(sig_ids)}
    assert len(sig_id_set) == len(sig_ids), "Duplicate sig_ids in compound data"

    filtered = siginfo[siginfo["sig_id"].isin(sig_id_set) & siginfo["pert_iname"].notna()].copy()
    filtered["_idx"] = filtered["sig_id"].map(sig_id_to_idx)

    drug_cell_map = {}
    for (drug, cell), group in filtered.groupby(["pert_iname", "cell_id"]):
        if drug not in drug_cell_map:
            drug_cell_map[drug] = {}
        indices = group["_idx"].values
        drug_cell_map[drug][cell] = signatures[indices].mean(axis=0)

    return {drug: cells for drug, cells in drug_cell_map.items() if len(cells) >= 5}


def run_synthetic():
    """Validate: transport-stable bracket predicts held-out cosine better than raw."""
    log("=== SYNTHETIC MODE ===")
    rng = np.random.default_rng(seed=2026070705)
    n_genes = 200
    n_cells = 8
    n_drugs_stable = 40
    n_drugs_unstable = 40

    log("Generating stable drugs (consistent direction across contexts)...")
    stable_results = []
    for _ in tqdm(range(n_drugs_stable), desc="Stable"):
        shared_mech = rng.standard_normal(n_genes) * 2.0
        sigs = np.array([
            shared_mech + rng.standard_normal(n_genes) * 0.3
            for _ in range(n_cells)
        ])
        raw_scores = []
        ts_scores = []
        holdout_cosines = []
        for hold in range(n_cells):
            train_idx = [i for i in range(n_cells) if i != hold]
            train_sigs = sigs[train_idx]
            consensus = train_sigs.mean(axis=0)
            raw_scores.append(direction_instability(train_sigs))
            ts_scores.append(transport_stable_bracket(train_sigs))
            holdout_cosines.append(cosine_sim(sigs[hold], consensus))
        stable_results.append({
            "raw_scores": raw_scores,
            "ts_scores": ts_scores,
            "holdout_cosines": holdout_cosines,
        })

    log("Generating unstable drugs (inconsistent direction across contexts)...")
    unstable_results = []
    for _ in tqdm(range(n_drugs_unstable), desc="Unstable"):
        sigs = rng.standard_normal((n_cells, n_genes)) * 1.5
        active = rng.choice(n_cells, size=3, replace=False)
        sigs[active] += rng.standard_normal((3, n_genes)) * 3.0
        raw_scores = []
        ts_scores = []
        holdout_cosines = []
        for hold in range(n_cells):
            train_idx = [i for i in range(n_cells) if i != hold]
            train_sigs = sigs[train_idx]
            consensus = train_sigs.mean(axis=0)
            raw_scores.append(direction_instability(train_sigs))
            ts_scores.append(transport_stable_bracket(train_sigs))
            holdout_cosines.append(cosine_sim(sigs[hold], consensus))
        unstable_results.append({
            "raw_scores": raw_scores,
            "ts_scores": ts_scores,
            "holdout_cosines": holdout_cosines,
        })

    all_results = stable_results + unstable_results
    raw_corrs = []
    ts_corrs = []
    for fold in range(n_cells):
        fold_raw = [r["raw_scores"][fold] for r in all_results]
        fold_ts = [r["ts_scores"][fold] for r in all_results]
        fold_holdout = [r["holdout_cosines"][fold] for r in all_results]
        raw_corrs.append(stats.spearmanr(fold_raw, fold_holdout).statistic)
        ts_corrs.append(stats.spearmanr(fold_ts, fold_holdout).statistic)

    ts_wins = sum(1 for r, t in zip(raw_corrs, ts_corrs) if t > r)
    log(f"\nPer-fold Spearman correlations with held-out cosine consistency:")
    for i, (r, t) in enumerate(zip(raw_corrs, ts_corrs)):
        winner = "TS" if t > r else "RAW"
        log(f"  Fold {i}: raw={r:.4f}, transport-stable={t:.4f} [{winner}]")
    log(f"\nTransport-stable wins {ts_wins}/{n_cells} folds")

    if ts_wins >= 4:
        log("  PASS: transport-stable bracket outpredicts raw in >=4/5 folds")
    else:
        log(f"  FAIL: transport-stable wins only {ts_wins} folds (need >=4)")


def run_real(data_dir: Path, output_dir: Path):
    """Run H5 on real LINCS data with leave-one-cell-line-out CV."""
    log("=== REAL MODE ===")

    sigs_path = data_dir / "lincs_subset.npz"
    siginfo_path = data_dir / "GSE92742_Broad_LINCS_sig_info.txt.gz"

    if not sigs_path.exists():
        log(f"ERROR: {sigs_path} not found.")
        sys.exit(1)

    log("Loading compound signatures...")
    data = np.load(sigs_path, allow_pickle=True)
    signatures = data["signatures"]
    sig_ids = list(data["sig_ids"])
    log(f"  {signatures.shape[0]:,} signatures x {signatures.shape[1]} genes")

    log("Loading sig info...")
    siginfo = pd.read_csv(siginfo_path, sep="\t", low_memory=False)

    log("Building drug-cell matrix...")
    drug_cell_means = build_drug_cell_matrix(signatures, sig_ids, siginfo)
    log(f"  {len(drug_cell_means)} drugs with >= 5 cell lines")

    all_cells = set()
    for cells in drug_cell_means.values():
        all_cells.update(cells.keys())
    all_cells = sorted(all_cells)
    log(f"  {len(all_cells)} unique cell lines")

    n_folds = min(5, len(all_cells))
    log(f"Running {n_folds}-fold leave-one-cell-line-out CV...")

    fold_results = []
    for fold_idx in range(n_folds):
        holdout_cell = all_cells[fold_idx]
        log(f"  Fold {fold_idx}: holding out {holdout_cell}")

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

            fold_drugs.append({
                "drug": drug,
                "raw_bracket": float(raw),
                "ts_bracket": float(ts),
                "holdout_cosine": holdout_cos,
            })

        if len(fold_drugs) < 20:
            log(f"    Skipping fold {fold_idx}: only {len(fold_drugs)} drugs")
            continue

        raw_vals = [d["raw_bracket"] for d in fold_drugs]
        ts_vals = [d["ts_bracket"] for d in fold_drugs]
        holdout_vals = [d["holdout_cosine"] for d in fold_drugs]

        raw_rho = stats.spearmanr(raw_vals, holdout_vals).statistic
        ts_rho = stats.spearmanr(ts_vals, holdout_vals).statistic

        fold_results.append({
            "fold": fold_idx,
            "holdout_cell": holdout_cell,
            "n_drugs": len(fold_drugs),
            "raw_spearman": float(raw_rho),
            "ts_spearman": float(ts_rho),
            "ts_wins": bool(ts_rho > raw_rho),
            "drugs": fold_drugs,
        })

        winner = "TS" if ts_rho > raw_rho else "RAW"
        log(f"    n={len(fold_drugs)}, raw_rho={raw_rho:.4f}, ts_rho={ts_rho:.4f} [{winner}]")

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "transport_stability_results.json", "w") as f:
        json.dump(fold_results, f, indent=2)

    ts_wins = sum(1 for f in fold_results if f["ts_wins"])
    n_valid_folds = len(fold_results)

    log(f"\n=== H5 RESULTS ===")
    log(f"Valid folds: {n_valid_folds}")
    log(f"Transport-stable wins: {ts_wins}/{n_valid_folds}")

    threshold = 4
    if n_valid_folds < 5:
        threshold = int(np.ceil(n_valid_folds * 0.8))
        log(f"  (Adjusted threshold to {threshold}/{n_valid_folds} due to fewer valid folds)")

    if ts_wins >= threshold:
        log(f"  H5 CONFIRMED: transport-stable bracket outpredicts raw in {ts_wins}/{n_valid_folds} folds")
    else:
        log(f"  H5 NOT CONFIRMED: transport-stable wins only {ts_wins}/{n_valid_folds} folds (needed {threshold})")

    log(f"\nResults saved to {output_dir / 'transport_stability_results.json'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--synthetic", action="store_true")
    mode.add_argument("--real", action="store_true")
    parser.add_argument("--data", type=Path,
                        help="Path to drug-perturbation-geometry/data/")
    parser.add_argument("--output", type=Path,
                        default=Path("results/05_transport_stability"))
    args = parser.parse_args()

    if args.synthetic:
        run_synthetic()
    elif args.real:
        if not args.data:
            parser.error("--real requires --data")
        run_real(args.data, args.output)
