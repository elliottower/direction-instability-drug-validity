"""Robustness check R2: Lambda sensitivity for transport-stable bracket.

Tests whether the H5-full result (TS outpredicts raw in all 66 folds) is
robust to the choice of frechet_penalty (lambda) in:

    TS(sigs) = DI(sigs) - lambda * frechet_var(sigs)

The main result uses lambda=1.0. This sweeps [0.25, 0.5, 0.75, 1.0, 1.5,
2.0, 3.0, 5.0] and records fold wins, Wilcoxon p-value, and mean delta_rho
at each value.

MUST be pre-registered before running on real data.

Usage:
    uv run python experiments/r2_lambda_sensitivity.py --real --data ../drug-perturbation-geometry/data/
    uv run python experiments/r2_lambda_sensitivity.py --synthetic
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


LAMBDA_VALUES = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]


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


def compute_per_drug_scores(drug_cell_means, holdout_cell, lam):
    """For one fold and one lambda, compute per-drug raw/TS/holdout scores."""
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
        ts = transport_stable_bracket(train_sigs, frechet_penalty=lam)
        holdout_cos = cosine_sim(holdout_sig, consensus)

        fold_drugs.append({
            "raw_bracket": float(raw),
            "ts_bracket": float(ts),
            "holdout_cosine": holdout_cos,
        })
    return fold_drugs


def run_real(data_dir: Path, output_dir: Path):
    """Run lambda sensitivity on real LINCS data."""
    log("=== REAL MODE: LAMBDA SENSITIVITY (R2) ===")

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

    # Pre-compute raw DI and holdout cosines (lambda-independent)
    # to avoid redundant computation across lambda values.
    log("Pre-computing raw DI and holdout cosines for all folds...")
    fold_data = {}  # holdout_cell -> list of {raw, holdout_cos, train_sigs}
    valid_cells = []
    for holdout_cell in tqdm(all_cells, desc="Pre-computing folds"):
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
            holdout_cos = cosine_sim(holdout_sig, consensus)

            # Pre-compute the Frechet variance (lambda-independent part of TS)
            norms = np.linalg.norm(train_sigs, axis=1, keepdims=True)
            norms = np.maximum(norms, 1e-10)
            unit = train_sigs / norms
            frechet_mean = unit.mean(axis=0)
            frechet_mean = frechet_mean / (np.linalg.norm(frechet_mean) + 1e-10)
            deviations = np.arccos(np.clip(unit @ frechet_mean, -1, 1))
            frechet_var = float(np.mean(deviations**2))

            fold_drugs.append({
                "raw_bracket": float(raw),
                "frechet_var": frechet_var,
                "holdout_cosine": holdout_cos,
            })

        if len(fold_drugs) >= 20:
            fold_data[holdout_cell] = fold_drugs
            valid_cells.append(holdout_cell)
        else:
            log(f"    Skipping {holdout_cell}: only {len(fold_drugs)} drugs")

    log(f"  {len(valid_cells)} valid folds (>= 20 drugs)")

    # Sweep lambda values using pre-computed components
    log(f"\nSweeping {len(LAMBDA_VALUES)} lambda values...")
    lambda_results = []

    for lam in tqdm(LAMBDA_VALUES, desc="Lambda sweep"):
        fold_results = []
        for holdout_cell in valid_cells:
            drugs = fold_data[holdout_cell]

            raw_vals = [d["raw_bracket"] for d in drugs]
            # TS = raw - lambda * frechet_var
            ts_vals = [d["raw_bracket"] - lam * d["frechet_var"] for d in drugs]
            holdout_vals = [d["holdout_cosine"] for d in drugs]

            raw_rho = stats.spearmanr(raw_vals, holdout_vals).statistic
            ts_rho = stats.spearmanr(ts_vals, holdout_vals).statistic

            fold_results.append({
                "holdout_cell": holdout_cell,
                "n_drugs": len(drugs),
                "raw_spearman": float(raw_rho),
                "ts_spearman": float(ts_rho),
                "ts_wins": bool(ts_rho > raw_rho),
            })

        ts_wins = sum(1 for f in fold_results if f["ts_wins"])
        n_valid = len(fold_results)
        frac = ts_wins / n_valid if n_valid > 0 else 0

        raw_rhos = [f["raw_spearman"] for f in fold_results]
        ts_rhos = [f["ts_spearman"] for f in fold_results]
        improvements = np.array([t - r for r, t in zip(raw_rhos, ts_rhos)])

        wilcoxon_result = stats.wilcoxon(improvements, alternative="greater")

        result = {
            "lambda": lam,
            "n_valid_folds": n_valid,
            "ts_wins": ts_wins,
            "win_fraction": float(frac),
            "wilcoxon_statistic": float(wilcoxon_result.statistic),
            "wilcoxon_p": float(wilcoxon_result.pvalue),
            "mean_delta_rho": float(np.mean(improvements)),
            "std_delta_rho": float(np.std(improvements)),
            "mean_raw_rho": float(np.mean(raw_rhos)),
            "mean_ts_rho": float(np.mean(ts_rhos)),
            "per_fold": fold_results,
        }
        lambda_results.append(result)

        log(f"  lambda={lam:.2f}: wins={ts_wins}/{n_valid} ({frac:.0%}), "
            f"p={wilcoxon_result.pvalue:.2e}, mean_delta={np.mean(improvements):.4f}")

    # Sanity check: at lambda=0, TS ≡ raw, so win fraction should be ~50%
    lam0_result = next((r for r in lambda_results if r["lambda"] == 0.0), None)
    if lam0_result is not None:
        lam0_frac = lam0_result["win_fraction"]
        if lam0_frac >= 0.80:
            log(f"WARNING: lambda=0 shows TS winning {lam0_frac:.0%} of folds. "
                f"At lambda=0, TS≡raw, so this indicates a bug (ties broken in TS's favor).")
        else:
            log(f"Sanity check PASSED: lambda=0 win fraction = {lam0_frac:.0%} (expected ~50%)")

    # Save results
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "lambda_sensitivity.json"
    with open(output_file, "w") as f:
        json.dump(lambda_results, f, indent=2)
    log(f"\nResults saved to {output_file}")

    # Print summary table
    log("\n" + "=" * 85)
    log("LAMBDA SENSITIVITY SUMMARY")
    log("=" * 85)
    log(f"{'lambda':>8}  {'wins':>6}  {'frac':>6}  {'Wilcoxon p':>12}  "
        f"{'mean_delta':>12}  {'std_delta':>10}  {'mean_TS_rho':>12}")
    log("-" * 85)
    for r in lambda_results:
        log(f"{r['lambda']:>8.2f}  "
            f"{r['ts_wins']:>3}/{r['n_valid_folds']:<3}  "
            f"{r['win_fraction']:>5.0%}  "
            f"{r['wilcoxon_p']:>12.2e}  "
            f"{r['mean_delta_rho']:>12.4f}  "
            f"{r['std_delta_rho']:>10.4f}  "
            f"{r['mean_ts_rho']:>12.4f}")
    log("=" * 85)


def run_synthetic():
    """Quick validation on synthetic data to verify the sweep logic."""
    log("=== SYNTHETIC MODE: LAMBDA SENSITIVITY ===")
    rng = np.random.default_rng(seed=2026070805)
    n_genes = 200
    n_cells = 15
    n_drugs = 80

    log(f"Generating {n_drugs} drugs across {n_cells} cell lines...")
    all_sigs = {}
    for d in range(n_drugs):
        shared = rng.standard_normal(n_genes) * (2.0 if d < 40 else 0.5)
        noise_scale = 0.3 if d < 40 else 1.5
        sigs = {}
        for c in range(n_cells):
            sigs[f"cell_{c}"] = shared + rng.standard_normal(n_genes) * noise_scale
        all_sigs[f"drug_{d}"] = sigs

    cell_names = [f"cell_{c}" for c in range(n_cells)]

    log(f"\nSweeping {len(LAMBDA_VALUES)} lambda values...")
    for lam in LAMBDA_VALUES:
        fold_results = []
        for holdout_cell in cell_names:
            fold_drugs = []
            for drug, cells in all_sigs.items():
                if holdout_cell not in cells:
                    continue
                train_cells = {c: s for c, s in cells.items() if c != holdout_cell}
                if len(train_cells) < 4:
                    continue
                train_sigs = np.array(list(train_cells.values()))
                holdout_sig = cells[holdout_cell]
                consensus = train_sigs.mean(axis=0)
                fold_drugs.append({
                    "raw_bracket": float(direction_instability(train_sigs)),
                    "ts_bracket": float(transport_stable_bracket(train_sigs, frechet_penalty=lam)),
                    "holdout_cosine": cosine_sim(holdout_sig, consensus),
                })

            if len(fold_drugs) < 20:
                continue

            raw_vals = [d["raw_bracket"] for d in fold_drugs]
            ts_vals = [d["ts_bracket"] for d in fold_drugs]
            holdout_vals = [d["holdout_cosine"] for d in fold_drugs]
            raw_rho = stats.spearmanr(raw_vals, holdout_vals).statistic
            ts_rho = stats.spearmanr(ts_vals, holdout_vals).statistic
            fold_results.append({"ts_wins": bool(ts_rho > raw_rho),
                                 "delta": ts_rho - raw_rho})

        ts_wins = sum(1 for f in fold_results if f["ts_wins"])
        n_valid = len(fold_results)
        frac = ts_wins / n_valid if n_valid > 0 else 0
        deltas = np.array([f["delta"] for f in fold_results])
        mean_d = float(np.mean(deltas)) if len(deltas) > 0 else 0.0
        log(f"  lambda={lam:.2f}: wins={ts_wins}/{n_valid} ({frac:.0%}), mean_delta={mean_d:.4f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--synthetic", action="store_true")
    mode.add_argument("--real", action="store_true")
    parser.add_argument("--data", type=Path,
                        help="Path to drug-perturbation-geometry/data/")
    parser.add_argument("--output", type=Path,
                        default=Path("results/r2_lambda_sensitivity"))
    args = parser.parse_args()

    if args.synthetic:
        run_synthetic()
    elif args.real:
        if not args.data:
            parser.error("--real requires --data")
        run_real(args.data, args.output)
