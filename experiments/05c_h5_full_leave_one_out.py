"""Experiment 05c: Full leave-one-cell-line-out cross-validation (H5-full).

Extends the original H5 (5 alphabetical folds) to ALL cell lines with >= 20
eligible drugs. Tests whether the transport-stable advantage generalizes
beyond the alphabetical selection and whether the rare/common cell-line
regime distinction holds across the full set.

Pre-registered in PREREGISTRATION_EXTENDED.md before running on real data.

Decision criterion: TS outpredicts raw in >= 80% of valid folds.

Usage:
    uv run python experiments/05c_h5_full_leave_one_out.py --synthetic
    uv run python experiments/05c_h5_full_leave_one_out.py --real --data PATH_TO_LINCS
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
    """Validate: full leave-one-out should show TS advantage in most folds."""
    log("=== SYNTHETIC MODE ===")
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

    fold_results = []
    cell_names = [f"cell_{c}" for c in range(n_cells)]
    for fold_idx, holdout_cell in enumerate(tqdm(cell_names, desc="Folds")):
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
                "drug": drug,
                "raw_bracket": float(direction_instability(train_sigs)),
                "ts_bracket": float(transport_stable_bracket(train_sigs)),
                "holdout_cosine": cosine_sim(holdout_sig, consensus),
            })

        if len(fold_drugs) < 20:
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
        })
        winner = "TS" if ts_rho > raw_rho else "RAW"
        log(f"  {holdout_cell}: n={len(fold_drugs)}, raw={raw_rho:.4f}, ts={ts_rho:.4f} [{winner}]")

    ts_wins = sum(1 for f in fold_results if f["ts_wins"])
    n_valid = len(fold_results)
    frac = ts_wins / n_valid if n_valid > 0 else 0
    log(f"\nTS wins {ts_wins}/{n_valid} folds ({frac:.1%})")
    log(f"Criterion (>= 80%): {'PASS' if frac >= 0.8 else 'FAIL'}")


def run_real(data_dir: Path, output_dir: Path):
    """Run H5-full on real LINCS data: all cell lines, not just first 5."""
    log("=== REAL MODE: FULL LEAVE-ONE-OUT ===")

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
    log(f"  {len(all_cells)} unique cell lines (ALL will be tested)")

    fold_results = []
    for fold_idx, holdout_cell in enumerate(tqdm(all_cells, desc="Leave-one-out folds")):
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
            log(f"    Skipping {holdout_cell}: only {len(fold_drugs)} drugs")
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
        })

        winner = "TS" if ts_rho > raw_rho else "RAW"
        log(f"    {holdout_cell}: n={len(fold_drugs)}, raw={raw_rho:.4f}, ts={ts_rho:.4f} [{winner}]")

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "h5_full_leave_one_out_results.json"
    with open(output_file, "w") as f:
        json.dump(fold_results, f, indent=2)

    ts_wins = sum(1 for f in fold_results if f["ts_wins"])
    n_valid = len(fold_results)
    frac = ts_wins / n_valid if n_valid > 0 else 0

    n_drugs_per_fold = [f["n_drugs"] for f in fold_results]
    median_n = np.median(n_drugs_per_fold) if n_drugs_per_fold else 0

    raw_rhos = [f["raw_spearman"] for f in fold_results]
    ts_rhos = [f["ts_spearman"] for f in fold_results]
    improvements = [t - r for r, t in zip(raw_rhos, ts_rhos)]
    imp_arr = np.array(improvements)

    wilcoxon_result = stats.wilcoxon(imp_arr, alternative="greater")

    log(f"\n=== H5-FULL RESULTS ===")
    log(f"Valid folds: {n_valid} / {len(all_cells)} cell lines")
    log(f"Excluded folds (< 20 drugs): {len(all_cells) - n_valid}")
    log(f"Median drugs per fold: {median_n:.0f}")

    log(f"\n--- CO-PRIMARY 1: Fold win fraction ---")
    log(f"Transport-stable wins: {ts_wins}/{n_valid} ({frac:.1%})")
    log(f"Criterion (>= 80%): {'CONFIRMED' if frac >= 0.8 else 'NOT CONFIRMED'}")

    log(f"\n--- CO-PRIMARY 2: Wilcoxon signed-rank on Delta-rho ---")
    log(f"Wilcoxon statistic: {wilcoxon_result.statistic:.1f}")
    log(f"p-value (one-sided, H_a: median > 0): {wilcoxon_result.pvalue:.2e}")
    log(f"Criterion (p < 0.05): {'CONFIRMED' if wilcoxon_result.pvalue < 0.05 else 'NOT CONFIRMED'}")

    log(f"\n--- SECONDARY: Mean Delta-rho ---")
    log(f"Mean raw rho: {np.mean(raw_rhos):.4f} (std {np.std(raw_rhos):.4f})")
    log(f"Mean TS rho:  {np.mean(ts_rhos):.4f} (std {np.std(ts_rhos):.4f})")
    log(f"Mean Delta-rho (TS - raw): {np.mean(improvements):.4f}")
    boot_deltas = []
    boot_rng = np.random.default_rng(seed=2026070805)
    for _ in range(10000):
        boot_idx = boot_rng.choice(len(improvements), size=len(improvements), replace=True)
        boot_deltas.append(np.mean(imp_arr[boot_idx]))
    ci_lo, ci_hi = np.percentile(boot_deltas, [2.5, 97.5])
    log(f"95% bootstrap CI: [{ci_lo:.4f}, {ci_hi:.4f}]")

    log(f"\n--- EXPLORATORY: Frequency-improvement correlation ---")
    n_drugs_arr = np.array(n_drugs_per_fold)
    freq_corr = stats.spearmanr(n_drugs_arr, imp_arr)
    log(f"Spearman(drug_count, Delta-rho): rho={freq_corr.statistic:.4f} (p={freq_corr.pvalue:.2e})")
    log(f"  (Negative = rarer cell lines show larger TS advantage)")

    log(f"\nResults saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--synthetic", action="store_true")
    mode.add_argument("--real", action="store_true")
    parser.add_argument("--data", type=Path,
                        help="Path to drug-perturbation-geometry/data/")
    parser.add_argument("--output", type=Path,
                        default=Path("results/05c_h5_full"))
    args = parser.parse_args()

    if args.synthetic:
        run_synthetic()
    elif args.real:
        if not args.data:
            parser.error("--real requires --data")
        run_real(args.data, args.output)
