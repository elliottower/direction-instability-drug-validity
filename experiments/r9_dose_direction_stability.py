"""Assumption check R9-DAG: Dose-direction stability.

Tests whether signature DIRECTION is independent of dose (magnitude),
which is the separability assumption underlying the causal DAG claim
that cosine normalization blocks the magnitude path.

For each drug x cell-line pair profiled at >= 2 doses in LINCS L1000,
compute cosine similarity between dose-level mean signatures. If
directions are dose-invariant, cosines should be near 1.0. If they
shift at extreme doses, the DAG assumption fails in that regime.

MUST be pre-registered before running on real data.

Usage:
    uv run python experiments/r9_dose_direction_stability.py --data ../drug-perturbation-geometry/data/
"""
import argparse
import json
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def cosine_sim(a: np.ndarray, b: np.ndarray) -> float:
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-10 or nb < 1e-10:
        return float("nan")
    return float(a @ b / (na * nb))


def run(data_dir: Path, output_dir: Path):
    log("=== R9-DAG: Dose-direction stability ===")

    sigs_path = data_dir / "lincs_subset.npz"
    siginfo_path = data_dir / "GSE92742_Broad_LINCS_sig_info.txt.gz"

    if not sigs_path.exists():
        log(f"ERROR: {sigs_path} not found.")
        sys.exit(1)

    log("Loading signatures...")
    data = np.load(sigs_path, allow_pickle=True)
    signatures = data["signatures"]
    sig_ids = list(data["sig_ids"])
    sig_id_to_idx = {sid: i for i, sid in enumerate(sig_ids)}
    log(f"  {signatures.shape[0]:,} signatures x {signatures.shape[1]} genes")

    log("Loading sig info...")
    siginfo = pd.read_csv(siginfo_path, sep="\t", low_memory=False)

    sig_id_set = set(sig_ids)
    cmpd = siginfo[
        (siginfo["sig_id"].isin(sig_id_set))
        & (siginfo["pert_type"] == "trt_cp")
        & (siginfo["pert_iname"].notna())
    ].copy()
    log(f"  {len(cmpd):,} compound profiles in subset")

    cmpd["dose_num"] = pd.to_numeric(cmpd["pert_dose"], errors="coerce")
    cmpd = cmpd[cmpd["dose_num"].notna() & (cmpd["dose_num"] > 0)].copy()
    cmpd["_idx"] = cmpd["sig_id"].map(sig_id_to_idx)
    log(f"  {len(cmpd):,} with valid numeric dose")

    # Group by drug x cell x dose, compute mean signature per dose level
    log("Computing per-dose mean signatures...")
    pair_doses = defaultdict(dict)  # (drug, cell) -> {dose: mean_sig}
    for (drug, cell, dose), group in tqdm(
        cmpd.groupby(["pert_iname", "cell_id", "dose_num"]),
        desc="Grouping",
    ):
        indices = group["_idx"].values
        mean_sig = signatures[indices].mean(axis=0)
        pair_doses[(drug, cell)][dose] = mean_sig

    # Filter to pairs with >= 2 distinct doses
    multi_dose = {k: v for k, v in pair_doses.items() if len(v) >= 2}
    log(f"  {len(multi_dose):,} drug x cell pairs with >= 2 doses")
    log(f"  {len(set(k[0] for k in multi_dose)):,} unique drugs")

    # Compute pairwise cosine between dose levels within each pair
    log("Computing within-pair dose cosines...")
    all_cosines = []
    all_dose_ratios = []
    pair_results = []

    for (drug, cell), dose_sigs in tqdm(multi_dose.items(), desc="Cosines"):
        doses = sorted(dose_sigs.keys())
        pair_cosines = []
        pair_ratios = []
        for i in range(len(doses)):
            for j in range(i + 1, len(doses)):
                cos = cosine_sim(dose_sigs[doses[i]], dose_sigs[doses[j]])
                if not np.isnan(cos):
                    ratio = doses[j] / doses[i]
                    pair_cosines.append(cos)
                    pair_ratios.append(ratio)
                    all_cosines.append(cos)
                    all_dose_ratios.append(ratio)

        if pair_cosines:
            pair_results.append({
                "drug": drug,
                "cell": cell,
                "n_doses": len(doses),
                "dose_range": [float(doses[0]), float(doses[-1])],
                "mean_cosine": float(np.mean(pair_cosines)),
                "min_cosine": float(np.min(pair_cosines)),
                "n_comparisons": len(pair_cosines),
            })

    cosines = np.array(all_cosines)
    ratios = np.array(all_dose_ratios)
    log(f"  {len(cosines):,} total dose-pair comparisons")

    # Overall statistics
    median_cos = float(np.median(cosines))
    mean_cos = float(np.mean(cosines))
    frac_below_05 = float(np.mean(cosines < 0.5))
    frac_below_03 = float(np.mean(cosines < 0.3))
    frac_above_09 = float(np.mean(cosines > 0.9))
    frac_above_07 = float(np.mean(cosines > 0.7))

    log(f"\n=== OVERALL RESULTS ===")
    log(f"  Median cosine: {median_cos:.4f}")
    log(f"  Mean cosine:   {mean_cos:.4f}")
    log(f"  Fraction > 0.9: {frac_above_09:.3f} ({int(frac_above_09*len(cosines))}/{len(cosines)})")
    log(f"  Fraction > 0.7: {frac_above_07:.3f}")
    log(f"  Fraction < 0.5: {frac_below_05:.3f} ({int(frac_below_05*len(cosines))}/{len(cosines)})")
    log(f"  Fraction < 0.3: {frac_below_03:.3f}")

    # Stratify by dose ratio
    log(f"\n=== STRATIFIED BY DOSE RATIO ===")
    ratio_bins = [(1, 2, "1-2x"), (2, 5, "2-5x"), (5, 20, "5-20x"),
                  (20, 100, "20-100x"), (100, 10000, "100x+")]
    stratified = {}
    for lo, hi, label in ratio_bins:
        mask = (ratios >= lo) & (ratios < hi)
        if mask.sum() > 0:
            bin_cos = cosines[mask]
            s = {
                "label": label,
                "n": int(mask.sum()),
                "median_cosine": float(np.median(bin_cos)),
                "mean_cosine": float(np.mean(bin_cos)),
                "frac_below_05": float(np.mean(bin_cos < 0.5)),
                "frac_above_09": float(np.mean(bin_cos > 0.9)),
            }
            stratified[label] = s
            log(f"  {label:>8s}: n={s['n']:>5d}, median={s['median_cosine']:.3f}, "
                f"<0.5={s['frac_below_05']:.3f}, >0.9={s['frac_above_09']:.3f}")

    # Pre-registered success criterion
    passes_median = median_cos > 0.7
    passes_fraction = frac_below_05 < 0.20
    log(f"\n=== PRE-REGISTERED CRITERIA ===")
    log(f"  Median > 0.7: {'PASS' if passes_median else 'FAIL'} ({median_cos:.4f})")
    log(f"  <20% below 0.5: {'PASS' if passes_fraction else 'FAIL'} ({frac_below_05:.3f})")

    if passes_median and passes_fraction:
        log("  OVERALL: PASS — separability assumption holds for majority of dataset")
    else:
        log("  OVERALL: FAIL — separability assumption violated; report as limitation")

    # Find worst-offender drugs (lowest mean cosine across cells)
    drug_cosines = defaultdict(list)
    for pr in pair_results:
        drug_cosines[pr["drug"]].append(pr["mean_cosine"])
    drug_means = {d: np.mean(cs) for d, cs in drug_cosines.items()}
    worst = sorted(drug_means.items(), key=lambda x: x[1])[:15]
    log(f"\n=== WORST-OFFENDER DRUGS (lowest mean dose-cosine) ===")
    for drug, mc in worst:
        n_cells = len(drug_cosines[drug])
        log(f"  {drug:30s} mean_cos={mc:.3f} ({n_cells} cell lines)")

    # Save results
    output_dir.mkdir(parents=True, exist_ok=True)
    results = {
        "n_pairs": len(multi_dose),
        "n_drugs": len(set(k[0] for k in multi_dose)),
        "n_comparisons": len(cosines),
        "median_cosine": median_cos,
        "mean_cosine": mean_cos,
        "frac_above_09": frac_above_09,
        "frac_above_07": frac_above_07,
        "frac_below_05": frac_below_05,
        "frac_below_03": frac_below_03,
        "passes_median_criterion": passes_median,
        "passes_fraction_criterion": passes_fraction,
        "stratified_by_dose_ratio": stratified,
        "worst_drugs": [{"drug": d, "mean_cosine": float(c)} for d, c in worst],
        "percentiles": {
            f"p{p}": float(np.percentile(cosines, p))
            for p in [5, 10, 25, 50, 75, 90, 95]
        },
    }

    output_file = output_dir / "dose_direction_stability.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nResults saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--data", type=Path, required=True,
        help="Path to drug-perturbation-geometry/data/",
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/r9_dose_stability"),
    )
    args = parser.parse_args()
    run(args.data, args.output)
