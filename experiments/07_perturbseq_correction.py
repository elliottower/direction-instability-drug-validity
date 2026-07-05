"""Experiment 07: Essential-gene correction on Perturb-seq direction instability.

Tests whether projecting out the universal-essential response subspace
changes the Grassmannian geodesic distance ranking for genetic knockdowns,
analogous to the toxicity correction tested on drug signatures in Exp 01.

Uses pre-computed subspaces from genetic-perturbation-holonomy repo
(Replogle et al. 2022 Perturb-seq: K562 + RPE1, 1,676 shared gene targets,
k=5 PCA subspaces, 2,000 HVGs).

Three modes:
  --synthetic   Validate pipeline on data with known essential-response structure
  --real --data DIR  Full pre-registered analysis on extracted subspaces

Usage:
    uv run python experiments/07_perturbseq_correction.py --synthetic
    uv run python experiments/07_perturbseq_correction.py --real \
        --data ../genetic-perturbation-holonomy/data/extracted/
"""
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from geometry.grassmannian import (
    frechet_mean_subspace,
    geodesic_distance,
    project_out_subspace,
    subspace_overlap,
)


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


ESSENTIALITY_THRESHOLD = -0.5


def run_synthetic():
    """Validate: projecting out a shared confound axis changes rankings.

    Plants a confound where essential knockdowns share a response subspace
    in BOTH cell types, making them appear to "transport" (low geodesic
    distance) when the true gene-specific responses are diverse. After
    projecting out the shared confound, essential genes' false transport
    is exposed and rankings change.
    """
    log("=== SYNTHETIC MODE ===")
    rng = np.random.default_rng()
    d = 200
    k = 5
    k_confound = 3
    n_essential = 80
    n_non_essential = 120
    n_total = n_essential + n_non_essential

    confound_basis, _ = np.linalg.qr(rng.standard_normal((d, k_confound)))

    subs_a = []
    subs_b = []
    labels = []

    log("Generating essential knockdowns (shared confound = false transport)...")
    shared_a = confound_basis @ rng.standard_normal((k_confound, k)) * 3.0
    shared_b = shared_a + confound_basis @ rng.standard_normal((k_confound, k)) * 0.3
    for _ in range(n_essential):
        gene_a = rng.standard_normal((d, k)) * rng.uniform(0.3, 1.5)
        gene_b = rng.standard_normal((d, k)) * rng.uniform(0.3, 1.5)
        Qa, _ = np.linalg.qr(shared_a + gene_a)
        Qb, _ = np.linalg.qr(shared_b + gene_b)
        subs_a.append(Qa)
        subs_b.append(Qb)
        labels.append("essential")

    log("Generating non-essential knockdowns (gene-specific, varying distances)...")
    for _ in range(n_non_essential):
        base = rng.standard_normal((d, k))
        Qa, _ = np.linalg.qr(base + rng.standard_normal((d, k)) * 0.1)
        Qb, _ = np.linalg.qr(base + rng.standard_normal((d, k)) * rng.uniform(0.3, 1.5))
        subs_a.append(Qa)
        subs_b.append(Qb)
        labels.append("non_essential")

    log("Computing raw geodesic distances...")
    raw_dists = np.array([geodesic_distance(subs_a[i], subs_b[i])
                          for i in tqdm(range(n_total), desc="Raw")])

    ess_mask = np.array([l == "essential" for l in labels])
    log(f"Raw distances: essential={raw_dists[ess_mask].mean():.4f}, "
        f"non-essential={raw_dists[~ess_mask].mean():.4f}")

    log("Computing Frechet mean of essential-response subspaces...")
    ess_subs_a = [subs_a[i] for i in range(n_total) if ess_mask[i]]
    ess_subs_b = [subs_b[i] for i in range(n_total) if ess_mask[i]]
    mean_a = frechet_mean_subspace(ess_subs_a)
    mean_b = frechet_mean_subspace(ess_subs_b)

    log("Projecting out essential-response subspace...")
    corr_subs_a = [project_out_subspace(subs_a[i], mean_a) for i in range(n_total)]
    corr_subs_b = [project_out_subspace(subs_b[i], mean_b) for i in range(n_total)]

    corr_dists = np.array([geodesic_distance(corr_subs_a[i], corr_subs_b[i])
                           for i in tqdm(range(n_total), desc="Corrected")])

    rho, p = stats.spearmanr(raw_dists, corr_dists)
    log(f"\nRaw vs corrected ranking: Spearman rho={rho:.4f}, p={p:.2e}")
    log(f"Corrected distances: essential={corr_dists[ess_mask].mean():.4f}, "
        f"non-essential={corr_dists[~ess_mask].mean():.4f}")

    ess_raw_mean = raw_dists[ess_mask].mean()
    ess_corr_mean = corr_dists[ess_mask].mean()
    non_raw_mean = raw_dists[~ess_mask].mean()
    non_corr_mean = corr_dists[~ess_mask].mean()

    log(f"\nCorrection effect:")
    log(f"  Essential: {ess_raw_mean:.4f} -> {ess_corr_mean:.4f} (delta={ess_corr_mean - ess_raw_mean:+.4f})")
    log(f"  Non-essential: {non_raw_mean:.4f} -> {non_corr_mean:.4f} (delta={non_corr_mean - non_raw_mean:+.4f})")

    if rho < 0.70:
        log(f"PASS: correction substantially disrupts ranking (rho={rho:.3f} < 0.70)")
    elif rho < 0.83:
        log(f"MARGINAL: correction has moderate effect (0.70 <= rho={rho:.3f} < 0.83)")
    else:
        log(f"MINIMAL: correction barely changes ranking (rho={rho:.3f} >= 0.83)")


def run_real(data_dir: Path, output_dir: Path):
    """Run HP1-HP3 on real Perturb-seq subspaces."""
    log("=== REAL MODE ===")

    subspaces_path = data_dir / "subspaces.npz"
    metadata_path = data_dir / "metadata.json"
    depmap_path = data_dir / "frozen_depmap_labels.json"

    for p in [subspaces_path, metadata_path, depmap_path]:
        if not p.exists():
            log(f"ERROR: {p} not found")
            sys.exit(1)

    log("Loading subspaces...")
    data = np.load(subspaces_path)
    subs_k562 = data["subspaces_k562"]
    subs_rpe1 = data["subspaces_rpe1"]
    raw_dists = data["distances"]
    raw_overlaps = data["overlaps"]
    n_genes, d, k = subs_k562.shape
    log(f"  {n_genes} genes, d={d} features, k={k} components")

    log("Loading metadata...")
    with open(metadata_path) as f:
        meta = json.load(f)
    gene_names = meta["gene_names"]

    log("Loading DepMap labels...")
    with open(depmap_path) as f:
        depmap = json.load(f)

    ess_indices = []
    ess_genes = []
    for i, g in enumerate(gene_names):
        if g in depmap:
            k562_eff = depmap[g]["k562_effect"]
            rpe1_eff = depmap[g]["rpe1_effect"]
            if k562_eff < ESSENTIALITY_THRESHOLD and rpe1_eff < ESSENTIALITY_THRESHOLD:
                ess_indices.append(i)
                ess_genes.append(g)
    ess_indices = np.array(ess_indices)
    n_ess = len(ess_indices)
    log(f"  Universal essential (Chronos < {ESSENTIALITY_THRESHOLD} in both): {n_ess}")

    # --- HP3: Characterize the failure-mode inversion ---
    log("\n=== HP3: Failure-mode inversion ===")
    ess_mask = np.zeros(n_genes, dtype=bool)
    ess_mask[ess_indices] = True

    ess_dists = raw_dists[ess_mask]
    non_dists = raw_dists[~ess_mask]
    d_cohen = (ess_dists.mean() - non_dists.mean()) / np.sqrt(
        (ess_dists.var() + non_dists.var()) / 2)
    u, p_mw = stats.mannwhitneyu(ess_dists, non_dists, alternative="greater")

    log(f"  Essential: mean={ess_dists.mean():.4f}, median={np.median(ess_dists):.4f}")
    log(f"  Non-essential: mean={non_dists.mean():.4f}, median={np.median(non_dists):.4f}")
    log(f"  Cohen's d (essential > non): {d_cohen:.4f}")
    log(f"  Mann-Whitney (essential > non): p={p_mw:.4e}")

    log("  Bootstrap (200 iterations)...")
    boot_ds = []
    rng = np.random.default_rng()
    for _ in range(200):
        boot_ess = rng.choice(ess_dists, size=n_ess, replace=True)
        boot_non = rng.choice(non_dists, size=len(non_dists), replace=True)
        boot_d = (boot_ess.mean() - boot_non.mean()) / np.sqrt(
            (boot_ess.var() + boot_non.var()) / 2)
        boot_ds.append(boot_d)
    boot_ds = np.array(boot_ds)
    log(f"  Bootstrap d: {boot_ds.mean():.4f} [{np.percentile(boot_ds, 2.5):.4f}, "
        f"{np.percentile(boot_ds, 97.5):.4f}]")

    log("  Severity stratification:")
    for thresh in [-0.5, -0.7, -1.0]:
        severe = [i for i, g in enumerate(gene_names)
                  if g in depmap and depmap[g]["k562_effect"] < thresh
                  and depmap[g]["rpe1_effect"] < thresh]
        if severe:
            sev_dists = raw_dists[severe]
            sev_d = (sev_dists.mean() - non_dists.mean()) / np.sqrt(
                (sev_dists.var() + non_dists.var()) / 2)
            log(f"    Chronos < {thresh}: n={len(severe)}, mean_dist={sev_dists.mean():.4f}, "
                f"d={sev_d:.4f}")

    # --- HP1: Essential-response subspace correction ---
    log("\n=== HP1: Essential-response subspace correction ===")

    log("Computing Frechet mean of essential-response subspaces...")
    ess_subs_k562 = [subs_k562[i] for i in ess_indices]
    ess_subs_rpe1 = [subs_rpe1[i] for i in ess_indices]
    mean_k562 = frechet_mean_subspace(ess_subs_k562)
    mean_rpe1 = frechet_mean_subspace(ess_subs_rpe1)
    log(f"  Frechet mean computed (K562: {mean_k562.shape}, RPE1: {mean_rpe1.shape})")

    mean_overlap = subspace_overlap(mean_k562, mean_rpe1)
    log(f"  K562-RPE1 essential mean subspace overlap: {mean_overlap:.4f}")

    log("Projecting out essential-response subspace from all genes...")
    corr_dists = np.zeros(n_genes)
    for i in tqdm(range(n_genes), desc="Correcting"):
        corr_k562 = project_out_subspace(subs_k562[i], mean_k562)
        corr_rpe1 = project_out_subspace(subs_rpe1[i], mean_rpe1)
        corr_dists[i] = geodesic_distance(corr_k562, corr_rpe1)

    rho_hp1, p_hp1 = stats.spearmanr(raw_dists, corr_dists)
    log(f"\n  Raw vs corrected ranking: Spearman rho={rho_hp1:.4f}, p={p_hp1:.2e}")

    if rho_hp1 < 0.70:
        log(f"  HP1 BITES: correction substantially disrupts ranking (rho < 0.70)")
    elif rho_hp1 < 0.83:
        log(f"  HP1 MARGINAL: moderate disruption (0.70 <= rho < 0.83)")
    else:
        log(f"  HP1 DOESN'T BITE: minimal disruption (rho >= 0.83, cf. drug paper 0.83)")

    corr_ess = corr_dists[ess_mask]
    corr_non = corr_dists[~ess_mask]
    log(f"\n  Corrected distances:")
    log(f"    Essential: {corr_ess.mean():.4f} (raw: {ess_dists.mean():.4f}, "
        f"delta={corr_ess.mean() - ess_dists.mean():+.4f})")
    log(f"    Non-essential: {corr_non.mean():.4f} (raw: {non_dists.mean():.4f}, "
        f"delta={corr_non.mean() - non_dists.mean():+.4f})")

    d_corr = (corr_ess.mean() - corr_non.mean()) / np.sqrt(
        (corr_ess.var() + corr_non.var()) / 2)
    log(f"    Cohen's d (ess > non) after correction: {d_corr:.4f} (raw: {d_cohen:.4f})")

    # Top movers: genes whose rank changes most after correction
    raw_ranks = stats.rankdata(raw_dists)
    corr_ranks = stats.rankdata(corr_dists)
    rank_changes = np.abs(corr_ranks - raw_ranks)
    top_movers = np.argsort(rank_changes)[-20:][::-1]
    log(f"\n  Top 20 rank movers after correction:")
    for idx in top_movers:
        g = gene_names[idx]
        is_ess = "ESS" if ess_mask[idx] else "   "
        log(f"    {is_ess} {g:15s} raw_rank={int(raw_ranks[idx]):4d} -> "
            f"corr_rank={int(corr_ranks[idx]):4d} (delta={int(rank_changes[idx]):+4d})")

    # --- HP2: Pathway function prediction ---
    log("\n=== HP2: Not tested (GO annotations not yet curated) ===")

    # --- Save results ---
    output_dir.mkdir(parents=True, exist_ok=True)
    results = {
        "n_genes": n_genes,
        "n_essential": n_ess,
        "essentiality_threshold": ESSENTIALITY_THRESHOLD,
        "hp3_failure_mode_inversion": {
            "essential_mean_dist": float(ess_dists.mean()),
            "non_essential_mean_dist": float(non_dists.mean()),
            "cohens_d": float(d_cohen),
            "mann_whitney_p": float(p_mw),
            "bootstrap_d_mean": float(boot_ds.mean()),
            "bootstrap_d_ci_lower": float(np.percentile(boot_ds, 2.5)),
            "bootstrap_d_ci_upper": float(np.percentile(boot_ds, 97.5)),
        },
        "hp1_correction": {
            "raw_vs_corrected_rho": float(rho_hp1),
            "raw_vs_corrected_p": float(p_hp1),
            "bites": rho_hp1 < 0.70,
            "marginal": 0.70 <= rho_hp1 < 0.83,
            "essential_mean_subspace_overlap": float(mean_overlap),
            "essential_raw_mean": float(ess_dists.mean()),
            "essential_corr_mean": float(corr_ess.mean()),
            "non_essential_raw_mean": float(non_dists.mean()),
            "non_essential_corr_mean": float(corr_non.mean()),
            "cohens_d_after_correction": float(d_corr),
        },
    }

    np.savez(output_dir / "corrected_distances.npz",
             raw_dists=raw_dists, corr_dists=corr_dists,
             ess_mask=ess_mask)

    def make_serializable(obj):
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        return obj

    results = json.loads(json.dumps(results, default=make_serializable))
    with open(output_dir / "perturbseq_correction_results.json", "w") as f:
        json.dump(results, f, indent=2)
    log(f"\nSaved to {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--synthetic", action="store_true")
    mode.add_argument("--real", action="store_true")
    parser.add_argument("--data", type=Path,
                        help="Path to genetic-perturbation-holonomy/data/extracted/")
    parser.add_argument("--output", type=Path,
                        default=Path("results/07_perturbseq_correction"))
    args = parser.parse_args()

    if args.synthetic:
        run_synthetic()
    elif args.real:
        if not args.data:
            parser.error("--real requires --data")
        run_real(args.data, args.output)
