"""Modal wrapper: JUMP-CP direction instability analysis.

Downloads ~2.7 GB of JUMP-CP compound profiles on Modal, computes
direction instability across source contexts, applies cell-health
feature correction, and saves results to a Modal volume.

Usage:
    modal run experiments/modal_08_jump_cp.py --detach
"""
import modal

app = modal.App("jump-cp-instability")

vol = modal.Volume.from_name("bracket-norm-results")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "numpy",
        "pandas",
        "pyarrow",
        "scipy",
        "tqdm",
        "matplotlib",
    )
)


@app.function(
    image=image,
    volumes={"/results": vol},
    timeout=86400,
    memory=65536,
    cpu=8,
)
def run_jump_cp_analysis():
    """Download JUMP-CP profiles and compute direction instability."""
    import json
    from datetime import datetime
    from pathlib import Path

    import numpy as np
    import pandas as pd
    from scipy import stats
    from tqdm import tqdm

    def log(msg: str):
        print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

    output_dir = Path("/results/08_jump_cp")
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Cell-health feature substrings (Cell Painting features) ---
    CELL_HEALTH_SUBSTRINGS = [
        "AreaShape_Area",
        "AreaShape_Compactness",
        "AreaShape_Eccentricity",
        "Intensity_MeanIntensity_DNA",
        "Intensity_IntegratedIntensity_DNA",
        "Intensity_MeanIntensity_Mito",
        "RadialDistribution_",
        "Granularity_",
        "Texture_Entropy",
    ]

    # --- Use the INTERPRETABLE (pre-harmony) profiles ---
    # Named Cell Painting features, not PCA components
    log("=== JUMP-CP DIRECTION INSTABILITY ===")
    url = (
        "https://cellpainting-gallery.s3.amazonaws.com/cpg0016-jump-assembled/"
        "source_all/workspace/profiles_assembled/COMPOUND/v1.0/"
        "profiles_var_mad_int.parquet"
    )
    local_parquet = Path("/tmp/jump_profiles.parquet")

    # --- Download once to local disk ---
    if not local_parquet.exists():
        import urllib.request
        log("Step 0: Downloading 11.5 GB parquet to local disk...")
        urllib.request.urlretrieve(url, local_parquet)
        log(f"  Downloaded to {local_parquet}")
    else:
        log(f"Step 0: Using cached {local_parquet}")

    # --- Read full file from disk (much faster than streaming) ---
    log("\nStep 1: Loading parquet from disk...")
    df = pd.read_parquet(local_parquet)
    log(f"  Shape: {df.shape}")

    meta_cols = [c for c in df.columns if c.startswith("Metadata_")]
    feat_cols = [c for c in df.columns if not c.startswith("Metadata_")]
    log(f"  Metadata: {meta_cols}")
    log(f"  Feature columns: {len(feat_cols)}")
    log(f"  Sample features: {feat_cols[:10]}")

    n_sources = df["Metadata_Source"].nunique()
    n_compounds = df["Metadata_JCP2022"].nunique()
    log(f"  Sources: {n_sources}")
    log(f"  Compounds: {n_compounds}")
    log(f"  Source distribution:\n{df['Metadata_Source'].value_counts().to_string()}")

    # Find compounds in >= 5 sources
    ctx_counts = df.groupby("Metadata_JCP2022")["Metadata_Source"].nunique()
    min_contexts = 5
    valid_compounds = ctx_counts[ctx_counts >= min_contexts].index.tolist()
    log(f"  Compounds in >= {min_contexts} sources: {len(valid_compounds)}")

    if len(valid_compounds) < 200:
        min_contexts = 3
        valid_compounds = ctx_counts[ctx_counts >= min_contexts].index.tolist()
        log(f"  Fallback: compounds in >= {min_contexts} sources: {len(valid_compounds)}")

    # Save metadata summary
    with open(output_dir / "metadata_summary.json", "w") as f:
        json.dump({
            "total_rows": len(df),
            "n_sources": int(n_sources),
            "n_compounds_total": int(n_compounds),
            "n_valid_compounds": len(valid_compounds),
            "min_contexts": min_contexts,
            "source_counts": df["Metadata_Source"].value_counts().to_dict(),
        }, f, indent=2)
    vol.commit()
    log("  Metadata summary saved.")

    # Identify cell-health features by name
    health_cols = []
    for col in feat_cols:
        for substr in CELL_HEALTH_SUBSTRINGS:
            if substr in col:
                health_cols.append(col)
                break
    non_health_cols = [c for c in feat_cols if c not in health_cols]
    log(f"  Cell-health features: {len(health_cols)}")
    log(f"  Non-health features: {len(non_health_cols)}")
    log(f"  Health feature examples: {health_cols[:10]}")

    # --- Aggregate per compound × source ---
    log("\nStep 2: Aggregating per compound × source...")
    valid_set = set(valid_compounds)
    df_valid = df[df["Metadata_JCP2022"].isin(valid_set)]
    del df  # free ~11 GB
    log(f"  Valid compound rows: {len(df_valid)}")

    grouped = df_valid.groupby(
        ["Metadata_JCP2022", "Metadata_Source"]
    )[feat_cols].mean()
    del df_valid
    log(f"  Aggregated shape: {grouped.shape}")

    # --- Compute direction instability ---
    log("\nStep 4: Computing direction instability per compound...")

    def direction_instability(signatures: np.ndarray) -> float:
        """D = 1 - mean(pairwise cosines of unit-normalized signatures)."""
        norms = np.linalg.norm(signatures, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-10)
        unit = signatures / norms
        K = unit.shape[0]
        if K < 2:
            return np.nan
        cosines = []
        for i in range(K):
            for j in range(i + 1, K):
                cosines.append(np.dot(unit[i], unit[j]))
        return 1.0 - float(np.mean(cosines))

    raw_instabilities = {}
    corrected_instabilities = {}

    for cpd in tqdm(valid_compounds, desc="Direction instability"):
        if cpd not in grouped.index.get_level_values(0):
            continue
        cpd_data = grouped.loc[cpd]
        if len(cpd_data) < min_contexts:
            continue

        # Raw: all features
        sigs_raw = cpd_data[feat_cols].values.astype(np.float64)
        valid_rows = ~np.any(np.isnan(sigs_raw), axis=1)
        sigs_raw = sigs_raw[valid_rows]
        if sigs_raw.shape[0] >= min_contexts:
            raw_instabilities[cpd] = direction_instability(sigs_raw)

        # Corrected: non-health features only
        sigs_corr = cpd_data[non_health_cols].values.astype(np.float64)
        sigs_corr = sigs_corr[valid_rows]
        if sigs_corr.shape[0] >= min_contexts:
            corrected_instabilities[cpd] = direction_instability(sigs_corr)

    # Filter NaN
    raw_instabilities = {k: v for k, v in raw_instabilities.items()
                         if not np.isnan(v)}
    corrected_instabilities = {k: v for k, v in corrected_instabilities.items()
                               if not np.isnan(v)}

    log(f"\n  Computed raw D for {len(raw_instabilities)} compounds")
    log(f"  Computed corrected D for {len(corrected_instabilities)} compounds")

    # --- Correlation ---
    log("\nStep 5: Computing raw vs corrected correlation...")
    shared_cpds = sorted(set(raw_instabilities) & set(corrected_instabilities))
    raw_vals = np.array([raw_instabilities[c] for c in shared_cpds])
    corr_vals = np.array([corrected_instabilities[c] for c in shared_cpds])

    rho, p = stats.spearmanr(raw_vals, corr_vals)
    log(f"  Spearman rho (raw vs corrected): {rho:.4f}, p={p:.2e}")
    log(f"  N compounds: {len(shared_cpds)}")

    all_raw = np.array(list(raw_instabilities.values()))
    all_corr = np.array(list(corrected_instabilities.values()))
    log(f"  Raw D: mean={all_raw.mean():.4f}, median={np.median(all_raw):.4f}, "
        f"std={all_raw.std():.4f}")
    log(f"  Corrected D: mean={all_corr.mean():.4f}, median={np.median(all_corr):.4f}, "
        f"std={all_corr.std():.4f}")

    if rho < 0.70:
        verdict = "BITES"
    elif rho < 0.83:
        verdict = "MARGINAL"
    else:
        verdict = "DOESN'T BITE"
    log(f"\n  Verdict: {verdict}")
    log(f"  (cf. LINCS drugs rho=0.83, Perturb-seq rho=0.91)")

    # Top movers
    raw_ranks = stats.rankdata(raw_vals)
    corr_ranks = stats.rankdata(corr_vals)
    rank_changes = np.abs(corr_ranks - raw_ranks)
    top_movers_idx = np.argsort(rank_changes)[-20:][::-1]
    log(f"\n  Top 20 rank movers after cell-health correction:")
    for idx in top_movers_idx:
        cpd = shared_cpds[idx]
        log(f"    {cpd:20s} raw_rank={int(raw_ranks[idx]):5d} -> "
            f"corr_rank={int(corr_ranks[idx]):5d} "
            f"(delta={int(rank_changes[idx]):+5d})")

    # --- Save results ---
    log("\nStep 6: Saving results...")
    results = {
        "n_compounds": len(shared_cpds),
        "min_contexts": min_contexts,
        "n_features_total": len(feat_cols),
        "n_features_health": len(health_cols),
        "n_features_non_health": len(non_health_cols),
        "health_features": health_cols[:50],
        "raw_instability": {
            "mean": float(all_raw.mean()),
            "median": float(np.median(all_raw)),
            "std": float(all_raw.std()),
            "min": float(all_raw.min()),
            "max": float(all_raw.max()),
        },
        "corrected_instability": {
            "mean": float(all_corr.mean()),
            "median": float(np.median(all_corr)),
            "std": float(all_corr.std()),
        },
        "correction": {
            "raw_vs_corrected_rho": float(rho),
            "raw_vs_corrected_p": float(p),
            "verdict": verdict,
        },
        "top_movers": [
            {
                "compound": shared_cpds[idx],
                "raw_rank": int(raw_ranks[idx]),
                "corrected_rank": int(corr_ranks[idx]),
                "rank_change": int(rank_changes[idx]),
            }
            for idx in top_movers_idx
        ],
    }

    with open(output_dir / "jump_cp_results.json", "w") as f:
        json.dump(results, f, indent=2)

    np.savez(
        output_dir / "jump_cp_instabilities.npz",
        compounds=np.array(shared_cpds, dtype=object),
        raw=raw_vals,
        corrected=corr_vals,
    )

    vol.commit()
    log(f"\nResults saved to /results/08_jump_cp/")
    log("=== DONE ===")

    return results


@app.local_entrypoint()
def main():
    call = run_jump_cp_analysis.spawn()
    print(f"Spawned: {call.object_id}")
    print("Function will run to completion independently.")
    print("Check results in Modal volume 'bracket-norm-results' at /08_jump_cp/")
