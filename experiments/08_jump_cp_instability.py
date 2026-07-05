"""Experiment 08: Direction instability in JUMP-CP morphological profiles.

Tests whether the validity ladder generalizes to a third perturbation domain:
Cell Painting morphological features across cell-type/batch contexts.

JUMP-CP (Joint Undertaking in Morphological Profiling) provides ~116,000
compound profiles across U2OS cells in multiple batches/plates. We compute
direction instability of morphological signatures across plate contexts and
test whether toxicity correction (cell-health features) changes the ranking.

Data source: JUMP-CP consortium pre-computed profiles
  https://github.com/jump-cellpainting/datasets
  S3: cellpainting-gallery/cpg0016-jump

This script:
  1. Downloads consensus compound profiles (parquet) from JUMP
  2. Computes direction instability across plate replicates
  3. Tests toxicity correction (cell-health morphological features)
  4. Reports the same rho(raw, corrected) metric as Exp 01 and 07

Usage:
    uv run python experiments/08_jump_cp_instability.py --download
    uv run python experiments/08_jump_cp_instability.py --analyze
"""
import argparse
import json
import sys
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

DATA_DIR = Path("data/jump_cp")

JUMP_PROFILES_URL = (
    "https://github.com/jump-cellpainting/datasets/raw/main/"
    "profile_index.csv"
)

# JUMP-CP well-level profiles from cellpainting-gallery
# Using the cpg0016 compound set, source_4 (U2OS)
JUMP_PARQUET_BASE = (
    "https://cellpainting-gallery.s3.amazonaws.com/cpg0016-jump/"
    "{source}/workspace/profiles/{batch}/{plate}/"
    "{plate}.parquet"
)

# Cell-health / toxicity-related CellPainting features
# These capture generic stress: cell area shrinkage, granularity changes,
# intensity changes in DNA/mito channels that indicate apoptosis
CELL_HEALTH_PREFIXES = [
    "Cells_AreaShape_Area",
    "Cells_AreaShape_Compactness",
    "Cells_AreaShape_Eccentricity",
    "Nuclei_AreaShape_Area",
    "Nuclei_Intensity_MeanIntensity_DNA",
    "Nuclei_Intensity_IntegratedIntensity_DNA",
    "Cytoplasm_AreaShape_Area",
    "Cells_Intensity_MeanIntensity_Mito",
    "Nuclei_RadialDistribution_",
    "Cells_Granularity_",
    "Nuclei_Texture_Entropy",
]


def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def download_profile_index():
    """Download the JUMP-CP profile index listing available plates."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    index_path = DATA_DIR / "profile_index.csv"
    if index_path.exists():
        log(f"  Profile index already cached: {index_path}")
        return index_path
    log("  Downloading JUMP-CP profile index...")
    urllib.request.urlretrieve(JUMP_PROFILES_URL, index_path)
    log(f"  Saved to {index_path}")
    return index_path


def download_jump_profiles():
    """Download a subset of JUMP-CP compound profiles for analysis.

    Strategy: Use the pre-aggregated 'consensus' profiles from the
    JUMP datasets repo, which provide compound-level aggregated features.
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # The JUMP consortium provides pre-computed consensus profiles
    # We use the compound profiles aggregated per source (batch context)
    consensus_url = (
        "https://github.com/jump-cellpainting/datasets/raw/main/"
        "profiles/jump-consensus-compound-profiles.parquet"
    )
    consensus_path = DATA_DIR / "jump-consensus-compound-profiles.parquet"

    if consensus_path.exists():
        log(f"  Consensus profiles already cached: {consensus_path}")
        return consensus_path

    log("  Downloading JUMP-CP consensus compound profiles...")
    log(f"  URL: {consensus_url}")
    try:
        urllib.request.urlretrieve(consensus_url, consensus_path)
        log(f"  Saved to {consensus_path}")
    except Exception as e:
        log(f"  Direct download failed: {e}")
        log("  Trying alternative: plate-level profiles...")
        return download_plate_profiles()
    return consensus_path


def download_plate_profiles():
    """Fallback: download individual plate-level profiles."""
    # Use the JUMP datasets metadata to find available profiles
    meta_url = (
        "https://raw.githubusercontent.com/jump-cellpainting/datasets/"
        "main/metadata/compound.csv.gz"
    )
    meta_path = DATA_DIR / "compound_metadata.csv.gz"
    if not meta_path.exists():
        log("  Downloading compound metadata...")
        urllib.request.urlretrieve(meta_url, meta_path)

    # Also get the plate-compound mapping
    plate_url = (
        "https://raw.githubusercontent.com/jump-cellpainting/datasets/"
        "main/metadata/plate.csv.gz"
    )
    plate_path = DATA_DIR / "plate_metadata.csv.gz"
    if not plate_path.exists():
        log("  Downloading plate metadata...")
        urllib.request.urlretrieve(plate_url, plate_path)

    return meta_path


def compute_direction_instability(signatures: np.ndarray) -> float:
    """Compute direction instability from K signature vectors.

    Args:
        signatures: (K, d) array of signature vectors across contexts.

    Returns:
        D = 1 - mean(pairwise cosines of unit-normalized signatures).
    """
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
    return 1.0 - np.mean(cosines)


def identify_cell_health_features(columns: list[str]) -> list[str]:
    """Identify cell-health/toxicity-related features."""
    health_cols = []
    for col in columns:
        for prefix in CELL_HEALTH_PREFIXES:
            if prefix in col:
                health_cols.append(col)
                break
    return health_cols


def run_download():
    """Download JUMP-CP data."""
    log("=== DOWNLOAD JUMP-CP DATA ===")
    result = download_jump_profiles()
    log(f"\nDownload complete: {result}")

    if result.suffix == ".parquet":
        df = pd.read_parquet(result)
        log(f"  Shape: {df.shape}")
        log(f"  Columns (first 20): {list(df.columns[:20])}")
        if "Metadata_InChIKey" in df.columns:
            n_compounds = df["Metadata_InChIKey"].nunique()
            log(f"  Unique compounds (InChIKey): {n_compounds}")
        if "Metadata_Source" in df.columns:
            log(f"  Sources: {df['Metadata_Source'].value_counts().to_dict()}")
        elif "Metadata_Plate" in df.columns:
            n_plates = df["Metadata_Plate"].nunique()
            log(f"  Unique plates: {n_plates}")


def run_analyze(output_dir: Path):
    """Compute direction instability on JUMP-CP profiles."""
    log("=== ANALYZE JUMP-CP DIRECTION INSTABILITY ===")

    consensus_path = DATA_DIR / "jump-consensus-compound-profiles.parquet"
    if not consensus_path.exists():
        log("ERROR: Run --download first")
        sys.exit(1)

    log("Loading consensus profiles...")
    df = pd.read_parquet(consensus_path)
    log(f"  Shape: {df.shape}")
    log(f"  Columns (first 30): {list(df.columns[:30])}")

    # Identify metadata vs feature columns
    meta_cols = [c for c in df.columns if c.startswith("Metadata_")]
    feat_cols = [c for c in df.columns if not c.startswith("Metadata_")]
    log(f"  Metadata columns: {len(meta_cols)}")
    log(f"  Feature columns: {len(feat_cols)}")
    log(f"  Metadata: {meta_cols}")

    # Determine the context variable (plate, batch, or source)
    context_col = None
    for candidate in ["Metadata_Source", "Metadata_Batch", "Metadata_Plate"]:
        if candidate in df.columns:
            context_col = candidate
            break
    if context_col is None:
        log("ERROR: No context column found")
        sys.exit(1)
    log(f"  Context column: {context_col}")
    log(f"  Context values: {df[context_col].nunique()}")

    # Compound identifier
    compound_col = None
    for candidate in ["Metadata_InChIKey", "Metadata_JCP2022",
                      "Metadata_pert_id", "Metadata_broad_sample"]:
        if candidate in df.columns:
            compound_col = candidate
            break
    if compound_col is None:
        log("ERROR: No compound identifier found")
        sys.exit(1)
    log(f"  Compound column: {compound_col}")

    # Filter to compounds tested in >= 3 contexts
    context_counts = df.groupby(compound_col)[context_col].nunique()
    min_contexts = 3
    valid_compounds = context_counts[context_counts >= min_contexts].index
    log(f"  Compounds in >= {min_contexts} contexts: {len(valid_compounds)}")

    if len(valid_compounds) < 100:
        log(f"  Trying min_contexts=2...")
        min_contexts = 2
        valid_compounds = context_counts[context_counts >= min_contexts].index
        log(f"  Compounds in >= {min_contexts} contexts: {len(valid_compounds)}")

    df_valid = df[df[compound_col].isin(valid_compounds)].copy()
    log(f"  Working dataset: {df_valid.shape}")

    # Drop features with NaN
    feat_matrix = df_valid[feat_cols].values
    nan_cols = np.any(np.isnan(feat_matrix), axis=0)
    clean_feat_cols = [feat_cols[i] for i in range(len(feat_cols)) if not nan_cols[i]]
    log(f"  Features after dropping NaN columns: {len(clean_feat_cols)}")

    # Compute direction instability per compound
    log("\nComputing direction instability...")
    compounds = list(valid_compounds)
    instabilities = {}

    for cpd in tqdm(compounds, desc="Computing D"):
        cpd_rows = df_valid[df_valid[compound_col] == cpd]
        sigs = cpd_rows[clean_feat_cols].values.astype(np.float64)
        # Remove any rows with NaN
        valid_rows = ~np.any(np.isnan(sigs), axis=1)
        sigs = sigs[valid_rows]
        if sigs.shape[0] >= min_contexts:
            instabilities[cpd] = compute_direction_instability(sigs)

    instabilities = {k: v for k, v in instabilities.items() if not np.isnan(v)}
    log(f"  Computed D for {len(instabilities)} compounds")

    if len(instabilities) == 0:
        log("ERROR: No valid instability scores computed")
        sys.exit(1)

    all_D = np.array(list(instabilities.values()))
    log(f"  Direction instability: mean={all_D.mean():.4f}, "
        f"median={np.median(all_D):.4f}, std={all_D.std():.4f}")
    log(f"  Range: [{all_D.min():.4f}, {all_D.max():.4f}]")

    # --- Toxicity correction ---
    log("\nApplying cell-health feature correction...")
    health_features = identify_cell_health_features(clean_feat_cols)
    log(f"  Cell-health features identified: {len(health_features)}")

    non_health_features = [f for f in clean_feat_cols if f not in health_features]
    log(f"  Non-health features: {len(non_health_features)}")

    corrected_instabilities = {}
    for cpd in tqdm(compounds, desc="Corrected D"):
        if cpd not in instabilities:
            continue
        cpd_rows = df_valid[df_valid[compound_col] == cpd]
        sigs = cpd_rows[non_health_features].values.astype(np.float64)
        valid_rows = ~np.any(np.isnan(sigs), axis=1)
        sigs = sigs[valid_rows]
        if sigs.shape[0] >= min_contexts:
            corrected_instabilities[cpd] = compute_direction_instability(sigs)

    corrected_instabilities = {k: v for k, v in corrected_instabilities.items()
                               if not np.isnan(v)}

    # Match compounds present in both
    shared_cpds = sorted(set(instabilities) & set(corrected_instabilities))
    log(f"  Compounds with both raw and corrected: {len(shared_cpds)}")

    raw_vals = np.array([instabilities[c] for c in shared_cpds])
    corr_vals = np.array([corrected_instabilities[c] for c in shared_cpds])

    rho, p = stats.spearmanr(raw_vals, corr_vals)
    log(f"\n  Raw vs corrected ranking: Spearman rho={rho:.4f}, p={p:.2e}")

    corr_D = np.array(list(corrected_instabilities.values()))
    log(f"  Corrected D: mean={corr_D.mean():.4f}, "
        f"median={np.median(corr_D):.4f}, std={corr_D.std():.4f}")

    if rho < 0.70:
        verdict = "BITES"
    elif rho < 0.83:
        verdict = "MARGINAL"
    else:
        verdict = "DOESN'T BITE"
    log(f"\n  Verdict: {verdict} (rho={rho:.4f}, cf. drug=0.83, Perturb-seq=0.91)")

    # Top movers
    raw_ranks = stats.rankdata(raw_vals)
    corr_ranks = stats.rankdata(corr_vals)
    rank_changes = np.abs(corr_ranks - raw_ranks)
    top_movers_idx = np.argsort(rank_changes)[-20:][::-1]
    log(f"\n  Top 20 rank movers after cell-health correction:")
    for idx in top_movers_idx:
        cpd = shared_cpds[idx]
        log(f"    {cpd[:20]:20s} raw_rank={int(raw_ranks[idx]):5d} -> "
            f"corr_rank={int(corr_ranks[idx]):5d} (delta={int(rank_changes[idx]):+5d})")

    # Save results
    output_dir.mkdir(parents=True, exist_ok=True)
    results = {
        "n_compounds": len(shared_cpds),
        "n_features_total": len(clean_feat_cols),
        "n_features_health": len(health_features),
        "n_features_non_health": len(non_health_features),
        "min_contexts": min_contexts,
        "context_column": context_col,
        "compound_column": compound_col,
        "raw_instability": {
            "mean": float(all_D.mean()),
            "median": float(np.median(all_D)),
            "std": float(all_D.std()),
        },
        "corrected_instability": {
            "mean": float(corr_D.mean()),
            "median": float(np.median(corr_D)),
            "std": float(corr_D.std()),
        },
        "correction": {
            "raw_vs_corrected_rho": float(rho),
            "raw_vs_corrected_p": float(p),
            "verdict": verdict,
        },
    }

    with open(output_dir / "jump_cp_results.json", "w") as f:
        json.dump(results, f, indent=2)
    np.savez(
        output_dir / "jump_cp_instabilities.npz",
        compounds=np.array(shared_cpds, dtype=object),
        raw=raw_vals,
        corrected=corr_vals,
    )
    log(f"\nSaved to {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--download", action="store_true",
                      help="Download JUMP-CP profiles")
    mode.add_argument("--analyze", action="store_true",
                      help="Compute direction instability and correction")
    parser.add_argument("--output", type=Path,
                        default=Path("results/08_jump_cp"))
    args = parser.parse_args()

    if args.download:
        run_download()
    elif args.analyze:
        run_analyze(args.output)
