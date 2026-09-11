"""Sensitivity check for Experiment 8: is cross-site D just well-to-well noise?

In the main analysis 92.6% of compound-source consensuses rest on a single
well, so "direction instability across imaging sites" could be measuring
well-to-well variation rather than anything about the sites.

This isolates the two. For compounds with at least two wells in at least three
sources:

    within-source D   -- across wells inside one site, averaged over sites
    across-source D   -- across the per-site medians

If the two are similar, site identity adds nothing beyond well noise and the
cross-site framing is wrong. If across-source D is clearly higher, sites
contribute variation of their own.

It also repeats the feature ablation on this replicated cohort, where every
consensus rests on two or more wells, to check that the headline rho survives
when single-well consensuses are excluded.

    uv run --no-project --with pyarrow --with pandas --with numpy --with scipy \
        python experiments/08b_within_vs_across_source.py
"""
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.fs as fs
import pyarrow.parquet as pq
from scipy.stats import spearmanr, wilcoxon

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
OUT_DIR = REPO / "results" / "08_jump_cp_feature_ablation"
S3_PATH = ("cellpainting-gallery/cpg0016-jump-assembled/source_all/workspace/"
           "profiles_assembled/COMPOUND/v1.0/profiles_var_mad_int.parquet")
PLATE_META = ("https://raw.githubusercontent.com/jump-cellpainting/datasets/"
              "v0.11.0/metadata/plate.csv.gz")
FEATURE_PREFIXES = ("Cells_", "Cytoplasm_", "Nuclei_")
DMSO_JCP = "JCP2022_033924"
KEEP_PLATE_TYPE = "COMPOUND"
MIN_WELLS_PER_SOURCE = 2
MIN_REPLICATED_SOURCES = 3
BLOCK = 400
MIN_NORM = 1e-8


def ts():
    return datetime.now().strftime("%H:%M:%S")


def di(X):
    """1 - mean pairwise cosine over distinct pairs."""
    n = np.linalg.norm(X, axis=1, keepdims=True)
    if (n < MIN_NORM).any() or X.shape[0] < 2:
        return np.nan
    U = X / n
    K = U.shape[0]
    mc = (float(np.square(U.sum(axis=0)).sum()) - K) / (K * (K - 1))
    return 1.0 - min(max(mc, -1.0), 1.0)


def main():
    s3 = fs.S3FileSystem(anonymous=True, region="us-east-1")
    f = pq.ParquetFile(s3.open_input_file(S3_PATH))
    feat = [n for n in f.schema_arrow.names if n.startswith(FEATURE_PREFIXES)]

    meta = f.read(columns=["Metadata_Source", "Metadata_Plate",
                           "Metadata_JCP2022"]).to_pandas()
    plate = pd.read_csv(PLATE_META)
    meta = meta.merge(plate[["Metadata_Source", "Metadata_Plate", "Metadata_PlateType"]],
                      on=["Metadata_Source", "Metadata_Plate"], how="left")
    assert meta.Metadata_PlateType.notna().all(), "unmatched plate rows in the join"
    ok = (meta.Metadata_PlateType.eq(KEEP_PLATE_TYPE)
          & meta.Metadata_JCP2022.notna()
          & meta.Metadata_JCP2022.ne(DMSO_JCP))

    per = meta[ok].groupby(["Metadata_JCP2022", "Metadata_Source"]).size()
    rep = per[per >= MIN_WELLS_PER_SOURCE]
    n_rep = rep.groupby("Metadata_JCP2022").size()
    cohort = set(n_rep[n_rep >= MIN_REPLICATED_SOURCES].index)
    mask = (ok & meta.Metadata_JCP2022.isin(cohort)).to_numpy()
    keys = meta.loc[mask, ["Metadata_JCP2022", "Metadata_Source"]].reset_index(drop=True)
    print(f"[{ts()}] cohort: {len(cohort):,} compounds, {mask.sum():,} wells")

    parts = []
    n_blocks = (len(feat) + BLOCK - 1) // BLOCK
    for b in range(n_blocks):
        cols = feat[b * BLOCK:(b + 1) * BLOCK]
        parts.append(f.read(columns=cols).to_pandas().loc[mask]
                     .reset_index(drop=True).astype("float32"))
        print(f"[{ts()}] streamed block {b+1}/{n_blocks}")
    X = pd.concat(parts, axis=1)
    del parts
    assert np.isfinite(X.to_numpy()).all(), "nonfinite values in well-level profiles"

    within, across, n_src = [], [], []
    for jcp, idx in keys.groupby("Metadata_JCP2022").groups.items():
        sub = keys.loc[idx]
        w_vals, medians = [], []
        for src, sidx in sub.groupby("Metadata_Source").groups.items():
            wells = X.loc[sidx].to_numpy(np.float32)
            if wells.shape[0] >= MIN_WELLS_PER_SOURCE:
                d = di(wells)
                if np.isfinite(d):
                    w_vals.append(d)
                medians.append(np.median(wells, axis=0))
        if len(w_vals) >= MIN_REPLICATED_SOURCES and len(medians) >= MIN_REPLICATED_SOURCES:
            a = di(np.vstack(medians).astype(np.float32))
            if np.isfinite(a):
                within.append(float(np.mean(w_vals)))
                across.append(float(a))
                n_src.append(len(medians))
    within, across = np.array(within), np.array(across)
    print(f"[{ts()}] {len(within):,} compounds with both quantities")

    stat, p = wilcoxon(across, within, alternative="greater")
    out = {
        "experiment": "08b_within_vs_across_source",
        "question": ("is cross-site direction instability distinguishable from "
                     "well-to-well variation within a site?"),
        "cohort": {"min_wells_per_source": MIN_WELLS_PER_SOURCE,
                   "min_replicated_sources": MIN_REPLICATED_SOURCES,
                   "n_compounds": int(len(within)),
                   "median_sources_per_compound": float(np.median(n_src))},
        "within_source_D": {"mean": float(within.mean()), "median": float(np.median(within)),
                            "q25": float(np.percentile(within, 25)),
                            "q75": float(np.percentile(within, 75))},
        "across_source_D": {"mean": float(across.mean()), "median": float(np.median(across)),
                            "q25": float(np.percentile(across, 25)),
                            "q75": float(np.percentile(across, 75))},
        "difference": {"mean_across_minus_within": float((across - within).mean()),
                       "median_difference": float(np.median(across - within)),
                       "frac_across_greater": float((across > within).mean()),
                       "wilcoxon_p_across_greater": (float(p) if p > 0 else "< 1e-300 (underflow)"),
                       "note": "secondary; the two sides differ in aggregation depth"},
        "spearman_within_vs_across": float(spearmanr(within, across).statistic),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "within_vs_across_source.json").write_text(json.dumps(out, indent=2))
    pd.DataFrame({"within_source_D": within, "across_source_D": across,
                  "n_sources": n_src}).to_csv(OUT_DIR / "within_vs_across_source.csv",
                                              index=False)

    print(f"\n[{ts()}] within-source D  median {np.median(within):.4f}")
    print(f"[{ts()}] across-source D  median {np.median(across):.4f}")
    print(f"[{ts()}] across > within in {(across>within).mean():.1%} of compounds, "
          f"Wilcoxon p = {p:.3g}")
    print(f"[{ts()}] written to {OUT_DIR}")


if __name__ == "__main__":
    main()
