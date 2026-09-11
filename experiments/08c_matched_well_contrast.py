"""Experiment 8, sensitivity analyses 1-3: is the cross-site signal real?

08b compared within-source D (raw wells) against across-source D (source
medians) and found +0.23. Those are different estimands at unequal aggregation
depth, so the comparison supports a direction but not a magnitude. These three
tests put both sides on equal footing.

1. Matched raw-well contrast. For each compound, split its well pairs into
   same-source and different-source and compare mean cosine distance. Both
   sides are raw wells, so neither benefits from consensus aggregation.

       Delta_c = D_different_c - D_same_c

2. Source-label permutation. Within each compound, shuffle source labels across
   its wells preserving the number of plates per source label, and recompute Delta_c. This
   null absorbs source counts, well counts and the observed vectors, so it is
   the right reference for test 1. The cosine matrix is computed once per
   compound and reused, making permutations nearly free.

3. Replicated-cohort ablation. Repeat the frozen feature ablation using
   consensuses that all rest on two or more wells, with its size-matched random
   null, to check that rho = 0.9977 was not an artifact of single-well
   consensuses.

Cohort: compounds with at least 2 wells in at least 3 sources, COMPOUND plates
only, DMSO excluded.

    uv run --no-project --with pyarrow --with pandas --with numpy --with scipy \
        python experiments/08c_matched_well_contrast.py
"""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.fs as fs
import pyarrow.parquet as pq
from scipy.stats import spearmanr

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
OUT_DIR = REPO / "results" / "08_jump_cp_feature_ablation"
CACHE = REPO / "data" / "jump_cp" / "cpg0016"   # overridden by the Modal wrapper
S3_PATH = ("cellpainting-gallery/cpg0016-jump-assembled/source_all/workspace/"
           "profiles_assembled/COMPOUND/v1.0/profiles_var_mad_int.parquet")
PLATE_META = ("https://raw.githubusercontent.com/jump-cellpainting/datasets/"
              "v0.11.0/metadata/plate.csv.gz")
FEATURE_PREFIXES = ("Cells_", "Cytoplasm_", "Nuclei_")
DMSO_JCP = "JCP2022_033924"
KEEP_PLATE_TYPE = "COMPOUND"
EXPECTED_ROWS = 803_853        # pinned release, cross-checked in the main script
EXPECTED_FEATURES = 3_180   # row-level, not whole-compound exclusion
MIN_WELLS_PER_SOURCE = 2
MIN_REPLICATED_SOURCES = 3
BLOCK = 150                # smaller blocks: peak is one block, not the file
MIN_NORM = 1e-8
N_PERM = 1000
N_RANDOM_ABLATIONS = 200
SEED = 20260908


CKPT = OUT_DIR / "08c_checkpoints"


def ts():
    return datetime.now().strftime("%H:%M:%S")


def log(msg):
    print(f"[{ts()}] {msg}", flush=True)


def fingerprint(**extra):
    """Bind a cached artifact to the configuration that produced it."""
    import hashlib
    fp = {"plate_filter": f"row-level Metadata_PlateType == {KEEP_PLATE_TYPE}",
          "min_wells_per_source": MIN_WELLS_PER_SOURCE,
          "min_replicated_sources": MIN_REPLICATED_SOURCES,
          "seed": SEED, "n_perm": N_PERM, **extra}
    return hashlib.sha256(json.dumps(fp, sort_keys=True, default=str).encode()).hexdigest()[:16]


def check_fingerprint(path, fp):
    """Refuse to resume from a checkpoint written under a different configuration.

    Inspects the stamp independently of `path`, so it still works when `path` is
    a sentinel that is never created (e.g. a directory-stage marker).
    """
    stamp = path.with_suffix(path.suffix + ".fingerprint")
    if stamp.exists() and stamp.read_text().strip() != fp:
        raise SystemExit(f"stale checkpoint for {path.name}: "
                         f"{stamp.read_text().strip()} != {fp}. "
                         "Delete the checkpoint directory and rerun.")
    if path.exists() and not stamp.exists():
        raise SystemExit(f"{path.name} exists without a fingerprint; delete and rerun.")
    stamp.parent.mkdir(parents=True, exist_ok=True)
    return path.exists()


def stamp_fingerprint(path, fp):
    """Write the stamp once the artifact exists, never before.

    Stamping ahead of the computation leaves a stamp with no data behind it when
    a run dies mid-stage, and the next run reads that as a valid configuration.
    """
    assert path.exists(), f"refusing to stamp {path.name}: artifact absent"
    path.with_suffix(path.suffix + ".fingerprint").write_text(fp)


def unit(X):
    n = np.linalg.norm(X, axis=1, keepdims=True)
    return X / n, n.ravel()


def di_from_unit(U):
    K = U.shape[0]
    mc = (float(np.square(U.sum(axis=0)).sum()) - K) / (K * (K - 1))
    return 1.0 - min(max(mc, -1.0), 1.0)


def load_cohort():
    s3 = fs.S3FileSystem(anonymous=True, region="us-east-1")
    f = pq.ParquetFile(s3.open_input_file(S3_PATH))
    feat = [n for n in f.schema_arrow.names if n.startswith(FEATURE_PREFIXES)]

    assert f.metadata.num_rows == EXPECTED_ROWS, "row count differs from the pinned release"
    assert len(feat) == EXPECTED_FEATURES, f"{len(feat)} features, expected {EXPECTED_FEATURES}"
    wm = CACHE / "well_metadata.parquet"
    assert wm.exists(), "run the main script's --eligible stage first"
    meta = pd.read_parquet(wm)          # already carries Metadata_PlateType; do not re-join
    required = {"Metadata_Source", "Metadata_Plate", "Metadata_JCP2022", "Metadata_PlateType"}
    missing = required - set(meta.columns)
    assert not missing, f"cached metadata missing columns: {missing}"
    ok = (meta.Metadata_PlateType.eq(KEEP_PLATE_TYPE)
          & meta.Metadata_JCP2022.notna()
          & meta.Metadata_JCP2022.ne(DMSO_JCP))

    per = meta[ok].groupby(["Metadata_JCP2022", "Metadata_Source"]).size()
    rep = per[per >= MIN_WELLS_PER_SOURCE]
    n_rep = rep.groupby("Metadata_JCP2022").size()
    cohort = set(n_rep[n_rep >= MIN_REPLICATED_SOURCES].index)
    keep_pairs = set(rep.index)                       # only replicated sources

    sel = ok & meta.Metadata_JCP2022.isin(cohort)
    sel &= pd.MultiIndex.from_frame(
        meta[["Metadata_JCP2022", "Metadata_Source"]]).isin(keep_pairs)
    mask = sel.to_numpy()
    keys = meta.loc[mask, ["Metadata_JCP2022", "Metadata_Source",
                           "Metadata_Plate"]].reset_index(drop=True)
    log(f"cohort {len(cohort):,} compounds, {mask.sum():,} wells")

    # One preallocated array, filled block by block. pd.concat would hold the
    # blocks, the concatenated frame and a to_numpy copy simultaneously -- three
    # copies of a 1.2 GB matrix, which OOMs a 16 GB machine.
    n_rows = int(mask.sum())
    XV = np.empty((n_rows, len(feat)), dtype=np.float32)
    take_idx = np.flatnonzero(mask)
    n_blocks = (len(feat) + BLOCK - 1) // BLOCK
    for b in range(n_blocks):
        lo, hi = b * BLOCK, min((b + 1) * BLOCK, len(feat))
        # subset rows in Arrow before pandas sees them: .to_pandas() on the full
        # 803,853 rows costs ~2.6 GB a block and OOMs a 16 GB machine
        tbl = f.read(columns=feat[lo:hi]).take(take_idx)
        arr = np.column_stack([c.to_numpy(zero_copy_only=False).astype(np.float32,
                                                                      copy=False)
                               for c in tbl.columns])
        assert np.isfinite(arr).all(), f"nonfinite values in block {b}"
        XV[:, lo:hi] = arr
        del tbl, arr
        log(f"streamed block {b+1}/{n_blocks}")
    assert len(keys) == n_rows
    return keys, XV, feat



def main():
    keys, XV, feat = load_cohort()
    pos = keys.groupby("Metadata_JCP2022", sort=True).indices   # positional
    SRC_ALL = keys.Metadata_Source.to_numpy()      # hoisted: materialising these
    PLATE_ALL = keys.Metadata_Plate.to_numpy()     # per compound copied 93k rows each time
    perm_rng = np.random.default_rng(SEED + 1)
    abl_rng = np.random.default_rng(SEED + 2)

    # ---- tests 1 and 2: matched raw-well contrast and its permutation null
    CKPT.mkdir(parents=True, exist_ok=True)
    ck = CKPT / "contrast.npz"
    if check_fingerprint(ck, fingerprint(stage="contrast")):
        z = np.load(ck, allow_pickle=True)
        same, diff, delta, n_src, perm_mean = (z["same"], z["diff"], z["delta"],
                                               z["n_src"], z["perm_mean"])
        perm_delta, singleton = z["perm_delta"], z["singleton"]
        log(f"resumed contrast from checkpoint: {len(delta):,} compounds")
    else:
        same, diff, delta, n_src, singleton = [], [], [], [], []
        perm_rows = []
        done = 0
        for gi, (jcp, idx) in enumerate(pos.items()):
            U, norms = unit(XV[idx].astype(np.float64))
            if (norms < MIN_NORM).any():
                continue
            C = 1.0 - np.clip(U @ U.T, -1.0, 1.0)
            src = SRC_ALL[idx]
            iu = np.triu_indices(len(idx), k=1)
            d = C[iu]
            same_mask = src[iu[0]] == src[iu[1]]
            if not same_mask.any() or same_mask.all():
                continue
            plates = PLATE_ALL[idx]
            pl_codes, pl_uniq = pd.factorize(plates)
            s_, d_ = float(d[same_mask].mean()), float(d[~same_mask].mean())
            same.append(s_); diff.append(d_); delta.append(d_ - s_)
            singleton.append(bool(np.bincount(pl_codes).max() == 1))
            n_src.append(len(set(src)))

            # all N_PERM relabelings at once, permuting plate -> source so that
            # plate blocks stay intact (plates per source label are preserved)
            src_codes = pd.factorize(src)[0]
            # one source label per plate. searchsorted returns positions in the
            # SORTED view, so map back through the sorter to get a real row.
            order = np.argsort(pl_codes, kind="stable")
            rep = order[np.searchsorted(pl_codes, np.arange(len(pl_uniq)),
                                        sorter=order)]
            pl_src = src_codes[rep]
            assert (src_codes == pl_src[pl_codes]).all(), \
                "a plate maps to more than one source"
            P = np.argsort(perm_rng.random((N_PERM, len(pl_uniq))), axis=1)
            sp = pl_src[P][:, pl_codes]                      # (N_PERM, n_wells)
            m = sp[:, iu[0]] == sp[:, iu[1]]                 # (N_PERM, n_pairs)
            ns = m.sum(axis=1)
            assert ((ns > 0) & (ns < m.shape[1])).all(), "degenerate permutation split"
            tot = d.sum(); npair = len(d)
            s_sum = (m * d).sum(axis=1)
            perm_rows.append(((tot - s_sum) / (npair - ns) - s_sum / ns).astype(np.float32))
            done += 1
            if done % 2000 == 0:
                log(f"contrast {done:,}/{len(pos):,}")
                np.savez(CKPT / "contrast_partial.npz", same=np.array(same),
                         diff=np.array(diff), delta=np.array(delta),
                         n_src=np.array(n_src), perm=np.vstack(perm_rows), done=done)
        same, diff, delta, n_src = map(np.array, (same, diff, delta, n_src))
        singleton = np.array(singleton)
        perm_delta = np.vstack(perm_rows)                  # (n_compounds, N_PERM)
        perm_mean = np.median(perm_delta, axis=0)           # per-permutation MEDIAN
        np.savez(ck, same=same, diff=diff, delta=delta, n_src=n_src,
                 perm_mean=perm_mean, perm_delta=perm_delta, singleton=singleton)
        stamp_fingerprint(ck, fingerprint(stage="contrast"))
        log(f"contrast checkpointed: {len(delta):,} compounds")

    obs = float(np.median(delta))
    p_perm = (1 + int((perm_mean >= obs).sum())) / (1 + N_PERM)   # median vs median

    # Self-check: for compounds whose every plate holds exactly one well, plate
    # permutation is equivalent to well permutation, so this null must sit near
    # zero. A displaced value means the permutation is still wrong.
    sing_null = (np.median(perm_delta[singleton], axis=0) if singleton.any()
                 else np.array([np.nan]))
    log(f"singleton-plate compounds: {int(singleton.sum()):,}/{len(singleton):,}; "
        f"their permutation null median = {float(np.median(sing_null)):+.5f}")
    log(f"matched contrast: median Delta = {obs:+.4f}, permutation p = {p_perm:.4g}")

    # ---- test 3: ablation on the replicated cohort
    frozen = (OUT_DIR / "ablated_features.txt").read_text().split()
    missing = sorted(set(frozen) - set(feat))
    assert not missing, f"frozen ablation features absent from schema: {missing[:5]}"
    ab = set(frozen)
    keep_mask = np.fromiter((c not in ab for c in feat), bool, len(feat))
    assert int((~keep_mask).sum()) == len(frozen)
    # Per-compound source medians in numpy. pd.DataFrame(XV).groupby().median()
    # over 3,180 columns promotes the result to float64 and allocates far more
    # than the result itself; this holds one (n_sources, 3180) block per compound.
    by_cmpd = {}
    for (jcp, _src), rows in keys.groupby(["Metadata_JCP2022", "Metadata_Source"],
                                          sort=True).indices.items():
        by_cmpd.setdefault(jcp, []).append(np.median(XV[rows], axis=0))
    arrs = [np.vstack(v) for v in by_cmpd.values()
            if len(v) >= MIN_REPLICATED_SOURCES]
    assert all(a.dtype == np.float32 and a.shape[1] == len(feat) for a in arrs)

    def rank_scores(mask=None):
        out = np.empty(len(arrs))
        for i, A in enumerate(arrs):
            B = A[:, mask] if mask is not None else A
            U, _ = unit(B.astype(np.float64))
            out[i] = di_from_unit(U)
        return out

    raw, abl = rank_scores(), rank_scores(keep_mask)
    rho = float(spearmanr(raw, abl).statistic)
    nck = CKPT / "ablation_null.npy"
    null_fp = fingerprint(stage="ablation_null", n_random=N_RANDOM_ABLATIONS)
    check_fingerprint(nck, null_fp)
    null = np.load(nck) if nck.exists() else np.array([])
    if len(null):
        log(f"resumed {len(null)} nulls from checkpoint")
    while len(null) < N_RANDOM_ABLATIONS:
        m = np.ones(len(feat), bool)
        m[abl_rng.choice(len(feat), size=len(frozen), replace=False)] = False
        null = np.append(null, float(spearmanr(raw, rank_scores(m)).statistic))
        if len(null) % 25 == 0:
            np.save(nck, null)
            stamp_fingerprint(nck, null_fp)
            log(f"replicated-cohort null {len(null)}/{N_RANDOM_ABLATIONS}")
    np.save(nck, null)
    stamp_fingerprint(nck, null_fp)
    p_rand = (1 + int((null <= rho).sum())) / (1 + N_RANDOM_ABLATIONS)

    out = {
        "experiment": "08c_matched_well_contrast",
        "cohort": {"min_wells_per_source": MIN_WELLS_PER_SOURCE,
                   "min_replicated_sources": MIN_REPLICATED_SOURCES,
                   "n_compounds_contrast": int(len(delta)),
                   "n_compounds_ablation": len(arrs),
                   "median_sources": float(np.median(n_src))},
        "matched_raw_well_contrast": {
            "same_source_D_median": float(np.median(same)),
            "different_source_D_median": float(np.median(diff)),
            "delta_median": obs,
            "delta_mean": float(delta.mean()),
            "delta_q25": float(np.percentile(delta, 25)),
            "delta_q75": float(np.percentile(delta, 75)),
            "frac_positive": float((delta > 0).mean()),
            "permutation": {"n": N_PERM, "statistic": "median of per-compound Delta",
                            "null_median": float(np.median(perm_mean)),
                            "null_q975": float(np.percentile(perm_mean, 97.5)),
                            "empirical_p": round(p_perm, 5),
                            "singleton_plate_null_median": float(np.median(sing_null)),
                            "n_singleton_plate_compounds": int(singleton.sum()),
                            "singleton_check": ("plate permutation reduces to well "
                                                "permutation for these; the null must "
                                                "sit near zero")},
        },
        "replicated_cohort_ablation": {
            "n_compounds": len(arrs),
            "spearman_raw_vs_ablated": rho,
            "n_features_ablated": len(frozen),
            "random_null_median": float(np.median(null)),
            "random_null_q025": float(np.percentile(null, 2.5)),
            "empirical_p_lower_tail": round(p_rand, 5),
        },
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "matched_well_contrast.json").write_text(json.dumps(out, indent=2))
    pd.DataFrame({"same_source_D": same, "different_source_D": diff,
                  "delta": delta, "n_sources": n_src}).to_csv(
        OUT_DIR / "matched_well_contrast.csv", index=False)

    print(f"\n[{ts()}] same-source D median      {np.median(same):.4f}")
    print(f"[{ts()}] different-source D median {np.median(diff):.4f}")
    print(f"[{ts()}] Delta median {obs:+.4f}, positive in {(delta>0).mean():.1%}, "
          f"permutation p = {p_perm:.4g}")
    log(f"replicated-cohort ablation rho = {rho:.4f} "
          f"(null median {np.median(null):.4f}, p = {p_rand:.4g})")
    print(f"[{ts()}] written to {OUT_DIR}")


if __name__ == "__main__":
    main()
