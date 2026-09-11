"""Experiment 8: morphological feature-family ablation on JUMP-CP compounds.

Does removing a structured block of morphological features reorganize a
direction-instability ranking of compounds across imaging sources?

Naming. The manuscript previously called this a "cell-health correction" and
reported rho = 0.987 on 25,254 compounds. That value was hardcoded in a figure
script; no computation produced it. There is also no published 743-feature
"cell health" column set to remove -- Way et al. (2021) define 70 phenotypic
readouts and train models to predict them, not a list of CellProfiler columns.
This script therefore tests what can actually be tested: whether deleting a
frozen, explicitly enumerated block of size/intensity/granularity/radial
features changes the ranking more than size-matched random deletions do.

What it measures. A JUMP-CP source is a partner imaging site, not a cell type
(cpg0016 compounds are U2OS). D across sources indexes cross-site technical
reproducibility, not context-dependence across cell types.

Why the interpretable profiles. The manifest's `compound` subset is Harmony-
aligned, which both destroys feature identity and removes the between-source
variation this experiment measures. `compound_interpretable`
(profiles_var_mad_int.parquet, 12.1 GB) retains original CellProfiler names and
negative-control MAD normalization only.

    # stage 1 -- eligibility from metadata only (seconds, ~1 MB)
    uv run --no-project --with pyarrow --with pandas python experiments/08_jump_cp_feature_ablation.py --eligible

    # stage 2 -- stream feature columns in blocks, write consensus (~7 min, no download)
    uv run --no-project --with pyarrow --with pandas --with numpy python experiments/08_jump_cp_feature_ablation.py --consensus

    # stage 3 -- ablation, random controls (local, fast)
    uv run --no-project --with pyarrow --with pandas --with numpy --with scipy python experiments/08_jump_cp_feature_ablation.py --analyze

Nothing is downloaded. The 12.1 GB parquet has a single row group, so row
filtering cannot be pushed down, but column pushdown works and medians are
independent per feature: the file is streamed once in column blocks, peaking
around 2 GB of memory, and only the consensus is written to disk. Each block is
committed before the next is read, so an interrupted run resumes.
"""
import argparse
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
CACHE = REPO / "data" / "jump_cp" / "cpg0016"
OUT_DIR = REPO / "results" / "08_jump_cp_feature_ablation"

S3_PATH = ("cellpainting-gallery/cpg0016-jump-assembled/source_all/workspace/"
           "profiles_assembled/COMPOUND/v1.0/profiles_var_mad_int.parquet")
MANIFEST = ("https://raw.githubusercontent.com/jump-cellpainting/datasets/"
            "v0.11.0/manifests/profile_index.json")
SUBSET = "compound_interpretable"
EXPECTED_ETAG = "67212e3cbdaa25de511f318cdd0503dc-3"
EXPECTED_ROWS = 803_853
EXPECTED_FEATURES = 3_180

FEATURE_PREFIXES = ("Cells_", "Cytoplasm_", "Nuclei_")
META_COLS = ["Metadata_Source", "Metadata_JCP2022"]
DMSO_JCP = "JCP2022_033924"          # negative control, present on every plate
PLATE_META = ("https://raw.githubusercontent.com/jump-cellpainting/datasets/"
              "v0.11.0/metadata/plate.csv.gz")
KEEP_PLATE_TYPE = "COMPOUND"         # the assembled table also carries TARGET2 wells,
                                     # whose compounds are deliberately over-replicated
BLOCK = 300                          # feature columns streamed per pass
MIN_SOURCES = 5
CONTROL_PLATE_FRACTION = 0.50   # a compound on >= half of all plates is a control
# Frozen 2026-09-10 from results/qc_control_audit/control_audit.json. Each sits on
# 93.6-100% of eligible COMPOUND plates; the most frequent retained compound sits
# at 6.9%, so the threshold is nowhere near a compound. Recomputed and asserted
# below, so refreshed metadata cannot silently change the cohort.
POSITIVE_CONTROLS = (
    "JCP2022_012818", "JCP2022_025848", "JCP2022_035095", "JCP2022_037716",
    "JCP2022_046054", "JCP2022_050797", "JCP2022_064022", "JCP2022_085227",
)
EXPECTED_CONTROL_EXCLUDED_COHORT = 24_992   # metadata-derived; the run must reproduce it
MIN_WELLS_PER_SOURCE = 1   # JUMP compounds are typically one well per source;
                           # replication is across sources, not within one
MIN_NORM = 1e-8
SEED = 20260908
N_RANDOM_ABLATIONS = 200
PREREG_BITES = 0.70
PREREG_PRESERVES = 0.83

ABLATION_FAMILIES = [
    r"^(Cells|Cytoplasm|Nuclei)_AreaShape_(Area|Compactness|Eccentricity|Extent|FormFactor|Solidity)$",
    r"^(Cells|Cytoplasm|Nuclei)_Intensity_\w+_(DNA|Mito)$",
    r"^(Cells|Cytoplasm|Nuclei)_Granularity_\d+_(DNA|Mito)$",
    r"^(Cells|Cytoplasm|Nuclei)_RadialDistribution_\w+_(DNA|Mito)_\w+$",
]


def ts():
    return datetime.now().strftime("%H:%M:%S")


def fingerprint(**extra):
    """Bind a cached artifact to the configuration that produced it."""
    import hashlib
    fp = {"etag": EXPECTED_ETAG, "plate_filter": "row-level Metadata_PlateType == COMPOUND",
          "min_sources": MIN_SOURCES, "seed": SEED if "SEED" in globals() else None, **extra}
    return hashlib.sha256(json.dumps(fp, sort_keys=True, default=str).encode()).hexdigest()[:16]


def check_fingerprint(path, fp):
    """Refuse to resume from a checkpoint written under a different configuration.

    Inspects the stamp independently of `path`, so it still works when `path` is
    a sentinel that is never created (e.g. a directory-stage marker).
    """
    stamp = path.with_suffix(path.suffix + ".fingerprint")
    if stamp.exists():
        observed = stamp.read_text().strip()
        if observed != fp:
            raise SystemExit(f"stale checkpoint for {path.name}: {observed} != {fp}. "
                             "Delete the checkpoint directory and rerun.")
        return True
    if path.exists():
        raise SystemExit(f"{path.name} exists without a fingerprint; delete and rerun.")
    stamp.parent.mkdir(parents=True, exist_ok=True)
    stamp.write_text(fp)
    return False


def _open(attempts=5):
    """Remote handle on the pinned parquet. Reads nothing until asked.

    S3 HeadObject occasionally times out; retry with backoff rather than
    losing a stage to a transient network error.
    """
    import time
    import pyarrow.fs as fs
    import pyarrow.parquet as pq
    last = None
    for a in range(attempts):
        try:
            s3 = fs.S3FileSystem(anonymous=True, region="us-east-1")
            # verify the object against the pinned manifest before reading it:
            # row count alone cannot distinguish two objects of the same shape
            import json as _json, urllib.request as _u
            with _u.urlopen(MANIFEST) as r:
                entry = next(e for e in _json.load(r) if e["subset"] == SUBSET)
            assert entry["etag"] == EXPECTED_ETAG, (
                f"manifest ETag moved: {entry['etag']} != {EXPECTED_ETAG}")
            req = _u.Request("https://cellpainting-gallery.s3.amazonaws.com/"
                             + S3_PATH.split("/", 1)[1], method="HEAD")
            live = _u.urlopen(req).headers.get("ETag", "").strip('"')
            assert live == EXPECTED_ETAG, (
                f"object ETag {live} != pinned {EXPECTED_ETAG}")
            f = pq.ParquetFile(s3.open_input_file(S3_PATH))
            break
        except OSError as e:                       # noqa: PERF203 - transient network
            last = e
            print(f"[{ts()}] S3 open failed ({a+1}/{attempts}): {e}", flush=True)
            time.sleep(5 * (a + 1))
    else:
        raise last
    assert f.metadata.num_rows == EXPECTED_ROWS, (
        f"row count changed: {f.metadata.num_rows} != {EXPECTED_ROWS}; the pinned "
        "release moved, so re-pin deliberately rather than proceeding")
    return f


def eligible():
    """Stage 1. Metadata columns only -- seconds, about a megabyte."""
    import pandas as pd
    f = _open()
    feat = [n for n in f.schema_arrow.names if n.startswith(FEATURE_PREFIXES)]
    assert len(feat) == EXPECTED_FEATURES, f"{len(feat)} features, expected {EXPECTED_FEATURES}"

    well_meta = f.read(columns=["Metadata_Source", "Metadata_Plate",
                               "Metadata_JCP2022"]).to_pandas()
    plate = pd.read_csv(PLATE_META)
    well_meta = well_meta.merge(
        plate[["Metadata_Source", "Metadata_Plate", "Metadata_PlateType"]],
        on=["Metadata_Source", "Metadata_Plate"], how="left", validate="many_to_one")
    assert well_meta.Metadata_PlateType.notna().all(), "unmatched plate rows in the join"
    assert len(well_meta) == EXPECTED_ROWS, "merge changed row count"

    # row-level restriction, applied before eligibility and before any median
    well_meta["row_ok"] = (well_meta.Metadata_PlateType.eq(KEEP_PLATE_TYPE)
                           & well_meta.Metadata_JCP2022.notna()
                           & well_meta.Metadata_JCP2022.ne(DMSO_JCP))
    CACHE.mkdir(parents=True, exist_ok=True)
    well_meta.to_parquet(CACHE / "well_metadata.parquet")
    print(f"[{ts()}] row filter keeps {int(well_meta.row_ok.sum()):,} of "
          f"{len(well_meta):,} wells", flush=True)
    df = well_meta.loc[well_meta.row_ok, META_COLS].copy()
    top = df.Metadata_JCP2022.value_counts().head(10)

    per_source = df.groupby(["Metadata_JCP2022", "Metadata_Source"]).size()
    per_source = per_source[per_source >= MIN_WELLS_PER_SOURCE]
    n_src = per_source.groupby("Metadata_JCP2022").size()
    keep = sorted(n_src[n_src >= MIN_SOURCES].index)
    single_well = float((per_source[per_source.index.get_level_values(0).isin(set(keep))] == 1).mean())

    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "eligible_compounds.json").write_text(json.dumps({
        "s3_path": S3_PATH, "manifest": MANIFEST, "etag": EXPECTED_ETAG,
        "n_wells": int(len(df)), "n_features": len(feat),
        "n_sources": int(df.Metadata_Source.nunique()),
        "min_sources": MIN_SOURCES, "min_wells_per_source": MIN_WELLS_PER_SOURCE,
        "excluded_control": DMSO_JCP,
        "top_10_by_well_count_before_positive_control_exclusion": {k: int(v) for k, v in top.items()},
        "n_eligible": len(keep),
        "frac_source_consensuses_from_one_well": round(single_well, 4),
        "single_well_note": "computed on the row-filtered, DMSO-excluded, eligible cohort",
        "eligible": keep,
    }, indent=2))
    print(f"[{ts()}] {df.Metadata_Source.nunique()} sources, {len(feat)} features")
    print(f"[{ts()}] {len(keep):,} compounds in >= {MIN_SOURCES} sources")
    print(f"[{ts()}] {single_well:.1%} of source consensuses rest on a single well")
    print(f"[{ts()}] top well counts before control exclusion: {dict(list(top.items())[:3])}")


def consensus():
    """Stage 2. Stream feature columns in blocks; median per compound-source."""
    import pandas as pd
    elig_path = CACHE / "eligible_compounds.json"
    assert elig_path.exists(), "run --eligible first"
    keep = set(json.loads(elig_path.read_text())["eligible"])

    f = _open()
    feat = [n for n in f.schema_arrow.names if n.startswith(FEATURE_PREFIXES)]
    wm = pd.read_parquet(CACHE / "well_metadata.parquet")
    meta = wm[META_COLS]
    mask = (wm.row_ok & wm.Metadata_JCP2022.isin(keep)).to_numpy()
    keys = meta.loc[mask, META_COLS].reset_index(drop=True)
    print(f"[{ts()}] {mask.sum():,} of {len(meta):,} wells belong to eligible compounds")

    BLOCKS = CACHE / "consensus_blocks"; BLOCKS.mkdir(parents=True, exist_ok=True)
    eligible_sha = hashlib.sha256("\n".join(sorted(keep)).encode()).hexdigest()
    fp = fingerprint(stage="consensus", eligible_sha256=eligible_sha, block=BLOCK)
    if any(BLOCKS.glob("block_*.parquet")) and not (BLOCKS / "_stage.fingerprint").exists():
        raise SystemExit("consensus blocks exist without a fingerprint; delete and rerun")
    check_fingerprint(BLOCKS / "_stage", fp)
    n_blocks = (len(feat) + BLOCK - 1) // BLOCK
    for b in range(n_blocks):
        out = BLOCKS / f"block_{b:03d}.parquet"
        if out.exists():
            print(f"[{ts()}] block {b+1}/{n_blocks}: cached"); continue
        cols = feat[b * BLOCK:(b + 1) * BLOCK]
        vals = f.read(columns=cols).to_pandas().loc[mask].reset_index(drop=True)
        assert np.isfinite(vals.to_numpy()).all(), f'nonfinite values in block {b}'
        block = pd.concat([keys, vals], axis=1)
        med = block.groupby(META_COLS, sort=True)[cols].median().astype("float32")
        if b == 0:
            med["n_wells"] = block.groupby(META_COLS, sort=True).size()
        med.reset_index().to_parquet(out)          # committed before the next read
        print(f"[{ts()}] block {b+1}/{n_blocks}: {len(cols)} features, {len(med):,} consensuses")
    print(f"[{ts()}] consensus complete in {BLOCKS}")


def direction_instability(U):
    """1 - mean pairwise cosine over distinct pairs, from unit rows.

    Uses the identity sum_{i!=j} u_i.u_j = ||sum_i u_i||^2 - K, so no K x K
    matrix is built. Rows must already be unit-normalized and finite.
    """
    K = U.shape[0]
    mean_cos = (float(np.square(U.sum(axis=0)).sum()) - K) / (K * (K - 1))
    return 1.0 - min(max(mean_cos, -1.0), 1.0)


def _spearman(a, b):
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).statistic)


def _positive_controls():
    """Plate-level positive controls, recomputed and checked against the frozen list.

    The assembled profile table carries no control-role column, so controls are
    identified by plate coverage. Returning the frozen tuple rather than the
    recomputed set means the cohort is fixed by the code under review, and the
    assert is what catches metadata drift.
    """
    wm = pd.read_parquet(CACHE / "well_metadata.parquet")
    elig = wm[wm.row_ok]
    ppc = elig.groupby("Metadata_JCP2022").Metadata_Plate.nunique()
    frac = ppc / elig.Metadata_Plate.nunique()
    recomputed = tuple(sorted(frac[frac >= CONTROL_PLATE_FRACTION].index))
    assert recomputed == POSITIVE_CONTROLS, (
        f"control list drifted from the frozen one:\n  recomputed {recomputed}\n"
        f"  frozen     {POSITIVE_CONTROLS}")
    return set(POSITIVE_CONTROLS)


def analyze():
    """Stage 3. Frozen ablation, size-matched random null, leave-one-source-out."""
    import pandas as pd

    BLOCKS = CACHE / "consensus_blocks"
    files = sorted(BLOCKS.glob("block_*.parquet"))
    assert files, "run --consensus first"
    cons = pd.concat([pd.read_parquet(f).set_index(META_COLS) for f in files],
                     axis=1).reset_index()
    feature_cols = [c for c in cons.columns if c.startswith(FEATURE_PREFIXES)]

    # --- plate-type restriction: the assembled table also carries TARGET2 wells
    # wells were restricted to KEEP_PLATE_TYPE before aggregation in eligible()/consensus();
    # an aggregate-level filter here could not undo a contaminated median.
    expected_blocks = (EXPECTED_FEATURES + BLOCK - 1) // BLOCK
    assert len(files) == expected_blocks, f"{len(files)} blocks, expected {expected_blocks}"
    assert len(feature_cols) == EXPECTED_FEATURES, f"{len(feature_cols)} features"
    assert len(feature_cols) == len(set(feature_cols)), "duplicate feature columns"
    print(f"[{ts()}] {cons.Metadata_JCP2022.nunique():,} compounds, "
          f"{len(cons):,} compound-source consensuses (wells were filtered before aggregation)")

    # --- plate-level positive controls, removed before any ranking is formed.
    # These sit on ~every plate and are experimental controls, not screened
    # compounds; their replicate counts also dwarf the cohort median.
    controls = _positive_controls()
    n_before = cons.Metadata_JCP2022.nunique()
    cons = cons[~cons.Metadata_JCP2022.isin(controls)]
    n_dropped = n_before - cons.Metadata_JCP2022.nunique()
    assert n_dropped == len(POSITIVE_CONTROLS), (
        f"expected to drop {len(POSITIVE_CONTROLS)} controls, dropped {n_dropped}")
    print(f"[{ts()}] excluded {n_dropped} plate-level positive controls")

    n_src = cons.groupby("Metadata_JCP2022").size()
    cons = cons[cons.Metadata_JCP2022.isin(n_src[n_src >= MIN_SOURCES].index)]
    n_after_source_filter = int(cons.Metadata_JCP2022.nunique())
    assert n_after_source_filter == EXPECTED_CONTROL_EXCLUDED_COHORT, (
        f"cohort is {n_after_source_filter:,}, expected "
        f"{EXPECTED_CONTROL_EXCLUDED_COHORT:,}")

    # --- ablation block, frozen to a committed file before any outcome is computed
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    frozen = OUT_DIR / "ablated_features.txt"
    pat = re.compile("|".join(ABLATION_FAMILIES))
    resolved = sorted(c for c in feature_cols if pat.match(c))
    assert frozen.exists(), ("frozen ablation list missing; refusing to create it "
                             "during an outcome run")
    committed = frozen.read_text().split()
    assert committed == resolved, ("the schema no longer resolves to the committed "
                                   "ablation list; re-freeze deliberately")
    ablated = committed
    fhash = hashlib.sha256("\n".join(ablated).encode()).hexdigest()
    assert ablated and len(ablated) < len(feature_cols)

    ablated_set = set(ablated)
    keep_mask = np.fromiter((c not in ablated_set for c in feature_cols),
                            bool, len(feature_cols))
    print(f"[{ts()}] {len(feature_cols)} features; ablating {len(ablated)} "
          f"(sha256 {fhash[:12]})")

    # --- precompute unit-normalized float32 arrays once
    groups, dropped = [], {"nonfinite": 0, "zero_norm": 0, "too_few": 0}
    for jcp, grp in cons.groupby("Metadata_JCP2022", sort=True):
        X = grp[feature_cols].to_numpy(np.float32)
        if X.shape[0] < 2:
            dropped["too_few"] += 1; continue
        if not np.isfinite(X).all():
            dropped["nonfinite"] += 1; continue
        norms = np.linalg.norm(X, axis=1)
        if (norms < MIN_NORM).any():
            dropped["zero_norm"] += 1; continue
        groups.append((jcp, list(grp.Metadata_Source), X))
    del cons
    print(f"[{ts()}] {len(groups):,} compounds retained; dropped {dropped}")
    assert not any(dropped.values()), dropped
    assert len(groups) == EXPECTED_CONTROL_EXCLUDED_COHORT, (
        f"scored {len(groups):,}, expected {EXPECTED_CONTROL_EXCLUDED_COHORT:,}")

    def score(mask=None):
        out = np.empty(len(groups))
        for i, (_, _, X) in enumerate(groups):
            Y = (X[:, mask] if mask is not None else X).astype(np.float64, copy=False)
            U = Y / np.linalg.norm(Y, axis=1, keepdims=True)
            out[i] = direction_instability(U)
        return out

    raw = score()
    abl = score(keep_mask)
    rho = _spearman(raw, abl)
    print(f"[{ts()}] observed rho = {rho:.4f}")

    # --- size-matched random ablations
    abl_rng = np.random.default_rng(20260908 + 2)
    null = np.empty(N_RANDOM_ABLATIONS)
    for r in range(N_RANDOM_ABLATIONS):
        m = np.ones(len(feature_cols), bool)
        m[abl_rng.choice(len(feature_cols), size=len(ablated), replace=False)] = False
        null[r] = _spearman(raw, score(m))
        if (r + 1) % 25 == 0:
            print(f"[{ts()}] null {r+1}/{N_RANDOM_ABLATIONS}")
    p_random = (1 + int((null <= rho).sum())) / (1 + N_RANDOM_ABLATIONS)

    # --- leave one source out
    all_sources = sorted({s for _, srcs, _ in groups for s in srcs})
    loo = {}
    for s in all_sources:
        idx = [(np.array([j for j, v in enumerate(srcs) if v != s]), i)
               for i, (_, srcs, _) in enumerate(groups)]
        keep_i = [(rows, i) for rows, i in idx if len(rows) >= MIN_SOURCES - 1]
        if len(keep_i) < 100:
            loo[s] = {"n": len(keep_i), "rho": None}; continue
        r_, a_ = np.empty(len(keep_i)), np.empty(len(keep_i))
        for k, (rows, i) in enumerate(keep_i):
            X = groups[i][2][rows].astype(np.float64, copy=False)
            U = X / np.linalg.norm(X, axis=1, keepdims=True)
            r_[k] = direction_instability(U)
            Y = X[:, keep_mask]
            V = Y / np.linalg.norm(Y, axis=1, keepdims=True)
            a_[k] = direction_instability(V)
        loo[s] = {"n": len(keep_i), "rho": round(_spearman(r_, a_), 4),
                  "note": "fixed cohort; >= MIN_SOURCES-1 remaining sources"}
        assert len(keep_i) == EXPECTED_CONTROL_EXCLUDED_COHORT, (
            f"LOO fold {s} holds {len(keep_i):,}, expected "
            f"{EXPECTED_CONTROL_EXCLUDED_COHORT:,}")
        print(f"[{ts()}] LOO {s}: n={len(keep_i):,} rho={loo[s]['rho']}")

    pr_raw = pd.Series(raw).rank(pct=True)
    pr_abl = pd.Series(abl).rank(pct=True)
    shift = (pr_raw - pr_abl).abs()

    out = {
        "experiment": "08_jump_cp_feature_ablation",
        "what_it_measures": ("cross-site technical reproducibility of a compound ranking "
                             "under deletion of a frozen morphological feature block; "
                             "JUMP-CP sources are imaging sites, not cell types"),
        "analysis_status": {
            "registration": "unregistered sensitivity analysis",
            "role": ("descriptive assessment of ranking sensitivity to deletion of "
                     "a frozen morphological feature block"),
        },
        "historical_decision_rule": {
            "bites_below": PREREG_BITES, "preserves_at_or_above": PREREG_PRESERVES,
            "interpretive_role": ("reported for continuity only; not treated as "
                                  "inferential evidence"),
        },
        "cohort": {"min_sources": MIN_SOURCES, "plate_type": KEEP_PLATE_TYPE},
        "positive_controls_excluded": {
            "n": len(POSITIVE_CONTROLS),
            "ids": list(POSITIVE_CONTROLS),
            "rule": (f"compound present on >= {CONTROL_PLATE_FRACTION:.0%} of "
                     "eligible COMPOUND plates"),
            "note": ("not a registered exclusion; a correctness repair, since "
                     "these are experimental controls rather than screened "
                     "compounds"),
        },
        "cohort_flow": {
            "before_positive_control_exclusion": n_before,
            "positive_controls_excluded": n_dropped,
            "after_positive_control_exclusion": n_before - n_dropped,
            "after_source_filter": n_after_source_filter,
            "scored": len(groups),
        },
        "n_compounds": len(groups),
        "n_features_total": len(feature_cols),
        "n_features_ablated": len(ablated),
        "ablated_feature_sha256": fhash,
        "dropped": dropped,
        "spearman_raw_vs_ablated": float(rho),
        "spearman_display": round(float(rho), 4),
        "random_ablation_null": {
            "n": N_RANDOM_ABLATIONS,
            "median": float(np.median(null)),
            "q025": float(np.percentile(null, 2.5)),
            "q05": float(np.percentile(null, 5)),
            "q975": float(np.percentile(null, 97.5)),
            "n_as_or_more_disruptive": int((null <= rho).sum()),
            "monte_carlo_denominator": N_RANDOM_ABLATIONS + 1,
            "seed": 20260908 + 2,
            "numpy_version": np.__version__,
            "empirical_p_lower_tail": round(p_random, 4),
        },
        "leave_one_source_out": loo,
        "percentile_rank_shift": {
            "median": round(float(shift.median()), 5),
            "p90": round(float(shift.quantile(0.90)), 5),
            "p99": round(float(shift.quantile(0.99)), 5),
            "frac_over_1pct": round(float((shift > 0.01).mean()), 4),
            "frac_over_5pct": round(float((shift > 0.05).mean()), 4),
        },
        "historical_threshold_classification": (
            "bites" if rho < PREREG_BITES else
            "preserves" if rho >= PREREG_PRESERVES else "indeterminate"),
    }
    (OUT_DIR / "feature_ablation_results.json").write_text(json.dumps(out, indent=2))
    pd.DataFrame({"compound": [g[0] for g in groups],
                  "n_sources": [len(g[1]) for g in groups],
                  "raw_instability": raw,
                  "ablated_instability": abl}).to_csv(
        OUT_DIR / "compound_instability.csv", index=False)

    print(f"\n[{ts()}] n = {out['n_compounds']:,}   "
          f"raw-vs-ablated Spearman rho = {rho:.4f}")
    print(f"[{ts()}] random null: median {np.median(null):.4f}, "
          f"2.5-97.5 pct [{np.percentile(null,2.5):.4f}, {np.percentile(null,97.5):.4f}], "
          f"empirical p = {p_random:.4f}")
    rr = [v['rho'] for v in loo.values() if v['rho'] is not None]
    if rr:
        print(f"[{ts()}] leave-one-source-out rho: {min(rr):.4f} to {max(rr):.4f}")
    print(f"[{ts()}] written to {OUT_DIR}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--eligible", action="store_true")
    g.add_argument("--consensus", action="store_true")
    g.add_argument("--analyze", action="store_true")
    a = ap.parse_args()
    eligible() if a.eligible else consensus() if a.consensus else analyze()
