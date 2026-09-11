"""Pre-run check for 08c: do the rewritten paths return what they replaced?

Two rewrites were made for memory, and a memory rewrite that changes an answer is
worse than the OOM it fixes. Each is checked against the code it replaced.

1. Streaming. `f.read(cols).take(idx)` replaced `f.read(cols).to_pandas().loc[mask]`.
   Checked on real columns from the pinned S3 release.
2. Consensus. A numpy per-group median replaced `pd.DataFrame(XV).groupby().median()`.
   Checked on a synthetic frame with the same group structure.
3. The fingerprint state machine, which previously stamped before computing and so
   marked a dead run's checkpoint as valid.

    uv run --no-project --with pyarrow --with pandas --with numpy \
        python experiments/checks/verify_08c_paths.py
"""
import sys
import numpy as np
import pandas as pd
import pyarrow.fs as fs
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
import importlib
m = importlib.import_module("08c_matched_well_contrast")

fails = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}", flush=True)
    if not ok:
        fails.append(name)


# ---- 1. streaming: arrow take vs pandas loc, on real data
s3 = fs.S3FileSystem(anonymous=True, region="us-east-1")
f = pq.ParquetFile(s3.open_input_file(m.S3_PATH))
feat = [n for n in f.schema_arrow.names if n.startswith(m.FEATURE_PREFIXES)]
check("schema unchanged", len(feat) == m.EXPECTED_FEATURES and
      f.metadata.num_rows == m.EXPECTED_ROWS, f"{len(feat)} feat, {f.metadata.num_rows:,} rows")

rng = np.random.default_rng(0)
mask = np.zeros(f.metadata.num_rows, bool)
mask[rng.choice(f.metadata.num_rows, size=4000, replace=False)] = True   # unsorted draw
cols = feat[:12] + feat[-12:]
tbl = f.read(columns=cols)
old = tbl.to_pandas().loc[mask].to_numpy(np.float32)
new = np.column_stack([c.to_numpy(zero_copy_only=False).astype(np.float32, copy=False)
                       for c in tbl.take(np.flatnonzero(mask)).columns])
check("arrow take == pandas loc", old.shape == new.shape and np.array_equal(old, new),
      f"{new.shape}, max|diff| = {np.abs(old - new).max():.3g}")
del tbl, old, new

# ---- 2. consensus: numpy per-group median vs pandas groupby.median
n, p = 6000, 300
XV = rng.normal(size=(n, p)).astype(np.float32)
keys = pd.DataFrame({
    "Metadata_JCP2022": rng.choice([f"JCP{i:05d}" for i in range(400)], n),
    "Metadata_Source": rng.choice([f"source_{i}" for i in range(10)], n)})

cons = pd.DataFrame(XV).groupby([keys.Metadata_JCP2022, keys.Metadata_Source]).median()
ck = cons.index.get_level_values(0)
old = {j: cons.to_numpy(np.float32)[ck == j] for j in sorted(set(ck))}

by = {}
for (jcp, _src), rows in keys.groupby(["Metadata_JCP2022", "Metadata_Source"],
                                      sort=True).indices.items():
    by.setdefault(jcp, []).append(np.median(XV[rows], axis=0))
new = {j: np.vstack(v) for j, v in by.items()}

check("same compounds", set(old) == set(new), f"{len(new)} compounds")
shapes = all(old[j].shape == new[j].shape for j in old)
worst = max(np.abs(old[j] - new[j]).max() for j in old) if shapes else float("nan")
check("consensus values agree", shapes and worst < 1e-6, f"max|diff| = {worst:.3g}")
check("consensus dtype preserved", all(v.dtype == np.float32 for v in new.values()))

# source ordering within a compound must line up, or DI is computed on a
# different set even when the values match elementwise
order_ok = all(np.array_equal(np.argsort(old[j].sum(1)), np.argsort(new[j].sum(1)))
               for j in old)
check("source order within compound matches", order_ok)

# DI is what actually gets ranked, so compare that rather than the medians alone
def di(A):
    U = A.astype(np.float64)
    U /= np.linalg.norm(U, axis=1, keepdims=True)
    K = U.shape[0]
    return 1.0 - (float(np.square(U.sum(0)).sum()) - K) / (K * (K - 1))

js = [j for j in old if old[j].shape[0] >= m.MIN_REPLICATED_SOURCES]
dmax = max(abs(di(old[j]) - di(new[j])) for j in js)
check("DI agrees on the ranked cohort", dmax < 1e-9, f"n = {len(js)}, max|diff| = {dmax:.3g}")

# ---- 3. fingerprint state machine
import tempfile, pathlib
with tempfile.TemporaryDirectory() as td:
    d = pathlib.Path(td)
    art = d / "x.npy"
    check("absent artifact does not resume", m.check_fingerprint(art, "aaa") is False)
    check("a dead run leaves no stamp",
          not (d / "x.npy.fingerprint").exists() and m.check_fingerprint(art, "aaa") is False)
    np.save(art, np.arange(3)); m.stamp_fingerprint(art, "aaa")
    check("stamped artifact resumes", m.check_fingerprint(art, "aaa") is True)
    try:
        m.check_fingerprint(art, "bbb"); check("config change refuses resume", False)
    except SystemExit:
        check("config change refuses resume", True)
    try:
        m.stamp_fingerprint(d / "missing.npy", "ccc")
        check("refuses to stamp an absent artifact", False)
    except AssertionError:
        check("refuses to stamp an absent artifact", True)

print(f"\n{len(fails)} failure(s)" + (f": {fails}" if fails else ""))
raise SystemExit(1 if fails else 0)
