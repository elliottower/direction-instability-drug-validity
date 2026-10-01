"""Modal wrapper: rebuild the LINCS matrices and emit the H3 cohort bundle.

Implements the rebuild specified in PREREG_H3_MAGNITUDE_AND_SHARED_AXIS.md
(frozen at f288507). Logic that computes anything lives in the pinned loader and
in 03c_h3_sensitivity.py; this file arranges workers, the volume, checkpoints and
the fail-closed checks.

Both matrices come from one pinned GCTX, restricted to the same frozen order of
978 landmark genes. The drug consensus uses the pinned loader's selection rule
rather than a reimplementation of it. Every shard carries the stage fingerprint
that produced it and is refused if the fingerprint moves.

    modal run --detach experiments/modal_h3_rebuild.py --stage all
"""
import modal

app = modal.App("di-h3-sensitivity")
vol = modal.Volume.from_name("di-h3", create_if_missing=True)

LOADER = "/Users/elliottower/Documents/GitHub/drug-perturbation-geometry/data/lincs_loader.py"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.1.3", "pandas==2.2.3", "scipy==1.14.1", "cmapPy==4.0.1",
                 "h5py==3.12.1")
    .add_local_dir("experiments", remote_path="/app/experiments")
    .add_local_file(LOADER, remote_path="/app/lincs_loader.py")
)

COMMON = dict(image=image, volumes={"/vol": vol}, timeout=86400, memory=65536, cpu=8.0)

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl"
GCTX_GZ = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz"
GCTX = GCTX_GZ.removesuffix(".gz")
N_LANDMARK = 978
MIN_HAIRPINS = 3
MIN_COMMON_RECORDS = 700
# the loader at the registered commit drug-perturbation-geometry@1dc20a2
EXPECTED_LOADER_SHA256 = "b7ec2cd46be4a800dc546594e82491f785fdf504d4bd8bc783be293b343a1611"
SHARD = 250

# pinned in the registration, verified rather than recorded. The GCTX is the one
# input whose checksum cannot exist before first retrieval; it is recorded once
# and required to match thereafter.
PINS = {
    "GSE92742_Broad_LINCS_sig_info.txt.gz":
        "19da29c0ee12ddf27f9698cd0da40beaff58657dcde9d382aae068737e831299",
    "GSE92742_Broad_LINCS_gene_info.txt.gz": None,     # stamped on first run
    "lincs_shrna_siginfo.csv.gz":
        "bd396fa0e1a2f00c1b5f2c8d2b35f9a056f5e5353382475655869038037ec014",
    # the corrected artifact, which the amendment frozen at 7f57136 makes the
    # reference; it replaces the superseded 65e5d10e2720... of f288507
    "phenotype_projection_results.json":
        "fd69e26fc9a3917323065b631688baeab8b283f735c8bf5b16210ba67bd21425",
}


def _sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _sha256_text(text):
    import hashlib
    return hashlib.sha256(text.encode()).hexdigest()


def _ts():
    from datetime import datetime
    return datetime.now().strftime("%H:%M:%S")


def _landmark_order(raw):
    """The 978 landmark gene ids, in one deterministic frozen order."""
    import pandas as pd
    g = pd.read_csv(raw / "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="\t", low_memory=False)
    lm = g[g.pr_is_lm == 1].sort_values("pr_gene_id")
    ids = [str(i) for i in lm.pr_gene_id]
    assert len(ids) == N_LANDMARK, f"{len(ids)} landmark genes, expected {N_LANDMARK}"
    return ids, _sha256_text("\n".join(ids))


def _fingerprint(raw):
    """Everything a shard depends on. A shard from a different one is refused."""
    import json
    _, order_hash = _landmark_order(raw)
    parts = {
        "gctx": (raw / "gctx.sha256").read_text().strip(),
        "sig_info": _sha256(raw / "GSE92742_Broad_LINCS_sig_info.txt.gz"),
        "shrna_sig_info": _sha256(raw / "lincs_shrna_siginfo.csv.gz"),
        "deposited": _sha256(raw / "phenotype_projection_results.json"),
        "landmark_order": order_hash,
        "loader": _sha256("/app/lincs_loader.py"),
        "extract_code": _sha256("/app/experiments/modal_h3_rebuild.py"),
        "shard_size": SHARD,
        "min_hairpins": MIN_HAIRPINS,
    }
    return parts, _sha256_text(json.dumps(parts, sort_keys=True))


@app.function(**COMMON)
def stage_fetch():
    """Download, decompress and verify. Partial files never become inputs."""
    import gzip
    import json
    import shutil
    import urllib.request
    from pathlib import Path

    raw = Path("/vol/raw"); raw.mkdir(parents=True, exist_ok=True)

    gz = raw / GCTX_GZ
    if not gz.exists():
        part = gz.with_suffix(gz.suffix + ".part")
        print(f"[{_ts()}] downloading {GCTX_GZ}")
        urllib.request.urlretrieve(f"{GEO}/{GCTX_GZ}", part)
        part.replace(gz)                                  # only a complete file is named
        vol.commit()

    gctx = raw / GCTX
    if not gctx.exists():
        part = gctx.with_suffix(gctx.suffix + ".part")
        print(f"[{_ts()}] decompressing to {GCTX}")
        with gzip.open(gz, "rb") as src, open(part, "wb") as dst:
            shutil.copyfileobj(src, dst)
        part.replace(gctx)
        vol.commit()
    with open(gctx, "rb") as fh:
        assert fh.read(8) == b"\x89HDF\r\n\x1a\n", "decompressed GCTX is not HDF5"

    digest = _sha256(gctx)
    stamp = raw / "gctx.sha256"
    if stamp.exists():
        assert stamp.read_text().strip() == digest, (
            f"GCTX checksum moved: {digest} against recorded {stamp.read_text().strip()}")
    else:
        stamp.write_text(digest)

    missing = [n for n in PINS if not (raw / n).exists()]
    assert not missing, f"upload these inputs to /vol/raw first: {missing}"
    recorded = {}
    for name, want in PINS.items():
        got = _sha256(raw / name)
        if want is not None:
            assert got == want, f"{name} sha256 {got}, pinned {want}"
            continue
        # not available locally at freeze: stamped on first retrieval, then required
        stamp_i = raw / f"{name}.sha256"
        if stamp_i.exists():
            assert stamp_i.read_text().strip() == got, (
                f"{name} sha256 moved: {got} against stamped {stamp_i.read_text().strip()}")
        else:
            stamp_i.write_text(got)
        recorded[name] = got

    assert _sha256("/app/lincs_loader.py") == EXPECTED_LOADER_SHA256, (
        "the loader in the image is not the one pinned at 1dc20a2")

    ids, order_hash = _landmark_order(raw)
    parts, fp = _fingerprint(raw)
    (raw / "stage_fingerprint.json").write_text(json.dumps({"parts": parts, "fingerprint": fp},
                                                           indent=2))
    vol.commit()
    print(json.dumps({"gctx_sha256": digest, "landmark_order_sha256": order_hash,
                      "recorded_on_first_run": recorded, "stage_fingerprint": fp}, indent=2))
    return fp


@app.function(**COMMON)
def stage_extract():
    """Per-drug consensus signatures, restricted to the landmarks, checkpointed per shard."""
    import json
    import sys
    from pathlib import Path

    import numpy as np
    from cmapPy.pandasGEXpress import parse

    sys.path.insert(0, "/app")
    import lincs_loader                                    # the pinned loader

    raw, shards = Path("/vol/raw"), Path("/vol/shards")
    shards.mkdir(parents=True, exist_ok=True)
    ids, _ = _landmark_order(raw)
    _, fp = _fingerprint(raw)

    deposited = json.loads((raw / "phenotype_projection_results.json").read_text())
    wanted = [r["drug"] for r in deposited]
    siginfo = lincs_loader.load_siginfo(raw / "GSE92742_Broad_LINCS_sig_info.txt.gz",
                                        compound_only=True)
    print(f"[{_ts()}] {len(wanted):,} deposited drugs; fingerprint {fp[:12]}")

    (shards / "landmark_gene_ids.json").write_text(json.dumps(ids))
    vol.commit()

    n_shards = (len(wanted) + SHARD - 1) // SHARD
    for s in range(n_shards):
        out = shards / f"shard_{s:03d}.npz"
        if out.exists():
            with np.load(out, allow_pickle=True) as z:
                seen = str(z["fingerprint"]) if "fingerprint" in z.files else ""
            if seen == fp:
                print(f"[{_ts()}] shard {s+1}/{n_shards}: cached")
                continue
            print(f"[{_ts()}] shard {s+1}/{n_shards}: fingerprint moved, rebuilding")
        names = wanted[s * SHARD:(s + 1) * SHARD]
        sub = siginfo[siginfo.pert_iname.isin(names) & siginfo.pert_iname.notna()]
        cells_by_drug = {}
        gct = parse.parse(str(raw / GCTX), cid=sorted(set(sub.sig_id.astype(str))), rid=ids)
        assert gct.data_df.shape[0] == N_LANDMARK, f"parsed {gct.data_df.shape[0]} rows"
        assert gct.data_df.columns.is_unique, "parsed matrix has duplicate signature ids"
        df = gct.data_df
        df.index = df.index.astype(str)
        df = df.reindex(index=ids)                 # rid= selects, it does not order
        assert not df.isna().any().any(), "a landmark gene is missing from the parse"
        assert list(df.index) == ids, "landmark order is not the frozen order"
        mat = df.T
        payload = {}
        # The deposited artifact averages every signature for a drug-cell pair
        # (03_phenotype_projection.py), rather than selecting one per cell line
        # as lincs_loader.get_consensus_signatures would. This file follows the
        # artifact it must reproduce; Deviation 8 records the difference. The
        # sig-info column that rule needs, distil_ss, is absent here in any case.
        for drug, grp in sub.groupby("pert_iname"):
            per_cell, cell_ids = [], []
            for cell_id, cell_rows in grp.groupby("cell_id"):
                sids = [x for x in sorted(set(cell_rows.sig_id.astype(str))) if x in mat.index]
                if sids:
                    per_cell.append(np.vstack([mat.loc[x].to_numpy(np.float64)
                                               for x in sids]).mean(axis=0))
                    cell_ids.append(str(cell_id))
            if per_cell:
                payload[drug] = np.vstack(per_cell)
                cells_by_drug[drug] = cell_ids
        # The cell lines travel in one JSON blob rather than one array per drug.
        # A suffixed key does not survive the npz round trip intact, and a key
        # that collides with a drug name silently replaces that drug's matrix.
        np.savez_compressed(out, fingerprint=np.array(fp),
                            __cells__=np.array(json.dumps(cells_by_drug)), **payload)
        vol.commit()                                       # checkpoint inside the unit
        print(f"[{_ts()}] shard {s+1}/{n_shards}: {len(payload)} drugs")
    return f"{n_shards} shards at {fp[:12]}"


@app.function(**COMMON)
def stage_shrna():
    """Rebuild the shRNA consensus directions from the same GCTX and landmark order."""
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd
    from cmapPy.pandasGEXpress import parse

    raw = Path("/vol/raw")
    ids, order_hash = _landmark_order(raw)
    info = pd.read_csv(raw / "lincs_shrna_siginfo.csv.gz")
    deposited = json.loads((raw / "phenotype_projection_results.json").read_text())
    targets = sorted({r["target"] for r in deposited})
    info = info[info.pert_iname.isin(targets)]

    counts = info.groupby("pert_iname").sig_id.nunique()
    usable = sorted(counts[counts >= MIN_HAIRPINS].index)
    print(f"[{_ts()}] {len(usable):,} of {len(targets):,} targets carry >= {MIN_HAIRPINS} hairpins")

    all_sids = sorted(set(info["sig_id"].astype(str)))
    gct = parse.parse(str(raw / GCTX), cid=all_sids, rid=ids)
    assert gct.data_df.shape[0] == N_LANDMARK, f"parsed {gct.data_df.shape[0]} rows"
    df = gct.data_df
    df.index = df.index.astype(str)
    df = df.reindex(index=ids)                     # identical ordering to the drug parse
    assert not df.isna().any().any(), "a landmark gene is missing from the parse"
    assert list(df.index) == ids, "landmark order is not the frozen order"
    mat = df.T

    assert mat.index.is_unique, "parsed shRNA matrix has duplicate signature ids"
    names, dirs, hairpins = [], [], {}
    for gene in usable:
        # eligibility above counted distinct ids, so the mean must too: a
        # duplicated metadata row would otherwise weight one hairpin twice
        sids = sorted(set(info.loc[info.pert_iname == gene, "sig_id"].astype(str)))
        used = [s for s in sids if s in mat.index]
        if len(used) < MIN_HAIRPINS:
            continue
        rows = [mat.loc[s].to_numpy(np.float64) for s in used]
        v = np.vstack(rows).mean(axis=0)
        norm = np.linalg.norm(v)
        assert np.isfinite(v).all() and norm > 0, f"{gene}: degenerate consensus"
        names.append(gene); dirs.append(v / norm)
        hairpins[gene] = {"n": len(used), "sig_ids": used}

    _, fp = _fingerprint(raw)
    out = Path("/vol/results"); out.mkdir(parents=True, exist_ok=True)
    tmp = out / "shrna_consensus.npz.part"
    with open(tmp, "wb") as fh:
        # the directions carry the gene axis they were built on, so a comparison
        # never has to assume their coordinate system
        np.savez_compressed(fh, genes=np.array(names), directions=np.array(dirs),
                            gene_ids=np.array(ids), fingerprint=np.array(fp))
    tmp.replace(out / "shrna_consensus.npz")
    (out / "shrna_hairpins.json").write_text(json.dumps(hairpins, indent=2, sort_keys=True))
    vol.commit()
    print(json.dumps({"targets_with_consensus": len(names),
                      "landmark_order_sha256": order_hash,
                      "consensus_sha256": _sha256(out / "shrna_consensus.npz"),
                      "hairpins_sha256": _sha256(out / "shrna_hairpins.json")}, indent=2))
    return f"{len(names)} consensus directions"


@app.function(**COMMON)
def stage_bundle():
    """Assemble the cohort bundle, failing closed before the statistics ever see it."""
    import json
    from pathlib import Path

    import numpy as np

    raw, shards, out = Path("/vol/raw"), Path("/vol/shards"), Path("/vol/results")
    _, fp = _fingerprint(raw)

    deposited = json.loads((raw / "phenotype_projection_results.json").read_text())

    # An absent shard is missing computation, not reconstruction attrition. The
    # 700 floor must never be satisfied by work that did not run.
    n_shards = (len(deposited) + SHARD - 1) // SHARD
    expected = {shards / f"shard_{s:03d}.npz" for s in range(n_shards)}
    observed = set(shards.glob("shard_*.npz"))
    assert observed == expected, {
        "missing": sorted(p.name for p in expected - observed),
        "unexpected": sorted(p.name for p in observed - expected)}

    sigs, cells = {}, {}
    for f in sorted(expected):
        with np.load(f, allow_pickle=True) as z:
            assert str(z["fingerprint"]) == fp, f"{f.name} was built under a different fingerprint"
            keys = {k for k in z.files if k not in ("fingerprint", "__cells__")}
            overlap = set(sigs) & keys
            assert not overlap, f"drug identifiers appear in more than one shard: {sorted(overlap)}"
            sigs.update({k: z[k] for k in keys})
            shard_cells = json.loads(str(z["__cells__"]))
            assert set(shard_cells) == keys, (
                f"{f.name}: {len(keys)} drugs but cell lines for {len(shard_cells)}")
            cells.update({drug: [str(c) for c in ids] for drug, ids in shard_cells.items()})
    for drug, matrix in sigs.items():
        assert matrix.shape[0] == len(cells[drug]), (
            f"{drug}: {matrix.shape[0]} rows against {len(cells[drug])} cell lines")

    with np.load(out / "shrna_consensus.npz", allow_pickle=True) as z:
        assert str(z["fingerprint"]) == fp, "the shRNA consensus predates the current fingerprint"
        consensus = {str(g): d for g, d in zip(z["genes"], z["directions"])}
    drugs, targets, mats, dirs = [], [], [], []
    dropped = {"no_signatures": [], "no_consensus": []}
    for rec in deposited:
        d, t = rec["drug"], rec["target"]
        if d not in sigs:
            dropped["no_signatures"].append(d); continue
        if t not in consensus:
            dropped["no_consensus"].append(d); continue
        drugs.append(d); targets.append(t)
        mats.append(np.asarray(sigs[d], dtype=np.float64)); dirs.append(consensus[t])

    assert len(drugs) == len(set(drugs)), "drug identifiers are not unique"
    assert len(drugs) == len(targets) == len(mats) == len(dirs), "ragged bundle"
    assert len(drugs) >= MIN_COMMON_RECORDS, (
        f"{len(drugs)} common records, registration requires {MIN_COMMON_RECORDS}")
    for d, m, u in zip(drugs, mats, dirs):
        assert m.ndim == 2 and m.shape[1] == N_LANDMARK, f"{d}: signatures are {m.shape}"
        assert m.shape[0] >= 5, f"{d}: {m.shape[0]} contexts"
        assert np.isfinite(m).all() and np.isfinite(u).all(), f"{d}: nonfinite"
        assert u.shape == (N_LANDMARK,), f"{d}: direction is {u.shape}"
        assert np.isclose(np.linalg.norm(u), 1.0, atol=1e-10), f"{d}: direction is not unit"
        assert np.all(np.linalg.norm(m, axis=1) > 0), f"{d}: zero-norm signature"

    np.savez_compressed(out / "cohort_bundle.npz", drugs=np.array(drugs),
                        targets=np.array(targets),
                        signatures=np.array(mats, dtype=object),
                        directions=np.array(dirs),
                        cell_lines=np.array([cells.get(d, []) for d in drugs], dtype=object))
    from importlib.metadata import version
    parts, _ = _fingerprint(raw)
    manifest = {
        "runtime_packages": {n: version(n) for n in
                             ("numpy", "pandas", "scipy", "cmapPy", "h5py")},
        "stage_fingerprint": fp, "fingerprint_parts": parts,
        "n_bundled": len(drugs), "n_deposited": len(deposited),
        "dropped_counts": {k: len(v) for k, v in dropped.items()},
        "dropped_identifiers": dropped,
        "cohort_identifier_sha256": _sha256_text("\n".join(sorted(drugs))),
        "bundle_sha256": _sha256(out / "cohort_bundle.npz"),
        "shrna_consensus_sha256": _sha256(out / "shrna_consensus.npz"),
        "shrna_hairpins_sha256": _sha256(out / "shrna_hairpins.json"),
    }
    (out / "rebuild_manifest.json").write_text(json.dumps(manifest, indent=2))
    vol.commit()
    print(json.dumps({k: manifest[k] for k in
                      ("n_bundled", "n_deposited", "dropped_counts", "stage_fingerprint")}, indent=2))
    return f"bundle: {len(drugs)} drugs"


@app.local_entrypoint()
def main(stage: str = "all"):
    """stage: all | fetch | extract | shrna | bundle"""
    if stage in ("all", "fetch"):
        print(stage_fetch.remote())
    if stage in ("all", "extract"):
        print(stage_extract.remote())
    if stage in ("all", "shrna"):
        print(stage_shrna.remote())
    if stage in ("all", "bundle"):
        print(stage_bundle.remote())
    print("\nretrieve, then run the registered analysis:")
    print("    modal volume get di-h3 results ./results/03c_h3_sensitivity/")
    print("    uv run --no-project --with numpy --with scipy python \\")
    print("      experiments/03c_h3_sensitivity.py \\")
    print("      --bundle results/03c_h3_sensitivity/cohort_bundle.npz \\")
    print("      --manifest results/03c_h3_sensitivity/rebuild_manifest.json")
