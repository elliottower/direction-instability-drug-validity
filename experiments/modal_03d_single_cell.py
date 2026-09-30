"""Modal wrapper: R0.1, R0.3, R0.4 and R0.6 on the released single-cell data.

Registered in `experiments/PREREG_H3_REFERENCE_DISCORDANCE.md` (frozen 7f57136).
Everything that computes a number lives in `geometry/single_cell.py`; this file
arranges workers, the volume, checkpoints and the fail-closed checks.

The two normalized single-cell releases are 10.7 GB and 8.7 GB, so each stage
checkpoints inside its unit: targets are processed in batches, every batch is
written and the volume committed, and a rerun skips the targets already done.

    modal run --detach experiments/modal_03d_single_cell.py --stage all
    modal volume get di-h3-sc results/single_cell_audit.json \\
        ./results/03d_h3_reference_discordance/
"""
import modal

app = modal.App("di-h3-single-cell")
vol = modal.Volume.from_name("di-h3-sc", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.1.3", "pandas==2.2.3", "anndata==0.11.4", "h5py==3.12.1",
                 "scipy==1.14.1")
    .add_local_dir("geometry", remote_path="/app/geometry")
)

COMMON = dict(image=image, volumes={"/vol": vol}, timeout=86400, memory=131072, cpu=8.0)

FIGSHARE = "https://ndownloader.figshare.com/files"
RELEASES = {
    "K562_essential": {"single_cell": ("35773075", "f1e221fbf6eac774c21c4242ed440c3f"),
                       "bulk": ("35780870", "30496767641cd2e660ee6ecb5baee132")},
    "RPE1_essential": {"single_cell": ("35775554", "2c36a053960f3fae157adacdbccd4485"),
                       "bulk": ("35775512", "6f1e7d6a09e2f869759e3c4526b7f171")},
}
BATCH = 100                      # targets per checkpoint
SEED_SPLIT = 20260924            # registered
SEED_SAMPLE = 20260925           # registered
SAMPLE_CELLS = 10_000
CONTROLS = ("non-targeting", "control", "")


def _ts():
    from datetime import datetime
    return datetime.now().strftime("%H:%M:%S")


def _sha256(path):
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _md5(path):
    import hashlib
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _gem_groups(adata):
    """The batch label, under whichever name the release uses."""
    for column in ("gem_group", "gemgroup", "batch", "sample"):
        if column in adata.obs:
            return adata.obs[column].to_numpy(), column
    raise AssertionError(f"no batch column among {list(adata.obs.columns)}")


def _targets(adata):
    for column in ("gene", "gene_symbol", "perturbation"):
        if column in adata.obs:
            return adata.obs[column].astype(str).to_numpy(), column
    raise AssertionError(f"no target column among {list(adata.obs.columns)}")


@app.function(**COMMON)
def stage_fetch():
    """Download and verify. A partial file never becomes an input."""
    import json
    import urllib.request
    from pathlib import Path

    raw = Path("/vol/raw")
    raw.mkdir(parents=True, exist_ok=True)
    recorded = {}
    for release, files in RELEASES.items():
        for kind, (file_id, md5) in files.items():
            target = raw / f"{release}_{kind}.h5ad"
            if not target.exists():
                partial = target.with_suffix(".part")
                print(f"[{_ts()}] downloading {release} {kind}")
                urllib.request.urlretrieve(f"{FIGSHARE}/{file_id}", partial)
                partial.replace(target)
                vol.commit()
            got = _md5(target)
            assert got == md5, f"{target.name}: md5 {got}, published {md5}"
            recorded[target.name] = {"md5": got, "sha256": _sha256(target)}
    (raw / "checksums.json").write_text(json.dumps(recorded, indent=2))
    vol.commit()
    return recorded


@app.function(**COMMON)
def stage_audit():
    """R0.1 representation and R0.3 design, for both releases."""
    import json
    import sys
    from pathlib import Path

    import anndata as ad
    import numpy as np

    sys.path.insert(0, "/app")
    from geometry.single_cell import group_design, representation_audit

    out = {}
    for release in RELEASES:
        path = Path("/vol/raw") / f"{release}_single_cell.h5ad"
        adata = ad.read_h5ad(path, backed="r")
        rng = np.random.default_rng(SEED_SAMPLE)
        rows = np.sort(rng.choice(adata.shape[0], size=min(SAMPLE_CELLS, adata.shape[0]),
                                  replace=False))
        sample = adata[rows].to_memory().X
        sample = sample.toarray() if hasattr(sample, "toarray") else np.asarray(sample)
        labels, label_column = _targets(adata)
        groups, group_column = _gem_groups(adata)
        out[release] = {
            "shape": list(adata.shape),
            "representation": representation_audit(sample),
            "columns": {"target": label_column, "batch": group_column},
            "layers": list(adata.layers.keys()), "has_raw": adata.raw is not None,
            "design": group_design(labels, groups, CONTROLS)}
        adata.file.close()
        print(f"[{_ts()}] audited {release}")
    path = Path("/vol/results")
    path.mkdir(parents=True, exist_ok=True)
    (path / "audit.json").write_text(json.dumps(out, indent=2))
    vol.commit()
    return {release: out[release]["representation"] for release in out}


@app.function(**COMMON)
def stage_targets():
    """R0.4 verification and R0.6 reliability, one batch of targets at a time."""
    import json
    import sys
    from pathlib import Path

    import anndata as ad
    import numpy as np

    sys.path.insert(0, "/app")
    from geometry.single_cell import cosine, split_half_reliability

    results = Path("/vol/results")
    results.mkdir(parents=True, exist_ok=True)
    for release in RELEASES:
        shard_dir = results / release
        shard_dir.mkdir(parents=True, exist_ok=True)
        single_cell = ad.read_h5ad(Path("/vol/raw") / f"{release}_single_cell.h5ad", backed="r")
        bulk = ad.read_h5ad(Path("/vol/raw") / f"{release}_bulk.h5ad")
        bulk_gene = {str(name).split("_")[1]: i for i, name in enumerate(bulk.obs_names)
                     if str(name).count("_") >= 2}

        labels, _ = _targets(single_cell)
        groups, _ = _gem_groups(single_cell)
        targets = sorted({t for t in np.unique(labels) if t not in CONTROLS})
        done = {path.stem for path in shard_dir.glob("*.json")}
        print(f"[{_ts()}] {release}: {len(targets)} targets, {len(done)} batches already done")

        for start in range(0, len(targets), BATCH):
            name = f"batch_{start // BATCH:04d}"
            if name in done:
                continue
            batch = {}
            for target in targets[start:start + BATCH]:
                rows = np.flatnonzero(labels == target)
                if len(rows) < 2:
                    continue
                cells = single_cell[rows].to_memory().X
                cells = cells.toarray() if hasattr(cells, "toarray") else np.asarray(cells)
                entry = {"n_cells": int(len(rows)),
                         "reliability": split_half_reliability(cells, groups[rows], SEED_SPLIT)}
                if target in bulk_gene:
                    released = np.asarray(bulk.X[bulk_gene[target]], dtype=np.float64).ravel()
                    mean = cells.mean(axis=0)
                    if len(released) == len(mean):
                        entry["cosine_with_released_bulk"] = cosine(mean, released)
                        entry["max_abs_difference"] = float(np.abs(mean - released).max())
                batch[target] = entry
            (shard_dir / f"{name}.json").write_text(json.dumps(batch, indent=2))
            vol.commit()                      # checkpoint inside the unit
            print(f"[{_ts()}] {release}: {name} written, {len(batch)} targets")
        single_cell.file.close()
    return "targets done"


@app.function(**COMMON)
def stage_merge():
    """One audit file, in the shape `03d_h3_reference_discordance.py` expects."""
    import json
    from pathlib import Path

    import numpy as np

    results = Path("/vol/results")
    audit = json.loads((results / "audit.json").read_text())
    reliability, verification = {}, {}
    naming = {"K562_essential": "C1-K562", "RPE1_essential": "C1-RPE1"}
    for release, name in naming.items():
        merged = {}
        for shard in sorted((results / release).glob("*.json")):
            merged.update(json.loads(shard.read_text()))
        reliability[name] = {target: entry["reliability"] for target, entry in merged.items()
                             if entry.get("reliability") is not None}
        cosines = [entry["cosine_with_released_bulk"] for entry in merged.values()
                   if "cosine_with_released_bulk" in entry]
        verification[name] = {
            "n_targets_compared": len(cosines),
            "median_cosine": float(np.median(cosines)) if cosines else None,
            "min_cosine": float(np.min(cosines)) if cosines else None,
            "threshold": 0.99,
            "describes": ("C1" if cosines and float(np.median(cosines)) >= 0.99
                          else "the single-cell construction, not C1"),
            "per_target_cosine": {t: e["cosine_with_released_bulk"] for t, e in merged.items()
                                  if "cosine_with_released_bulk" in e}}
    out = {"R0.1_representation": {r: audit[r]["representation"] for r in audit},
           "R0.3_design": {r: audit[r]["design"] for r in audit},
           "R0.4_verification": verification,
           "R0.6_reliability": reliability,
           "seeds": {"split_half": SEED_SPLIT, "representation_sample": SEED_SAMPLE},
           "checksums": json.loads((Path("/vol/raw") / "checksums.json").read_text())}
    (results / "single_cell_audit.json").write_text(json.dumps(out, indent=2))
    vol.commit()
    return {name: len(reliability[name]) for name in reliability}


@app.local_entrypoint()
def main(stage: str = "all"):
    """stage: all | fetch | audit | targets | merge"""
    if stage in ("all", "fetch"):
        print(stage_fetch.remote())
    if stage in ("all", "audit"):
        print(stage_audit.remote())
    if stage in ("all", "targets"):
        print(stage_targets.remote())
    if stage in ("all", "merge"):
        print(stage_merge.remote())
    print("\nretrieve with:")
    print("    modal volume get di-h3-sc results/single_cell_audit.json \\")
    print("        ./results/03d_h3_reference_discordance/")
