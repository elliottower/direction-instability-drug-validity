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
    # the frozen analysis bases, which R0.4 and R0.6 are evaluated on
    .add_local_dir("registry", remote_path="/app/registry")
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
# Amendment 4: R0.4 and R0.6 are evaluated on the construction each release
# defines, not on the release's full gene axis
CONSTRUCTION = {"K562_essential": "C1-K562", "RPE1_essential": "C1-RPE1"}
BASES = "/app/registry/frozen/analysis_bases.json"
BASES_PIN = "/app/registry/frozen/analysis_bases_pin.json"
GENE_INFO = "/app/registry/frozen/landmark_gene_ids.json"
GENERATOR = "/app/experiments/03i_freeze_analysis_bases.py"


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


def _validated_bases():
    """The bases, through the same validator the R0-R7 driver uses.

    The single-cell stages produce R0.4 and R0.6 outside the driver, so a basis
    trusted here and checked there is checked in one of the two places it is used.
    """
    import json
    from pathlib import Path

    _app_importable()
    from geometry.references import validated_analysis_bases

    gene_ids = json.loads(Path(GENE_INFO).read_text())
    pin = json.loads(Path(BASES_PIN).read_text())
    bases = validated_analysis_bases(BASES, gene_ids, pin, generator=GENERATOR)
    return bases, {"analysis_bases_sha256": pin["file_sha256"],
                   "analysis_bases_pin_sha256": _sha256(BASES_PIN),
                   "generator_sha256": pin["generator_sha256"],
                   "landmark_axis_sha256": pin["landmark_axis_sha256"]}


def _shard_dir(results, release, basis):
    """Where this basis's shards live, so a superseded run is never resumed into.

    R0.4 and R0.6 ran once on each release's full gene axis, and those shards
    exist. Keying the directory by the basis digest means the amended run starts
    empty and the diagnostic run that exposed the defect stays where it is.
    """
    return results / release / f"basis_{basis['basis_sha256'][:16]}"


def _or_none(value):
    """A float, or null where it is not a number, so a shard stays valid JSON."""
    import math
    return None if value is None or not math.isfinite(value) else float(value)


def _app_importable():
    """Make `/app` importable, which is where the image carries the repository."""
    import site
    site.addsitedir("/app")


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
            if not (got == md5):
                raise AssertionError(f"{target.name}: md5 {got}, published {md5}")
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
    """R0.4 verification and R0.6 reliability, one batch of targets at a time.

    Both run on the construction's frozen analysis basis, which is what the
    registration defines C1 on: "restricted to landmarks present in the file,
    matched by symbol through the pinned gene-info file". Running them on the
    release's full gene axis made R0.4 verify one vector space while R0.6
    assessed another, and let a gene that is not a landmark into the check.
    Deviation 13 records that; Amendment 4 fixes the basis it is repaired on.
    """
    import json
    from pathlib import Path

    import anndata as ad
    import numpy as np

    _app_importable()
    from geometry.references import require_finite
    from geometry.single_cell import (basis_columns, batches, cosine, gene_names,
                                      split_half_reliability, target_seed)

    bases, provenance = _validated_bases()
    preflight = Path("/vol/results/basis_preflight.json")
    if not preflight.exists():
        raise AssertionError(
            "basis_preflight.json is absent: Amendment 4 confirms every operand is finite on "
            "the frozen basis before R0.4 and R0.6 run. Run --stage preflight first.")
    confirmed = json.loads(preflight.read_text())
    if not confirmed.get("every_operand_is_finite") or \
            confirmed.get("analysis_bases_sha256") != provenance["analysis_bases_sha256"]:
        raise AssertionError(
            "the preflight did not pass on the bases this run uses: "
            f"finite={confirmed.get('every_operand_is_finite')}, preflight bases "
            f"{confirmed.get('analysis_bases_sha256')}, this run's "
            f"{provenance['analysis_bases_sha256']}")

    results = Path("/vol/results")
    results.mkdir(parents=True, exist_ok=True)
    for release in RELEASES:
        construction = CONSTRUCTION[release]
        basis = bases[construction]
        symbols = basis["landmark_symbols"]
        shard_dir = _shard_dir(results, release, basis)
        shard_dir.mkdir(parents=True, exist_ok=True)
        single_cell = ad.read_h5ad(Path("/vol/raw") / f"{release}_single_cell.h5ad", backed="r")
        bulk = ad.read_h5ad(Path("/vol/raw") / f"{release}_bulk.h5ad")
        bulk_gene = {str(name).split("_")[1]: i for i, name in enumerate(bulk.obs_names)
                     if str(name).count("_") >= 2}

        # the same ordered coordinates in both files, so the two routes are
        # compared on one space rather than on two of the same length
        cell_columns = basis_columns(gene_names(single_cell), symbols)
        bulk_columns = basis_columns(gene_names(bulk), symbols)
        print(f"[{_ts()}] {release}: {construction} on {len(symbols)} landmarks, "
              f"basis {basis['basis_sha256'][:16]}")

        labels, _ = _targets(single_cell)
        groups, _ = _gem_groups(single_cell)
        targets = sorted({t for t in np.unique(labels) if t not in CONTROLS})
        done = {path.stem for path in shard_dir.glob("*.json")}
        print(f"[{_ts()}] {release}: {len(targets)} targets, {len(done)} batches already done")

        for name, members in batches(targets, BATCH):
            if name in done:
                continue
            batch = {}
            for target in members:
                rows = np.flatnonzero(labels == target)
                if len(rows) < 2:
                    batch[target] = {"n_cells": int(len(rows)), "status": "not computed",
                                     "reason": "fewer than two cells in the release"}
                    continue
                cells = single_cell[rows].to_memory().X
                cells = cells.toarray() if hasattr(cells, "toarray") else np.asarray(cells)
                # a non-finite value on a frozen basis coordinate stops the run; it
                # is not a second automatic exclusion
                cells = require_finite(
                    cells[:, cell_columns],
                    f"{construction}: {target}'s cells on the frozen basis")
                entry = {"n_cells": int(len(rows)),
                         "n_landmarks": len(symbols),
                         "basis_sha256": basis["basis_sha256"],
                         "analysis_bases_sha256": provenance["analysis_bases_sha256"],
                         "reliability": _or_none(split_half_reliability(
                             cells, groups[rows], target_seed(SEED_SPLIT, str(target))))}
                if target not in bulk_gene:
                    entry["cosine_status"] = "not computed"
                    entry["cosine_reason"] = ("the release holds no pseudobulk row for "
                                              "this target")
                else:
                    released = require_finite(
                        np.asarray(bulk.X[bulk_gene[target]]).ravel()[bulk_columns],
                        f"{construction}: {target}'s released pseudobulk on the frozen basis")
                    mean = cells.mean(axis=0)
                    entry["cosine_with_released_bulk"] = cosine(mean, released)
                    entry["max_abs_difference"] = float(np.abs(mean - released).max())
                batch[target] = entry
            # a shard is named only once it is complete, so an interrupted write
            # cannot be mistaken for a finished batch on resume
            partial = shard_dir / f"{name}.json.part"
            partial.write_text(json.dumps(batch, indent=2, sort_keys=True, allow_nan=False))
            partial.replace(shard_dir / f"{name}.json")
            vol.commit()                      # checkpoint inside the unit
            print(f"[{_ts()}] {release}: {name} written, {len(batch)} targets")
        single_cell.file.close()
    return "targets done"


@app.function(**COMMON)
def stage_merge():
    """One audit file, in the shape `03d_h3_reference_discordance.py` expects."""
    import json
    import sys
    from pathlib import Path

    import numpy as np

    sys.path.insert(0, "/app")
    from geometry.single_cell import merge_batches, verification_summary

    bases, provenance = _validated_bases()
    results = Path("/vol/results")
    audit = json.loads((results / "audit.json").read_text())
    reliability, verification = {}, {}
    naming = {"K562_essential": "C1-K562", "RPE1_essential": "C1-RPE1"}
    for release, name in naming.items():
        basis = bases[name]
        shard_dir = _shard_dir(results, release, basis)
        merged = merge_batches({shard.stem: json.loads(shard.read_text())
                                for shard in sorted(shard_dir.glob("*.json"))})
        stale = {target: entry.get("basis_sha256") for target, entry in merged.items()
                 if "basis_sha256" in entry and entry["basis_sha256"] != basis["basis_sha256"]}
        if stale:
            raise AssertionError(
                f"{release}: {len(stale)} shard entries were computed on another basis, first "
                f"{sorted(stale.items())[:3]}; they are not merged into this one")
        reliability[name] = {target: entry["reliability"] for target, entry in merged.items()
                             if entry.get("reliability") is not None
                             and np.isfinite(entry["reliability"])}
        verification[name] = verification_summary(name, merged)
    out = {"R0.1_representation": {r: audit[r]["representation"] for r in audit},
           "R0.3_design": {r: audit[r]["design"] for r in audit},
           "R0.4_verification": verification,
           "R0.6_reliability": reliability,
           "seeds": {"split_half": SEED_SPLIT, "representation_sample": SEED_SAMPLE},
           "checksums": json.loads((Path("/vol/raw") / "checksums.json").read_text()),
           "analysis_bases": {**provenance,
                              "per_construction": {name: bases[name]["basis_sha256"]
                                                   for name in naming.values()}},
           "preflight": json.loads((results / "basis_preflight.json").read_text())}
    (results / "single_cell_audit.json").write_text(json.dumps(out, indent=2, allow_nan=False))
    vol.commit()
    return {name: len(reliability[name]) for name in reliability}


@app.local_entrypoint()
def main(stage: str = "all"):
    """stage: all | fetch | audit | preflight | targets | merge"""
    if stage in ("all", "fetch"):
        print(stage_fetch.remote())
    if stage in ("all", "audit"):
        print(stage_audit.remote())
    if stage in ("all", "preflight"):
        print(stage_basis_preflight.remote())
    if stage in ("all", "targets"):
        print(stage_targets.remote())
    if stage in ("all", "merge"):
        print(stage_merge.remote())
    print("\nretrieve with:")
    print("    modal volume get di-h3-sc results/single_cell_audit.json \\")
    print("        ./results/03d_h3_reference_discordance/")


@app.function(**COMMON)
def probe_nonfinite_values():
    """How many non-finite values each release holds, and how many targets they reach.

    R0.4 returned a NaN median cosine for RPE1 and the code routed that into the same
    branch as a cosine below 0.99, so C1-RPE1 was labelled as describing the
    single-cell construction rather than C1. 1,943 of 2,393 per-target cosines are
    NaN while the 450 that compute have a median of 1.000000, so the label rests on
    undefined arithmetic rather than on disagreement. This finds where the
    non-finite values are.
    """
    import json
    from pathlib import Path

    import anndata as ad
    import numpy as np

    found = {}
    for release in RELEASES:
        single_cell = ad.read_h5ad(Path("/vol/raw") / f"{release}_single_cell.h5ad",
                                   backed="r")
        labels, _ = _targets(single_cell)
        n_cells, n_genes = single_cell.shape
        bad_cells, bad_values, bad_genes = 0, 0, set()
        for start in range(0, n_cells, 20000):          # streamed, never all in memory
            block = single_cell[start:start + 20000].to_memory().X
            block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
            mask = ~np.isfinite(block)
            if mask.any():
                bad_values += int(mask.sum())
                rows = np.flatnonzero(mask.any(axis=1))
                bad_cells += len(rows)
                bad_genes.update(int(g) for g in np.flatnonzero(mask.any(axis=0)))
        bad_row_labels = set()
        if bad_values:
            for start in range(0, n_cells, 20000):
                block = single_cell[start:start + 20000].to_memory().X
                block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
                rows = np.flatnonzero((~np.isfinite(block)).any(axis=1))
                bad_row_labels.update(str(labels[start + r]) for r in rows)
        found[release] = {
            "shape": [int(n_cells), int(n_genes)],
            "non_finite_values": bad_values,
            "cells_holding_one": bad_cells,
            "genes_holding_one": len(bad_genes),
            "targets_reached": len(bad_row_labels),
            "example_targets": sorted(bad_row_labels)[:8],
        }
        print(json.dumps({release: found[release]}, indent=2))
    out = Path("/vol/results/nonfinite_probe.json")
    out.write_text(json.dumps(found, indent=2))
    vol.commit()
    return found


#: The raw counterparts of the normalized releases, by Figshare file id and published
#: md5. Figshare publishes no md5 for rpe1_raw_singlecell_01.h5ad, so that one is
#: recorded by the digest computed on download and cannot be checked against the
#: publisher. That limitation is reported rather than papered over.
RAW_FILES = {
    "rpe1_raw_bulk_01.h5ad": ("35775581", "74765fa87635467a869ea972356ae0e7"),
    "rpe1_raw_singlecell_01.h5ad": ("35775606", None),
    "K562_essential_raw_bulk_01.h5ad": ("35773070", "8321d5d3ffc99db2a5c71edca4189735"),
}


@app.function(**COMMON)
def stage_scan_raw_rpe1():
    """Are ATF3 and CCL2 finite before normalization?

    Perplexity's round-19 ruling: the normalized RPE1 files cannot be repaired by
    clipping, zero-substitution, nanmean or dropping the offending coordinate, and the
    raw files are the next diagnostic layer. If the raw values are finite then the
    infinities arise in the published normalization and a deterministic reconstruction
    may be justifiable; if they are not, the C1-RPE1 arm is held as not computable.

    Downloads to the same volume the normalized files sit on, verifies each md5 that
    Figshare publishes, streams the finite-value scan, and reports ATF3 and CCL2
    specifically alongside the whole-matrix counts.
    """
    import json
    import urllib.request
    from pathlib import Path

    import anndata as ad
    import numpy as np

    raw = Path("/vol/raw")
    raw.mkdir(parents=True, exist_ok=True)
    out = Path("/vol/results")
    out.mkdir(parents=True, exist_ok=True)
    wanted = {"ENSG00000162772": "ATF3", "ENSG00000108691": "CCL2"}
    found = {}

    for name, (file_id, published) in RAW_FILES.items():
        target = raw / name
        if not target.exists():
            partial = target.with_suffix(".part")
            print(f"[{_ts()}] downloading {name}")
            urllib.request.urlretrieve(f"{FIGSHARE}/{file_id}", partial)
            partial.replace(target)
            vol.commit()
        got = _md5(target)
        if published is None:
            checksum = {"md5": got, "published": None,
                        "note": "Figshare publishes no md5 for this file"}
        else:
            if got != published:
                raise AssertionError(f"{name}: md5 {got}, published {published}")
            checksum = {"md5": got, "published": published, "matches": True}

        adata = ad.read_h5ad(target, backed="r")
        n_obs, n_var = adata.shape
        names = [str(v) for v in adata.var_names]
        columns = {ens: names.index(ens) for ens in wanted if ens in names}
        counts = {"nan": 0, "posinf": 0, "neginf": 0}
        obs_hit, gene_hit = 0, set()
        per_gene = {ens: {"nan": 0, "posinf": 0, "neginf": 0, "min": None, "max": None}
                    for ens in columns}
        for start in range(0, n_obs, 20000):
            block = adata[start:start + 20000].to_memory().X
            block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
            block = np.asarray(block, dtype=np.float64)
            bad = ~np.isfinite(block)
            if bad.any():
                counts["nan"] += int(np.isnan(block).sum())
                counts["posinf"] += int(np.isposinf(block).sum())
                counts["neginf"] += int(np.isneginf(block).sum())
                obs_hit += int(bad.any(axis=1).sum())
                gene_hit.update(int(g) for g in np.flatnonzero(bad.any(axis=0)))
            for ens, column in columns.items():
                values = block[:, column]
                finite = values[np.isfinite(values)]
                per_gene[ens]["nan"] += int(np.isnan(values).sum())
                per_gene[ens]["posinf"] += int(np.isposinf(values).sum())
                per_gene[ens]["neginf"] += int(np.isneginf(values).sum())
                if finite.size:
                    lo, hi = float(finite.min()), float(finite.max())
                    cur = per_gene[ens]
                    cur["min"] = lo if cur["min"] is None else min(cur["min"], lo)
                    cur["max"] = hi if cur["max"] is None else max(cur["max"], hi)
            json.dump({"file": name, "rows_scanned": min(start + 20000, n_obs),
                       "of": n_obs, "non_finite_so_far": sum(counts.values())},
                      open(out / "raw_scan_progress.json", "w"))   # RULE ONE
            vol.commit()

        found[name] = {
            "checksum": checksum,
            "shape": [int(n_obs), int(n_var)],
            "non_finite_total": sum(counts.values()),
            "non_finite_by_kind": counts,
            "observations_affected": obs_hit,
            "gene_columns_affected": len(gene_hit),
            "the_two_genes": {wanted[ens]: {"ensembl": ens, "present": ens in columns,
                                            **per_gene.get(ens, {})} for ens in wanted},
        }
        print(json.dumps({name: found[name]}, indent=2))
        (out / "raw_rpe1_scan.json").write_text(json.dumps(found, indent=2))
        vol.commit()

    return found


@app.function(**COMMON)
def stage_zero_variance_mechanism_r2():
    """Whether ATF3 and CCL2 have zero within-gem-group variance in the raw RPE1 data.

    The raw files hold no non-finite value and both genes reach a minimum of 0.0, so
    the infinities are introduced by the published gem-group Z-normalization. A gene
    whose counts are constant within a gem group has zero standard deviation there, and
    dividing by it yields +/-inf for every cell in that group. This checks that
    directly, per gem group, and reports which groups are degenerate for each gene.
    """
    import json
    from pathlib import Path

    import anndata as ad
    import numpy as np

    raw = Path("/vol/raw")
    out = Path("/vol/results")
    out.mkdir(parents=True, exist_ok=True)
    wanted = {"ENSG00000162772": "ATF3", "ENSG00000108691": "CCL2"}

    counts = ad.read_h5ad(raw / "rpe1_raw_singlecell_01.h5ad", backed="r")
    normalized = ad.read_h5ad(raw / "RPE1_essential_single_cell.h5ad", backed="r")
    names = [str(v) for v in counts.var_names]
    columns = {ens: names.index(ens) for ens in wanted if ens in names}
    groups, _ = _gem_groups(counts)
    groups = np.asarray(groups)

    found = {}
    for ens, column in columns.items():
        raw_values = np.empty(counts.n_obs, dtype=np.float64)
        norm_values = np.empty(counts.n_obs, dtype=np.float64)
        for start in range(0, counts.n_obs, 20000):
            stop = min(start + 20000, counts.n_obs)
            for source, into in ((counts, raw_values), (normalized, norm_values)):
                block = source[start:stop].to_memory().X
                block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
                into[start:stop] = np.asarray(block[:, column], dtype=np.float64)
            json.dump({"gene": wanted[ens], "rows": stop, "of": int(counts.n_obs)},
                      open(out / "zero_variance_progress.json", "w"))    # RULE ONE
            vol.commit()

        per_group = {}
        for group in np.unique(groups):
            rows = groups == group
            block = raw_values[rows]
            sd = float(block.std())
            nonfinite = int((~np.isfinite(norm_values[rows])).sum())
            per_group[str(group)] = {"n_cells": int(rows.sum()), "raw_sd": sd,
                                     "raw_all_equal": bool(block.min() == block.max()),
                                     "raw_value_if_constant": (float(block[0]) if
                                                               block.min() == block.max()
                                                               else None),
                                     "non_finite_after_normalization": nonfinite}
        degenerate = sorted(g for g, v in per_group.items() if v["raw_all_equal"])
        poisoned = sorted(g for g, v in per_group.items()
                          if v["non_finite_after_normalization"])
        found[wanted[ens]] = {
            "ensembl": ens,
            "n_gem_groups": int(len(per_group)),
            "groups_with_constant_raw_values": degenerate,
            "groups_with_non_finite_normalized_values": poisoned,
            "the_two_sets_agree": degenerate == poisoned,
            "non_finite_cells_total": int((~np.isfinite(norm_values)).sum()),
            "per_group": per_group,
        }
        print(json.dumps({wanted[ens]: {k: v for k, v in found[wanted[ens]].items()
                                        if k != "per_group"}}, indent=2))
        (out / "zero_variance_mechanism.json").write_text(json.dumps(found, indent=2))
        vol.commit()
    return {g: {k: v for k, v in d.items() if k != "per_group"} for g, d in found.items()}


@app.function(**COMMON)
def stage_control_variance_in_the_two_groups():
    """Why exactly one gem group per gene goes non-finite.

    The first guess was wrong: no gem group has constant raw counts for either gene,
    yet ATF3 is non-finite only in group 46 and CCL2 only in group 9, together the
    7,913 cells. Z-normalization is against the control cells of the group, so the
    denominator is the control standard deviation, not the all-cell one. This reports
    the control-cell counts for each gene in its own affected group beside an
    unaffected group, so the divisor is visible rather than inferred.
    """
    import json
    from pathlib import Path

    import anndata as ad
    import numpy as np

    raw = Path("/vol/raw")
    out = Path("/vol/results")
    counts = ad.read_h5ad(raw / "rpe1_raw_singlecell_01.h5ad", backed="r")
    names = [str(v) for v in counts.var_names]
    labels, _ = _targets(counts)
    groups, _ = _gem_groups(counts)
    labels, groups = np.asarray(labels), np.asarray(groups)
    is_control = np.isin(labels, list(CONTROLS))

    cases = {"ATF3": ("ENSG00000162772", "46"), "CCL2": ("ENSG00000108691", "9")}
    found = {"controls_label_values": sorted(set(map(str, labels[is_control])))[:5],
             "n_control_cells": int(is_control.sum())}
    for gene, (ens, affected) in cases.items():
        column = names.index(ens)
        values = np.empty(counts.n_obs, dtype=np.float64)
        for start in range(0, counts.n_obs, 20000):
            stop = min(start + 20000, counts.n_obs)
            block = counts[start:stop].to_memory().X
            block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
            values[start:stop] = np.asarray(block[:, column], dtype=np.float64)
            json.dump({"gene": gene, "rows": stop}, open(out / "cv_progress.json", "w"))
            vol.commit()
        per_group = {}
        for group in np.unique(groups):
            rows = (groups == group) & is_control
            block = values[rows]
            per_group[str(group)] = {
                "n_control_cells": int(rows.sum()),
                "control_sd": float(block.std()) if block.size else None,
                "control_mean": float(block.mean()) if block.size else None,
                "control_all_equal": bool(block.size and block.min() == block.max()),
                "control_unique_values": sorted(set(block.tolist()))[:6] if block.size else [],
            }
        degenerate = sorted(g for g, v in per_group.items() if v["control_all_equal"])
        found[gene] = {
            "ensembl": ens,
            "group_that_went_non_finite": affected,
            "groups_whose_controls_are_constant": degenerate,
            "the_affected_group_is_explained": affected in degenerate,
            "affected_group_controls": per_group.get(affected),
            "per_group": per_group,
        }
        print(json.dumps({gene: {k: v for k, v in found[gene].items() if k != "per_group"}},
                         indent=2))
        (out / "control_variance_mechanism.json").write_text(json.dumps(found, indent=2))
        vol.commit()
    return {g: {k: v for k, v in d.items() if k != "per_group"}
            for g, d in found.items() if isinstance(d, dict) and "ensembl" in d}


@app.function(**COMMON)
def stage_can_the_normalization_be_reproduced():
    """Does per-gem-group control Z-normalization reproduce the published RPE1 file?

    Perplexity's round-19 hierarchy permits reconstructing the normalized source from
    raw only if the published normalization can be reproduced and the reconstruction is
    verified against all unaffected values, not merely the two defective genes. This is
    that test, and it is the gate on whether the C1-RPE1 arm survives.

    The candidate transformation, from the mechanism already established: for each gem
    group, subtract the mean of that group's non-targeting control cells and divide by
    their standard deviation, per gene. Agreement is measured on the finite positions
    only, because the published file has no recoverable value where it holds +inf.

    A sample of genes is used rather than all 8,749, because the question is whether the
    transformation is the right one, and a wrong transformation disagrees everywhere.
    The two defective genes are always included.
    """
    import json
    from pathlib import Path

    import anndata as ad
    import numpy as np

    raw_dir = Path("/vol/raw")
    out = Path("/vol/results")
    counts = ad.read_h5ad(raw_dir / "rpe1_raw_singlecell_01.h5ad", backed="r")
    published = ad.read_h5ad(raw_dir / "RPE1_essential_single_cell.h5ad", backed="r")
    names = [str(v) for v in counts.var_names]
    if [str(v) for v in published.var_names] != names:
        raise AssertionError("the raw and published files do not share a gene axis")

    labels, _ = _targets(counts)
    groups, _ = _gem_groups(counts)
    labels, groups = np.asarray(labels), np.asarray(groups)
    is_control = np.isin(labels, list(CONTROLS))

    rng = np.random.default_rng(20261006)
    sample = sorted(set(rng.choice(len(names), size=40, replace=False).tolist())
                    | {names.index("ENSG00000162772"), names.index("ENSG00000108691")})

    worst, compared, nonfinite_published, per_gene = 0.0, 0, 0, {}
    for column in sample:
        raw_values = np.empty(counts.n_obs, dtype=np.float64)
        pub_values = np.empty(counts.n_obs, dtype=np.float64)
        for start in range(0, counts.n_obs, 20000):
            stop = min(start + 20000, counts.n_obs)
            for source, into in ((counts, raw_values), (published, pub_values)):
                block = source[start:stop].to_memory().X
                block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
                into[start:stop] = np.asarray(block[:, column], dtype=np.float64)

        rebuilt = np.full(counts.n_obs, np.nan)
        for group in np.unique(groups):
            rows = groups == group
            control = raw_values[rows & is_control]
            if not control.size:
                continue
            sd = control.std()
            rebuilt[rows] = ((raw_values[rows] - control.mean()) / sd if sd > 0
                             else np.inf * np.sign(raw_values[rows] - control.mean()))

        usable = np.isfinite(pub_values) & np.isfinite(rebuilt)
        if usable.any():
            error = float(np.abs(rebuilt[usable] - pub_values[usable]).max())
            worst = max(worst, error)
            compared += int(usable.sum())
        else:
            error = None
        nonfinite_published += int((~np.isfinite(pub_values)).sum())
        per_gene[names[column]] = {
            "n_finite_in_both": int(usable.sum()),
            "max_abs_difference": error,
            "n_non_finite_published": int((~np.isfinite(pub_values)).sum()),
            "n_non_finite_rebuilt": int((~np.isfinite(rebuilt)).sum()),
        }
        json.dump({"genes_done": len(per_gene), "of": len(sample),
                   "worst_so_far": worst}, open(out / "repro_progress.json", "w"))
        vol.commit()

    result = {
        "candidate": ("per gem group, subtract the mean of that group's non-targeting "
                      "control cells and divide by their standard deviation, per gene"),
        "n_genes_sampled": len(sample),
        "n_values_compared": compared,
        "max_abs_difference_on_finite_positions": worst,
        "n_non_finite_in_published_sample": nonfinite_published,
        "reproduces": worst < 1e-4,
        "per_gene": per_gene,
    }
    result["reading"] = (
        "the published normalization is reproduced from the raw counts by this "
        f"transformation to {worst:.3e} on {compared} finite positions, so a "
        "deterministic reconstruction is justifiable"
        if result["reproduces"] else
        f"this transformation does not reproduce the published file: {worst:.3e}. The "
        "reconstruction route is not justified on it and the C1-RPE1 arm stays not "
        "computable unless the exact published procedure is obtained.")
    (out / "normalization_reproduction.json").write_text(json.dumps(result, indent=2))
    vol.commit()
    print(json.dumps({k: v for k, v in result.items() if k != "per_gene"}, indent=2))
    return {k: v for k, v in result.items() if k != "per_gene"}


@app.function(**COMMON)
def probe_raw_obs_fields():
    """What per-cell fields the raw RPE1 file carries, for the UMI-total step."""
    import json
    from pathlib import Path

    import anndata as ad
    import numpy as np

    counts = ad.read_h5ad(Path("/vol/raw") / "rpe1_raw_singlecell_01.h5ad", backed="r")
    found = {"shape": list(map(int, counts.shape)),
             "obs_columns": {c: str(counts.obs[c].dtype) for c in counts.obs.columns},
             "var_columns": {c: str(counts.var[c].dtype) for c in counts.var.columns},
             "layers": list(counts.layers.keys()) if counts.layers else []}
    for column in counts.obs.columns:
        values = counts.obs[column]
        if values.dtype.kind in "if":
            found[f"obs::{column}"] = {"min": float(values.min()), "max": float(values.max()),
                                       "median": float(values.median())}
    # the total over the genes this file holds, for the first block of cells
    block = counts[:2000].to_memory().X
    block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
    totals = np.asarray(block, dtype=np.float64).sum(axis=1)
    found["summed_over_this_file_genes"] = {"min": float(totals.min()),
                                            "max": float(totals.max()),
                                            "median": float(np.median(totals))}
    print(json.dumps(found, indent=2))
    return found


@app.function(**COMMON)
def stage_reproduce_with_the_published_scale_factor():
    """The paper's two documented steps, using the scale factor the raw file carries.

    STAR Methods: "(i) UMI count normalization: We scale expression within all cells so
    that their total UMI counts equal the median UMI count of core control cells within
    the experiment. (ii) Relative z-normalization: Within each gemgroup, for each gene,
    we compute the mean and standard deviation of expression within control cells and
    use these to z-normalize expression."

    The earlier attempt omitted step (i) and disagreed by 1.8e+03. The raw file carries
    `core_scale_factor` per cell, which is step (i)'s factor as the authors computed it,
    so neither the target nor the core-control set has to be inferred.
    """
    import json
    from pathlib import Path

    import anndata as ad
    import numpy as np

    raw_dir = Path("/vol/raw")
    out = Path("/vol/results")
    counts = ad.read_h5ad(raw_dir / "rpe1_raw_singlecell_01.h5ad", backed="r")
    published = ad.read_h5ad(raw_dir / "RPE1_essential_single_cell.h5ad", backed="r")
    names = [str(v) for v in counts.var_names]
    if [str(v) for v in published.var_names] != names:
        raise AssertionError("the raw and published files do not share a gene axis")

    scale = np.asarray(counts.obs["core_scale_factor"], dtype=np.float64)
    labels, _ = _targets(counts)
    groups, _ = _gem_groups(counts)
    labels, groups = np.asarray(labels), np.asarray(groups)
    is_control = np.isin(labels, list(CONTROLS))

    rng = np.random.default_rng(20261006)
    sample = sorted(set(rng.choice(len(names), size=40, replace=False).tolist())
                    | {names.index("ENSG00000162772"), names.index("ENSG00000108691")})

    worst, compared, per_gene = 0.0, 0, {}
    for column in sample:
        raw_values = np.empty(counts.n_obs, dtype=np.float64)
        pub_values = np.empty(counts.n_obs, dtype=np.float64)
        for start in range(0, counts.n_obs, 20000):
            stop = min(start + 20000, counts.n_obs)
            for source, into in ((counts, raw_values), (published, pub_values)):
                block = source[start:stop].to_memory().X
                block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
                into[start:stop] = np.asarray(block[:, column], dtype=np.float64)

        scaled = raw_values * scale                       # step (i)
        rebuilt = np.full(counts.n_obs, np.nan)
        for group in np.unique(groups):                   # step (ii)
            rows = groups == group
            control = scaled[rows & is_control]
            if not control.size:
                continue
            sd = control.std()
            centered = scaled[rows] - control.mean()
            rebuilt[rows] = centered / sd if sd > 0 else np.inf * np.sign(centered)

        usable = np.isfinite(pub_values) & np.isfinite(rebuilt)
        error = float(np.abs(rebuilt[usable] - pub_values[usable]).max()) if usable.any() else None
        if error is not None:
            worst = max(worst, error)
            compared += int(usable.sum())
        per_gene[names[column]] = {
            "n_finite_in_both": int(usable.sum()), "max_abs_difference": error,
            "n_non_finite_published": int((~np.isfinite(pub_values)).sum()),
            "n_non_finite_rebuilt": int((~np.isfinite(rebuilt)).sum())}
        json.dump({"genes_done": len(per_gene), "of": len(sample), "worst": worst},
                  open(out / "repro2_progress.json", "w"))
        vol.commit()

    result = {
        "procedure": ("raw counts times the file's own core_scale_factor, then per gem "
                      "group a z against the non-targeting control cells of that group"),
        "documented_in": "Replogle et al. 2022 Cell, STAR Methods, internal normalization",
        "n_genes_sampled": len(sample), "n_values_compared": compared,
        "max_abs_difference_on_finite_positions": worst,
        "reproduces": worst < 1e-3,
        "non_finite_positions_agree": all(
            v["n_non_finite_published"] == v["n_non_finite_rebuilt"] for v in per_gene.values()),
        "per_gene": per_gene,
    }
    result["reading"] = (
        f"the published normalization reproduces from the raw counts to {worst:.3e} over "
        f"{compared} finite positions, by the authors' two documented steps and their own "
        "per-cell scale factor, so a deterministic reconstruction is justifiable"
        if result["reproduces"] else
        f"this does not reproduce the published file either: {worst:.3e}")
    (out / "normalization_reproduction_v2.json").write_text(json.dumps(result, indent=2))
    vol.commit()
    print(json.dumps({k: v for k, v in result.items() if k != "per_gene"}, indent=2))
    return {k: v for k, v in result.items() if k != "per_gene"}


@app.function(**COMMON)
def stage_basis_preflight():
    """Is every value R0.4 and R0.6 will read finite on the frozen basis?

    This computes no registered statistic. It confirms that each release lands the
    frozen basis exactly once per symbol, and that the values on those coordinates
    are finite in both the single-cell and the pseudobulk file, so a defect
    surfaces before 2,393 targets of work rather than part way through it.

    Amendment 4 excludes a coordinate only through the frozen basis, so anything
    non-finite found here is a refusal and not a second exclusion.
    """
    import json
    from pathlib import Path

    import anndata as ad
    import numpy as np

    _app_importable()
    from geometry.single_cell import basis_columns, gene_names

    bases, provenance = _validated_bases()
    report, clean = {}, True
    for release in RELEASES:
        construction = CONSTRUCTION[release]
        basis = bases[construction]
        symbols = basis["landmark_symbols"]
        entry = {"construction": construction, "n_landmarks": len(symbols),
                 "basis_sha256": basis["basis_sha256"],
                 "n_landmarks_matched_before_the_finite_filter": basis["n_landmarks_matched"],
                 "excluded_as_undefined": [gene["symbol"]
                                           for gene in basis["excluded_as_undefined"]]}
        for kind in ("single_cell", "bulk"):
            path = Path("/vol/raw") / f"{release}_{kind}.h5ad"
            adata = ad.read_h5ad(path, backed="r") if kind == "single_cell" \
                else ad.read_h5ad(path)
            columns = basis_columns(gene_names(adata), symbols)
            n_rows = adata.shape[0]
            non_finite, rows_affected, scanned = 0, 0, 0
            for start in range(0, n_rows, 20_000):
                block = adata[start:min(start + 20_000, n_rows)]
                values = block.to_memory().X if kind == "single_cell" else block.X
                values = values.toarray() if hasattr(values, "toarray") else np.asarray(values)
                offending = ~np.isfinite(values[:, columns])
                non_finite += int(offending.sum())
                rows_affected += int(offending.any(axis=1).sum())
                scanned += int(offending.size)
                print(f"[{_ts()}] {release} {kind}: {start + len(values)}/{n_rows} rows, "
                      f"{non_finite} non-finite so far", flush=True)
            if kind == "single_cell":
                adata.file.close()
            entry[kind] = {"file": path.name, "shape": [int(n_rows), len(symbols)],
                           "n_values_scanned": scanned,
                           "n_non_finite_on_the_basis": non_finite,
                           "n_rows_affected": rows_affected,
                           "every_basis_value_is_finite": non_finite == 0}
            clean = clean and non_finite == 0
        report[release] = entry

    report["computed_no_registered_statistic"] = True
    report["every_operand_is_finite"] = clean
    report.update(provenance)
    report["reading"] = (
        "every value R0.4 and R0.6 read on each construction's frozen basis is finite"
        if clean else
        "a non-finite value remains on a frozen basis coordinate; R0.4 and R0.6 refuse "
        "rather than exclude it, and the basis is not widened or narrowed to suit")
    out = Path("/vol/results")
    out.mkdir(parents=True, exist_ok=True)
    (out / "basis_preflight.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    vol.commit()
    # the report is written and committed first, so a refusal leaves a record of
    # what it refused rather than only a traceback
    if not clean:
        raise AssertionError(
            "a frozen basis coordinate holds a non-finite value in a source R0.4 or R0.6 "
            "reads. Amendment 4 refuses rather than excluding it a second time; see "
            "basis_preflight.json for which release, which file and how many values.")
    print(json.dumps(report, indent=2))
    return f"every operand finite on {len(bases)} frozen bases"
