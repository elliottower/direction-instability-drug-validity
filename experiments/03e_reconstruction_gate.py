"""The reconstruction gate of the corrected-baseline registration (frozen 7f57136).

An independent rebuild from the GEO GCTX is compared with the extraction the
corrected artifact was computed from. The comparison joins on identifiers, never
on row or column order: the defect of Deviation 9 was an identifier-versus-position
error, and a gate that compares by position could not see it.

Both halves of the registered gate are here. Compound matrices are joined on
(drug, cell line) and genes on Entrez id. shRNA signatures are joined on signature
id, and each target consensus is then rebuilt independently on both sides and
compared, because the defect this gate exists to catch was in the target
directions rather than in the compound data. Identifier sequences are hashed
separately from the numeric arrays. The tolerance
is the float32 one, because the GCTX holds float32 values whose magnitude reaches
about ten.

    PYTHONPATH=. uv run --no-project --with numpy --with pandas python \\
        experiments/03e_reconstruction_gate.py \\
        --rebuilt results/03c_h3_sensitivity/gctx_rebuild \\
        --extraction ../drug-perturbation-geometry/data
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from cmapPy.pandasGEXpress import parse

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
OUT = REPO / "results" / "03c_h3_sensitivity"
REFERENCE_ARTIFACT = (REPO / "results" / "03_phenotype_projection"
                      / "phenotype_projection_results.json")
RTOL, ATOL = 1e-5, 1e-5
N_LANDMARK = 978
MIN_SIGNATURES = 3        # distinct signature ids, not reagents


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def align(rebuilt_rows, rebuilt_ids, extraction_rows, extraction_ids):
    """Reorder both matrices into one identifier order, or say why it cannot be done."""
    rebuilt_ids, extraction_ids = list(map(str, rebuilt_ids)), list(map(str, extraction_ids))
    assert len(set(rebuilt_ids)) == len(rebuilt_ids), "the rebuild repeats an identifier"
    assert len(set(extraction_ids)) == len(extraction_ids), "the extraction repeats an identifier"
    assert set(rebuilt_ids) == set(extraction_ids), (
        f"identifier sets differ: {len(set(rebuilt_ids) - set(extraction_ids))} only in the "
        f"rebuild, {len(set(extraction_ids) - set(rebuilt_ids))} only in the extraction")
    order = {identifier: i for i, identifier in enumerate(extraction_ids)}
    permutation = np.array([order[identifier] for identifier in rebuilt_ids])
    return np.asarray(rebuilt_rows), np.asarray(extraction_rows)[permutation], rebuilt_ids


def compare(rebuilt, extraction, cohort=None, rtol=RTOL, atol=ATOL):
    """Per-drug comparison after aligning rows by cell line and columns by gene id.

    `cohort` is the set of drugs that must be present on both sides. Comparing the
    observed intersection instead would let a rebuild holding one drug pass, so
    the cohort is required rather than inferred.
    """
    rebuilt_matrices, rebuilt_cells, rebuilt_genes = rebuilt
    extraction_matrices, extraction_cells, extraction_genes = extraction
    cohort = set(rebuilt_matrices) if cohort is None else set(cohort)
    missing_rebuilt = cohort - set(rebuilt_matrices)
    missing_extraction = cohort - set(extraction_matrices)
    assert not missing_rebuilt, (
        f"{len(missing_rebuilt)} cohort drugs absent from the rebuild: "
        f"{sorted(missing_rebuilt)[:5]}")
    assert not missing_extraction, (
        f"{len(missing_extraction)} cohort drugs absent from the extraction: "
        f"{sorted(missing_extraction)[:5]}")
    rebuilt_matrices = {d: rebuilt_matrices[d] for d in cohort}
    extraction_matrices = {d: extraction_matrices[d] for d in cohort}

    gene_order = {gene: i for i, gene in enumerate(map(str, extraction_genes))}
    assert set(map(str, rebuilt_genes)) == set(gene_order), "the gene axes differ"
    gene_permutation = np.array([gene_order[gene] for gene in map(str, rebuilt_genes)])

    worst, per_drug = 0.0, {}
    for drug in sorted(cohort):
        left, right, cells = align(rebuilt_matrices[drug], rebuilt_cells[drug],
                                   extraction_matrices[drug], extraction_cells[drug])
        right = right[:, gene_permutation]     # the extraction's genes, in the rebuild's order
        difference = np.abs(left - right)
        allowed = atol + rtol * np.abs(right)
        per_drug[drug] = {"n_cell_lines": len(cells),
                          "max_abs_difference": float(difference.max()),
                          "within_tolerance": bool((difference <= allowed).all())}
        worst = max(worst, float(difference.max()))
    return {"n_drugs": len(per_drug), "max_abs_difference": worst,
            "all_within_tolerance": all(entry["within_tolerance"] for entry in per_drug.values()),
            "rtol": rtol, "atol": atol,
            "n_outside_cohort": {"rebuild": len(set(rebuilt[0]) - cohort),
                                 "extraction": len(set(extraction[0]) - cohort)},
            "per_drug": per_drug}


def load_rebuild(directory: Path):
    """Per-drug matrices, their cell lines, and the gene axis, from the Modal rebuild."""
    matrices, cells = {}, {}
    genes = None
    for shard in sorted(Path(directory).glob("shard_*.npz")):
        with np.load(shard, allow_pickle=True) as data:
            assert "__cells__" in data.files, (
                f"{shard.name} carries no cell-line identifiers; the join would fall back to "
                "row order, which is what this gate exists to avoid")
            shard_cells = json.loads(str(data["__cells__"]))
            for key in data.files:
                if key in ("fingerprint", "__cells__"):
                    continue
                matrices[key] = data[key]
            cells.update({drug: [str(c) for c in ids] for drug, ids in shard_cells.items()})
    gene_file = Path(directory) / "landmark_gene_ids.json"
    if gene_file.exists():
        genes = [str(g) for g in json.loads(gene_file.read_text())]
    assert genes is not None, f"{gene_file} is missing; the rebuild must record its gene axis"
    missing = set(matrices) - set(cells)
    assert not missing, f"{len(missing)} drugs arrived without cell-line identifiers"
    return matrices, cells, genes


def load_extraction(data_dir: Path, cohort=None):
    """The same objects from the extraction the corrected artifact was computed from.

    `cohort` restricts the build to the drugs the gate compares. Without it this
    builds a matrix for every compound in the metadata, which is tens of thousands
    of drugs the gate never looks at.
    """
    compounds = np.load(data_dir / "lincs_subset.npz", allow_pickle=True)
    signatures = compounds["signatures"]      # decompressed once
    position = {str(sig_id): i for i, sig_id in enumerate(compounds["sig_ids"])}
    siginfo = pd.read_csv(data_dir / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(position) & siginfo.pert_iname.notna()].copy()
    if cohort is not None:
        siginfo = siginfo[siginfo.pert_iname.isin(set(cohort))]
    siginfo["_row"] = siginfo.sig_id.astype(str).map(position)

    matrices, cells = {}, {}
    for drug, group in siginfo.groupby("pert_iname"):
        per_cell, identifiers = [], []
        for cell_id, rows in group.groupby("cell_id"):
            per_cell.append(signatures[rows._row.values].mean(axis=0))
            identifiers.append(str(cell_id))
        matrices[drug] = np.vstack(per_cell)
        cells[drug] = identifiers
    return matrices, cells, [str(g) for g in compounds["gene_ids"]]


def compare_shrna(rebuilt, extraction, stored=None,
                  min_signatures=MIN_SIGNATURES, rtol=RTOL, atol=ATOL):
    """The half of the gate that guards the object Deviation 9 corrupted.

    A perfect compound reconstruction cannot detect a recurrence of that defect,
    because the defect was in building the target directions. This compares the
    shRNA signatures themselves, joined by signature id, then rebuilds each target
    consensus independently on both sides and compares the consensus and the unit
    direction it becomes.

    Each side is (signatures by signature id, target -> its signature ids, genes).

    `stored` is (directions by target, gene axis) as the rebuild wrote them.
    Comparing those against the directions recomputed here from verified signatures
    is the check that bites: Deviation 9 was a rebuild whose signatures were all
    correct and whose consensus construction was not. The stored axis is aligned on
    its own identifiers, never on the extraction's permutation, because assuming a
    coordinate system is the error this gate exists to catch.

    A value that moved is counted in the result rather than raised, so the record
    survives the failure. A structural violation — a missing identifier, a repeated
    gene, a direction whose shape does not match its axis — raises, because there is
    no number to report about it.
    """
    rebuilt_signatures, rebuilt_membership, rebuilt_genes = rebuilt
    extraction_signatures, extraction_membership, extraction_genes = extraction

    for name, axis in (("rebuild", rebuilt_genes), ("extraction", extraction_genes)):
        assert len(set(map(str, axis))) == len(axis), (
            f"the {name} gene axis repeats an identifier, which makes its coordinate map "
            "ambiguous")
    assert set(rebuilt_signatures) == set(extraction_signatures), (
        f"signature-id sets differ: {len(set(rebuilt_signatures) - set(extraction_signatures))} "
        f"only in the rebuild, "
        f"{len(set(extraction_signatures) - set(rebuilt_signatures))} only in the extraction")
    gene_order = {gene: i for i, gene in enumerate(map(str, extraction_genes))}
    assert set(map(str, rebuilt_genes)) == set(gene_order), "the shRNA gene axes differ"
    permutation = np.array([gene_order[gene] for gene in map(str, rebuilt_genes)])

    # a value mismatch is recorded rather than raised: the gate's own result file
    # is written by its caller, and an exception here left no record of the one
    # failure the gate exists to document
    worst_signature, signature_failures = 0.0, []
    for sig_id, vector in rebuilt_signatures.items():
        other = np.asarray(extraction_signatures[sig_id])[permutation]
        difference = np.abs(np.asarray(vector) - other)
        worst_signature = max(worst_signature, float(difference.max()))
        if not (difference <= atol + rtol * np.abs(other)).all():
            signature_failures.append({"sig_id": sig_id,
                                       "max_abs_difference": float(difference.max())})

    assert set(rebuilt_membership) == set(extraction_membership), (
        f"target sets differ: {sorted(set(rebuilt_membership) ^ set(extraction_membership))[:5]}")

    stored_permutation = None
    if stored is not None:
        stored_directions, stored_genes = stored
        stored_genes = [str(g) for g in stored_genes]
        assert len(set(stored_genes)) == len(stored_genes), (
            "the stored consensus repeats a gene identifier")
        assert set(stored_genes) == set(map(str, rebuilt_genes)), (
            "the stored consensus does not share the rebuild's gene axis")
        assert set(stored_directions) == set(rebuilt_membership), (
            "the stored consensus covers a different set of targets than the rebuild: "
            f"{sorted(set(stored_directions) ^ set(rebuilt_membership))[:5]}")
        stored_index = {gene: i for i, gene in enumerate(stored_genes)}
        stored_permutation = np.array([stored_index[gene] for gene in map(str, rebuilt_genes)])

    worst = {"consensus": 0.0, "direction": 0.0, "stored": 0.0}
    worst_target = {"consensus": None, "direction": None, "stored": None}
    per_target, target_failures = {}, []
    for target in sorted(rebuilt_membership):
        left_ids, right_ids = sorted(rebuilt_membership[target]), sorted(extraction_membership[target])
        assert left_ids == right_ids, (
            f"{target}: the two sides group different signatures into it")
        assert len(left_ids) >= min_signatures, (
            f"{target}: {len(left_ids)} signatures, below the eligibility floor")
        left = np.vstack([rebuilt_signatures[i] for i in left_ids]).mean(axis=0)
        right = np.vstack([np.asarray(extraction_signatures[i])[permutation]
                           for i in right_ids]).mean(axis=0)
        # a zero or nonfinite consensus must say so, not surface later as a
        # division warning and an obscure tolerance failure
        for side, vector in (("rebuild", left), ("extraction", right)):
            assert np.isfinite(vector).all(), f"{target}: nonfinite consensus on the {side}"
            assert np.linalg.norm(vector) > 0, f"{target}: zero-norm consensus on the {side}"
        # the registered tolerance is the elementwise float32 rule. Reducing it to
        # one scalar threshold compares unnormalized consensus entries, whose
        # magnitude reaches about ten, against a unit-scale bound: too lax below
        # one and too strict above it.
        consensus_difference = np.abs(left - right)
        consensus_error = float(consensus_difference.max())
        consensus_ok = bool((consensus_difference <= atol + rtol * np.abs(right)).all())
        # a normalized direction can agree while the consensus behind it does not,
        # so both are compared
        left_unit = left / np.linalg.norm(left)
        right_unit = right / np.linalg.norm(right)
        direction_difference = np.abs(left_unit - right_unit)
        direction_error = float(direction_difference.max())
        direction_ok = bool((direction_difference <= atol + rtol * np.abs(right_unit)).all())
        entry = {"n_signatures": len(left_ids),
                 "max_abs_consensus_difference": consensus_error,
                 "max_abs_direction_difference": direction_error,
                 "within_tolerance": consensus_ok and direction_ok}
        if stored is not None:
            assert target in stored_directions, f"{target}: the rebuild stored no direction"
            raw_stored = np.asarray(stored_directions[target])
            assert raw_stored.shape == (len(stored_genes),), (
                f"{target}: the stored direction is {raw_stored.shape} on an axis of "
                f"{len(stored_genes)} genes")
            assert np.isfinite(raw_stored).all(), f"{target}: nonfinite stored direction"
            norm = float(np.linalg.norm(raw_stored))
            assert abs(norm - 1.0) <= 1e-5, (
                f"{target}: the stored direction has norm {norm:.6g}; the artifact's contract "
                "is a unit vector")
            stored_vector = raw_stored[stored_permutation]
            stored_difference = np.abs(stored_vector - left_unit)
            stored_error = float(stored_difference.max())
            stored_ok = bool((stored_difference <= atol + rtol * np.abs(left_unit)).all())
            entry["max_abs_stored_direction_difference"] = stored_error
            entry["within_tolerance"] = entry["within_tolerance"] and stored_ok
            if stored_error > worst["stored"]:
                worst["stored"], worst_target["stored"] = stored_error, target
        if consensus_error > worst["consensus"]:
            worst["consensus"], worst_target["consensus"] = consensus_error, target
        if direction_error > worst["direction"]:
            worst["direction"], worst_target["direction"] = direction_error, target
        per_target[target] = entry
        if not entry["within_tolerance"]:
            target_failures.append(target)

    return {"n_signatures": len(rebuilt_signatures), "n_targets": len(per_target),
            "max_abs_signature_difference": worst_signature,
            "n_signature_failures": len(signature_failures),
            "worst_signature_failures": sorted(
                signature_failures, key=lambda record: -record["max_abs_difference"])[:10],
            "signatures_within_tolerance": not signature_failures,
            "max_abs_consensus_difference": worst["consensus"],
            "max_abs_direction_difference": worst["direction"],
            "max_abs_stored_direction_difference": (worst["stored"] if stored is not None
                                                    else None),
            "stored_directions_compared": stored is not None,
            "worst_target": worst_target,
            "n_target_failures": len(target_failures),
            "targets_outside_tolerance": target_failures[:20],
            "all_within_tolerance": bool(not signature_failures and not target_failures),
            "rtol": rtol, "atol": atol, "per_target": per_target}


def load_shrna_rebuild_strict(directory: Path):
    """The loader the registered gate uses: a missing consensus artifact is fatal.

    Comparing source signatures and independently recomputed consensuses cannot
    show that the directions carried downstream are the right ones. Deviation 9
    was exactly a correct extraction with a wrong persisted direction, so the gate
    refuses to run without that artifact rather than quietly checking less.
    """
    signatures, membership, genes, stored = load_shrna_rebuild(directory)
    assert stored is not None, (
        f"{Path(directory) / 'shrna_consensus.npz'} is missing; the registered gate must "
        "compare the stored target directions against independently recomputed directions")
    return signatures, membership, genes, stored


def load_shrna_rebuild(directory: Path):
    """Independently rebuilt shRNA signatures, their target grouping, and the gene axis."""
    path = Path(directory) / "shrna_signatures.npz"
    assert path.exists(), f"{path} is missing; the rebuild must emit shRNA signatures with ids"
    with np.load(path, allow_pickle=True) as data:
        ids = [str(i) for i in data["sig_ids"]]
        matrix = data["signatures"]          # decompressed once, not once per id
        signatures = {sig_id: np.asarray(matrix[i]) for i, sig_id in enumerate(ids)}
        genes = [str(g) for g in data["gene_ids"]]
        membership = {str(target): [str(i) for i in members]
                      for target, members in json.loads(str(data["membership"])).items()}
    assert len(set(ids)) == len(ids), "the rebuild repeats a signature id"

    stored = None
    consensus_path = Path(directory) / "shrna_consensus.npz"
    if consensus_path.exists():
        with np.load(consensus_path, allow_pickle=True) as data:
            assert "gene_ids" in data.files, (
                f"{consensus_path} carries no gene axis; the gate cannot tell which "
                "coordinate system its directions use")
            directions = data["directions"]
            stored = ({str(target): np.asarray(directions[i])
                       for i, target in enumerate(data["genes"])},
                      [str(g) for g in data["gene_ids"]])
    return signatures, membership, genes, stored


def shrna_from_parsed_frame(frame, membership, landmark_ids):
    """Signatures keyed by id, carrying the gene labels the frame itself reports.

    The pairing of a value with a gene comes from the parsed frame's own index and
    is never reordered here. Reordering by labels read from the source is correct
    and is what coordinate harmonization requires; assigning expected labels
    positionally, before the source's own identifier-to-value pairing has been
    established, is what cannot be undone. Deviation 11 was a matrix whose columns
    did not follow the labels it declared, and reindexing an artifact that is
    already falsely labeled neither repairs nor reveals it.
    """
    labels = [str(gene) for gene in frame.index]
    assert len(set(labels)) == len(labels), "the parsed frame repeats a gene identifier"
    wanted = set(map(str, landmark_ids))
    keep = [i for i, gene in enumerate(labels) if gene in wanted]
    genes = [labels[i] for i in keep]
    assert len(genes) == len(wanted), (
        f"the parse covers {len(genes)} of {len(wanted)} landmark genes")

    columns = {str(sig_id): i for i, sig_id in enumerate(frame.columns)}
    assert len(columns) == frame.shape[1], "the parsed frame repeats a signature id"
    values = frame.to_numpy()
    signatures = {}
    for members in membership.values():
        for sig_id in members:
            if sig_id not in signatures:
                signatures[sig_id] = values[keep, columns[sig_id]]
    return signatures, genes


def load_shrna_from_gctx(gctx_path: Path, siginfo_path: Path, gene_info_path: Path,
                         targets_wanted=None, min_signatures=MIN_SIGNATURES):
    """The gate's own parse of the pinned GCTX, by a different route from production.

    The production extractor asks cmapPy for a fixed row order and reindexes the
    result onto it. This asks only for the landmark rows, reads the identifiers the
    file reports for them, and leaves them in that order: `compare_shrna` aligns
    the two axes on identifiers, so no order has to be assumed anywhere.
    """
    gene_info = pd.read_csv(gene_info_path, sep="\t", low_memory=False)
    landmark_ids = [str(gene) for gene in gene_info[gene_info.pr_is_lm == 1].pr_gene_id]
    assert len(landmark_ids) == N_LANDMARK, f"{len(landmark_ids)} landmark genes"

    siginfo = pd.read_csv(siginfo_path)
    membership = {}
    for gene, group in siginfo.groupby("pert_iname"):
        members = sorted(set(group.sig_id.astype(str)))
        if len(members) < min_signatures:
            continue
        if targets_wanted is None or gene in targets_wanted:
            membership[str(gene)] = members

    wanted = sorted({sig_id for members in membership.values() for sig_id in members})
    gct = parse.parse(str(gctx_path), cid=wanted, rid=landmark_ids)
    signatures, genes = shrna_from_parsed_frame(gct.data_df, membership, landmark_ids)
    return signatures, membership, genes


def _decoded(values) -> list:
    """Identifiers as text, without letting two of them collapse into one."""
    out = [value.decode() if isinstance(value, bytes) else str(value) for value in values]
    assert len(set(out)) == len(out), "the source metadata repeats an identifier"
    return out


def read_gctx_slice(gctx_path, wanted_signatures, wanted_genes, block=2000):
    """Read the registered signatures and genes straight out of the GCTX's HDF5.

    This is the gate's second source, and it is independent of production in the
    way that matters: the identifiers come from the file's own row and column
    metadata, the integer positions are built here from those identifiers, only the
    needed slices are read, and no expected label is ever attached to a value. The
    production route asks a parser for an order and trusts what comes back, so a
    parser whose identifier lookup is wrong would be wrong in both places.

    Returns (values, genes, signature ids, hashes). `genes` and the signature ids
    are read from the file, never from the request.
    """
    with h5py.File(gctx_path, "r") as handle:
        row_ids = _decoded(handle["/0/META/ROW/id"][:])
        column_ids = _decoded(handle["/0/META/COL/id"][:])
        matrix = handle["/0/DATA/0/matrix"]

        # orientation is decided by the metadata lengths, never by assuming one:
        # a transposed matrix must not pass on a shape coincidence
        fits = [name for name, shape in (("genes_then_signatures",
                                          (len(row_ids), len(column_ids))),
                                         ("signatures_then_genes",
                                          (len(column_ids), len(row_ids))))
                if tuple(matrix.shape) == shape]
        assert len(fits) == 1, (
            f"{Path(gctx_path).name} is {tuple(matrix.shape)} against {len(row_ids)} row and "
            f"{len(column_ids)} column identifiers, which fits {len(fits)} orientations; "
            "a shape coincidence must not decide which axis is which")
        orientation = fits[0]

        gene_position = {gene: i for i, gene in enumerate(row_ids)}
        signature_position = {sig: i for i, sig in enumerate(column_ids)}
        missing_genes = sorted(set(map(str, wanted_genes)) - set(gene_position))
        missing_signatures = sorted(set(map(str, wanted_signatures)) - set(signature_position))
        assert not missing_genes, (
            f"{len(missing_genes)} requested genes are absent from the source: "
            f"{missing_genes[:5]}")
        assert not missing_signatures, (
            f"{len(missing_signatures)} requested signatures are absent from the source: "
            f"{missing_signatures[:5]}")

        gene_rows = sorted(gene_position[gene] for gene in set(map(str, wanted_genes)))
        signature_columns = sorted(signature_position[sig]
                                   for sig in set(map(str, wanted_signatures)))
        genes = [row_ids[i] for i in gene_rows]
        signatures = [column_ids[i] for i in signature_columns]

        blocks = []
        for start in range(0, len(signature_columns), block):
            columns = signature_columns[start:start + block]
            if orientation == "signatures_then_genes":
                blocks.append(np.asarray(matrix[columns, :])[:, gene_rows])
            else:
                blocks.append(np.asarray(matrix[:, columns])[gene_rows, :].T)
        values = np.vstack(blocks) if blocks else np.empty((0, len(gene_rows)))

        hashes = {"source_row_axis_sha256": sha256_text("\n".join(row_ids)),
                  "source_column_axis_sha256": sha256_text("\n".join(column_ids)),
                  "selected_gene_axis_sha256": sha256_text("\n".join(genes)),
                  "selected_signature_axis_sha256": sha256_text("\n".join(signatures)),
                  "orientation": orientation}

    assert values.shape == (len(signatures), len(genes)), f"read {values.shape}"
    assert np.isfinite(values).all(), f"{Path(gctx_path).name} holds nonfinite values"
    return values, genes, signatures, hashes


def shrna_from_source(gctx_path, siginfo_path, gene_info_path, targets_wanted=None,
                      min_signatures=MIN_SIGNATURES):
    """The shRNA half's second source: the pinned GCTX, read independently."""
    gene_info = pd.read_csv(gene_info_path, sep="\t", low_memory=False)
    landmark_ids = [str(gene) for gene in gene_info[gene_info.pr_is_lm == 1].pr_gene_id]
    assert len(landmark_ids) == N_LANDMARK, f"{len(landmark_ids)} landmark genes"

    siginfo = pd.read_csv(siginfo_path)
    membership = {}
    for gene, group in siginfo.groupby("pert_iname"):
        members = sorted(set(group.sig_id.astype(str)))
        if len(members) < min_signatures:
            continue
        if targets_wanted is None or gene in targets_wanted:
            membership[str(gene)] = members

    wanted = sorted({sig_id for members in membership.values() for sig_id in members})
    values, genes, signatures, hashes = read_gctx_slice(gctx_path, wanted, landmark_ids)
    assert set(genes) == set(landmark_ids), "the source returned a different gene set"
    assert set(signatures) == set(wanted), "the source returned a different signature set"
    rows = {sig_id: i for i, sig_id in enumerate(signatures)}
    return ({sig_id: values[rows[sig_id]] for sig_id in signatures}, membership, genes), hashes


def compounds_from_source(gctx_path, siginfo_path, gene_info_path, cohort):
    """The compound half's second source, by the same reader.

    Amendment 2 retains `lincs_subset.npz` on the strength of its earlier pass.
    Retaining it is a claim about its values, so the claim is checked against the
    source rather than carried forward.
    """
    gene_info = pd.read_csv(gene_info_path, sep="\t", low_memory=False)
    landmark_ids = [str(gene) for gene in gene_info[gene_info.pr_is_lm == 1].pr_gene_id]
    siginfo = pd.read_csv(siginfo_path, sep="\t", low_memory=False)
    siginfo = siginfo[siginfo.pert_iname.isin(set(cohort))].copy()
    siginfo["sig_id"] = siginfo.sig_id.astype(str)

    wanted = sorted(set(siginfo.sig_id))
    values, genes, signatures, hashes = read_gctx_slice(gctx_path, wanted, landmark_ids)
    assert set(genes) == set(landmark_ids), "the source returned a different gene set"
    rows = {sig_id: i for i, sig_id in enumerate(signatures)}

    matrices, cells = {}, {}
    for drug, group in siginfo.groupby("pert_iname"):
        per_cell, identifiers = [], []
        for cell_id, members in group.groupby("cell_id"):
            positions = [rows[sig_id] for sig_id in sorted(set(members.sig_id))]
            per_cell.append(values[positions].mean(axis=0))
            identifiers.append(str(cell_id))
        matrices[drug] = np.vstack(per_cell)
        cells[drug] = identifiers
    return (matrices, cells, genes), hashes


def load_shrna_extraction(data_dir: Path, targets_wanted=None):
    """The same objects from the pinned extraction, grouped by the pinned metadata rule.

    The signature matrix is decompressed once. Indexing an NpzFile inside a loop
    re-reads and re-inflates the whole array on every access, which on 154,993
    signatures is hundreds of gigabytes of allocation rather than one.
    """
    shrna = np.load(data_dir / "lincs_shrna.npz", allow_pickle=True)
    ids = [str(sig_id) for sig_id in shrna["sig_ids"]]
    assert len(set(ids)) == len(ids), "the extraction repeats a signature id"
    row_of = {sig_id: i for i, sig_id in enumerate(ids)}
    genes = [str(g) for g in shrna["gene_ids"]]

    siginfo = pd.read_csv(data_dir / "lincs_shrna_siginfo.csv.gz")
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(row_of)]
    membership = {}
    for gene, group in siginfo.groupby("pert_iname"):
        members = sorted(set(group.sig_id.astype(str)))
        if len(members) < MIN_SIGNATURES:
            continue
        if targets_wanted is None or gene in targets_wanted:
            membership[str(gene)] = members

    # only the signatures the compared targets use are held in memory
    wanted = sorted({sig_id for members in membership.values() for sig_id in members})
    matrix = shrna["signatures"]
    signatures = {sig_id: np.asarray(matrix[row_of[sig_id]]) for sig_id in wanted}
    del matrix
    return signatures, membership, genes


def compare_both(compound_rebuilt, compound_extraction, shrna_rebuilt, shrna_extraction,
                 cohort, stored=None):
    """Both registered halves, with a structural failure named rather than thrown.

    An assertion inside either comparison used to propagate before anything was
    written, so the single failure the gate exists to document left no artifact of
    itself. It is caught here and recorded; the caller writes the result and then
    raises.
    """
    halves = ("compound_signatures", "shrna_signatures_and_consensuses")
    result = {}
    try:
        result["compound_signatures"] = compare(compound_rebuilt, compound_extraction,
                                                cohort=cohort)
        result["shrna_signatures_and_consensuses"] = compare_shrna(
            shrna_rebuilt, shrna_extraction, stored=stored)
    except AssertionError as failure:
        result["structural_failure"] = str(failure)
    result["all_within_tolerance"] = bool(
        not result.get("structural_failure")
        and all(result.get(half, {}).get("all_within_tolerance") for half in halves))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rebuilt", type=Path, required=True,
                        help="directory holding the GCTX rebuild's shards")
    parser.add_argument("--extraction", type=Path, required=True)
    parser.add_argument("--cohort", type=Path, default=REFERENCE_ARTIFACT,
                        help="the artifact whose drugs and targets both sides must hold")
    parser.add_argument("--gctx", type=Path,
                        help="the pinned GCTX. Given, the shRNA half compares the canonical "
                             "artifact against this gate's own parse of the source, which is "
                             "what Amendment 2 requires; `lincs_shrna.npz` is then read by "
                             "nothing. Omitted, the retired extraction is the second source.")
    parser.add_argument("--shrna-siginfo", type=Path,
                        help="the pinned shRNA metadata, with --gctx")
    parser.add_argument("--gene-info", type=Path,
                        help="the pinned gene info, with --gctx")
    parser.add_argument("--compound-siginfo", type=Path,
                        help="the pinned compound metadata, with --gctx. Both halves are then "
                             "verified against the source, which is what retaining "
                             "`lincs_subset.npz` claims.")
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    assert (bool(args.gctx) == bool(args.shrna_siginfo) == bool(args.gene_info)
            == bool(args.compound_siginfo)), (
        "--gctx needs --shrna-siginfo, --compound-siginfo and --gene-info, which pin the "
        "membership rules and the landmark set it reads against")

    records = json.loads(args.cohort.read_text())
    drugs = sorted({record["drug"] for record in records})
    assert len(drugs) == len(records), (
        f"{len(records)} records hold {len(drugs)} distinct drugs; the cohort unit is the drug")
    targets = sorted({str(record["target"]) for record in records})

    rebuilt_signatures, rebuilt_membership, rebuilt_genes, stored = load_shrna_rebuild_strict(
        args.rebuilt)
    source_hashes = {}
    if args.gctx:
        shrna_second, source_hashes["shrna"] = shrna_from_source(
            args.gctx, args.shrna_siginfo, args.gene_info, targets_wanted=set(targets))
        compound_second, source_hashes["compound"] = compounds_from_source(
            args.gctx, args.compound_siginfo, args.gene_info, cohort=drugs)
    else:
        shrna_second = load_shrna_extraction(args.extraction, targets_wanted=set(targets))
        compound_second = load_extraction(args.extraction, cohort=drugs)
    result = compare_both(
        load_rebuild(args.rebuilt),
        compound_second,
        (rebuilt_signatures, rebuilt_membership, rebuilt_genes),
        shrna_second,
        cohort=drugs, stored=stored)
    if args.gctx:
        result["second_source"] = {
            "route": "the gate's own HDF5 read of " + str(args.gctx),
            "gctx_sha256": sha256_file(args.gctx), "axes": source_hashes,
            "halves_verified_against_the_source": ["compound", "shrna"]}
        result["canonical_artifact_sha256"] = {
            name: sha256_file(Path(args.rebuilt) / name)
            for name in ("shrna_signatures.npz", "shrna_consensus.npz")
            if (Path(args.rebuilt) / name).exists()}
    else:
        result["second_source"] = {"route": "lincs_shrna.npz, retired by Amendment 2"}

    result["cohort"] = {"file": str(args.cohort), "sha256": sha256_file(args.cohort),
                        "n_records": len(records), "n_unique_drugs": len(drugs),
                        "n_unique_targets": len(targets)}
    result["identifier_hashes"] = {"drugs": sha256_text("\n".join(drugs)),
                                   "targets": sha256_text("\n".join(targets))}
    result["gate_code_sha256"] = sha256_file(Path(__file__))
    passed = result["all_within_tolerance"]
    result["gate"] = ("reconstruction" if passed
                      else "failed: the analyses registered against this gate are void")

    args.output.mkdir(parents=True, exist_ok=True)
    path = args.output / "reconstruction_gate.json"
    path.write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("compound_signatures", "shrna_signatures_and_consensuses")},
                     indent=2))
    assert passed, f"the gate failed; the record is at {path}"


if __name__ == "__main__":
    main()
