"""Target-perturbation references, built as the registration frozen at 7f57136 fixes them.

A reference gives each target gene a direction in the LINCS landmark space. Two
families are built here:

C0   the pooled construction of `03b_h3_crispri_ground_truth.py`: the mean over a
     target's cells minus the mean over all non-targeting cells, placed in the 978
     landmarks with zeros where a landmark is absent from the file. This is the
     construction the known values were computed under.
C1   the official Replogle pseudobulk, z-normalized against controls within each
     GEM group, restricted to the landmarks the file measures. C1 is primary, and
     the same procedure serves the RPE1 and genome-wide releases.

A direction is unit-normalized after its rows are averaged and after it is
restricted to a gene space, never before.
"""
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

N_LANDMARK = 978
MIN_CELLS_PER_TARGET = 10     # as in the pooled construction
MIN_LANDMARKS = 500           # C1 analyses are void below this


@dataclass(frozen=True)
class Reference:
    """Directions for one reference, on one gene space."""

    name: str
    directions: dict            # gene symbol -> unit vector over `positions`
    positions: np.ndarray       # indices into the 978-landmark order
    audit: dict = field(default_factory=dict)

    def __post_init__(self):
        for gene, vector in self.directions.items():
            if vector.shape != self.positions.shape:
                raise AssertionError(
                    f"{gene}: {vector.shape} on {self.positions.shape}")

    def restricted_to(self, positions: np.ndarray) -> "Reference":
        """The same directions on a narrower gene space, renormalized after cutting."""
        keep = np.flatnonzero(np.isin(self.positions, positions))
        if len(keep) != len(positions):
            raise AssertionError(
                f"{self.name} does not cover the requested gene space: "
                f"{len(keep)} of {len(positions)} positions")
        cut = {gene: unit(vector[keep]) for gene, vector in self.directions.items()}
        return Reference(f"{self.name}|restricted", cut, np.asarray(positions),
                         {**self.audit, "restricted_from": len(self.positions)})

    def without_gene(self, gene: str, landmark_symbols) -> "Reference":
        """The reference with one gene's own coordinate dropped, where it is present."""
        symbols = np.asarray(landmark_symbols)[self.positions]
        keep = symbols != gene
        if keep.all():
            return self
        positions = self.positions[keep]
        cut = {g: unit(v[keep]) for g, v in self.directions.items()}
        return Reference(f"{self.name}|no-{gene}", cut, positions, dict(self.audit))


def require_finite(array, what) -> np.ndarray:
    """Refuse a non-finite operand instead of carrying it into a statistic.

    Amendment 4 removes a coordinate only through a frozen analysis basis. A
    non-finite value anywhere else is a defect that basis does not cover, so it
    stops the calculation where it is found rather than being dropped there.
    """
    array = np.asarray(array, dtype=np.float64)
    offending = ~np.isfinite(array)
    if offending.any():
        raise AssertionError(
            f"{what}: {int(offending.sum())} of {array.size} values are not finite, "
            "and a coordinate leaves an analysis only through its frozen basis")
    return array


def unit(vector: np.ndarray) -> np.ndarray:
    vector = require_finite(vector, "a direction before normalization")
    norm = np.linalg.norm(vector)
    if not norm > 0:
        raise AssertionError("a direction with zero norm cannot be normalized")
    return vector / norm


def sha256_file(path) -> str:
    """The digest of a file, read in blocks so a large release does not land in memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def basis_sha256(gene_ids) -> str:
    """The canonical hash of an ordered basis, over its landmark gene ids."""
    return hashlib.sha256("\n".join(str(gene) for gene in gene_ids).encode()).hexdigest()


def frozen_landmark_order(gene_info_path) -> list:
    """The 978 landmark gene ids, ascending by Entrez id, as the rebuild freezes them."""
    gene_info = pd.read_csv(gene_info_path, sep="\t", low_memory=False)
    landmark = gene_info[gene_info["pr_is_lm"] == 1].sort_values("pr_gene_id")
    ids = [str(gene) for gene in landmark["pr_gene_id"]]
    if len(ids) != N_LANDMARK:
        raise AssertionError(f"{len(ids)} landmark genes, expected {N_LANDMARK}")
    if len(set(ids)) != len(ids):
        raise AssertionError("the frozen landmark order repeats a gene id")
    return ids


def check_declared_axis(declared, frozen, source) -> None:
    """Refuse an artifact whose declared gene axis is not the frozen landmark order.

    This is a statement about labels, and labels can be wrong about the matrix they
    sit on: Deviation 11 was a matrix whose columns did not follow the identifiers
    it declared, and both extractions declared the same identifiers, so this check
    passes on it. The coordinates themselves are verified against a parse of the
    pinned source by the reconstruction gate. Amendment 2 registers both layers and
    neither replaces the other.
    """
    declared, frozen = [str(gene) for gene in declared], [str(gene) for gene in frozen]
    if len(set(declared)) != len(declared):
        raise AssertionError(f"{source} repeats a gene identifier")
    if declared == frozen:
        return
    position = next((i for i, (left, right) in enumerate(zip(declared, frozen))
                     if left != right), None)
    detail = (f"; first difference at position {position}: {declared[position]} against "
              f"{frozen[position]}" if position is not None else "")
    raise AssertionError(
        f"{source} declares {len(declared)} genes against the frozen order's "
        f"{len(frozen)}{detail}" if len(declared) != len(frozen) else
        f"{source} declares a gene axis that is not the frozen landmark order{detail}")


def landmark_symbols(gene_info_path, gene_ids) -> list:
    """Symbols for the landmark gene ids, in the order the extraction holds them."""
    gene_info = pd.read_csv(gene_info_path, sep="\t", low_memory=False)
    landmark = gene_info[gene_info["pr_is_lm"] == 1]
    entrez_to_symbol = dict(zip(landmark["pr_gene_id"].astype(str), landmark["pr_gene_symbol"]))
    return [entrez_to_symbol.get(str(gene_id), str(gene_id)) for gene_id in gene_ids]


def _shared_positions(symbols, file_genes):
    """Landmark positions whose symbol appears exactly once on both sides.

    A symbol duplicated in either gene list cannot be matched without choosing
    between rows, so it is dropped and counted.
    """
    symbols = np.asarray(symbols)
    landmark_counts = pd.Series(symbols).value_counts()
    file_counts = pd.Series(file_genes).value_counts()
    positions, matched, duplicated, unmatched = [], [], [], []
    for position, symbol in enumerate(symbols):
        if landmark_counts.get(symbol, 0) > 1 or file_counts.get(symbol, 0) > 1:
            duplicated.append(symbol)
            continue
        if file_counts.get(symbol, 0) == 1:
            positions.append(position)
            matched.append(symbol)
        else:
            unmatched.append(symbol)
    audit = {"n_matched": len(matched), "n_unmatched_landmarks": len(unmatched),
             "n_duplicated_symbols": len(set(duplicated)),
             "duplicated_symbols": sorted(set(duplicated))[:50]}
    return np.array(positions, dtype=int), matched, audit


def _dense(block) -> np.ndarray:
    """A float64 array from a released matrix block, sparse or not."""
    return np.asarray(block.toarray() if hasattr(block, "toarray") else block, dtype=np.float64)


def _file_gene_names(adata) -> np.ndarray:
    """The gene names a released file carries, from `gene_name` where it has one."""
    return (adata.var["gene_name"].astype(str).to_numpy() if "gene_name" in adata.var
            else adata.var_names.astype(str).to_numpy())


def _bulk_columns(adata, symbols):
    """Matched landmark positions, their symbols, and the file columns holding them."""
    file_genes = _file_gene_names(adata)
    positions, matched, audit = _shared_positions(symbols, file_genes)
    column_of = {gene: i for i, gene in enumerate(file_genes)}
    columns = np.array([column_of[symbol] for symbol in matched], dtype=int)
    return positions, matched, columns, audit


def finite_landmark_basis(h5ad_path, symbols, gene_ids, name) -> dict:
    """One construction's analysis basis, and why each coordinate is in it or not.

    The basis is the ordered intersection of the frozen landmark axis with the
    genes this source represents exactly once, narrowed to the coordinates the
    source defines with finite values throughout. A coordinate leaves the basis
    because the source does not define it, never because of a quantity computed
    from it, so the rule is outcome-independent and generalizes past any named
    gene. Amendment 4 registers it.

    This is the one reader allowed to see a non-finite value: it measures where
    they are. Every other loader refuses them.
    """
    adata = ad.read_h5ad(h5ad_path)
    positions, matched, columns, audit = _bulk_columns(adata, symbols)
    matrix = _dense(adata.X[:, columns])
    offending = ~np.isfinite(matrix)
    defined = ~offending.any(axis=0)

    excluded = []
    for index in np.flatnonzero(~defined):
        column = matrix[:, index]
        excluded.append({
            "symbol": matched[index],
            "landmark_gene_id": str(gene_ids[positions[index]]),
            "landmark_position": int(positions[index]),
            "n_non_finite_values": int((~np.isfinite(column)).sum()),
            "n_rows": int(len(column)),
            "n_nan": int(np.isnan(column).sum()),
            "n_posinf": int(np.isposinf(column).sum()),
            "n_neginf": int(np.isneginf(column).sum()),
            "reason": "the source leaves this coordinate undefined"})

    kept = np.flatnonzero(defined)
    basis_positions = positions[kept]
    basis_ids = [str(gene_ids[position]) for position in basis_positions]
    return {
        "construction": name,
        "file": Path(h5ad_path).name,
        "file_sha256": _sha256(h5ad_path),
        "n_rows_in_source": int(matrix.shape[0]),
        "n_landmarks_in_frozen_axis": int(len(symbols)),
        "n_landmarks_matched": int(len(positions)),
        "n_landmarks_in_basis": int(len(basis_positions)),
        "n_excluded_as_undefined": len(excluded),
        "excluded_as_undefined": excluded,
        "landmark_positions": [int(position) for position in basis_positions],
        "landmark_gene_ids": basis_ids,
        "landmark_symbols": [matched[index] for index in kept],
        "basis_sha256": basis_sha256(basis_ids),
        "registered_floor": MIN_LANDMARKS,
        "meets_registered_floor": bool(len(basis_positions) >= MIN_LANDMARKS),
    }


# every field of the genome-wide record that its phenotype-filtered twin must
# repeat, since they read one file and take one basis
SHARED_BASIS_FIELDS = ("file", "file_sha256", "n_rows_in_source", "n_landmarks_matched",
                       "n_landmarks_in_basis", "landmark_positions", "landmark_gene_ids",
                       "landmark_symbols", "basis_sha256", "excluded_as_undefined",
                       "n_excluded_as_undefined", "registered_floor",
                       "meets_registered_floor")


def validated_analysis_bases(path, gene_ids, expected, generator=None) -> dict:
    """The committed bases, checked against an external pin before anything reads them.

    One validator for every program that consumes the bases, because a basis
    checked in the driver and trusted in the single-cell stage is checked in one
    of the two places it is used.

    A basis file that carried only its own hash would attest to itself: an edited
    list recomputes its own digest and passes. `expected` therefore comes from the
    amendment-specific external pin, `registry/frozen/analysis_bases_pin.json`,
    which lives outside the artifact it describes and pins the whole file, the
    generator, the schema, the construction names, each source and its digest,
    each basis digest and size, and the landmark axis. Each amendment writes its
    own pin artifact and pins that artifact's digest in its own text; an earlier
    amendment's pin file is never rewritten.

    Args:
        path: the bases artifact.
        gene_ids: the frozen landmark order this run is using.
        expected: the amendment's external pin, as a dict.
        generator: the generator script, where its digest is to be checked.
    """
    path = Path(path)
    if not path.exists():
        raise AssertionError(
            f"{path} is absent, and Amendment 4 evaluates every C1 construction on a "
            "frozen analysis basis. Generate it with "
            "experiments/03i_freeze_analysis_bases.py; a basis is not recomputed by the "
            "run it governs, because a basis computed inside that run is not frozen.")
    digest = _sha256(path)
    if digest != expected["file_sha256"]:
        raise AssertionError(
            f"{path.name} hashes to {digest} and the frozen identities pin "
            f"{expected['file_sha256']}; the bases artifact has changed since it was pinned")
    if generator is not None:
        generated_by = _sha256(generator)
        if generated_by != expected["generator_sha256"]:
            raise AssertionError(
                f"{Path(generator).name} hashes to {generated_by} and the frozen identities "
                f"pin {expected['generator_sha256']}; the code that measures the bases has "
                "changed since they were pinned")

    document = json.loads(path.read_text())
    if document.get("schema_version") != expected["schema_version"]:
        raise AssertionError(
            f"the bases artifact declares schema version {document.get('schema_version')} "
            f"and the frozen identities pin {expected['schema_version']}")
    axis = document.get("gene_info", {}).get("landmark_axis_sha256")
    if axis != expected["landmark_axis_sha256"] or axis != basis_sha256(gene_ids):
        raise AssertionError(
            f"the bases were frozen against landmark axis {axis}, the frozen identities pin "
            f"{expected['landmark_axis_sha256']}, and this run's axis hashes to "
            f"{basis_sha256(gene_ids)}")

    bases = document["constructions"]
    if sorted(bases) != sorted(expected["constructions"]):
        raise AssertionError(
            f"the bases artifact holds {sorted(bases)} and the frozen identities pin "
            f"{sorted(expected['constructions'])}")

    for name, record in bases.items():
        pinned = expected["constructions"][name]
        ids, positions = record["landmark_gene_ids"], record["landmark_positions"]
        for field in ("file", "file_sha256", "basis_sha256", "n_landmarks_in_basis",
                      "n_landmarks_matched"):
            if record[field] != pinned[field]:
                raise AssertionError(
                    f"{name}: {field} is {record[field]!r} and the frozen identities pin "
                    f"{pinned[field]!r}")
        if not (len(ids) == len(positions) == len(record["landmark_symbols"])
                == record["n_landmarks_in_basis"]):
            raise AssertionError(
                f"{name}: {len(ids)} gene ids, {len(positions)} positions, "
                f"{len(record['landmark_symbols'])} symbols and a declared "
                f"{record['n_landmarks_in_basis']} do not agree")
        if basis_sha256(ids) != record["basis_sha256"]:
            raise AssertionError(
                f"{name}: the basis does not hash to its recorded {record['basis_sha256']}")
        if sorted(set(positions)) != list(positions):
            raise AssertionError(
                f"{name}: the basis is not in the frozen landmark order, or repeats a position")
        if [str(gene_ids[position]) for position in positions] != [str(gene) for gene in ids]:
            raise AssertionError(
                f"{name}: the basis's gene ids are not this run's axis at its own positions, "
                "so it was frozen against a different landmark order")
        if len(ids) < MIN_LANDMARKS or not record["meets_registered_floor"]:
            raise AssertionError(
                f"{name}: {len(ids)} landmarks, below the registered floor of {MIN_LANDMARKS}")

    # C0 is not one of the file-backed bases: it declares the whole axis with zeros
    # where its file measures nothing. Its source and both spaces are pinned, so the
    # intersection R0.5 compares it on is not taken on trust.
    pooled, pinned_pool = document["pooled_construction"], expected["pooled_construction"]
    for field in ("file", "file_sha256", "n_landmarks_declared", "n_landmarks_measured",
                  "declared_sha256", "measured_sha256"):
        if pooled[field] != pinned_pool[field]:
            raise AssertionError(
                f"C0: {field} is {pooled[field]!r} and the pin names {pinned_pool[field]!r}")
    if basis_sha256(pooled["measured_gene_ids"]) != pooled["measured_sha256"]:
        raise AssertionError("C0: the measured coordinates do not hash to their record")
    if len(pooled["measured_positions"]) != pooled["n_landmarks_measured"]:
        raise AssertionError(
            f"C0: {len(pooled['measured_positions'])} measured positions against a declared "
            f"{pooled['n_landmarks_measured']}")

    # the intersections the registered comparisons run on, pinned by count
    for key, n_shared in expected["registered_comparisons"].items():
        if document["registered_comparisons"][key]["n_shared"] != n_shared:
            raise AssertionError(
                f"{key}: the artifact records "
                f"{document['registered_comparisons'][key]['n_shared']} shared landmarks and "
                f"the pin names {n_shared}")
    if sorted(document["registered_comparisons"]) != sorted(expected["registered_comparisons"]):
        raise AssertionError(
            f"the artifact records comparisons {sorted(document['registered_comparisons'])} "
            f"and the pin names {sorted(expected['registered_comparisons'])}")

    # the phenotype-filtered arm reads the genome-wide file and takes its basis, so
    # the two records agree on everything but their own names
    for name, source_of in document.get("shares_a_basis", {}).items():
        for field in SHARED_BASIS_FIELDS:
            if bases[name][field] != bases[source_of][field]:
                raise AssertionError(
                    f"{name} and {source_of} read one file and must share one basis, and "
                    f"their {field} differ")
    return bases


def _restricted_to_basis(basis, positions, matched, columns, audit, source):
    """The file's matched landmarks narrowed to a frozen basis, refusing a gap.

    A basis coordinate the file does not represent is not dropped: the basis
    describes this source, so a gap means the file is not the file the basis was
    frozen on.
    """
    basis = np.asarray(basis, dtype=int)
    if basis.size != len(set(basis.tolist())):
        raise AssertionError(f"{source}: the basis repeats a position")
    if not np.all(np.diff(basis) > 0):
        raise AssertionError(f"{source}: the basis is not in the frozen landmark order")
    missing = np.setdiff1d(basis, positions)
    if len(missing):
        raise AssertionError(
            f"{source} represents {len(basis) - len(missing)} of the frozen basis's "
            f"{len(basis)} landmark positions; {len(missing)} are absent, first "
            f"{sorted(missing.tolist())[:10]}. The basis was frozen on another file.")
    index_of = {int(position): index for index, position in enumerate(positions)}
    keep = np.array([index_of[int(position)] for position in basis], dtype=int)
    dropped = sorted(set(positions.tolist()) - set(basis.tolist()))
    narrowed = {**audit,
                "n_landmarks_in_basis": int(len(basis)),
                "n_excluded_by_basis": len(dropped),
                "excluded_by_basis": [matched[index_of[position]] for position in dropped]}
    return basis, [matched[index] for index in keep], columns[keep], narrowed


def load_replogle_bulk(h5ad_path, symbols, name, qualifying_rows=None, basis=None) -> Reference:
    """A reference from an official Replogle pseudobulk file.

    Args:
        h5ad_path: one of the released `*_normalized_bulk_01.h5ad` files.
        symbols: landmark symbols in extraction order.
        name: reference name, recorded in the output.
        qualifying_rows: optional predicate on the observation table. A gene is
            eligible when at least one row qualifies, and its direction is built
            from the qualifying rows alone. This is R4's phenotype-positive rule.
        basis: the construction's record from `registry/frozen/analysis_bases.json`.
            Directions are built on its coordinates alone, and its `file_sha256`
            must be this file's, so a basis cannot be applied to a source other
            than the one it was frozen on. Amendment 4 requires it wherever a
            source leaves a matched landmark undefined; where every matched
            landmark is defined, passing it changes nothing and omitting it still
            refuses a non-finite value.
    """
    adata = ad.read_h5ad(h5ad_path)
    obs = adata.obs.copy()
    obs["gene"] = [str(name_).split("_")[1] if str(name_).count("_") >= 2 else str(name_)
                   for name_ in adata.obs_names]
    positions, matched, columns, audit = _bulk_columns(adata, symbols)
    source_sha256 = _sha256(h5ad_path)
    if basis is not None:
        # the basis names the file it was frozen on, and this is that file or the
        # load stops: a basis applied to a different release is not a basis
        if basis["file_sha256"] != source_sha256:
            raise AssertionError(
                f"{name}: the frozen basis was measured on a file with sha256 "
                f"{basis['file_sha256']} and {Path(h5ad_path).name} hashes to "
                f"{source_sha256}, so the basis does not describe this source")
        positions, matched, columns, audit = _restricted_to_basis(
            basis["landmark_positions"], positions, matched, columns, audit, str(h5ad_path))
        audit["basis_sha256"] = basis["basis_sha256"]

    # on the basis, the source is finite or the load stops: no statistic is
    # computed from an undefined coordinate, and none is dropped quietly either
    matrix = require_finite(_dense(adata.X[:, columns]),
                            f"{name} read from {Path(h5ad_path).name} on its analysis basis")
    qualifies = np.ones(len(obs), dtype=bool) if qualifying_rows is None else \
        np.asarray(qualifying_rows(obs), dtype=bool)

    directions, rows_total, rows_used, dropped = {}, {}, {}, []
    for gene, rows in obs.groupby("gene", observed=True).indices.items():
        rows_total[gene] = int(len(rows))
        used = rows[qualifies[rows]]
        if len(used) == 0:
            dropped.append(gene)
            continue
        averaged = matrix[used].mean(axis=0)     # rows averaged before normalization
        if np.linalg.norm(averaged) == 0:
            dropped.append(gene)
            continue
        directions[gene] = unit(averaged)
        rows_used[gene] = int(len(used))

    audit.update({"file": str(h5ad_path), "file_sha256": source_sha256,
                  "positions_sha256": basis_sha256(positions.tolist()),
                  "n_landmarks_used": int(len(positions)),
                  "n_genes_in_file": int(obs["gene"].nunique()),
                  "n_genes_with_direction": len(directions),
                  "n_genes_dropped": len(dropped),
                  "rows_per_gene": rows_total, "rows_used_per_gene": rows_used,
                  "row_filter": qualifying_rows is not None})
    return Reference(name, directions, positions, audit)


def pooled_crispri_reference(h5ad_path, symbols, name="C0") -> Reference:
    """The pooled construction the known CRISPRi values were computed under.

    Mean over a target's cells minus the mean over all non-targeting cells, on the
    landmarks the file measures, placed back into the full 978 with zeros
    elsewhere. Reproduced here so that C1 can be compared against it.
    """
    adata = ad.read_h5ad(h5ad_path)
    file_genes = adata.var_names.astype(str).to_numpy()
    symbols = np.asarray(symbols)
    shared = [gene for gene in file_genes if gene in set(symbols)]
    file_column = {gene: i for i, gene in enumerate(file_genes)}
    columns = np.array([file_column[gene] for gene in shared])
    landmark_position = {symbol: i for i, symbol in enumerate(symbols)}
    targets = np.array([landmark_position[gene] for gene in shared])

    labels = adata.obs["gene"].astype(str).to_numpy()
    control = np.isin(labels, ["non-targeting", "control", ""])
    if control.sum() == 0:
        raise AssertionError("no control cells in the pooled construction")

    matrix = adata.X
    dense = (lambda rows: np.asarray(matrix[rows][:, columns].todense(), dtype=np.float64)) \
        if hasattr(matrix, "todense") else \
        (lambda rows: np.asarray(matrix[rows][:, columns], dtype=np.float64))

    control_mean = dense(control).mean(axis=0)
    directions, n_cells = {}, {}
    for gene in np.unique(labels):
        if gene in ("non-targeting", "control", ""):
            continue
        rows = labels == gene
        if rows.sum() < MIN_CELLS_PER_TARGET:
            continue
        difference = dense(rows).mean(axis=0) - control_mean
        full = np.zeros(N_LANDMARK)
        full[targets] = difference
        if np.linalg.norm(full) == 0:
            continue
        directions[gene] = unit(full)
        n_cells[gene] = int(rows.sum())

    audit = {"file": str(h5ad_path), "n_shared_landmarks": len(shared),
             "n_control_cells": int(control.sum()), "n_genes_with_direction": len(directions),
             "cells_per_gene": n_cells, "min_cells": MIN_CELLS_PER_TARGET,
             "note": "zeros outside the shared landmarks, as the construction that produced "
                     "the known values does"}
    return Reference(name, directions, np.arange(N_LANDMARK), audit)


def shared_space(*references) -> np.ndarray:
    """Landmark positions every reference measures."""
    positions = references[0].positions
    for reference in references[1:]:
        positions = np.intersect1d(positions, reference.positions)
    return positions


def alignment_matrix(mean_signatures: np.ndarray, reference: Reference, genes) -> np.ndarray:
    """(n_drugs, n_genes) squared cosine between each drug mean and each direction.

    `mean_signatures` is already restricted to `reference.positions`.
    """
    directions = np.vstack([reference.directions[gene] for gene in genes])
    means = mean_signatures / np.linalg.norm(mean_signatures, axis=1, keepdims=True)
    return (means @ directions.T) ** 2


def projected_dispersion(pair_differences: np.ndarray, direction: np.ndarray) -> float:
    """Mean |(s_i - s_j) . u| over the distinct context pairs of one drug."""
    return float(np.abs(pair_differences @ direction).mean())
