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
from dataclasses import dataclass, field

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
            assert vector.shape == self.positions.shape, f"{gene}: {vector.shape} on {self.positions.shape}"

    def restricted_to(self, positions: np.ndarray) -> "Reference":
        """The same directions on a narrower gene space, renormalized after cutting."""
        keep = np.flatnonzero(np.isin(self.positions, positions))
        assert len(keep) == len(positions), (
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


def unit(vector: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vector)
    assert norm > 0, "a direction with zero norm cannot be normalized"
    return vector / norm


def frozen_landmark_order(gene_info_path) -> list:
    """The 978 landmark gene ids, ascending by Entrez id, as the rebuild freezes them."""
    gene_info = pd.read_csv(gene_info_path, sep="\t", low_memory=False)
    landmark = gene_info[gene_info["pr_is_lm"] == 1].sort_values("pr_gene_id")
    ids = [str(gene) for gene in landmark["pr_gene_id"]]
    assert len(ids) == N_LANDMARK, f"{len(ids)} landmark genes, expected {N_LANDMARK}"
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
    assert len(set(declared)) == len(declared), f"{source} repeats a gene identifier"
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


def load_replogle_bulk(h5ad_path, symbols, name, qualifying_rows=None) -> Reference:
    """A reference from an official Replogle pseudobulk file.

    Args:
        h5ad_path: one of the released `*_normalized_bulk_01.h5ad` files.
        symbols: landmark symbols in extraction order.
        name: reference name, recorded in the output.
        qualifying_rows: optional predicate on the observation table. A gene is
            eligible when at least one row qualifies, and its direction is built
            from the qualifying rows alone. This is R4's phenotype-positive rule.
    """
    adata = ad.read_h5ad(h5ad_path)
    obs = adata.obs.copy()
    obs["gene"] = [str(name_).split("_")[1] if str(name_).count("_") >= 2 else str(name_)
                   for name_ in adata.obs_names]
    file_genes = adata.var["gene_name"].astype(str).to_numpy() if "gene_name" in adata.var else \
        adata.var_names.astype(str).to_numpy()

    positions, matched, audit = _shared_positions(symbols, file_genes)
    column_of = {gene: i for i, gene in enumerate(file_genes)}
    columns = np.array([column_of[symbol] for symbol in matched])

    matrix = np.asarray(adata.X[:, columns], dtype=np.float64)
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

    audit.update({"file": str(h5ad_path), "n_genes_in_file": int(obs["gene"].nunique()),
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
    assert control.sum() > 0, "no control cells in the pooled construction"

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
