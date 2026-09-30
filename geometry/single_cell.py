"""Single-cell audit helpers for R0.1, R0.3, R0.4 and R0.6 of the frozen registration.

The Modal wrapper (`experiments/modal_03d_single_cell.py`) arranges workers,
volumes and checkpoints; everything that computes a number lives here, so it can
be tested without Modal and without the 19 GB of released single-cell data.
"""
import hashlib

import numpy as np

N_SPLITS = 50


def representation_audit(sample: np.ndarray) -> dict:
    """What a matrix actually holds: counts, normalized values, or something else."""
    values = np.asarray(sample, dtype=np.float64).ravel()
    finite = values[np.isfinite(values)]
    return {"dtype": str(np.asarray(sample).dtype),
            "min": float(finite.min()), "max": float(finite.max()),
            "mean": float(finite.mean()),
            "fraction_nonzero": float((finite != 0).mean()),
            "fraction_integer_valued": float((finite == np.round(finite)).mean()),
            "fraction_negative": float((finite < 0).mean()),
            "n_values_sampled": int(finite.size)}


def group_design(labels, groups, control_labels=("non-targeting", "control", "")) -> dict:
    """Cells per target, and the batches each target appears in."""
    labels, groups = np.asarray(labels), np.asarray(groups)
    control = np.isin(labels, list(control_labels))
    per_target, batches = {}, {}
    for target in np.unique(labels[~control]):
        rows = labels == target
        per_target[str(target)] = int(rows.sum())
        batches[str(target)] = int(len(np.unique(groups[rows])))
    return {"n_cells": int(len(labels)), "n_control_cells": int(control.sum()),
            "n_targets": len(per_target), "n_gem_groups": int(len(np.unique(groups))),
            "cells_per_target": per_target, "gem_groups_per_target": batches}


def pseudobulk(matrix: np.ndarray, labels) -> dict:
    """Mean over each target's cells, which is what the released bulk should be."""
    labels = np.asarray(labels)
    return {str(target): np.asarray(matrix[labels == target], dtype=np.float64).mean(axis=0)
            for target in np.unique(labels)}


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    denominator = np.linalg.norm(a) * np.linalg.norm(b)
    if denominator == 0:
        return float("nan")
    return float(a @ b / denominator)


def split_half_reliability(matrix: np.ndarray, groups, seed: int, n_splits: int = N_SPLITS) -> float:
    """Median cosine between the two half-means over random halvings of the batches.

    Batches, not cells, are split: a direction that only reproduces within one GEM
    group is not a reliable direction for the target.
    """
    groups = np.asarray(groups)
    unique = np.unique(groups)
    if len(unique) < 2:
        return float("nan")
    rng = np.random.default_rng(seed)
    cosines = []
    for _ in range(n_splits):
        shuffled = rng.permutation(unique)
        left = np.isin(groups, shuffled[: len(shuffled) // 2])
        if left.sum() == 0 or (~left).sum() == 0:
            continue
        first = np.asarray(matrix[left], dtype=np.float64).mean(axis=0)
        second = np.asarray(matrix[~left], dtype=np.float64).mean(axis=0)
        cosines.append(cosine(first, second))
    return float(np.median(cosines)) if cosines else float("nan")


def split_half_reliability_by_unit(matrix: np.ndarray, units, seed: int,
                                   n_splits: int = N_SPLITS) -> float:
    """Split-half reliability where the units are hairpins rather than batches.

    A target with one hairpin measured in several cell lines has no hairpin split
    to make, so its signatures are halved at random instead.
    """
    units = np.asarray(units)
    if len(np.unique(units)) >= 2:
        return split_half_reliability(matrix, units, seed, n_splits)
    n = len(matrix)
    if n < 2:
        return float("nan")
    rng = np.random.default_rng(seed)
    cosines = []
    for _ in range(n_splits):
        order = rng.permutation(n)
        left, right = order[: n // 2], order[n // 2:]
        cosines.append(cosine(np.asarray(matrix[left], dtype=np.float64).mean(axis=0),
                              np.asarray(matrix[right], dtype=np.float64).mean(axis=0)))
    return float(np.median(cosines))


def target_seed(master_seed: int, target: str) -> int:
    """A per-target seed fixed by the master seed and the target's name.

    Deriving it from the name rather than from a position or a batch index is what
    makes an interrupted run and a clean run produce the same numbers: a target
    keeps its seed however the batches happen to fall.
    """
    digest = hashlib.blake2b(f"{master_seed}:{target}".encode(), digest_size=8).digest()
    return int.from_bytes(digest, "big") % (2 ** 32)


def batches(targets, size: int):
    """The work split into named batches, in one deterministic order."""
    ordered = sorted(targets)
    return [(f"batch_{start // size:04d}", ordered[start:start + size])
            for start in range(0, len(ordered), size)]


def merge_batches(written: dict) -> dict:
    """One result per target, refusing duplicates and conflicts.

    `written` maps a batch name to its target dictionary. A target appearing twice
    is a sign that the batching changed between runs, which would make a resumed
    run disagree with a clean one, so assembly stops rather than picking a winner.
    """
    merged, source = {}, {}
    for batch in sorted(written):
        for target, entry in written[batch].items():
            if target in merged:
                agreement = "identical" if merged[target] == entry else "conflicting"
                raise AssertionError(
                    f"{target} appears in {source[target]} and {batch} with {agreement} "
                    "values; the batching changed between runs, so assembly stops")
            merged[target] = entry
            source[target] = batch
    return merged
