"""Single-cell audit helpers for R0.1, R0.3, R0.4 and R0.6 of the frozen registration.

The Modal wrapper (`experiments/modal_03d_single_cell.py`) arranges workers,
volumes and checkpoints; everything that computes a number lives here, so it can
be tested without Modal and without the 19 GB of released single-cell data.
"""
import hashlib

import numpy as np

from .references import require_finite

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


def gene_names(adata) -> np.ndarray:
    """The symbols a released file carries, from `gene_name` where it has one."""
    return (adata.var["gene_name"].astype(str).to_numpy() if "gene_name" in adata.var
            else adata.var_names.astype(str).to_numpy())


def basis_columns(names, symbols) -> np.ndarray:
    """The columns of a released file holding each basis symbol, in basis order.

    A symbol the file does not carry exactly once is not skipped: the basis was
    frozen on this release, so a gap or a duplicate means the file is not the
    file the basis describes.
    """
    names = np.asarray([str(name) for name in names])
    position = {}
    for index, name in enumerate(names):
        position.setdefault(name, []).append(index)
    columns, missing, duplicated = [], [], []
    for symbol in symbols:
        found = position.get(str(symbol), [])
        if len(found) == 1:
            columns.append(found[0])
        elif found:
            duplicated.append(symbol)
        else:
            missing.append(symbol)
    if missing or duplicated:
        raise AssertionError(
            f"the frozen basis does not land on this file: {len(missing)} symbols absent "
            f"{missing[:10]}, {len(duplicated)} carried more than once {duplicated[:10]}")
    return np.asarray(columns, dtype=int)


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    """The cosine of two finite vectors, refusing a non-finite operand."""
    a = require_finite(a, "the first operand of a cosine")
    b = require_finite(b, "the second operand of a cosine")
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


def verification_summary(name, merged: dict) -> dict:
    """R0.4 for one construction: the counts, and a refusal where the median is undefined.

    A median that is not a number is not below a threshold, so it decides nothing.
    The registered rule fires when the median per-target cosine is below 0.99;
    routing an undefined median into that branch labels the construction on
    arithmetic rather than on disagreement, which is what Deviation 13 records.
    Here an incomplete set of cosines is reported as refused, with the counts that
    say how incomplete it is and a reason for every target that has none.
    """
    per_target = {target: entry["cosine_with_released_bulk"]
                  for target, entry in merged.items()
                  if "cosine_with_released_bulk" in entry}
    reasons = {target: entry.get("cosine_reason", entry.get("reason"))
               for target, entry in merged.items()
               if "cosine_with_released_bulk" not in entry}
    cosines = np.asarray(list(per_target.values()), dtype=np.float64)
    finite = cosines[np.isfinite(cosines)]
    complete = bool(cosines.size) and finite.size == cosines.size
    median = float(np.median(finite)) if finite.size else None
    # an undefined cosine is serialized as null beside the target it belongs to,
    # so the reader sees a target with no value rather than a value that is not a
    # number, and the file stays valid JSON for something other than Python
    not_finite = sorted(target for target, value in per_target.items()
                        if not np.isfinite(value))
    per_target = {target: (float(value) if np.isfinite(value) else None)
                  for target, value in per_target.items()}

    def first(key):
        return next((entry[key] for entry in merged.values() if key in entry), None)

    return {
        "construction": name,
        "basis_sha256": first("basis_sha256"),
        "n_landmarks": first("n_landmarks"),
        "n_targets_in_release": len(merged),
        "n_total": int(cosines.size),
        "n_finite": int(finite.size),
        "n_non_finite": int(cosines.size - finite.size),
        "n_not_computed": len(reasons),
        "reason_not_computed": reasons,
        "threshold": 0.99,
        "status": "computed" if complete else "refused",
        "median_cosine": median if complete else None,
        "min_cosine": float(np.min(finite)) if complete and finite.size else None,
        "n_below_threshold": int((finite < 0.99).sum()) if complete else None,
        "describes": None if not complete else
                     ("C1" if median >= 0.99 else "the single-cell construction, not C1"),
        "refusal": None if complete else
                   (f"{int(cosines.size - finite.size)} of {int(cosines.size)} per-target "
                    "cosines are not finite on the frozen basis and "
                    f"{len(reasons)} targets have none, so no median is reported and R0.6 "
                    "and R7c take no label from this module"),
        "per_target_cosine": per_target,
        "per_target_not_finite": not_finite,
    }
