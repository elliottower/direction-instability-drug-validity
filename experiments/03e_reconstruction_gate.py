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

import numpy as np
import pandas as pd

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
OUT = REPO / "results" / "03c_h3_sensitivity"
REFERENCE_ARTIFACT = (REPO / "results" / "03_phenotype_projection"
                      / "phenotype_projection_results.json")
RTOL, ATOL = 1e-5, 1e-5
MIN_HAIRPINS = 3


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


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
            for key in data.files:
                if key == "fingerprint":
                    continue
                if "\x00cells" in key:
                    cells[key.split("\x00")[0]] = [str(c) for c in data[key]]
                else:
                    matrices[key] = data[key]
    gene_file = Path(directory) / "landmark_gene_ids.json"
    if gene_file.exists():
        genes = [str(g) for g in json.loads(gene_file.read_text())]
    assert genes is not None, f"{gene_file} is missing; the rebuild must record its gene axis"
    missing = set(matrices) - set(cells)
    assert not missing, f"{len(missing)} drugs arrived without cell-line identifiers"
    return matrices, cells, genes


def load_extraction(data_dir: Path):
    """The same objects from the extraction the corrected artifact was computed from."""
    compounds = np.load(data_dir / "lincs_subset.npz", allow_pickle=True)
    signatures = compounds["signatures"]
    position = {str(sig_id): i for i, sig_id in enumerate(compounds["sig_ids"])}
    siginfo = pd.read_csv(data_dir / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(position) & siginfo.pert_iname.notna()].copy()
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


def compare_shrna(rebuilt, extraction, stored_directions=None,
                  min_signatures=MIN_HAIRPINS, rtol=RTOL, atol=ATOL):
    """The half of the gate that guards the object Deviation 9 corrupted.

    A perfect compound reconstruction cannot detect a recurrence of that defect,
    because the defect was in building the target directions. This compares the
    shRNA signatures themselves, joined by signature id, then rebuilds each target
    consensus independently on both sides and compares the consensus and the unit
    direction it becomes.

    Each side is (signatures by signature id, target -> its signature ids, genes).

    `stored_directions` is what the rebuild itself wrote out as each target's unit
    direction. Comparing it against the direction recomputed here from verified
    signatures is the check that bites: Deviation 9 was a rebuild whose signatures
    were all correct and whose consensus construction was not.
    """
    rebuilt_signatures, rebuilt_membership, rebuilt_genes = rebuilt
    extraction_signatures, extraction_membership, extraction_genes = extraction

    assert set(rebuilt_signatures) == set(extraction_signatures), (
        f"signature-id sets differ: {len(set(rebuilt_signatures) - set(extraction_signatures))} "
        f"only in the rebuild, "
        f"{len(set(extraction_signatures) - set(rebuilt_signatures))} only in the extraction")
    gene_order = {gene: i for i, gene in enumerate(map(str, extraction_genes))}
    assert set(map(str, rebuilt_genes)) == set(gene_order), "the shRNA gene axes differ"
    permutation = np.array([gene_order[gene] for gene in map(str, rebuilt_genes)])

    worst_signature = 0.0
    for sig_id, vector in rebuilt_signatures.items():
        other = np.asarray(extraction_signatures[sig_id])[permutation]
        difference = np.abs(np.asarray(vector) - other)
        assert (difference <= atol + rtol * np.abs(other)).all(), (
            f"signature {sig_id} differs by {difference.max():.3g}")
        worst_signature = max(worst_signature, float(difference.max()))

    assert set(rebuilt_membership) == set(extraction_membership), (
        f"target sets differ: {sorted(set(rebuilt_membership) ^ set(extraction_membership))[:5]}")

    worst = {"consensus": 0.0, "direction": 0.0, "stored": 0.0}
    per_target, worst_target = {}, None
    for target in sorted(rebuilt_membership):
        left_ids, right_ids = sorted(rebuilt_membership[target]), sorted(extraction_membership[target])
        assert left_ids == right_ids, (
            f"{target}: the two sides group different signatures into it")
        assert len(left_ids) >= min_signatures, (
            f"{target}: {len(left_ids)} signatures, below the eligibility floor")
        left = np.vstack([rebuilt_signatures[i] for i in left_ids]).mean(axis=0)
        right = np.vstack([np.asarray(extraction_signatures[i])[permutation]
                           for i in right_ids]).mean(axis=0)
        consensus_error = float(np.abs(left - right).max())
        # a normalized direction can agree while the consensus behind it does not,
        # so both are compared
        direction_error = float(np.abs(left / np.linalg.norm(left)
                                       - right / np.linalg.norm(right)).max())
        entry = {"n_signatures": len(left_ids),
                 "max_abs_consensus_difference": consensus_error,
                 "max_abs_direction_difference": direction_error}
        if stored_directions is not None:
            assert target in stored_directions, f"{target}: the rebuild stored no direction"
            stored = np.asarray(stored_directions[target])[permutation]
            recomputed = left / np.linalg.norm(left)
            stored_error = float(np.abs(stored - recomputed).max())
            entry["max_abs_stored_direction_difference"] = stored_error
            worst["stored"] = max(worst["stored"], stored_error)
        per_target[target] = entry
        if consensus_error > worst["consensus"]:
            worst_target = target
        worst["consensus"] = max(worst["consensus"], consensus_error)
        worst["direction"] = max(worst["direction"], direction_error)

    tolerance = atol + rtol * 1.0
    return {"n_signatures": len(rebuilt_signatures), "n_targets": len(per_target),
            "max_abs_signature_difference": worst_signature,
            "max_abs_consensus_difference": worst["consensus"],
            "max_abs_direction_difference": worst["direction"],
            "max_abs_stored_direction_difference": (worst["stored"] if stored_directions
                                                    is not None else None),
            "stored_directions_compared": stored_directions is not None,
            "worst_target": worst_target,
            "all_within_tolerance": bool(worst["consensus"] <= tolerance
                                         and worst["direction"] <= tolerance
                                         and worst["stored"] <= tolerance),
            "rtol": rtol, "atol": atol, "per_target": per_target}


def load_shrna_rebuild(directory: Path):
    """Independently rebuilt shRNA signatures, their target grouping, and the gene axis."""
    path = Path(directory) / "shrna_signatures.npz"
    assert path.exists(), f"{path} is missing; the rebuild must emit shRNA signatures with ids"
    with np.load(path, allow_pickle=True) as data:
        ids = [str(i) for i in data["sig_ids"]]
        signatures = {sig_id: data["signatures"][i] for i, sig_id in enumerate(ids)}
        genes = [str(g) for g in data["gene_ids"]]
        membership = {str(target): [str(i) for i in members]
                      for target, members in json.loads(str(data["membership"])).items()}
    assert len(set(ids)) == len(ids), "the rebuild repeats a signature id"

    stored = None
    consensus_path = Path(directory) / "shrna_consensus.npz"
    if consensus_path.exists():
        with np.load(consensus_path, allow_pickle=True) as data:
            stored = {str(gene): data["directions"][i]
                      for i, gene in enumerate(data["genes"])}
    return signatures, membership, genes, stored


def load_shrna_extraction(data_dir: Path, targets_wanted=None):
    """The same objects from the pinned extraction, grouped by the pinned metadata rule."""
    shrna = np.load(data_dir / "lincs_shrna.npz", allow_pickle=True)
    ids = [str(sig_id) for sig_id in shrna["sig_ids"]]
    assert len(set(ids)) == len(ids), "the extraction repeats a signature id"
    signatures = {sig_id: shrna["signatures"][i] for i, sig_id in enumerate(ids)}
    siginfo = pd.read_csv(data_dir / "lincs_shrna_siginfo.csv.gz")
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(signatures)]

    membership = {}
    for gene, group in siginfo.groupby("pert_iname"):
        members = sorted(set(group.sig_id.astype(str)))
        if len(members) < MIN_HAIRPINS:
            continue
        if targets_wanted is None or gene in targets_wanted:
            membership[str(gene)] = members
    return signatures, membership, [str(g) for g in shrna["gene_ids"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rebuilt", type=Path, required=True,
                        help="directory holding the GCTX rebuild's shards")
    parser.add_argument("--extraction", type=Path, required=True)
    parser.add_argument("--cohort", type=Path, default=REFERENCE_ARTIFACT,
                        help="the artifact whose drugs and targets both sides must hold")
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()

    records = json.loads(args.cohort.read_text())
    drugs = sorted({record["drug"] for record in records})
    targets = sorted({record["target"] for record in records})

    compounds = compare(load_rebuild(args.rebuilt), load_extraction(args.extraction), cohort=drugs)
    rebuilt_signatures, rebuilt_membership, rebuilt_genes, stored = load_shrna_rebuild(args.rebuilt)
    shrna = compare_shrna((rebuilt_signatures, rebuilt_membership, rebuilt_genes),
                          load_shrna_extraction(args.extraction, targets_wanted=set(targets)),
                          stored_directions=stored)

    result = {"cohort": {"file": str(args.cohort), "sha256": sha256_file(args.cohort),
                         "n_drugs": len(drugs), "n_targets": len(targets)},
              "compound_signatures": compounds, "shrna_signatures_and_consensuses": shrna}
    result["identifier_hashes"] = {"drugs": sha256_text("\n".join(drugs)),
                                   "targets": sha256_text("\n".join(targets))}
    passed = compounds["all_within_tolerance"] and shrna["all_within_tolerance"]
    result["gate"] = ("reconstruction" if passed
                      else "failed: the analyses registered against this gate are void")

    args.output.mkdir(parents=True, exist_ok=True)
    path = args.output / "reconstruction_gate.json"
    path.write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("compound_signatures", "shrna_signatures_and_consensuses")}
                     | {"compound_max_abs": compounds["max_abs_difference"],
                        "shrna_max_abs_consensus": shrna["max_abs_consensus_difference"]}, indent=2))
    assert passed, (
        f"the rebuild departs from the extraction: compounds by "
        f"{compounds['max_abs_difference']:.3g}, shRNA consensuses by "
        f"{shrna['max_abs_consensus_difference']:.3g}")


if __name__ == "__main__":
    main()
