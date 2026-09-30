"""The reconstruction gate of the corrected-baseline registration (frozen 7f57136).

An independent rebuild from the GEO GCTX is compared with the extraction the
corrected artifact was computed from. The comparison joins on identifiers, never
on row or column order: the defect of Deviation 9 was an identifier-versus-position
error, and a gate that compares by position could not see it.

Signature matrices are joined on (drug, cell line) and genes on Entrez id. Both
identifier sequences are hashed separately from the numeric arrays. The tolerance
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


def compare(rebuilt, extraction, rtol=RTOL, atol=ATOL):
    """Per-drug comparison after aligning rows by cell line and columns by gene id."""
    rebuilt_matrices, rebuilt_cells, rebuilt_genes = rebuilt
    extraction_matrices, extraction_cells, extraction_genes = extraction
    assert set(rebuilt_matrices) == set(extraction_matrices), (
        "the two sides hold different drugs")

    gene_order = {gene: i for i, gene in enumerate(map(str, extraction_genes))}
    assert set(map(str, rebuilt_genes)) == set(gene_order), "the gene axes differ"
    gene_permutation = np.array([gene_order[gene] for gene in map(str, rebuilt_genes)])

    worst, per_drug = 0.0, {}
    for drug in sorted(rebuilt_matrices):
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
            "rtol": rtol, "atol": atol, "per_drug": per_drug}


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


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rebuilt", type=Path, required=True,
                        help="directory holding the GCTX rebuild's shards")
    parser.add_argument("--extraction", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()

    rebuilt = load_rebuild(args.rebuilt)
    extraction = load_extraction(args.extraction)
    common = sorted(set(rebuilt[0]) & set(extraction[0]))
    rebuilt = ({d: rebuilt[0][d] for d in common}, {d: rebuilt[1][d] for d in common}, rebuilt[2])
    extraction = ({d: extraction[0][d] for d in common}, {d: extraction[1][d] for d in common},
                  extraction[2])

    result = compare(rebuilt, extraction)
    result["identifier_hashes"] = {
        "drugs": sha256_text("\n".join(common)),
        "rebuilt_genes": sha256_text("\n".join(rebuilt[2])),
        "extraction_genes": sha256_text("\n".join(extraction[2]))}
    result["gate"] = ("reconstruction" if result["all_within_tolerance"]
                      else "failed: the analyses registered against this gate are void")

    args.output.mkdir(parents=True, exist_ok=True)
    path = args.output / "reconstruction_gate.json"
    path.write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items() if k != "per_drug"}, indent=2))
    assert result["all_within_tolerance"], (
        f"the rebuild departs from the extraction by {result['max_abs_difference']:.3g}")


if __name__ == "__main__":
    main()
