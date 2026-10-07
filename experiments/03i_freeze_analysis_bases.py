#!/usr/bin/env python3
"""Freeze the analysis basis of every construction that reads an official Replogle file.

Amendment 4 registers the basis rule: a construction's basis is the ordered
intersection of the frozen LINCS landmark axis with the genes its source
represents exactly once and defines with finite values throughout. This script
measures that basis per source and writes it with its hash, so the loaders consume
a pinned list rather than recomputing one at run time, and so a coordinate that
leaves an analysis is named in a committed file before any statistic is computed.

It names no gene. Which coordinates are excluded is read off the released files.

    python experiments/03i_freeze_analysis_bases.py \
        --replogle data/replogle2022 \
        --gene-info data/lincs_extraction/GSE92742_Broad_LINCS_gene_info.txt.gz \
        --output registry/frozen/analysis_bases.json
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np

from geometry.references import (N_LANDMARK, Reference, basis_sha256, finite_landmark_basis,
                                 frozen_landmark_order, landmark_symbols, shared_space,
                                 sha256_file)

# one source file per construction, as the registration frozen at 7f57136 names them
SOURCES = {
    "C1-K562": "K562_essential_normalized_bulk_01.h5ad",
    "C1-RPE1": "rpe1_normalized_bulk_01.h5ad",
    "C1-GW": "K562_gwps_normalized_bulk_01.h5ad",
}

# R4's phenotype-positive arm is a row filter on the genome-wide file, so it reads
# the same source and takes the same basis. Deriving a basis from the qualifying
# rows alone would let a phenotype decide a coordinate, which is what the rule
# exists to prevent.
SHARES_A_BASIS = {"C1-GW-phenotype-positive": "C1-GW"}

# the comparisons the registration frozen at 7f57136 names, each on the landmarks
# its two constructions share. R0.5: C0 and C1; C1 and C1-RPE1; C1 and C1-GW.
# R6a: C1 and C1-RPE1, on "landmarks present in both files".
REGISTERED_COMPARISONS = [
    ("R0.5", "C0", "C1-K562"),
    ("R0.5", "C1-K562", "C1-RPE1"),
    ("R0.5", "C1-K562", "C1-GW"),
    ("R6a", "C1-K562", "C1-RPE1"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replogle", type=Path, required=True,
                        help="the directory holding the released pseudobulk files")
    parser.add_argument("--gene-info", type=Path, required=True,
                        help="the pinned LINCS gene-info file")
    parser.add_argument("--perturbseq", type=Path, required=True,
                        help="the scPerturb file C0 is built from, read for its gene axis "
                             "alone so C0's declared and measured spaces are both recorded")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    gene_ids = frozen_landmark_order(args.gene_info)
    symbols = landmark_symbols(args.gene_info, gene_ids)

    constructions = {}
    for name, filename in SOURCES.items():
        source = args.replogle / filename
        if not source.exists():
            raise SystemExit(f"{name}: {source} is absent, and a basis is not guessed")
        constructions[name] = finite_landmark_basis(source, symbols, gene_ids, name)
        record = constructions[name]
        print(f"{name:<8} {record['n_landmarks_matched']:>4} matched -> "
              f"{record['n_landmarks_in_basis']:>4} in basis   "
              f"{record['basis_sha256'][:16]}  "
              f"excluded: {[g['symbol'] for g in record['excluded_as_undefined']] or 'none'}")

    # C0 places its difference in the full 978 with zeros where the file measures
    # nothing, so its declared space is every landmark whatever the file holds.
    # The measured count is recorded beside it rather than inferred.
    perturbseq = ad.read_h5ad(args.perturbseq, backed="r")
    measured = sorted({str(gene) for gene in perturbseq.var_names} & set(symbols))
    perturbseq.file.close()
    measured_positions = [index for index, symbol in enumerate(symbols)
                          if str(symbol) in set(measured)]
    measured_ids = [str(gene_ids[index]) for index in measured_positions]
    pooled = {
        "construction": "C0",
        "file": args.perturbseq.name,
        "file_sha256": sha256_file(args.perturbseq),
        # declared and measured are different things and are named apart: C0's
        # declared space is the whole axis whatever the file holds, and a count of
        # 978 beside a count of 728 invites reading one as the other
        "n_landmarks_declared": N_LANDMARK,
        "n_landmarks_measured": len(measured_positions),
        "landmark_positions": list(range(N_LANDMARK)),
        "landmark_gene_ids": [str(gene) for gene in gene_ids],
        "landmark_symbols": [str(symbol) for symbol in symbols],
        "declared_sha256": basis_sha256(gene_ids),
        "measured_positions": measured_positions,
        "measured_gene_ids": measured_ids,
        "measured_sha256": basis_sha256(measured_ids),
        "why": "the pooled construction places its difference in the full 978 landmarks "
               "with zeros where the file measures nothing, so its declared space is the "
               "whole axis and no coordinate of it is undefined",
    }
    print(f"C0       {N_LANDMARK:>4} declared, {len(measured_positions):>4} measured in "
          f"{args.perturbseq.name}, measured basis {pooled['measured_sha256'][:16]}")

    for name, source_of in SHARES_A_BASIS.items():
        constructions[name] = {**constructions[source_of],
                               "construction": name,
                               "basis_taken_from": source_of,
                               "why": "a row filter on the same source does not change which "
                                      "coordinates that source defines"}

    # the intersections the registered comparisons run on, from the production
    # helper rather than from a count entered by hand
    spaces = {name: Reference(name, {}, np.asarray(record["landmark_positions"], dtype=int))
              for name, record in {**constructions, "C0": pooled}.items()}
    records = {**constructions, "C0": pooled}
    comparisons = {}
    for module, left, right in REGISTERED_COMPARISONS:
        key = f"{module}: {left} against {right}"
        unfiltered = {side: set(records[side]["landmark_positions"])
                      | {gene["landmark_position"]
                         for gene in records[side].get("excluded_as_undefined", [])}
                      for side in (left, right)}
        comparisons[key] = {
            "module": module, "left": left, "right": right,
            "n_left": len(spaces[left].positions), "n_right": len(spaces[right].positions),
            "n_shared": int(len(shared_space(spaces[left], spaces[right]))),
            "n_shared_if_nothing_were_excluded": len(unfiltered[left] & unfiltered[right]),
        }
        row = comparisons[key]
        row["coordinates_lost_to_the_rule"] = (row["n_shared_if_nothing_were_excluded"]
                                               - row["n_shared"])
        print(f"{key:<36} {row['n_shared']:>4} shared, "
              f"{row['coordinates_lost_to_the_rule']} lost to the rule")

    out = {
        "what_this_is":
            "the frozen analysis basis of each construction built from an official Replogle "
            "pseudobulk file: the ordered intersection of the frozen 978-landmark axis with "
            "the genes the source represents exactly once and defines with finite values. "
            "A coordinate is excluded by the state of the source, never by a quantity "
            "computed from it.",
        "schema_version": 1,
        "registered_by": "experiments/PREREG_H3_S1S3_AMENDMENT_4.md",
        "generated_by": "experiments/03i_freeze_analysis_bases.py",
        "gene_info": {"file": args.gene_info.name,
                      "n_landmarks": len(gene_ids),
                      "landmark_axis_sha256": basis_sha256(gene_ids)},
        "constructions": constructions,
        "pooled_construction": pooled,
        "shares_a_basis": SHARES_A_BASIS,
        "registered_comparisons": comparisons,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True))

    # the ordered axis, beside the bases, so a consumer that holds no gene-info
    # file can still check that a basis sits on the axis it claims
    axis = args.output.parent / "landmark_gene_ids.json"
    axis.write_text(json.dumps([str(gene) for gene in gene_ids], indent=2))
    print(f"\nwrote {args.output}")
    print(f"wrote {axis}")


if __name__ == "__main__":
    main()
