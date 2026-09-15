import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
build = __import__("03_phenotype_projection").build_shrna_consensus


def _fixture(n_genes=6):
    # three hairpins per gene; row i of the matrix carries the value i
    sig_ids = [f"S{i}" for i in range(9)]
    sigs = np.array([[float(i)] * n_genes for i in range(9)])
    info = pd.DataFrame({"sig_id": sig_ids,
                         "pert_iname": ["GENE_A"] * 3 + ["GENE_B"] * 3 + ["GENE_C"] * 3})
    return sigs, sig_ids, info


def test_consensus_is_invariant_to_metadata_row_order():
    sigs, sig_ids, info = _fixture()
    base = build(sigs, sig_ids, info)
    for shuffle in (info.iloc[::-1], info.sample(frac=1.0, random_state=None)):
        other = build(sigs, sig_ids, shuffle.reset_index(drop=True))
        for gene in base:
            assert other[gene] == pytest.approx(base[gene], abs=1e-12)


def test_consensus_follows_the_matrix_when_its_row_order_changes():
    sigs, sig_ids, info = _fixture()
    base = build(sigs, sig_ids, info)
    perm = np.array([8, 7, 6, 5, 4, 3, 2, 1, 0])
    other = build(sigs[perm], [sig_ids[i] for i in perm], info)
    for gene in base:
        assert other[gene] == pytest.approx(base[gene], abs=1e-12)


def test_gene_consensus_averages_only_its_own_hairpins():
    sigs, sig_ids, info = _fixture()
    c = build(sigs, sig_ids, info)
    assert c["GENE_A"] == pytest.approx([1.0] * 6, abs=1e-12)   # rows 0,1,2
    assert c["GENE_B"] == pytest.approx([4.0] * 6, abs=1e-12)   # rows 3,4,5
    assert c["GENE_C"] == pytest.approx([7.0] * 6, abs=1e-12)   # rows 6,7,8


def test_metadata_position_indexing_would_give_a_different_answer():
    sigs, sig_ids, info = _fixture()
    correct = build(sigs, sig_ids, info)
    scrambled_info = info.iloc[[8, 7, 6, 5, 4, 3, 2, 1, 0]].reset_index(drop=True)
    by_position = {}
    for gene, grp in scrambled_info.groupby("pert_iname"):
        by_position[gene] = sigs[grp.index.to_numpy()].mean(axis=0)
    assert by_position["GENE_A"] != pytest.approx(correct["GENE_A"], abs=1e-9)


def test_duplicate_metadata_rows_do_not_double_weight_a_hairpin():
    sigs, sig_ids, info = _fixture()
    dup = pd.concat([info, info.iloc[[0]]], ignore_index=True)
    assert build(sigs, sig_ids, dup)["GENE_A"] == pytest.approx(
        build(sigs, sig_ids, info)["GENE_A"], abs=1e-12)


def test_row_count_mismatch_fails_closed():
    sigs, sig_ids, info = _fixture()
    with pytest.raises(AssertionError, match="signature rows"):
        build(sigs[:-1], sig_ids, info)
