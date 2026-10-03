import anndata as ad
import numpy as np
import pandas as pd
import pytest

from geometry.references import (Reference, alignment_matrix, check_declared_axis,
                                 landmark_symbols, load_replogle_bulk,
                                 pooled_crispri_reference, shared_space, unit)


def _gene_info(tmp_path, symbols, extra_non_landmark=("ZZZ1",)):
    rows = [{"pr_gene_id": 1000 + i, "pr_gene_symbol": s, "pr_is_lm": 1}
            for i, s in enumerate(symbols)]
    rows += [{"pr_gene_id": 9000 + i, "pr_gene_symbol": s, "pr_is_lm": 0}
             for i, s in enumerate(extra_non_landmark)]
    path = tmp_path / "gene_info.txt.gz"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression="gzip")
    return path, [str(1000 + i) for i in range(len(symbols))]


def _bulk(tmp_path, perturbed_genes, measured_genes, values=None, rows_per_gene=1,
          energy_p=None, name="bulk.h5ad"):
    rng = np.random.default_rng()
    obs_names, obs_rows = [], []
    for gene in perturbed_genes:
        for r in range(rows_per_gene):
            obs_names.append(f"{r}_{gene}_P1P2_ENSG{r:08d}")
            obs_rows.append({"energy_test_p_value": 1.0 if energy_p is None else energy_p[gene][r]})
    x = rng.standard_normal((len(obs_names), len(measured_genes))) if values is None else values
    adata = ad.AnnData(X=np.asarray(x, dtype=np.float32),
                       obs=pd.DataFrame(obs_rows, index=obs_names),
                       var=pd.DataFrame({"gene_name": list(measured_genes)},
                                        index=list(measured_genes)))
    path = tmp_path / name
    adata.write_h5ad(path)
    return path


def test_landmark_symbols_follow_the_extraction_order(tmp_path):
    symbols = ["AAA", "BBB", "CCC"]
    path, gene_ids = _gene_info(tmp_path, symbols)

    assert landmark_symbols(path, gene_ids) == symbols
    assert landmark_symbols(path, list(reversed(gene_ids))) == list(reversed(symbols))
    # a gene id the file does not call a landmark keeps its id, rather than
    # silently borrowing a neighbour's symbol
    assert landmark_symbols(path, ["9000"]) == ["9000"]


def test_shared_space_drops_symbols_duplicated_on_either_side(tmp_path):
    symbols = ["AAA", "BBB", "CCC", "BBB"]          # BBB duplicated among landmarks
    _, gene_ids = _gene_info(tmp_path, symbols)
    path = _bulk(tmp_path, ["T1", "T2"], ["AAA", "BBB", "CCC", "DDD"])
    reference = load_replogle_bulk(path, symbols, "test")

    kept = [symbols[p] for p in reference.positions]
    assert kept == ["AAA", "CCC"]
    assert reference.audit["n_duplicated_symbols"] == 1
    assert reference.audit["n_unmatched_landmarks"] == 0


def test_rows_are_averaged_before_normalization(tmp_path):
    symbols = ["G1", "G2", "G3"]
    values = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])       # two rows of one gene
    path = _bulk(tmp_path, ["T1"], symbols, values=values, rows_per_gene=2)
    reference = load_replogle_bulk(path, symbols, "test")

    expected = unit(np.array([1.0, 1.0, 0.0]))
    assert reference.directions["T1"] == pytest.approx(expected, abs=1e-6)
    assert reference.audit["rows_per_gene"]["T1"] == 2
    assert reference.audit["rows_used_per_gene"]["T1"] == 2


def test_qualifying_rows_build_the_direction_from_those_rows_alone(tmp_path):
    symbols = ["G1", "G2", "G3"]
    values = np.array([[3.0, 0.0, 0.0],      # T1 qualifies
                       [0.0, 0.0, 3.0]])     # T2 does not
    path = _bulk(tmp_path, ["T1", "T2"], symbols, values=values, rows_per_gene=1,
                 energy_p={"T1": [0.01], "T2": [0.9]})
    # T1 has one row at p = 0.01; give it a second, nonqualifying row
    values2 = np.array([[3.0, 0.0, 0.0], [0.0, 3.0, 0.0]])
    path2 = _bulk(tmp_path, ["T1"], symbols, values=values2, rows_per_gene=2,
                  energy_p={"T1": [0.01, 0.9]}, name="bulk2.h5ad")

    eligible = load_replogle_bulk(path2, symbols, "test",
                                  qualifying_rows=lambda obs: obs["energy_test_p_value"] < 0.05)
    assert eligible.directions["T1"] == pytest.approx(unit(np.array([1.0, 0.0, 0.0])), abs=1e-6)
    assert eligible.audit["rows_per_gene"]["T1"] == 2
    assert eligible.audit["rows_used_per_gene"]["T1"] == 1

    all_rows = load_replogle_bulk(path2, symbols, "test")
    assert all_rows.directions["T1"] == pytest.approx(unit(np.array([1.0, 1.0, 0.0])), abs=1e-6)

    only_failing = load_replogle_bulk(path, symbols, "test",
                                      qualifying_rows=lambda obs: obs["energy_test_p_value"] < 0.05)
    assert "T2" not in only_failing.directions
    assert only_failing.audit["n_genes_dropped"] == 1


def test_restriction_renormalizes_on_the_narrower_space():
    positions = np.array([0, 1, 2, 3])
    direction = unit(np.array([1.0, 1.0, 1.0, 5.0]))
    reference = Reference("r", {"T": direction}, positions)
    narrowed = reference.restricted_to(np.array([0, 1, 2]))

    assert np.linalg.norm(narrowed.directions["T"]) == pytest.approx(1.0)
    assert narrowed.directions["T"] == pytest.approx(unit(np.ones(3)))


def test_dropping_the_target_gene_is_a_no_op_when_it_is_absent():
    symbols = ["A", "B", "C"]
    reference = Reference("r", {"TP53": unit(np.array([1.0, 2.0, 3.0]))}, np.arange(3))
    assert reference.without_gene("TP53", symbols) is reference

    with_gene = Reference("r", {"B": unit(np.array([1.0, 2.0, 3.0]))}, np.arange(3))
    cut = with_gene.without_gene("B", symbols)
    assert list(cut.positions) == [0, 2]
    assert cut.directions["B"] == pytest.approx(unit(np.array([1.0, 3.0])))


def test_pooled_construction_zero_fills_landmarks_the_file_does_not_measure(tmp_path):
    symbols = [f"G{i}" for i in range(12)]
    measured = symbols[:5]
    rng = np.random.default_rng()
    labels = ["non-targeting"] * 40 + ["T1"] * 20
    x = rng.standard_normal((60, len(measured)))
    x[40:] += 5.0                                  # a real perturbation effect
    adata = ad.AnnData(X=x.astype(np.float32),
                       obs=pd.DataFrame({"gene": labels}, index=[f"c{i}" for i in range(60)]),
                       var=pd.DataFrame(index=measured))
    path = tmp_path / "pooled.h5ad"
    adata.write_h5ad(path)

    reference = pooled_crispri_reference(path, symbols)
    direction = reference.directions["T1"]
    assert len(direction) == 978 or len(direction) == len(reference.positions)
    assert np.allclose(direction[5:], 0.0)
    assert not np.allclose(direction[:5], 0.0)
    assert reference.audit["n_shared_landmarks"] == 5


def test_pooled_construction_skips_targets_with_too_few_cells(tmp_path):
    symbols = [f"G{i}" for i in range(4)]
    labels = ["non-targeting"] * 30 + ["Rare"] * 3 + ["Common"] * 25
    rng = np.random.default_rng()
    adata = ad.AnnData(X=rng.standard_normal((58, 4)).astype(np.float32),
                       obs=pd.DataFrame({"gene": labels}, index=[f"c{i}" for i in range(58)]),
                       var=pd.DataFrame(index=symbols))
    path = tmp_path / "pooled2.h5ad"
    adata.write_h5ad(path)

    reference = pooled_crispri_reference(path, symbols)
    assert "Rare" not in reference.directions
    assert "Common" in reference.directions


def test_alignment_matrix_is_the_squared_cosine_and_ignores_magnitude():
    rng = np.random.default_rng()
    positions = np.arange(6)
    genes = ["T1", "T2"]
    reference = Reference("r", {g: unit(rng.standard_normal(6)) for g in genes}, positions)
    means = rng.standard_normal((9, 6))

    matrix = alignment_matrix(means, reference, genes)
    for i in range(9):
        for j, gene in enumerate(genes):
            cos = (means[i] @ reference.directions[gene]) / np.linalg.norm(means[i])
            assert matrix[i, j] == pytest.approx(cos ** 2)
    scaled = alignment_matrix(means * 17.0, reference, genes)
    assert scaled == pytest.approx(matrix)


def test_shared_space_is_the_intersection():
    a = Reference("a", {}, np.array([0, 1, 2, 5]))
    b = Reference("b", {}, np.array([1, 2, 3, 5]))
    assert list(shared_space(a, b)) == [1, 2, 5]


def test_the_declared_axis_check_passes_the_frozen_order_and_refuses_a_reordering():
    frozen = [str(100 + i) for i in range(12)]
    check_declared_axis(list(frozen), frozen, "artifact.npz")

    reordered = list(reversed(frozen))
    with pytest.raises(AssertionError, match="not the frozen landmark order"):
        check_declared_axis(reordered, frozen, "artifact.npz")

    with pytest.raises(AssertionError, match="declares 11 genes against the frozen order's 12"):
        check_declared_axis(frozen[:-1], frozen, "artifact.npz")

    with pytest.raises(AssertionError, match="repeats a gene identifier"):
        check_declared_axis(frozen[:-1] + [frozen[0]], frozen, "artifact.npz")


def test_the_declared_axis_check_cannot_see_a_matrix_that_lies_about_its_labels():
    # Deviation 11's shape, and the reason Amendment 2 registers a second layer:
    # the labels are right, the columns under them are not, and this check passes
    rng = np.random.default_rng()
    frozen = [str(100 + i) for i in range(12)]
    matrix = rng.standard_normal((5, 12))
    mislabeled = matrix[:, rng.permutation(12)]

    check_declared_axis(list(frozen), frozen, "truthful.npz")
    check_declared_axis(list(frozen), frozen, "mislabeled.npz")
    assert not np.allclose(matrix, mislabeled)
