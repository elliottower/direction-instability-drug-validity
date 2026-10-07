"""Amendment 4's finite-coordinate rule, on sources built to exercise it.

The rule names no gene. These build releases whose undefined coordinate is a
different gene each time, so a test passes because the rule reads the source and
not because a symbol was special-cased.
"""
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from geometry.references import (MIN_LANDMARKS, basis_sha256, finite_landmark_basis,
                                 load_replogle_bulk, require_finite, shared_space, unit)
from geometry.single_cell import basis_columns, cosine, verification_summary

SYMBOLS = [f"LM{i:03d}" for i in range(12)]
GENE_IDS = [str(100 + i) for i in range(12)]


def a_release(path, carried, n_rows=9, undefined=(), extra=()):
    """A pseudobulk release carrying `carried` symbols, with `undefined` set to +inf."""
    names = list(carried) + list(extra)
    values = np.arange(1, n_rows * len(names) + 1, dtype=np.float64)
    matrix = values.reshape(n_rows, len(names))
    for symbol in undefined:
        matrix[n_rows // 2, names.index(symbol)] = np.inf
    obs = pd.DataFrame(index=[f"row{i}_GENE{i % 3}_rep" for i in range(n_rows)])
    var = pd.DataFrame({"gene_name": names}, index=[f"ENSG{i:05d}" for i in range(len(names))])
    ad.AnnData(X=matrix, obs=obs, var=var).write_h5ad(path)
    return path


def test_a_coordinate_the_source_leaves_undefined_leaves_the_basis(tmp_path):
    source = a_release(tmp_path / "clean.h5ad", SYMBOLS, undefined=["LM004"])
    record = finite_landmark_basis(source, SYMBOLS, GENE_IDS, "C1-test")

    assert record["n_landmarks_matched"] == len(SYMBOLS)
    assert record["n_landmarks_in_basis"] == len(SYMBOLS) - 1
    assert [gene["symbol"] for gene in record["excluded_as_undefined"]] == ["LM004"]
    assert "LM004" not in record["landmark_symbols"]
    assert record["excluded_as_undefined"][0]["n_posinf"] == 1
    assert record["excluded_as_undefined"][0]["landmark_gene_id"] == GENE_IDS[4]


def test_the_rule_follows_the_source_rather_than_a_named_gene(tmp_path):
    first = finite_landmark_basis(
        a_release(tmp_path / "a.h5ad", SYMBOLS, undefined=["LM004"]),
        SYMBOLS, GENE_IDS, "A")
    second = finite_landmark_basis(
        a_release(tmp_path / "b.h5ad", SYMBOLS, undefined=["LM009", "LM001"]),
        SYMBOLS, GENE_IDS, "B")

    assert [gene["symbol"] for gene in second["excluded_as_undefined"]] == ["LM001", "LM009"]
    assert "LM004" in second["landmark_symbols"]
    assert "LM001" in first["landmark_symbols"] and "LM009" in first["landmark_symbols"]


def test_a_gene_outside_the_landmark_axis_never_enters_the_basis(tmp_path):
    source = a_release(tmp_path / "extra.h5ad", SYMBOLS, undefined=["NOTALANDMARK"],
                       extra=["NOTALANDMARK"])
    record = finite_landmark_basis(source, SYMBOLS, GENE_IDS, "C1-test")

    assert record["n_landmarks_in_basis"] == len(SYMBOLS)
    assert record["excluded_as_undefined"] == []
    assert "NOTALANDMARK" not in record["landmark_symbols"]


def test_one_source_excluding_a_coordinate_leaves_the_others_holding_it(tmp_path):
    affected = finite_landmark_basis(
        a_release(tmp_path / "affected.h5ad", SYMBOLS, undefined=["LM004"]),
        SYMBOLS, GENE_IDS, "affected")
    unaffected = finite_landmark_basis(
        a_release(tmp_path / "unaffected.h5ad", SYMBOLS),
        SYMBOLS, GENE_IDS, "unaffected")

    assert 4 not in affected["landmark_positions"]
    assert 4 in unaffected["landmark_positions"]
    assert affected["basis_sha256"] != unaffected["basis_sha256"]


def test_a_comparison_between_two_bases_is_their_intersection(tmp_path):
    affected = a_release(tmp_path / "affected.h5ad", SYMBOLS, undefined=["LM004"])
    unaffected = a_release(tmp_path / "unaffected.h5ad", SYMBOLS)
    narrow = finite_landmark_basis(affected, SYMBOLS, GENE_IDS, "narrow")
    wide = finite_landmark_basis(unaffected, SYMBOLS, GENE_IDS, "wide")

    left = load_replogle_bulk(affected, SYMBOLS, "narrow", basis=narrow)
    right = load_replogle_bulk(unaffected, SYMBOLS, "wide", basis=wide)
    shared = shared_space(left, right)

    assert 4 not in shared.tolist()
    assert sorted(shared.tolist()) == sorted(narrow["landmark_positions"])
    assert len(right.positions) == len(shared) + 1


def test_every_direction_built_on_a_basis_is_finite_with_a_nonzero_norm(tmp_path):
    source = a_release(tmp_path / "affected.h5ad", SYMBOLS, undefined=["LM004"])
    record = finite_landmark_basis(source, SYMBOLS, GENE_IDS, "C1-test")
    reference = load_replogle_bulk(source, SYMBOLS, "C1-test",
                                   basis=record)

    assert reference.directions
    for gene, direction in reference.directions.items():
        assert np.isfinite(direction).all(), gene
        assert np.linalg.norm(direction) == pytest.approx(1.0)


def test_an_undefined_coordinate_stops_a_load_that_does_not_exclude_it(tmp_path):
    source = a_release(tmp_path / "affected.h5ad", SYMBOLS, undefined=["LM004"])
    with pytest.raises(AssertionError, match="not finite"):
        load_replogle_bulk(source, SYMBOLS, "C1-test")


def test_an_undefined_coordinate_inside_a_basis_is_refused_not_excluded_again(tmp_path):
    source = a_release(tmp_path / "affected.h5ad", SYMBOLS, undefined=["LM004"])
    record = finite_landmark_basis(source, SYMBOLS, GENE_IDS, "C1-test")
    widened = {**record, "landmark_positions": sorted(record["landmark_positions"] + [4])}

    # the basis describes this file and holds a coordinate the file leaves
    # undefined, which is the one case the rule must not silently narrow
    with pytest.raises(AssertionError, match="not finite"):
        load_replogle_bulk(source, SYMBOLS, "C1-test", basis=widened)


def test_a_basis_frozen_on_another_file_is_refused(tmp_path):
    frozen = finite_landmark_basis(a_release(tmp_path / "clean.h5ad", SYMBOLS),
                                   SYMBOLS, GENE_IDS, "C1-test")
    later = a_release(tmp_path / "later.h5ad", SYMBOLS, undefined=["LM007"])

    assert 7 in frozen["landmark_positions"]
    with pytest.raises(AssertionError, match="does not describe this source"):
        load_replogle_bulk(later, SYMBOLS, "C1-test", basis=frozen)


def test_a_source_lacking_a_basis_coordinate_is_refused(tmp_path):
    narrower = a_release(tmp_path / "narrower.h5ad", [s for s in SYMBOLS if s != "LM002"])
    record = finite_landmark_basis(narrower, SYMBOLS, GENE_IDS, "C1-test")
    widened = {**record, "landmark_positions": sorted(record["landmark_positions"] + [2])}

    with pytest.raises(AssertionError, match="absent"):
        load_replogle_bulk(narrower, SYMBOLS, "C1-test", basis=widened)


def test_the_basis_hash_is_over_the_ordered_gene_ids(tmp_path):
    record = finite_landmark_basis(
        a_release(tmp_path / "affected.h5ad", SYMBOLS, undefined=["LM004"]),
        SYMBOLS, GENE_IDS, "C1-test")

    assert record["basis_sha256"] == basis_sha256(record["landmark_gene_ids"])
    assert record["basis_sha256"] != basis_sha256(list(reversed(record["landmark_gene_ids"])))
    assert record["basis_sha256"] != basis_sha256(GENE_IDS)


def test_the_registered_floor_is_reported_against_the_basis_not_the_match(tmp_path):
    record = finite_landmark_basis(
        a_release(tmp_path / "affected.h5ad", SYMBOLS, undefined=["LM004"]),
        SYMBOLS, GENE_IDS, "C1-test")

    assert record["registered_floor"] == MIN_LANDMARKS
    assert record["meets_registered_floor"] == (record["n_landmarks_in_basis"] >= MIN_LANDMARKS)


def test_basis_columns_refuses_a_file_the_basis_does_not_land_on():
    carried = ["LM000", "LM001", "LM002"]
    assert basis_columns(carried, ["LM002", "LM000"]).tolist() == [2, 0]
    with pytest.raises(AssertionError, match="absent"):
        basis_columns(carried, ["LM000", "LM009"])
    with pytest.raises(AssertionError, match="more than once"):
        basis_columns(carried + ["LM001"], ["LM001"])


def test_a_non_finite_operand_is_refused_rather_than_normalized():
    with pytest.raises(AssertionError, match="not finite"):
        unit(np.array([1.0, np.inf, 2.0]))
    with pytest.raises(AssertionError, match="not finite"):
        cosine(np.array([1.0, np.nan]), np.array([1.0, 1.0]))
    with pytest.raises(AssertionError, match="not finite"):
        require_finite(np.array([[0.0, -np.inf]]), "a test operand")


def test_an_undefined_median_is_refused_rather_than_read_as_below_the_threshold():
    merged = {"A": {"cosine_with_released_bulk": 1.0},
              "B": {"cosine_with_released_bulk": float("nan")},
              "C": {"status": "not computed", "reason": "fewer than two cells in the release"}}
    summary = verification_summary("C1-test", merged)

    assert summary["status"] == "refused"
    assert summary["median_cosine"] is None
    assert summary["describes"] is None
    assert summary["n_total"] == 2 and summary["n_finite"] == 1 and summary["n_non_finite"] == 1
    assert summary["n_not_computed"] == 1
    assert summary["reason_not_computed"]["C"] == "fewer than two cells in the release"


def test_a_complete_set_of_cosines_takes_the_registered_label():
    high = verification_summary("C1-test", {str(i): {"cosine_with_released_bulk": 0.999}
                                            for i in range(5)})
    low = verification_summary("C1-test", {str(i): {"cosine_with_released_bulk": 0.5}
                                           for i in range(5)})

    assert high["status"] == "computed" and high["describes"] == "C1"
    assert low["status"] == "computed"
    assert low["describes"] == "the single-cell construction, not C1"
    assert low["n_below_threshold"] == 5


def test_no_undefined_value_reaches_the_serialized_output():
    merged = {"A": {"cosine_with_released_bulk": float("nan"), "reliability": float("inf")},
              "B": {"cosine_with_released_bulk": 0.9}}
    rendered = json.dumps(verification_summary("C1-test", merged), allow_nan=False)

    assert "NaN" not in rendered and "Infinity" not in rendered
    assert json.loads(rendered)["median_cosine"] is None
