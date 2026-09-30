import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
_mod = __import__("03e_reconstruction_gate")


def _pair(n_drugs=6, n_cells=4, n_genes=20, shuffle_rows=False, shuffle_genes=False):
    """The same data on both sides, optionally stored in a different order."""
    rng = np.random.default_rng()
    genes = [str(1000 + i) for i in range(n_genes)]
    cells = [f"LINE{i}" for i in range(n_cells)]
    matrices = {f"drug{d}": rng.standard_normal((n_cells, n_genes)) for d in range(n_drugs)}

    left = ({d: m.copy() for d, m in matrices.items()},
            {d: list(cells) for d in matrices}, list(genes))
    right_matrices, right_cells, right_genes = {}, {}, list(genes)
    if shuffle_genes:
        order = rng.permutation(n_genes)
        right_genes = [genes[i] for i in order]
    for drug, matrix in matrices.items():
        row_order = rng.permutation(n_cells) if shuffle_rows else np.arange(n_cells)
        column_order = np.array([genes.index(g) for g in right_genes])
        right_matrices[drug] = matrix[row_order][:, column_order]
        right_cells[drug] = [cells[i] for i in row_order]
    return left, (right_matrices, right_cells, right_genes)


def test_identical_data_passes_whatever_order_it_is_stored_in():
    for rows, columns in ((False, False), (True, False), (False, True), (True, True)):
        left, right = _pair(shuffle_rows=rows, shuffle_genes=columns)
        # the rebuild is `right`, stored in its own order; the extraction is `left`
        result = _mod.compare(right, left)
        assert result["all_within_tolerance"]
        assert result["max_abs_difference"] == pytest.approx(0.0, abs=1e-12)


def test_a_value_that_moved_is_caught_even_when_the_order_also_moved():
    left, right = _pair(shuffle_rows=True, shuffle_genes=True)
    right[0]["drug2"][1, 3] += 0.5

    result = _mod.compare(right, left)
    assert not result["all_within_tolerance"]
    assert not result["per_drug"]["drug2"]["within_tolerance"]
    assert result["per_drug"]["drug0"]["within_tolerance"]


def test_float32_rounding_passes_at_the_registered_tolerance():
    left, right = _pair()
    rounded = {d: m.astype(np.float32).astype(np.float64) for d, m in right[0].items()}

    result = _mod.compare((rounded, right[1], right[2]), left)
    assert result["all_within_tolerance"]
    assert result["max_abs_difference"] > 0          # rounding really did move the values


def test_a_missing_cell_line_is_named_not_silently_dropped():
    left, right = _pair()
    right[1]["drug1"] = right[1]["drug1"][:-1]
    right[0]["drug1"] = right[0]["drug1"][:-1]

    with pytest.raises(AssertionError, match="identifier sets differ"):
        _mod.compare(right, left)


def test_a_repeated_identifier_is_refused():
    left, right = _pair()
    right[1]["drug0"][1] = right[1]["drug0"][0]

    with pytest.raises(AssertionError, match="repeats an identifier"):
        _mod.compare(right, left)


def test_a_different_gene_axis_is_refused():
    left, right = _pair()
    right = (right[0], right[1], right[2][:-1] + ["9999"])

    with pytest.raises(AssertionError, match="gene axes differ"):
        _mod.compare(right, left)


def test_align_returns_rows_in_a_shared_order():
    rng = np.random.default_rng()
    rows = rng.standard_normal((5, 3))
    ids = [f"C{i}" for i in range(5)]
    order = rng.permutation(5)

    left, right, used = _mod.align(rows[order], [ids[i] for i in order], rows, ids)
    assert left == pytest.approx(right)
    assert used == [ids[i] for i in order]


def _shrna_pair(n_targets=4, per_target=3, n_genes=15, shuffle_genes=False):
    """Signatures on both sides, grouped into targets, optionally in different gene order."""
    rng = np.random.default_rng()
    genes = [str(2000 + i) for i in range(n_genes)]
    signatures, membership = {}, {}
    for t in range(n_targets):
        ids = [f"SIG_{t}_{k}" for k in range(per_target)]
        membership[f"TARGET{t}"] = ids
        for sig_id in ids:
            signatures[sig_id] = rng.standard_normal(n_genes)
    left = ({k: v.copy() for k, v in signatures.items()},
            {k: list(v) for k, v in membership.items()}, list(genes))
    right_genes = list(genes)
    if shuffle_genes:
        right_genes = [genes[i] for i in rng.permutation(n_genes)]
    order = np.array([genes.index(g) for g in right_genes])
    right = ({k: v[order].copy() for k, v in signatures.items()},
             {k: list(v) for k, v in membership.items()}, right_genes)
    return left, right


def test_shrna_gate_passes_identical_data_in_a_different_gene_order():
    left, right = _shrna_pair(shuffle_genes=True)
    result = _mod.compare_shrna(right, left)
    assert result["all_within_tolerance"]
    assert result["max_abs_consensus_difference"] == pytest.approx(0.0, abs=1e-12)
    assert result["n_targets"] == 4


def test_shrna_gate_catches_a_changed_signature():
    left, right = _shrna_pair()
    right[0]["SIG_1_0"][2] += 0.5
    with pytest.raises(AssertionError, match="signature SIG_1_0 differs"):
        _mod.compare_shrna(right, left)


def test_shrna_gate_catches_signatures_assigned_to_the_wrong_target():
    left, right = _shrna_pair()
    # swap one signature between two targets: every vector is still correct, and
    # every target still has the same count, but the consensuses are not the same
    right[1]["TARGET0"][0], right[1]["TARGET1"][0] = right[1]["TARGET1"][0], right[1]["TARGET0"][0]
    with pytest.raises(AssertionError, match="group different signatures"):
        _mod.compare_shrna(right, left)


def test_shrna_gate_catches_a_missing_signature_and_a_missing_target():
    left, right = _shrna_pair()
    dropped = dict(right[0])
    dropped.pop("SIG_2_1")
    with pytest.raises(AssertionError, match="signature-id sets differ"):
        _mod.compare_shrna((dropped, right[1], right[2]), left)

    left, right = _shrna_pair()
    fewer = {k: v for k, v in right[1].items() if k != "TARGET3"}
    with pytest.raises(AssertionError, match="target sets differ"):
        _mod.compare_shrna((right[0], fewer, right[2]), left)


def test_a_scaled_target_is_caught_at_the_signature_stage():
    left, right = _shrna_pair(n_targets=2, per_target=3)
    for sig_id in right[1]["TARGET0"]:
        right[0][sig_id] = right[0][sig_id] * 3.0

    # the unit direction of that target is unchanged by scaling, so the signature
    # comparison is what catches it
    with pytest.raises(AssertionError, match="differs by"):
        _mod.compare_shrna(right, left)


def test_a_wrongly_built_stored_consensus_is_caught_even_when_signatures_agree():
    left, right = _shrna_pair(n_targets=3, per_target=3)
    signatures, membership, genes = right

    honest = {}
    for target, ids in membership.items():
        mean = np.vstack([signatures[i] for i in ids]).mean(axis=0)
        honest[target] = mean / np.linalg.norm(mean)
    assert _mod.compare_shrna(right, left, stored_directions=honest)["all_within_tolerance"]

    # Deviation 9's shape: every signature is right, and the consensus was built
    # from the wrong rows
    defective = dict(honest)
    wrong_rows = [signatures[i] for i in membership["TARGET1"][:1] + membership["TARGET2"][:2]]
    mean = np.vstack(wrong_rows).mean(axis=0)
    defective["TARGET0"] = mean / np.linalg.norm(mean)

    result = _mod.compare_shrna(right, left, stored_directions=defective)
    assert not result["all_within_tolerance"]
    assert result["max_abs_stored_direction_difference"] > 0.01
    assert result["per_target"]["TARGET0"]["max_abs_stored_direction_difference"] > 0.01
    assert result["per_target"]["TARGET1"]["max_abs_stored_direction_difference"] == pytest.approx(
        0.0, abs=1e-12)


def test_shrna_gate_refuses_a_target_below_the_eligibility_floor():
    left, right = _shrna_pair(per_target=2)
    with pytest.raises(AssertionError, match="below the eligibility floor"):
        _mod.compare_shrna(right, left)


def test_compound_gate_fails_when_a_cohort_drug_is_missing_from_either_side():
    left, right = _pair(n_drugs=5)
    cohort = sorted(left[0])

    short = ({d: m for d, m in right[0].items() if d != "drug3"},
             {d: c for d, c in right[1].items() if d != "drug3"}, right[2])
    with pytest.raises(AssertionError, match="absent from the rebuild"):
        _mod.compare(short, left, cohort=cohort)

    short_extraction = ({d: m for d, m in left[0].items() if d != "drug2"},
                        {d: c for d, c in left[1].items() if d != "drug2"}, left[2])
    with pytest.raises(AssertionError, match="absent from the extraction"):
        _mod.compare(right, short_extraction, cohort=cohort)


def test_a_rebuild_holding_one_drug_cannot_pass_the_cohort_check():
    left, right = _pair(n_drugs=6)
    cohort = sorted(left[0])
    single = ({"drug0": right[0]["drug0"]}, {"drug0": right[1]["drug0"]}, right[2])

    with pytest.raises(AssertionError, match="absent from the rebuild"):
        _mod.compare(single, left, cohort=cohort)


def test_extra_drugs_outside_the_cohort_are_counted_not_silently_dropped():
    left, right = _pair(n_drugs=5)
    cohort = sorted(left[0])[:4]

    result = _mod.compare(right, left, cohort=cohort)
    assert result["n_drugs"] == 4
    assert result["n_outside_cohort"] == {"rebuild": 1, "extraction": 1}
