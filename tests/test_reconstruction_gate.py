import json
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


def _shrna_pair(n_targets=4, per_target=3, n_genes=15, shuffle_genes=False, offset=0.0):
    """Signatures on both sides, grouped into targets, optionally in different gene order.

    `offset` lifts every coordinate away from zero, which is what the GCTX values
    do: the registered tolerance is relative, so a comparison at magnitude ten
    allows an order of magnitude more than one at magnitude one.
    """
    rng = np.random.default_rng()
    genes = [str(2000 + i) for i in range(n_genes)]
    signatures, membership = {}, {}
    for t in range(n_targets):
        ids = [f"SIG_{t}_{k}" for k in range(per_target)]
        membership[f"TARGET{t}"] = ids
        for sig_id in ids:
            signatures[sig_id] = offset + rng.standard_normal(n_genes)
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


def test_shrna_gate_reports_a_changed_signature_rather_than_raising():
    left, right = _shrna_pair()
    right[0]["SIG_1_0"][2] += 0.5

    # raising here wrote no record of the failure the gate exists to document
    result = _mod.compare_shrna(right, left)
    assert not result["all_within_tolerance"]
    assert not result["signatures_within_tolerance"]
    assert result["n_signature_failures"] == 1
    assert result["worst_signature_failures"][0]["sig_id"] == "SIG_1_0"
    assert result["max_abs_signature_difference"] == pytest.approx(0.5, abs=1e-9)


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
    # comparison is the only one that can catch it
    result = _mod.compare_shrna(right, left)
    assert not result["signatures_within_tolerance"]
    assert result["n_signature_failures"] == 3
    assert result["per_target"]["TARGET0"]["max_abs_direction_difference"] == pytest.approx(
        0.0, abs=1e-12)


def test_a_wrongly_built_stored_consensus_is_caught_even_when_signatures_agree():
    left, right = _shrna_pair(n_targets=3, per_target=3)
    signatures, membership, genes = right

    honest = {}
    for target, ids in membership.items():
        mean = np.vstack([signatures[i] for i in ids]).mean(axis=0)
        honest[target] = mean / np.linalg.norm(mean)
    assert _mod.compare_shrna(right, left,
                              stored=(honest, genes))["all_within_tolerance"]

    # Deviation 9's shape: every signature is right, and the consensus was built
    # from the wrong rows
    defective = dict(honest)
    wrong_rows = [signatures[i] for i in membership["TARGET1"][:1] + membership["TARGET2"][:2]]
    mean = np.vstack(wrong_rows).mean(axis=0)
    defective["TARGET0"] = mean / np.linalg.norm(mean)

    result = _mod.compare_shrna(right, left, stored=(defective, genes))
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


def test_stored_direction_is_aligned_on_its_own_axis_not_the_extraction_permutation():
    # the rebuild and the extraction disagree about gene order, which is the case
    # that made the previous version of this check wrong
    left, right = _shrna_pair(n_targets=3, per_target=3, shuffle_genes=True)
    signatures, membership, rebuilt_genes = right

    honest = {}
    for target, ids in membership.items():
        mean = np.vstack([signatures[i] for i in ids]).mean(axis=0)
        honest[target] = mean / np.linalg.norm(mean)

    result = _mod.compare_shrna(right, left, stored=(honest, rebuilt_genes))
    assert result["all_within_tolerance"]
    assert result["max_abs_stored_direction_difference"] == pytest.approx(0.0, abs=1e-12)

    # the same values stored against a consistently reordered axis still pass,
    # because the axis is what the alignment uses
    shuffled = list(reversed(rebuilt_genes))
    reordered = {t: v[[rebuilt_genes.index(g) for g in shuffled]] for t, v in honest.items()}
    assert _mod.compare_shrna(right, left,
                              stored=(reordered, shuffled))["all_within_tolerance"]

    # values in the rebuild's order, declared against a different axis, must not:
    # that is a mislabeled coordinate system, which is the error being guarded
    assert not _mod.compare_shrna(right, left,
                                  stored=(honest, shuffled))["all_within_tolerance"]


def test_a_stored_direction_with_two_coordinates_swapped_fails_and_names_the_target():
    left, right = _shrna_pair(n_targets=3, per_target=3, shuffle_genes=True)
    signatures, membership, rebuilt_genes = right
    honest = {}
    for target, ids in membership.items():
        mean = np.vstack([signatures[i] for i in ids]).mean(axis=0)
        honest[target] = mean / np.linalg.norm(mean)

    corrupted = {t: v.copy() for t, v in honest.items()}
    corrupted["TARGET2"][[0, 1]] = corrupted["TARGET2"][[1, 0]]

    result = _mod.compare_shrna(right, left, stored=(corrupted, rebuilt_genes))
    assert not result["all_within_tolerance"]
    assert result["per_target"]["TARGET2"]["max_abs_stored_direction_difference"] > 0
    assert result["per_target"]["TARGET0"]["max_abs_stored_direction_difference"] == pytest.approx(
        0.0, abs=1e-12)


def test_a_stored_axis_that_is_not_the_rebuild_axis_is_refused():
    left, right = _shrna_pair(n_targets=2, per_target=3)
    signatures, membership, genes = right
    honest = {t: np.ones(len(genes)) / np.sqrt(len(genes)) for t in membership}

    with pytest.raises(AssertionError, match="does not share the rebuild's gene axis"):
        _mod.compare_shrna(right, left, stored=(honest, genes[:-1] + ["9999"]))


def test_a_zero_norm_consensus_is_named_rather_than_dividing_by_zero():
    left, right = _shrna_pair(n_targets=2, per_target=2)
    for side in (left, right):
        for sig_id in side[1]["TARGET0"]:
            side[0][sig_id] = np.zeros_like(side[0][sig_id])

    with pytest.raises(AssertionError, match="zero-norm consensus"):
        _mod.compare_shrna(right, left, min_signatures=2)


def test_the_gate_refuses_a_rebuild_with_no_stored_consensus(tmp_path):
    rng = np.random.default_rng()
    np.savez_compressed(tmp_path / "shrna_signatures.npz",
                        sig_ids=np.array(["SIG1", "SIG2", "SIG3"]),
                        signatures=rng.standard_normal((3, 5)),
                        gene_ids=np.array([str(i) for i in range(5)]),
                        membership=json.dumps({"TARGET0": ["SIG1", "SIG2", "SIG3"]}))

    # the lenient loader tolerates the absence; the one the gate uses does not
    assert _mod.load_shrna_rebuild(tmp_path)[3] is None
    with pytest.raises(AssertionError, match="shrna_consensus.npz is missing"):
        _mod.load_shrna_rebuild_strict(tmp_path)


def test_the_gate_refuses_stored_directions_covering_other_targets():
    left, right = _shrna_pair(n_targets=3, per_target=3)
    signatures, membership, genes = right
    honest = {}
    for target, ids in membership.items():
        mean = np.vstack([signatures[i] for i in ids]).mean(axis=0)
        honest[target] = mean / np.linalg.norm(mean)

    extra = dict(honest)
    extra["TARGET_NOT_IN_THE_REBUILD"] = honest["TARGET0"]
    with pytest.raises(AssertionError, match="different set of targets"):
        _mod.compare_shrna(right, left, stored=(extra, genes))


def test_the_gate_refuses_a_nonfinite_or_unnormalized_stored_direction():
    left, right = _shrna_pair(n_targets=2, per_target=3)
    signatures, membership, genes = right
    honest = {}
    for target, ids in membership.items():
        mean = np.vstack([signatures[i] for i in ids]).mean(axis=0)
        honest[target] = mean / np.linalg.norm(mean)

    with_nan = {t: v.copy() for t, v in honest.items()}
    with_nan["TARGET0"][2] = np.nan
    with pytest.raises(AssertionError, match="nonfinite stored direction"):
        _mod.compare_shrna(right, left, stored=(with_nan, genes))

    unnormalized = {t: v * 4.0 for t, v in honest.items()}
    with pytest.raises(AssertionError, match="contract is a unit vector"):
        _mod.compare_shrna(right, left, stored=(unnormalized, genes))

    wrong_length = {t: v[:-1] for t, v in honest.items()}
    with pytest.raises(AssertionError, match="stored direction is"):
        _mod.compare_shrna(right, left, stored=(wrong_length, genes))


def test_the_gate_refuses_a_repeated_gene_identifier_on_any_axis():
    left, right = _shrna_pair(n_targets=2, per_target=3)
    duplicated = right[2][:-1] + [right[2][0]]
    with pytest.raises(AssertionError, match="repeats an identifier"):
        _mod.compare_shrna((right[0], right[1], duplicated), left)

    signatures, membership, genes = right
    honest = {t: np.ones(len(genes)) / np.sqrt(len(genes)) for t in membership}
    with pytest.raises(AssertionError, match="repeats a gene identifier"):
        _mod.compare_shrna(right, left, stored=(honest, genes[:-1] + [genes[0]]))


def _write_shard(path, matrices, cells, fingerprint="fp"):
    """A shard in the layout the Modal rebuild writes."""
    np.savez_compressed(path, fingerprint=np.array(fingerprint),
                        __cells__=np.array(json.dumps(cells)), **matrices)


def test_the_loader_keeps_matrices_and_cell_lines_apart(tmp_path):
    rng = np.random.default_rng()
    matrices = {"drugA": rng.standard_normal((3, 6)), "drugB": rng.standard_normal((2, 6))}
    cells = {"drugA": ["MCF7", "PC3", "A375"], "drugB": ["HT29", "A549"]}
    _write_shard(tmp_path / "shard_000.npz", matrices, cells)
    (tmp_path / "landmark_gene_ids.json").write_text(json.dumps([str(i) for i in range(6)]))

    loaded_matrices, loaded_cells, genes = _mod.load_rebuild(tmp_path)
    assert set(loaded_matrices) == {"drugA", "drugB"}
    for drug in matrices:
        # the failure this guards: a cell-line array arriving under a drug's own
        # name and replacing its matrix
        assert loaded_matrices[drug].dtype.kind == "f"
        assert loaded_matrices[drug] == pytest.approx(matrices[drug])
        assert loaded_cells[drug] == cells[drug]
    assert genes == [str(i) for i in range(6)]


def test_a_shard_without_cell_identifiers_is_refused(tmp_path):
    rng = np.random.default_rng()
    np.savez_compressed(tmp_path / "shard_000.npz", fingerprint=np.array("fp"),
                        drugA=rng.standard_normal((3, 6)))
    (tmp_path / "landmark_gene_ids.json").write_text(json.dumps([str(i) for i in range(6)]))

    with pytest.raises(AssertionError, match="carries no cell-line identifiers"):
        _mod.load_rebuild(tmp_path)


def test_the_extraction_loader_builds_only_the_cohort(tmp_path, monkeypatch):
    import pandas as pd

    rng = np.random.default_rng()
    sig_ids = [f"SIG{i}" for i in range(6)]
    np.savez_compressed(tmp_path / "lincs_subset.npz",
                        sig_ids=np.array(sig_ids),
                        signatures=rng.standard_normal((6, 4)),
                        gene_ids=np.array([str(i) for i in range(4)]))
    pd.DataFrame({"sig_id": sig_ids,
                  "pert_iname": ["wanted", "wanted", "other", "other", "third", "third"],
                  "cell_id": ["MCF7", "PC3"] * 3}).to_csv(
        tmp_path / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t", index=False,
        compression="gzip")

    everything, _, _ = _mod.load_extraction(tmp_path)
    assert set(everything) == {"wanted", "other", "third"}

    # the gate compares one cohort, so the loader must not build the rest
    restricted, cells, genes = _mod.load_extraction(tmp_path, cohort=["wanted"])
    assert set(restricted) == {"wanted"}
    assert restricted["wanted"].shape == (2, 4)
    assert cells["wanted"] == ["MCF7", "PC3"]
    assert genes == [str(i) for i in range(4)]


def test_consensus_tolerance_is_elementwise_not_a_unit_scale_threshold():
    left, right = _shrna_pair(n_targets=3, per_target=3, offset=10.0)
    # every element moves by half the registered relative tolerance: inside the
    # elementwise rule everywhere, and above the 2e-5 that one unit-scale
    # threshold allowed, which failed consensuses the tolerance admits
    for sig_id in right[1]["TARGET0"]:
        right[0][sig_id] = right[0][sig_id] * (1.0 + 0.5 * _mod.RTOL)

    result = _mod.compare_shrna(right, left)
    assert result["all_within_tolerance"]
    assert result["max_abs_consensus_difference"] > _mod.ATOL + _mod.RTOL

    left, right = _shrna_pair(n_targets=3, per_target=3, offset=10.0)
    for sig_id in right[1]["TARGET0"]:
        right[0][sig_id] = right[0][sig_id] * (1.0 + 20.0 * _mod.RTOL)

    outside = _mod.compare_shrna(right, left)
    assert not outside["all_within_tolerance"]
    assert outside["n_signature_failures"] == 3


def test_the_record_names_the_worst_target_for_each_comparison_separately():
    left, right = _shrna_pair(n_targets=3, per_target=3, shuffle_genes=True)
    signatures, membership, rebuilt_genes = right
    honest = {}
    for target, ids in membership.items():
        mean = np.vstack([signatures[i] for i in ids]).mean(axis=0)
        honest[target] = mean / np.linalg.norm(mean)

    # one target's extraction signatures moved, a different target's stored
    # direction is wrong: a single worst-target field can only name one of them
    for sig_id in membership["TARGET0"]:
        left[0][sig_id] = left[0][sig_id] + 0.5
    corrupted = {t: v.copy() for t, v in honest.items()}
    corrupted["TARGET2"][[0, 1]] = corrupted["TARGET2"][[1, 0]]

    result = _mod.compare_shrna(right, left, stored=(corrupted, rebuilt_genes))
    assert not result["all_within_tolerance"]
    assert result["worst_target"]["consensus"] == "TARGET0"
    assert result["worst_target"]["direction"] == "TARGET0"
    assert result["worst_target"]["stored"] == "TARGET2"
    assert result["per_target"]["TARGET2"]["max_abs_stored_direction_difference"] > 0
    assert result["per_target"]["TARGET0"]["max_abs_stored_direction_difference"] == pytest.approx(
        0.0, abs=1e-12)


def test_a_declared_axis_the_values_do_not_follow_is_caught():
    # Deviation 11's shape: correct values under a declared gene order they do not
    # use. Both sides agree on their gene_ids, so every check that compares labels
    # passes and only a comparison of values can see it.
    left, right = _shrna_pair(n_targets=3, per_target=3)
    rng = np.random.default_rng()
    scramble = rng.permutation(len(right[2]))
    while (scramble == np.arange(len(scramble))).all():
        scramble = rng.permutation(len(scramble))

    mislabeled = {sig_id: vector[scramble] for sig_id, vector in left[0].items()}
    assert left[2] == right[2]                       # the declared axes still match

    result = _mod.compare_shrna(right, (mislabeled, left[1], left[2]))
    assert not result["all_within_tolerance"]
    assert result["n_signature_failures"] == 9


def test_compare_both_records_a_structural_failure_instead_of_losing_it():
    compound_left, compound_right = _pair(n_drugs=4)
    shrna_left, shrna_right = _shrna_pair(n_targets=3, per_target=3)
    cohort = sorted(compound_left[0])

    passing = _mod.compare_both(compound_right, compound_left, shrna_right, shrna_left,
                                cohort=cohort)
    assert passing["all_within_tolerance"]
    assert "structural_failure" not in passing

    fewer = {k: v for k, v in shrna_right[1].items() if k != "TARGET2"}
    failed = _mod.compare_both(compound_right, compound_left,
                               (shrna_right[0], fewer, shrna_right[2]), shrna_left,
                               cohort=cohort)
    assert not failed["all_within_tolerance"]
    assert "target sets differ" in failed["structural_failure"]
    # the half that ran before the failure is still in the record
    assert failed["compound_signatures"]["all_within_tolerance"]


def _parsed_frame(n_genes=6, n_sigs=5, extra_rows=2, gene_order=None):
    """A frame shaped like a cmapPy parse: genes down the index, signatures across."""
    import pandas as pd

    rng = np.random.default_rng()
    genes = gene_order or [str(3000 + i) for i in range(n_genes)]
    rows = genes + [f"NOT_LANDMARK_{i}" for i in range(extra_rows)]
    sigs = [f"SIG_{i}" for i in range(n_sigs)]
    return pd.DataFrame(rng.standard_normal((len(rows), n_sigs)), index=rows, columns=sigs), genes


def test_the_parse_keeps_each_value_with_the_label_the_frame_gave_it():
    frame, genes = _parsed_frame()
    membership = {"TARGET0": ["SIG_0", "SIG_1", "SIG_2"]}

    signatures, parsed_genes = _mod.shrna_from_parsed_frame(frame, membership, genes)
    assert parsed_genes == genes                       # the frame's own order, not the request's
    assert set(signatures) == {"SIG_0", "SIG_1", "SIG_2"}
    for sig_id, vector in signatures.items():
        for gene, value in zip(parsed_genes, vector):
            assert value == pytest.approx(frame.loc[gene, sig_id])


def test_the_parse_reports_the_frames_order_not_the_requested_one():
    # the production route asks for a row order and reindexes onto it, which hides
    # a frame whose rows arrive in another order; this must not
    frame, genes = _parsed_frame()
    shuffled = list(reversed(genes))

    signatures, parsed_genes = _mod.shrna_from_parsed_frame(
        frame, {"TARGET0": ["SIG_0", "SIG_1", "SIG_2"]}, shuffled)
    assert parsed_genes == genes
    for gene, value in zip(parsed_genes, signatures["SIG_0"]):
        assert value == pytest.approx(frame.loc[gene, "SIG_0"])


def test_the_parse_refuses_a_frame_missing_a_landmark_or_repeating_an_identifier():
    frame, genes = _parsed_frame()
    with pytest.raises(AssertionError, match="of .* landmark genes"):
        _mod.shrna_from_parsed_frame(frame, {"T": ["SIG_0"]}, genes + ["9999"])

    duplicated = frame.copy()
    duplicated.index = [frame.index[0]] + list(frame.index[1:])
    duplicated = duplicated.set_axis([frame.index[0]] * 2 + list(frame.index[2:]), axis=0)
    with pytest.raises(AssertionError, match="repeats a gene identifier"):
        _mod.shrna_from_parsed_frame(duplicated, {"T": ["SIG_0"]}, genes)


def _write_gctx(path, gene_order, signature_order, orientation="signatures_then_genes"):
    """A GCTX in the real layout, with a sentinel per (gene, signature) pair.

    The value at a pair encodes both identifiers, so a value found under the wrong
    label is visible rather than merely different.
    """
    import h5py

    genes, signatures = list(gene_order), list(signature_order)
    sentinel = np.array([[1000.0 * (g + 1) + (s + 1) for g in range(len(genes))]
                         for s in range(len(signatures))])
    stored = sentinel if orientation == "signatures_then_genes" else sentinel.T
    with h5py.File(path, "w") as handle:
        handle.create_dataset("/0/DATA/0/matrix", data=stored.astype(np.float32))
        handle.create_dataset("/0/META/ROW/id",
                              data=np.array(genes, dtype=h5py.string_dtype()))
        handle.create_dataset("/0/META/COL/id",
                              data=np.array(signatures, dtype=h5py.string_dtype()))
    return {(genes[g], signatures[s]): sentinel[s, g]
            for g in range(len(genes)) for s in range(len(signatures))}


def test_the_source_reader_pairs_each_value_with_the_identifier_the_file_stores(tmp_path):
    # rows are physically C, A, B and are requested as B, C, A
    genes = ["C", "A", "B", "D"]
    signatures = [f"SIG{i}" for i in range(5)]
    expected = _write_gctx(tmp_path / "x.gctx", genes, signatures)

    values, read_genes, read_signatures, hashes = _mod.read_gctx_slice(
        tmp_path / "x.gctx", ["SIG3", "SIG0"], ["B", "C", "A"])
    assert set(read_genes) == {"A", "B", "C"}
    assert set(read_signatures) == {"SIG0", "SIG3"}
    for row, sig_id in enumerate(read_signatures):
        for column, gene in enumerate(read_genes):
            assert values[row, column] == pytest.approx(expected[(gene, sig_id)])
    assert hashes["orientation"] == "signatures_then_genes"
    assert hashes["source_row_axis_sha256"] != hashes["selected_gene_axis_sha256"]


def test_the_source_reader_reads_either_orientation_and_refuses_an_ambiguous_one(tmp_path):
    genes, signatures = ["G1", "G2", "G3"], [f"SIG{i}" for i in range(6)]
    expected = _write_gctx(tmp_path / "t.gctx", genes, signatures,
                           orientation="genes_then_signatures")

    values, read_genes, read_signatures, hashes = _mod.read_gctx_slice(
        tmp_path / "t.gctx", signatures, genes)
    assert hashes["orientation"] == "genes_then_signatures"
    for row, sig_id in enumerate(read_signatures):
        for column, gene in enumerate(read_genes):
            assert values[row, column] == pytest.approx(expected[(gene, sig_id)])

    # a square matrix fits both readings, and nothing in the file settles it
    square = ["G1", "G2", "G3"]
    _write_gctx(tmp_path / "square.gctx", square, ["SIG0", "SIG1", "SIG2"])
    with pytest.raises(AssertionError, match="fits 2 orientations"):
        _mod.read_gctx_slice(tmp_path / "square.gctx", ["SIG0"], square)


def test_the_source_reader_refuses_a_missing_or_repeated_identifier(tmp_path):
    genes, signatures = ["G1", "G2", "G3", "G4"], [f"SIG{i}" for i in range(5)]
    _write_gctx(tmp_path / "x.gctx", genes, signatures)

    with pytest.raises(AssertionError, match="requested genes are absent"):
        _mod.read_gctx_slice(tmp_path / "x.gctx", ["SIG0"], genes + ["G9"])
    with pytest.raises(AssertionError, match="requested signatures are absent"):
        _mod.read_gctx_slice(tmp_path / "x.gctx", ["SIG_NOT_THERE"], genes)

    import h5py
    with h5py.File(tmp_path / "dup.gctx", "w") as handle:
        handle.create_dataset("/0/DATA/0/matrix", data=np.zeros((5, 4), dtype=np.float32))
        handle.create_dataset("/0/META/ROW/id",
                              data=np.array(["G1", "G1", "G3", "G4"],
                                            dtype=h5py.string_dtype()))
        handle.create_dataset("/0/META/COL/id",
                              data=np.array(signatures, dtype=h5py.string_dtype()))
    with pytest.raises(AssertionError, match="repeats an identifier"):
        _mod.read_gctx_slice(tmp_path / "dup.gctx", ["SIG0"], ["G1"])


def test_a_globally_permuted_canonical_artifact_fails_against_the_source(tmp_path):
    # the defect, end to end: the canonical artifact declares the right axis and
    # its columns follow another one, and the source read is what catches it
    genes = [f"G{i}" for i in range(8)]
    signatures = [f"SIG{i}" for i in range(6)]
    _write_gctx(tmp_path / "x.gctx", genes, signatures)
    values, read_genes, read_signatures, _ = _mod.read_gctx_slice(
        tmp_path / "x.gctx", signatures, genes)

    membership = {"TARGET0": sorted(read_signatures[:3]),
                  "TARGET1": sorted(read_signatures[3:])}
    honest = {sig_id: values[i] for i, sig_id in enumerate(read_signatures)}
    source = (honest, membership, read_genes)

    assert _mod.compare_shrna((dict(honest), membership, list(read_genes)),
                              source)["all_within_tolerance"]

    rng = np.random.default_rng()
    scramble = rng.permutation(len(read_genes))
    while (scramble == np.arange(len(scramble))).all():
        scramble = rng.permutation(len(scramble))
    mislabeled = {sig_id: vector[scramble] for sig_id, vector in honest.items()}

    result = _mod.compare_shrna((mislabeled, membership, list(read_genes)), source)
    assert not result["all_within_tolerance"]
    assert result["n_signature_failures"] == len(read_signatures)
