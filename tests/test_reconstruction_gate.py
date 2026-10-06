import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
_mod = __import__("03e_reconstruction_gate")
_cov = __import__("03f_shrna_coverage_identity")
_reg = __import__("03g_build_cohort_manifests")
_freeze = __import__("03h_freeze_identities")


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
    # a target grouping a signature the side does not hold is refused by name,
    # before any comparison indexes it
    with pytest.raises(AssertionError, match=r"groups 1 signatures it does not hold"):
        _mod.compare_shrna((dropped, right[1], right[2]), left)

    # and a signature held by one side alone, inside nobody's membership, is
    # refused on the identifier sets
    extra = dict(right[0])
    extra["SIG_UNGROUPED"] = np.zeros(len(right[2]))
    with pytest.raises(AssertionError, match="signature-id sets differ"):
        _mod.compare_shrna((extra, right[1], right[2]), left)

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


def test_one_predicate_failing_does_not_stop_the_others_being_evaluated():
    # v1 wrapped every comparison in one `try`, so a structural failure in the shRNA
    # half meant the compound half it had already computed was never recorded
    order = []

    def ok(name):
        def run():
            order.append(name)
            return {"all_within_tolerance": True}
        return run

    def boom():
        order.append("legacy_artifact_integrity")
        _mod.require(False, "the map is not a total bijection")

    result = _mod.evaluate({
        "cohort_identity": ok("cohort_identity"),
        "compound_source_reconstruction": ok("compound_source_reconstruction"),
        "shrna_reconstruction_on_the_eligible_set": ok("shrna"),
        "shrna_eligibility_partition": boom,
        "legacy_artifact_integrity": ok("legacy_artifact_integrity"),
    })

    assert not result["release_authorized"]
    assert result["structural_failures"] == {
        "shrna_eligibility_partition": "GateError: the map is not a total bijection"}
    # the predicate after the failure still ran, and the ones before are still recorded
    assert "legacy_artifact_integrity" in order
    assert result["predicate_held"]["compound_source_reconstruction"]
    assert result["predicate_held"]["legacy_artifact_integrity"]
    assert not result["predicate_held"]["shrna_eligibility_partition"]
    # a production failure and a legacy failure are reported apart
    assert result["production_authorized"] is False
    assert result["legacy_custody_intact"] is True


def test_a_predicate_raising_something_other_than_an_assertion_is_still_recorded():
    # catching only AssertionError let a malformed report (KeyError) or bad JSON
    # (JSONDecodeError) escape the wrapper, so the predicates after it never ran
    ran = []

    def ok(name):
        def run():
            ran.append(name)
            return {"all_within_tolerance": True}
        return run

    for failure in (KeyError("inputs"), json.JSONDecodeError("bad", "{", 0),
                    TypeError("expected str, bytes or os.PathLike, not NoneType")):
        ran.clear()

        def boom(failure=failure):
            raise failure

        result = _mod.evaluate({
            "cohort_identity": boom,
            "compound_source_reconstruction": ok("compound"),
            "shrna_reconstruction_on_the_eligible_set": ok("shrna"),
            "shrna_eligibility_partition": ok("partition"),
            "legacy_artifact_integrity": ok("legacy"),
        })
        assert not result["release_authorized"]
        assert type(failure).__name__ in result["structural_failures"]["cohort_identity"]
        # every predicate after the one that raised still ran
        assert len(ran) == 4


def test_the_release_condition_cannot_be_narrowed_by_the_invocation():
    # a gate whose scope depends on which arguments were passed lets a narrower run
    # look like a passing one
    fewer = {name: (lambda: {"all_within_tolerance": True})
             for name in _mod.RELEASE_PREDICATES[:-1]}
    with pytest.raises(AssertionError, match="frozen release condition"):
        _mod.evaluate(fewer)

    extra = {name: (lambda: {"all_within_tolerance": True})
             for name in _mod.RELEASE_PREDICATES}
    extra["something_else"] = lambda: {"all_within_tolerance": True}
    with pytest.raises(AssertionError, match="frozen release condition"):
        _mod.evaluate(extra)


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


def test_an_empty_comparison_is_refused_rather_than_passing():
    # the gate's worst failure mode: nothing to compare, so nothing can fail
    genes = [str(100 + i) for i in range(5)]
    empty = ({}, {}, list(genes))

    with pytest.raises(AssertionError, match="registered targets absent"):
        _mod.compare_shrna(empty, empty, stored=({}, genes),
                           expected_targets=["TARGET0"])

    with pytest.raises(AssertionError, match="holds no signatures"):
        _mod.compare_shrna(empty, empty)

    # and the same for the compound side, where an empty cohort would compare nothing
    rng = np.random.default_rng()
    compound = ({"drugA": rng.standard_normal((2, 5))}, {"drugA": ["MCF7", "PC3"]},
                list(genes))
    with pytest.raises(AssertionError, match="cohort is empty"):
        _mod.compare(compound, compound, cohort=set())


def test_a_registered_target_missing_from_both_sides_is_refused():
    left, right = _shrna_pair(n_targets=3, per_target=3)
    # both sides agree, and agree about a universe smaller than the registered one
    with pytest.raises(AssertionError,
                       match="registered targets absent from the rebuild, 1: TARGET_GONE"):
        _mod.compare_shrna(right, left,
                           expected_targets=["TARGET0", "TARGET1", "TARGET2", "TARGET_GONE"])

    # every missing identifier is named, because the message is the whole record of a
    # structural failure: it raises before the report is written
    absent = [f"GONE{i}" for i in range(9)]
    with pytest.raises(AssertionError) as raised:
        _mod.compare_shrna(right, left,
                           expected_targets=["TARGET0", "TARGET1", "TARGET2", *absent])
    assert all(name in str(raised.value) for name in absent)

    honest = {}
    for target, ids in right[1].items():
        mean = np.vstack([right[0][i] for i in ids]).mean(axis=0)
        honest[target] = mean / np.linalg.norm(mean)
    short = {t: v for t, v in honest.items() if t != "TARGET2"}
    with pytest.raises(AssertionError, match="with no stored direction, 1: TARGET2"):
        _mod.compare_shrna(right, left, stored=(short, right[2]),
                           expected_targets=sorted(right[1]))


def test_a_matrix_whose_row_count_contradicts_its_cell_lines_is_refused():
    # one numeric row declaring two cell lines would broadcast against two equal
    # source rows and pass
    genes = [str(100 + i) for i in range(5)]
    one_row = ({"drugB": np.ones((1, 5))}, {"drugB": ["MCF7", "PC3"]}, list(genes))
    two_rows = ({"drugB": np.ones((2, 5))}, {"drugB": ["MCF7", "PC3"]}, list(genes))

    with pytest.raises(AssertionError, match=r"the rebuild matrix is \(1, 5\) against 2"):
        _mod.compare(one_row, two_rows, cohort=["drugB"])


def test_a_repeated_gene_identifier_on_a_compound_axis_is_refused():
    rng = np.random.default_rng()
    genes = [str(100 + i) for i in range(5)]
    duplicated = genes[:-1] + [genes[0]]
    matrix = rng.standard_normal((2, 5))
    left = ({"drugC": matrix}, {"drugC": ["MCF7", "PC3"]}, duplicated)

    with pytest.raises(AssertionError, match="compound gene axis repeats an identifier"):
        _mod.compare(left, left, cohort=["drugC"])


def test_an_empty_cohort_cannot_pass_the_compound_comparison():
    with pytest.raises(AssertionError, match="the cohort is empty"):
        _mod.compare(({}, {}, ["1"]), ({}, {}, ["1"]), cohort=[])


def test_a_consensus_repeating_a_target_is_refused_before_the_dictionary_keeps_one(tmp_path):
    rng = np.random.default_rng()
    genes = [str(i) for i in range(5)]
    np.savez_compressed(tmp_path / "shrna_signatures.npz",
                        sig_ids=np.array(["SIG1", "SIG2", "SIG3"]),
                        signatures=rng.standard_normal((3, 5)),
                        gene_ids=np.array(genes),
                        membership=json.dumps({"TARGET0": ["SIG1", "SIG2", "SIG3"]}))
    # the corrupted direction comes first and the correct one second, which a
    # dictionary comprehension would silently resolve in favour of the correct one
    np.savez_compressed(tmp_path / "shrna_consensus.npz",
                        genes=np.array(["TARGET0", "TARGET0"]),
                        directions=np.vstack([np.ones(5) / np.sqrt(5), rng.standard_normal(5)]),
                        gene_ids=np.array(genes))

    with pytest.raises(AssertionError, match="repeats a target"):
        _mod.load_shrna_rebuild(tmp_path)


def test_a_drug_arriving_from_two_shards_is_refused(tmp_path):
    rng = np.random.default_rng()
    cells = {"drugA": ["MCF7", "PC3"]}
    _write_shard(tmp_path / "shard_000.npz", {"drugA": rng.standard_normal((2, 6))}, cells)
    _write_shard(tmp_path / "shard_001.npz", {"drugA": rng.standard_normal((2, 6))}, cells)
    (tmp_path / "landmark_gene_ids.json").write_text(json.dumps([str(i) for i in range(6)]))

    with pytest.raises(AssertionError, match="arrives from more than one shard"):
        _mod.load_rebuild(tmp_path)


def test_the_gate_holds_no_assert_statements():
    # `assert` disappears under -O, so a contract written as one is not a contract.
    # An earlier pass converted some and left thirty-eight, which read as done.
    import ast

    source = (Path(__file__).resolve().parents[1] / "experiments"
              / "03e_reconstruction_gate.py").read_text()
    remaining = [node.lineno for node in ast.walk(ast.parse(source))
                 if isinstance(node, ast.Assert)]
    assert remaining == [], f"assert statements at lines {remaining}"


def test_the_gate_contracts_survive_python_minus_o():
    # `assert` disappears under -O, so a contract written as one is not a contract
    import subprocess

    probe = (
        "import importlib.util, numpy as np;"
        "spec = importlib.util.spec_from_file_location('g', "
        f"{str(Path(__file__).resolve().parents[1] / 'experiments' / '03e_reconstruction_gate.py')!r});"
        "g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g);"
        "r = g.compare(({}, {}, ['1']), ({}, {}, ['1']), cohort=[])"
    )
    finished = subprocess.run(["python", "-O", "-c", probe], capture_output=True, text=True)
    assert finished.returncode != 0
    assert "the cohort is empty" in finished.stderr


def _only_report(output):
    """The one report under a gate output directory, which is per run."""
    reports = sorted(output.glob("gate_v*/reconstruction_gate.json"))
    assert len(reports) == 1, f"expected one report, found {len(reports)}"
    return reports[0]


def _gate_fixture(tmp_path):
    """Every input the registered gate takes, small but structurally complete.

    Three drugs, because the shRNA eligibility boundary is part of the release condition:
    two have shRNA data and one does not, so the paired cohort, the excluded records and
    their partition are all exercised rather than being empty.

    The retained extraction is written permuted against its declared axis, which is what
    the deposited one is, so the legacy-integrity predicate has something to be right
    about in both directions.
    """
    import pandas as pd

    rng = np.random.default_rng()
    n_genes = _mod.N_LANDMARK
    genes = [str(10_000 + i) for i in range(n_genes)]
    targets = {"TARGETA": [f"SHRNA_A{i}" for i in range(3)],
               "TARGETB": [f"SHRNA_B{i}" for i in range(3)]}
    drug_cells = {"drugone": {"MCF7": ["CPD_1_1", "CPD_1_2"], "PC3": ["CPD_1_3"]},
                  "drugtwo": {"MCF7": ["CPD_2_1"], "A375": ["CPD_2_2"]},
                  "drugthree": {"MCF7": ["CPD_3_1"], "PC3": ["CPD_3_2"]}}
    drug_targets = {"drugone": "TARGETA", "drugtwo": "TARGETB", "drugthree": "TARGETC"}
    shrna_ids = [sig for ids in targets.values() for sig in ids]
    compound_ids = [sig for cells in drug_cells.values() for ids in cells.values()
                    for sig in ids]

    # the source, with its genes physically stored in an order of their own
    stored_order = list(rng.permutation(genes))
    values = {sig: rng.standard_normal(n_genes) for sig in shrna_ids + compound_ids}
    rebuild = tmp_path / "rebuilt"
    extraction = tmp_path / "extraction"
    registry = tmp_path / "registry"
    rebuild.mkdir(), extraction.mkdir(), registry.mkdir()

    gctx = tmp_path / "source.gctx"
    with h5py.File(gctx, "w") as handle:
        matrix = np.vstack([[values[sig][genes.index(gene)] for gene in stored_order]
                            for sig in shrna_ids + compound_ids])
        handle.create_dataset("/0/DATA/0/matrix", data=matrix.astype(np.float32))
        handle.create_dataset("/0/META/ROW/id",
                              data=np.array(stored_order, dtype=h5py.string_dtype()))
        handle.create_dataset("/0/META/COL/id",
                              data=np.array(shrna_ids + compound_ids,
                                            dtype=h5py.string_dtype()))

    gene_info = tmp_path / "gene_info.txt.gz"
    pd.DataFrame({"pr_gene_id": genes + ["999999"],
                  "pr_is_lm": [1] * n_genes + [0],
                  "pr_gene_symbol": [f"SYM{i}" for i in range(n_genes + 1)]}).to_csv(
        gene_info, sep="\t", index=False, compression="gzip")

    shrna_siginfo = tmp_path / "shrna_siginfo.csv.gz"
    pd.DataFrame({"sig_id": shrna_ids,
                  "pert_iname": [t for t, ids in targets.items() for _ in ids]}).to_csv(
        shrna_siginfo, index=False, compression="gzip")

    compound_siginfo = tmp_path / "sig_info.txt.gz"
    pd.DataFrame({"sig_id": compound_ids,
                  "pert_iname": [drug for drug, cells in drug_cells.items()
                                 for ids in cells.values() for _ in ids],
                  "cell_id": [cell for cells in drug_cells.values()
                              for cell, ids in cells.items() for _ in ids]}).to_csv(
        compound_siginfo, sep="\t", index=False, compression="gzip")

    # the canonical artifact and its derived consensus, in the frozen landmark order
    signatures = np.vstack([values[sig] for sig in shrna_ids])
    np.savez_compressed(rebuild / "shrna_signatures.npz", sig_ids=np.array(shrna_ids),
                        signatures=signatures, gene_ids=np.array(genes),
                        membership=json.dumps(targets))
    directions = []
    for ids in targets.values():
        mean = np.vstack([values[sig] for sig in ids]).mean(axis=0)
        directions.append(mean / np.linalg.norm(mean))
    np.savez_compressed(rebuild / "shrna_consensus.npz", genes=np.array(list(targets)),
                        directions=np.vstack(directions), gene_ids=np.array(genes))

    matrices, cells = {}, {}
    for drug, by_cell in drug_cells.items():
        matrices[drug] = np.vstack([np.vstack([values[s] for s in ids]).mean(axis=0)
                                    for ids in by_cell.values()])
        cells[drug] = list(by_cell)
    _write_shard(rebuild / "shard_000.npz", matrices, cells)
    (rebuild / "landmark_gene_ids.json").write_text(json.dumps(genes))
    (rebuild / "rebuild_manifest.json").write_text(json.dumps(
        {"n_drugs": len(drug_cells), "gene_axis": "frozen landmark order"}))

    # the retained extraction: its columns follow `actual`, its declared axis says `genes`
    actual = list(rng.permutation(genes))
    while actual == genes:
        actual = list(rng.permutation(genes))
    np.savez_compressed(
        extraction / "lincs_subset.npz", sig_ids=np.array(compound_ids),
        signatures=np.vstack([[values[sig][genes.index(gene)] for gene in actual]
                              for sig in compound_ids]),
        gene_ids=np.array(genes))
    # deliberately NOT the same bytes as the hashed metadata: a copy beside the
    # extraction that disagrees must not be the one the gate groups by
    pd.DataFrame({"sig_id": compound_ids,
                  "pert_iname": ["drugone"] * len(compound_ids),
                  "cell_id": ["MCF7"] * len(compound_ids)}).to_csv(
        extraction / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t", index=False,
        compression="gzip")

    axis_map = tmp_path / "legacy_axis_map.json"
    order = [genes.index(gene) for gene in actual]
    axis_map.write_text(json.dumps(
        {"name": "legacy_axis_map", "schema_version": 1,
         "declared_index_of_each_actual_column": order,
         "map_sha256": _mod.sha256_text(_mod.canonical(order)),
         "legacy_sha256_observed": _mod.sha256_file(extraction / "lincs_subset.npz"),
         # the stage that builds a real map refuses to write one unless this is 0.0
         "recovered_from": {"n_cohort_signatures_confirmed": len(compound_ids),
                            "max_abs_difference_under_the_map": 0.0}}))

    # the cohort manifests, built by the registered builder rather than hand-written
    paired_src = tmp_path / "paired_source.json"
    compound_src = tmp_path / "compound_source.json"
    paired_src.write_text(json.dumps(
        [{"drug": d, "target": drug_targets[d]} for d in ("drugone", "drugtwo")]))
    compound_src.write_text(json.dumps(
        [{"drug": d, "target": drug_targets[d]} for d in drug_cells]))
    (registry / "INDEX.json").write_text(json.dumps(
        _reg.build(paired_src, compound_src, registry), indent=2, sort_keys=True))

    # the coverage measurement, by the registered measurement rather than a stub
    deposited = tmp_path / "deposited.json"
    deposited.write_text(json.dumps(
        [{"drug": "drugone", "target": "TARGETA", "proj_shrna": 1.0, "enrich_shrna": 0.1},
         {"drug": "drugtwo", "target": "TARGETB", "proj_shrna": 2.0, "enrich_shrna": 0.2},
         {"drug": "drugthree", "target": "TARGETC"}]))
    coverage = tmp_path / "coverage.json"
    _cov.main(["--deposited", str(deposited),
               "--consensus", str(rebuild / "shrna_consensus.npz"),
               "--eligible", str(registry / "shrna_eligible_targets.json"),
               "--excluded", str(registry / "shrna_excluded_records.json"),
               "--compound-cohort", str(registry / "cohort_812_compound.json"),
               "--out", str(coverage)])

    # the expected hashes, generated from the files and held outside them
    frozen = tmp_path / "expected_identities.json"
    _freeze.main(["--repo", str(Path(_mod.__file__).resolve().parents[1]),
                  "--rebuilt", str(rebuild), "--registry", str(registry),
                  "--coverage", str(coverage), "--axis-map", str(axis_map),
                  "--deposited", str(deposited),
                  "--legacy", str(extraction / "lincs_subset.npz"),
                  "--frozen-out", str(frozen),
                  "--manifest-out", str(tmp_path / "implementation_manifest.json")])

    output = tmp_path / "out"
    return ["--rebuilt", str(rebuild), "--extraction", str(extraction),
            "--cohort-manifest", str(registry / "cohort_812_compound.json"),
            "--paired-manifest", str(registry / "cohort_795_shrna_paired.json"),
            "--eligible-targets", str(registry / "shrna_eligible_targets.json"),
            "--excluded-records", str(registry / "shrna_excluded_records.json"),
            "--axis-map", str(axis_map), "--coverage-report", str(coverage),
            "--deposited-records", str(deposited),
            "--frozen-identities", str(frozen),
            "--cohort-scope", "compound", "--gctx", str(gctx),
            "--shrna-siginfo", str(shrna_siginfo),
            "--compound-siginfo", str(compound_siginfo),
            "--gene-info", str(gene_info), "--output", str(output)], rebuild, output


def test_the_whole_gate_passes_on_consistent_inputs_and_records_its_provenance(tmp_path):
    argv, _, output = _gate_fixture(tmp_path)

    result = _mod.main(argv)
    assert result["gate"] == "reconstruction"
    assert list(result["predicates_evaluated"]) == list(_mod.RELEASE_PREDICATES)
    assert all(result["predicate_held"].values())
    written = json.loads(_only_report(output).read_text())
    assert written["release_authorized"]
    assert written["production_authorized"]
    assert written["legacy_custody_intact"]
    assert written["gate_version"] == _mod.GATE_VERSION
    # every input that decides membership, aggregation or coordinates is hashed
    assert set(written["input_sha256"]) == {
        "compound_cohort_manifest", "shrna_paired_manifest", "eligible_targets_manifest",
        "excluded_records_manifest", "axis_map", "coverage_report",
        "deposited_records", "frozen_identities",
        "gctx", "shrna_siginfo", "compound_siginfo", "gene_info",
        "landmark_gene_ids", "shrna_signatures.npz", "shrna_consensus.npz",
        "rebuild_manifest.json", "shard_000.npz", "lincs_subset.npz"}
    assert written["gate_code_sha256"]
    # the legacy artifact is right in both directions, which a bare disagreement is not
    legacy = written["legacy_artifact_integrity"]
    assert legacy["under_the_map"]["agrees_with_the_source"]
    assert not legacy["under_the_declared_labels"]["agrees_with_the_source"]
    assert legacy["map_is_a_total_bijection"]
    assert legacy["map_round_trips_to_the_identity"]


def test_the_whole_gate_fails_on_a_canonical_artifact_whose_columns_were_permuted(tmp_path):
    argv, rebuild, output = _gate_fixture(tmp_path)

    # Deviation 11, through the command line: the declared axis is untouched and
    # every column of the canonical artifact moves
    with np.load(rebuild / "shrna_signatures.npz", allow_pickle=True) as data:
        contents = {key: data[key] for key in data.files}
    rng = np.random.default_rng()
    scramble = rng.permutation(contents["signatures"].shape[1])
    while (scramble == np.arange(len(scramble))).all():
        scramble = rng.permutation(len(scramble))
    contents["signatures"] = contents["signatures"][:, scramble]
    np.savez_compressed(rebuild / "shrna_signatures.npz", **contents)

    with pytest.raises(AssertionError, match="the gate failed"):
        _mod.main(argv)

    written = json.loads(_only_report(output).read_text())
    assert written["gate"].startswith("failed")
    assert not written["release_authorized"]
    shrna = written["shrna_reconstruction_on_the_eligible_set"]
    assert shrna["n_signature_failures"] == 6
    # the predicates that do not depend on the corrupted artifact still ran and are recorded
    assert written["predicate_held"]["compound_source_reconstruction"]
    assert written["predicate_held"]["cohort_identity"]
    assert written["predicate_held"]["legacy_artifact_integrity"]


def test_a_legacy_artifact_that_agrees_under_its_declared_labels_fails_custody(tmp_path):
    # bare disagreement is too weak a predicate, and so is bare agreement: a file that
    # reads correctly under its declared axis is not the file the deposited numbers
    # came from, so custody of the deposited analysis has been lost
    argv, _, output = _gate_fixture(tmp_path)
    extraction = Path(argv[argv.index("--extraction") + 1])
    axis_map = Path(argv[argv.index("--axis-map") + 1])
    mapping = json.loads(axis_map.read_text())
    order = mapping["declared_index_of_each_actual_column"]

    with np.load(extraction / "lincs_subset.npz", allow_pickle=True) as data:
        contents = {key: data[key] for key in data.files}
    inverse = np.empty(len(order), dtype=int)
    for actual, declared in enumerate(order):
        inverse[declared] = actual
    contents["signatures"] = contents["signatures"][:, inverse]   # now honestly labeled
    np.savez_compressed(extraction / "lincs_subset.npz", **contents)
    # the relabeled file is registered as-is, so the byte check passes and the two-sided
    # comparison is what decides: a separate test covers the bytes changing underneath
    frozen_path = Path(argv[argv.index("--frozen-identities") + 1])
    frozen = json.loads(frozen_path.read_text())
    frozen["legacy_extraction_sha256"] = _mod.sha256_file(extraction / "lincs_subset.npz")
    frozen_path.write_text(json.dumps(frozen))

    with pytest.raises(AssertionError, match="the gate failed"):
        _mod.main(argv)
    written = json.loads(_only_report(output).read_text())
    assert not written["legacy_custody_intact"]
    assert written["production_authorized"]        # production is unaffected
    legacy = written["legacy_artifact_integrity"]
    assert legacy["under_the_declared_labels"]["agrees_with_the_source"]
    assert "no longer carries exactly the recorded defect" in legacy["reading"]


def test_a_legacy_map_that_is_not_a_bijection_is_refused(tmp_path):
    argv, _, _ = _gate_fixture(tmp_path)
    axis_map = Path(argv[argv.index("--axis-map") + 1])
    mapping = json.loads(axis_map.read_text())
    mapping["declared_index_of_each_actual_column"][1] = \
        mapping["declared_index_of_each_actual_column"][0]
    order = mapping["declared_index_of_each_actual_column"]
    mapping["map_sha256"] = _mod.sha256_text(_mod.canonical(order))
    mapping["pinned_map_sha256"] = mapping["map_sha256"]
    axis_map.write_text(json.dumps(mapping))

    with pytest.raises(AssertionError, match="the gate failed"):
        _mod.main(argv)


def test_a_map_recording_an_inexact_raw_reproduction_is_refused(tmp_path):
    # exactness is true of raw signatures and false of aggregated means, so the gate
    # checks the map artifact's own raw-level record rather than asserting exactness
    # where even the production rebuild sits at 3.12e-06
    argv, _, output = _gate_fixture(tmp_path)
    axis_map = Path(argv[argv.index("--axis-map") + 1])
    frozen_path = Path(argv[argv.index("--frozen-identities") + 1])
    mapping = json.loads(axis_map.read_text())
    mapping["recovered_from"]["max_abs_difference_under_the_map"] = 1e-9
    axis_map.write_text(json.dumps(mapping))
    frozen = json.loads(frozen_path.read_text())
    frozen["axis_map_file_sha256"] = _mod.sha256_file(axis_map)
    frozen_path.write_text(json.dumps(frozen))

    with pytest.raises(AssertionError, match="the gate failed"):
        _mod.main(argv)
    written = json.loads(_only_report(output).read_text())
    assert "exact raw-signature reproduction" in \
        written["structural_failures"]["legacy_artifact_integrity"]
    assert written["production_authorized"]


def test_a_legacy_map_whose_hash_is_not_the_pinned_one_is_refused(tmp_path):
    argv, _, _ = _gate_fixture(tmp_path)
    axis_map = Path(argv[argv.index("--axis-map") + 1])
    mapping = json.loads(axis_map.read_text())
    mapping["pinned_map_sha256"] = "0" * 64
    axis_map.write_text(json.dumps(mapping))

    with pytest.raises(AssertionError, match="the gate failed"):
        _mod.main(argv)


def test_the_shrna_universe_is_the_pinned_set_and_not_the_cohort_targets(tmp_path):
    # the v1 defect: the gate took the compound cohort's target symbols as the shRNA
    # universe, so a drug whose target has no shRNA data refused the whole half
    argv, _, output = _gate_fixture(tmp_path)
    eligible = Path(argv[argv.index("--eligible-targets") + 1])
    manifest = json.loads(eligible.read_text())

    # TARGETC belongs to the excluded record and must not be required of the rebuild
    assert "TARGETC" not in manifest["targets"]
    result = _mod.main(argv)
    assert result["predicate_held"]["shrna_reconstruction_on_the_eligible_set"]
    assert result["shrna_reconstruction_on_the_eligible_set"]["n_targets"] == 2


def _under_optimization(body):
    """Run a probe against the gate module with assertions removed."""
    import subprocess

    gate = Path(__file__).resolve().parents[1] / "experiments" / "03e_reconstruction_gate.py"
    probe = ("import importlib.util, numpy as np;"
             f"spec = importlib.util.spec_from_file_location('g', {str(gate)!r});"
             "g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g);"
             + body)
    return subprocess.run(["python", "-O", "-c", probe], capture_output=True, text=True)


def test_unequal_cell_sets_are_refused_even_under_optimization():
    # with `assert` removed, align() used to drop the unmatched cell line silently
    finished = _under_optimization(
        "genes = ['1', '2'];"
        "left = ({'d': np.ones((1, 2))}, {'d': ['C1']}, genes);"
        "right = ({'d': np.ones((2, 2))}, {'d': ['C1', 'C2']}, genes);"
        "g.compare(left, right, cohort=['d'])")
    assert finished.returncode != 0
    assert "declared cell lines" in finished.stderr or "identifier sets differ" in finished.stderr


def test_unequal_membership_is_refused_even_under_optimization():
    finished = _under_optimization(
        "genes = ['1', '2'];"
        "v = {s: np.ones(2) for s in ('A', 'B', 'C')};"
        "left = (dict(v), {'T': ['A', 'B']}, genes);"
        "right = (dict(v), {'T': ['A', 'B', 'C']}, genes);"
        "g.compare_shrna(left, right)")
    assert finished.returncode != 0
    assert "group different signatures" in finished.stderr


def test_a_signature_array_contradicting_its_identifiers_is_refused(tmp_path):
    rng = np.random.default_rng()
    genes = [str(i) for i in range(3)]
    ids = ["SIG1", "SIG2", "SIG3"]
    membership = json.dumps({"TARGET0": ids})

    # an unlabeled fourth row would never be read, so its contents never compared
    np.savez_compressed(tmp_path / "shrna_signatures.npz", sig_ids=np.array(ids),
                        signatures=rng.standard_normal((4, 3)), gene_ids=np.array(genes),
                        membership=membership)
    with pytest.raises(AssertionError, match=r"holds a \(4, 3\) array against 3 signature"):
        _mod.load_shrna_rebuild(tmp_path)

    # an extra dimension survives vstack by broadcasting
    np.savez_compressed(tmp_path / "shrna_signatures.npz", sig_ids=np.array(ids),
                        signatures=rng.standard_normal((3, 2, 3)), gene_ids=np.array(genes),
                        membership=membership)
    with pytest.raises(AssertionError, match=r"holds a \(3, 2, 3\) array"):
        _mod.load_shrna_rebuild(tmp_path)

    np.savez_compressed(tmp_path / "shrna_signatures.npz", sig_ids=np.array(ids),
                        signatures=np.full((3, 3), np.nan), gene_ids=np.array(genes),
                        membership=membership)
    with pytest.raises(AssertionError, match="nonfinite values"):
        _mod.load_shrna_rebuild(tmp_path)


def test_the_retained_extraction_is_refused_when_it_repeats_a_signature(tmp_path):
    import pandas as pd

    rng = np.random.default_rng()
    genes = [str(i) for i in range(4)]
    # the first occurrence is corrupted and the later one correct, which a position
    # lookup would resolve in favour of the later one
    np.savez_compressed(tmp_path / "lincs_subset.npz",
                        sig_ids=np.array(["CPD1", "CPD1", "CPD2"]),
                        signatures=rng.standard_normal((3, 4)),
                        gene_ids=np.array(genes))
    pd.DataFrame({"sig_id": ["CPD1", "CPD2"], "pert_iname": ["drug", "drug"],
                  "cell_id": ["MCF7", "PC3"]}).to_csv(
        tmp_path / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t", index=False,
        compression="gzip")

    with pytest.raises(AssertionError, match="repeats a signature"):
        _mod.load_extraction(tmp_path, cohort=["drug"])


def test_the_retained_extraction_uses_the_metadata_the_gate_hashed(tmp_path):
    # two copies of the metadata existed: one hashed on the command line, one sitting
    # beside the extraction and read instead
    import pandas as pd

    rng = np.random.default_rng()
    genes = [str(i) for i in range(4)]
    np.savez_compressed(tmp_path / "lincs_subset.npz", sig_ids=np.array(["CPD1", "CPD2"]),
                        signatures=rng.standard_normal((2, 4)), gene_ids=np.array(genes))
    beside = tmp_path / "GSE92742_Broad_LINCS_sig_info.txt.gz"
    pd.DataFrame({"sig_id": ["CPD1", "CPD2"], "pert_iname": ["drug", "drug"],
                  "cell_id": ["MCF7", "MCF7"]}).to_csv(beside, sep="\t", index=False,
                                                       compression="gzip")
    hashed = tmp_path / "hashed_sig_info.txt.gz"
    pd.DataFrame({"sig_id": ["CPD1", "CPD2"], "pert_iname": ["drug", "drug"],
                  "cell_id": ["MCF7", "PC3"]}).to_csv(hashed, sep="\t", index=False,
                                                      compression="gzip")

    _, beside_cells, _ = _mod.load_extraction(tmp_path, cohort=["drug"])
    _, hashed_cells, _ = _mod.load_extraction(tmp_path, cohort=["drug"], siginfo_path=hashed)
    assert beside_cells["drug"] == ["MCF7"]
    assert sorted(hashed_cells["drug"]) == ["MCF7", "PC3"]


def test_a_manifest_that_disagrees_with_its_own_declared_identity_is_refused(tmp_path):
    # B8 promises the counts and the drug and target hashes, so the gate enforces them
    # rather than carrying them as decoration
    argv, _, _ = _gate_fixture(tmp_path)
    manifest_path = Path(argv[argv.index("--cohort-manifest") + 1])
    frozen_path = Path(argv[argv.index("--frozen-identities") + 1])

    for field, value in (("n_records", 99), ("n_unique_drugs", 99),
                         ("n_unique_targets", 99), ("drug_sha256", "0" * 64),
                         ("target_sha256", "0" * 64), ("pair_sha256", "0" * 64)):
        manifest = json.loads(manifest_path.read_text())
        manifest[field] = value
        manifest_path.write_text(json.dumps(manifest))
        frozen = json.loads(frozen_path.read_text())
        frozen["manifest_sha256"]["cohort_812_compound.json"] = \
            _mod.sha256_file(manifest_path)
        frozen_path.write_text(json.dumps(frozen))

        with pytest.raises(AssertionError, match="cohort_812_compound"):
            _mod.main(argv)


def test_an_eligible_manifest_repeating_a_target_or_miscounting_is_refused(tmp_path):
    argv, _, _ = _gate_fixture(tmp_path)
    path = Path(argv[argv.index("--eligible-targets") + 1])
    frozen_path = Path(argv[argv.index("--frozen-identities") + 1])
    original = json.loads(path.read_text())

    for mutate in (lambda m: m["targets"].append(m["targets"][0]),
                   lambda m: m.update(n_targets=99)):
        manifest = json.loads(json.dumps(original))
        mutate(manifest)
        path.write_text(json.dumps(manifest))
        frozen = json.loads(frozen_path.read_text())
        frozen["manifest_sha256"]["shrna_eligible_targets.json"] = _mod.sha256_file(path)
        frozen_path.write_text(json.dumps(frozen))

        with pytest.raises(AssertionError, match="shrna_eligible_targets"):
            _mod.main(argv)


def test_a_manifest_whose_file_hash_is_not_the_frozen_one_is_refused(tmp_path):
    # the manifests are pinned outside themselves, so a regenerated manifest that is
    # internally consistent still does not pass
    argv, _, _ = _gate_fixture(tmp_path)
    frozen_path = Path(argv[argv.index("--frozen-identities") + 1])
    frozen = json.loads(frozen_path.read_text())
    frozen["manifest_sha256"]["cohort_795_shrna_paired.json"] = "0" * 64
    frozen_path.write_text(json.dumps(frozen))

    with pytest.raises(AssertionError, match="cohort_795_shrna_paired.json: sha256"):
        _mod.main(argv)


def test_a_coverage_report_or_deposited_file_that_is_not_the_frozen_one_is_refused(tmp_path):
    argv, _, output = _gate_fixture(tmp_path)
    frozen_path = Path(argv[argv.index("--frozen-identities") + 1])

    for field in ("coverage_report_sha256", "deposited_sha256"):
        frozen = json.loads(frozen_path.read_text())
        frozen[field] = "0" * 64
        frozen_path.write_text(json.dumps(frozen))
        with pytest.raises(AssertionError, match="the gate failed"):
            _mod.main(argv)
        written = sorted(output.glob("gate_v*/reconstruction_gate.json"))[-1]
        assert not json.loads(written.read_text())[
            "predicate_held"]["shrna_eligibility_partition"]


def test_a_coverage_report_missing_a_registered_equality_is_refused(tmp_path):
    argv, _, _ = _gate_fixture(tmp_path)
    coverage = Path(argv[argv.index("--coverage-report") + 1])
    frozen_path = Path(argv[argv.index("--frozen-identities") + 1])
    report = json.loads(coverage.read_text())
    del report["set_equalities"]["partition_is_exact"]
    coverage.write_text(json.dumps(report))
    frozen = json.loads(frozen_path.read_text())
    frozen["coverage_report_sha256"] = _mod.sha256_file(coverage)
    frozen_path.write_text(json.dumps(frozen))

    with pytest.raises(AssertionError, match="the gate failed"):
        _mod.main(argv)


def test_two_runs_write_separate_directories_and_neither_overwrites_the_other(tmp_path):
    # the report filename is fixed, so two scopes sharing a directory would leave only
    # the second report and no trace that the first existed
    argv, _, output = _gate_fixture(tmp_path)
    first = _mod.main(argv)

    scoped = list(argv)
    scoped[scoped.index("--cohort-scope") + 1] = "shrna_paired"
    second = _mod.main(scoped)

    reports = sorted(output.glob("gate_v*/reconstruction_gate.json"))
    assert len(reports) == 2
    scopes = {json.loads(r.read_text())["cohort_scope"] for r in reports}
    assert scopes == {"compound", "shrna_paired"}
    assert first["cohort"]["n_drugs"] == 3 and second["cohort"]["n_drugs"] == 2
    assert all(json.loads(r.read_text())["release_authorized"] for r in reports)
