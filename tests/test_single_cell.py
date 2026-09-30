import numpy as np
import pytest

from geometry.single_cell import (cosine, group_design, pseudobulk, representation_audit,
                                  split_half_reliability, split_half_reliability_by_unit)


def test_representation_audit_separates_counts_from_normalized_values():
    rng = np.random.default_rng()
    counts = rng.poisson(1.5, size=(500, 30))
    audit = representation_audit(counts)
    assert audit["fraction_integer_valued"] == pytest.approx(1.0)
    assert audit["fraction_negative"] == 0.0
    assert audit["min"] >= 0

    z_scores = rng.standard_normal((500, 30))
    audit = representation_audit(z_scores)
    assert audit["fraction_integer_valued"] < 0.01
    assert audit["fraction_negative"] == pytest.approx(0.5, abs=0.1)


def test_group_design_counts_cells_and_batches_per_target():
    labels = np.array(["non-targeting"] * 5 + ["TP53"] * 4 + ["MYC"] * 3)
    groups = np.array([1, 1, 2, 2, 3] + [1, 1, 2, 3] + [2, 2, 2])

    design = group_design(labels, groups)
    assert design["n_control_cells"] == 5
    assert design["n_targets"] == 2
    assert design["cells_per_target"] == {"TP53": 4, "MYC": 3}
    assert design["gem_groups_per_target"] == {"TP53": 3, "MYC": 1}


def test_pseudobulk_is_the_mean_over_a_target_cells():
    matrix = np.array([[1.0, 2.0], [3.0, 4.0], [10.0, 10.0]])
    labels = np.array(["A", "A", "B"])
    bulk = pseudobulk(matrix, labels)
    assert bulk["A"] == pytest.approx([2.0, 3.0])
    assert bulk["B"] == pytest.approx([10.0, 10.0])


def test_split_half_reliability_is_high_for_a_real_direction_and_low_for_noise():
    rng = np.random.default_rng()
    genes = 60
    truth = rng.standard_normal(genes)
    groups = np.repeat(np.arange(8), 25)

    signal = truth + rng.standard_normal((200, genes)) * 0.5
    noise = rng.standard_normal((200, genes))
    assert split_half_reliability(signal, groups, seed=int(rng.integers(1e6))) > 0.8
    assert abs(split_half_reliability(noise, groups, seed=int(rng.integers(1e6)))) < 0.5


def test_split_half_reliability_needs_two_batches():
    rng = np.random.default_rng()
    one_batch = np.zeros(40, dtype=int)
    assert np.isnan(split_half_reliability(rng.standard_normal((40, 10)), one_batch, seed=1))


def test_batch_confounded_direction_does_not_look_reliable():
    rng = np.random.default_rng()
    genes = 40
    # every batch carries its own direction, so no direction survives a batch split
    groups = np.repeat(np.arange(6), 30)
    matrix = np.vstack([np.tile(rng.standard_normal(genes), (30, 1)) for _ in range(6)])
    matrix += rng.standard_normal(matrix.shape) * 0.1

    scores = [split_half_reliability(matrix, groups, seed=s) for s in range(20)]
    assert np.median(scores) < 0.6


def test_hairpin_split_falls_back_to_halving_signatures():
    rng = np.random.default_rng()
    truth = rng.standard_normal(30)
    matrix = truth + rng.standard_normal((6, 30)) * 0.3

    one_hairpin = np.array(["TRCN1"] * 6)
    assert split_half_reliability_by_unit(matrix, one_hairpin, seed=2) > 0.7

    three_hairpins = np.array(["TRCN1", "TRCN1", "TRCN2", "TRCN2", "TRCN3", "TRCN3"])
    assert split_half_reliability_by_unit(matrix, three_hairpins, seed=2) > 0.7
    assert np.isnan(split_half_reliability_by_unit(matrix[:1], np.array(["TRCN1"]), seed=2))


def test_cosine_is_scale_free_and_signed():
    rng = np.random.default_rng()
    a = rng.standard_normal(20)
    assert cosine(a, 5 * a) == pytest.approx(1.0)
    assert cosine(a, -a) == pytest.approx(-1.0)
    assert np.isnan(cosine(a, np.zeros(20)))


def test_target_seed_follows_the_name_not_the_position():
    from geometry.single_cell import target_seed

    assert target_seed(7, "TP53") == target_seed(7, "TP53")
    assert target_seed(7, "TP53") != target_seed(8, "TP53")
    assert target_seed(7, "TP53") != target_seed(7, "MYC")
    # the same target in a different batch keeps its seed, which is what makes a
    # resumed run agree with a clean one
    first = [target_seed(7, t) for t in ["A", "B", "C"]]
    reordered = [target_seed(7, t) for t in ["C", "A", "B"]]
    assert sorted(first) == sorted(reordered)


def test_batches_are_deterministic_and_cover_every_target():
    from geometry.single_cell import batches

    targets = [f"T{i}" for i in range(250)]
    named = batches(targets, 100)
    assert [name for name, _ in named] == ["batch_0000", "batch_0001", "batch_0002"]
    assert sum(len(members) for _, members in named) == 250
    assert batches(list(reversed(targets)), 100) == named


def test_an_interrupted_run_resumes_to_the_same_result_as_a_clean_one():
    from geometry.single_cell import batches, merge_batches, split_half_reliability, target_seed

    rng = np.random.default_rng()
    genes, groups = 25, np.repeat(np.arange(6), 20)
    data = {f"T{i}": rng.standard_normal(genes) + rng.standard_normal((120, genes)) * 0.5
            for i in range(30)}

    def compute(targets):
        return {t: {"reliability": split_half_reliability(data[t], groups,
                                                          target_seed(99, t))} for t in targets}

    named = batches(list(data), 8)
    clean = merge_batches({name: compute(members) for name, members in named})

    # the first two batches survive an interruption; the rest are computed later
    written = {name: compute(members) for name, members in named[:2]}
    written.update({name: compute(members) for name, members in named[2:]})
    resumed = merge_batches(written)

    assert clean == resumed
    assert len(clean) == 30


def test_merge_refuses_a_target_written_by_two_batches():
    from geometry.single_cell import merge_batches

    with pytest.raises(AssertionError, match="conflicting"):
        merge_batches({"batch_0000": {"T1": {"reliability": 0.9}},
                       "batch_0001": {"T1": {"reliability": 0.2}}})
    with pytest.raises(AssertionError, match="identical"):
        merge_batches({"batch_0000": {"T1": {"reliability": 0.9}},
                       "batch_0001": {"T1": {"reliability": 0.9}}})
