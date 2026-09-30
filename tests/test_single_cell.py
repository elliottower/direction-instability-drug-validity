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
