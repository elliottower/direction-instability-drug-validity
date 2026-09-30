import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import spearmanr

from geometry.inference import (cluster_bootstrap, comparison_reading, is_practically_null,
                                own_target_percentile, percentile_interval, permutation_reading,
                                rank_partial_correlation, target_balanced_mean,
                                target_balanced_spearman, two_sided_permutation_p,
                                unique_target_permutations)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
_registered = __import__("03c_h3_sensitivity")


def test_partial_correlation_matches_the_registered_implementation():
    rng = np.random.default_rng()
    for _ in range(50):
        n = int(rng.integers(30, 200))
        y, x, c1, c2 = (rng.standard_normal(n) for _ in range(4))
        mine = rank_partial_correlation(y, x, (c1, c2))
        registered = _registered.partial_spearman(y, x, (c1, c2))
        assert mine == pytest.approx(registered, abs=1e-12)


def test_partial_correlation_refuses_a_collinear_covariate():
    rng = np.random.default_rng()
    x = rng.standard_normal(80)
    with pytest.raises(AssertionError):
        rank_partial_correlation(x, rng.standard_normal(80), (np.exp(x),))


def test_partial_correlation_removes_a_shared_covariate():
    rng = np.random.default_rng()
    inflated = []
    for _ in range(200):
        z = rng.standard_normal(150)
        y = z + 0.4 * rng.standard_normal(150)
        x = z + 0.4 * rng.standard_normal(150)
        raw = spearmanr(y, x).statistic
        adjusted = rank_partial_correlation(y, x, (z,))
        inflated.append(raw - adjusted)
    assert np.mean(inflated) > 0.4


def test_own_target_percentile_is_one_when_the_own_target_always_wins():
    rng = np.random.default_rng()
    alignment = rng.uniform(0, 1, size=(40, 12))
    own = rng.integers(0, 12, 40)
    alignment[np.arange(40), own] = 2.0
    assert own_target_percentile(alignment, own) == pytest.approx(np.ones(40))


def test_own_target_percentile_is_one_half_under_ties_and_under_noise():
    flat = np.ones((25, 8))
    own = np.zeros(25, dtype=int)
    assert own_target_percentile(flat, own) == pytest.approx(np.full(25, 0.5))

    rng = np.random.default_rng()
    q = []
    for _ in range(300):
        alignment = rng.standard_normal((30, 15))
        q.append(own_target_percentile(alignment, rng.integers(0, 15, 30)).mean())
    assert np.mean(q) == pytest.approx(0.5, abs=0.01)


def test_target_balanced_statistics_ignore_how_many_drugs_a_target_has():
    rng = np.random.default_rng()
    targets = np.array(["a"] * 20 + ["b"] + ["c"])
    values = np.concatenate([np.full(20, 1.0), [0.0], [2.0]])
    assert target_balanced_mean(values, targets) == pytest.approx(1.0)

    x, y = rng.standard_normal(3), rng.standard_normal(3)
    singletons = np.array(["p", "q", "r"])
    assert target_balanced_spearman(x, y, singletons) == pytest.approx(
        spearmanr(x, y).statistic)


def test_cluster_bootstrap_is_wider_than_resampling_drugs_when_targets_repeat():
    rng = np.random.default_rng()
    n_targets, per_target = 30, 8
    targets = np.repeat(np.arange(n_targets), per_target)
    # every drug of a target carries that target's value, which is what sharing a
    # reference direction does to the data
    values = np.repeat(rng.standard_normal(n_targets), per_target)

    cluster = cluster_bootstrap(targets, lambda idx: values[idx].mean(), 2000, seed=1)
    naive = np.array([values[rng.integers(0, len(values), len(values))].mean()
                      for _ in range(2000)])
    cluster_width = np.subtract(*reversed(percentile_interval(cluster[:, 0])))
    naive_width = np.subtract(*reversed(percentile_interval(naive)))
    assert cluster_width > 2 * naive_width


def test_cluster_bootstrap_carries_whole_targets_and_is_reproducible():
    targets = np.array([0, 0, 0, 1, 2, 2])
    sizes = {0: 3, 1: 1, 2: 2}
    seen = []

    def record(idx):
        seen.append(idx.copy())
        return 0.0

    cluster_bootstrap(targets, record, 200, seed=7)
    assert len(seen) == 200
    for idx in seen:
        counts = np.bincount(targets[idx], minlength=3)
        # a target is drawn whole or not at all, so its rows arrive in multiples
        # of its size, and the replicate holds three targets' worth of rows
        for target, size in sizes.items():
            assert counts[target] % size == 0
        assert sum(counts[t] // size for t, size in sizes.items()) == 3
    assert any(np.bincount(targets[idx], minlength=3)[0] == 0 for idx in seen)

    first = cluster_bootstrap(targets, lambda idx: float(idx.sum()), 50, seed=11)
    again = cluster_bootstrap(targets, lambda idx: float(idx.sum()), 50, seed=11)
    assert first == pytest.approx(again)


def test_unique_target_permutation_keeps_shared_targets_together():
    rng = np.random.default_rng()
    targets = np.array(["TOP2A"] * 12 + ["CDK1"] * 11 + ["PLK1"] * 2 + ["WEE1"])
    for assigned in unique_target_permutations(targets, 200, seed=int(rng.integers(1e6))):
        for target in np.unique(targets):
            rows = targets == target
            assert len(set(assigned[rows])) == 1
        # a permutation of the unique targets uses each exactly once
        _, first_rows = np.unique(targets, return_index=True)
        assert sorted(assigned[first_rows]) == list(range(len(first_rows)))


def test_unique_target_permutation_is_not_a_record_permutation():
    targets = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    splits = 0
    for assigned in unique_target_permutations(targets, 500, seed=3):
        splits += len(set(assigned[:4])) > 1
    assert splits == 0


def test_comparison_rule_needs_the_paired_interval_for_every_reading():
    wide = np.random.default_rng().normal(0, 0.5, 4000)          # spans zero
    tight_negative = np.full(4000, -0.2) + np.random.default_rng().normal(0, 0.01, 4000)
    alt_tight = np.full(4000, 0.05) + np.random.default_rng().normal(0, 0.01, 4000)

    # the failure the rule exists to prevent: a big drop with an interval that
    # spans zero must not read as attenuated
    assert comparison_reading(0.27, 0.01, wide, wide) == "inconclusive"
    assert comparison_reading(0.27, 0.05, tight_negative, alt_tight) == "attenuated"
    assert comparison_reading(0.27, -0.30, tight_negative, alt_tight) == "reversed"

    retained_alt = np.full(4000, 0.25) + np.random.default_rng().normal(0, 0.01, 4000)
    assert comparison_reading(0.27, 0.25, wide, retained_alt) == "retained"
    assert comparison_reading(0.10, 0.01, tight_negative, alt_tight) == (
        "reference too small to attenuate")


def test_practical_null_needs_a_tight_interval_around_zero():
    rng = np.random.default_rng()
    assert is_practically_null(rng.normal(0.0, 0.02, 5000))
    assert not is_practically_null(rng.normal(0.0, 0.5, 5000))
    assert not is_practically_null(rng.normal(0.30, 0.02, 5000))


def test_permutation_readings_separate_assignment_dependence_from_resistance():
    rng = np.random.default_rng()
    centred_null = rng.normal(0.0, 0.05, 5000)
    assert permutation_reading(0.40, centred_null) == "assignment-dependent"

    shifted_null = rng.normal(0.26, 0.05, 5000)
    assert permutation_reading(0.27, shifted_null) == "permutation-resistant"
    assert permutation_reading(0.12, rng.normal(0.05, 0.2, 5000)) == "indeterminate"


def test_permutation_p_counts_both_tails():
    null = np.linspace(-1, 1, 1001)
    assert two_sided_permutation_p(0.0, null) == pytest.approx(1.0, abs=1e-3)
    assert two_sided_permutation_p(2.0, null) == pytest.approx(1 / 1002, abs=1e-6)


def test_percentile_table_matches_the_direct_percentile_for_every_assignment():
    from geometry.inference import percentile_table

    rng = np.random.default_rng()
    for _ in range(30):
        alignment = rng.standard_normal((25, 9))
        if rng.random() < 0.3:                       # force ties
            alignment = np.round(alignment, 1)
        table = percentile_table(alignment)
        for column in range(9):
            own = np.full(25, column)
            assert table[np.arange(25), own] == pytest.approx(
                own_target_percentile(alignment, own))
