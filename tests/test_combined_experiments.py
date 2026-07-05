"""Tests for combined paper experiment utility functions.

Tests the core computation functions from experiments 12, 13, 14
using synthetic data with known answers.
"""
import numpy as np
import pytest
from importlib import import_module

exp12 = import_module("experiments.12_perturbseq_cosine_transport")
exp13 = import_module("experiments.13_lincs_metric_benchmark")
exp14 = import_module("experiments.14_jump_cp_moa_stratification")


class TestCosineDirectionInstability:
    def test_identical_vectors_give_zero(self):
        v = np.array([1.0, 2.0, 3.0])
        assert exp12.cosine_direction_instability(v, v) == pytest.approx(0.0, abs=1e-10)

    def test_orthogonal_vectors_give_one(self):
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([0.0, 1.0, 0.0])
        assert exp12.cosine_direction_instability(v1, v2) == pytest.approx(1.0, abs=1e-10)

    def test_antiparallel_vectors_give_zero(self):
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([-1.0, 0.0, 0.0])
        assert exp12.cosine_direction_instability(v1, v2) == pytest.approx(0.0, abs=1e-10)

    def test_magnitude_does_not_matter(self):
        v1 = np.array([1.0, 1.0, 0.0])
        v2 = np.array([1.0, 0.0, 0.0])
        d1 = exp12.cosine_direction_instability(v1, v2)
        d2 = exp12.cosine_direction_instability(v1 * 100, v2 * 0.001)
        assert d1 == pytest.approx(d2, abs=1e-10)

    def test_45_degree_angle(self):
        v1 = np.array([1.0, 0.0])
        v2 = np.array([1.0, 1.0])
        expected = 1.0 - np.cos(np.pi / 4)
        assert exp12.cosine_direction_instability(v1, v2) == pytest.approx(expected, abs=1e-10)

    def test_zero_vector_gives_nan(self):
        v1 = np.array([1.0, 0.0])
        v2 = np.array([0.0, 0.0])
        assert np.isnan(exp12.cosine_direction_instability(v1, v2))


class TestComputeSameGeneDi:
    def test_identical_pairs_give_zero(self):
        sigs = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        di = exp12.compute_same_gene_di(sigs, sigs)
        assert di[0] == pytest.approx(0.0, abs=1e-10)
        assert di[1] == pytest.approx(0.0, abs=1e-10)

    def test_orthogonal_pairs_give_one(self):
        k562 = np.array([[1.0, 0.0], [0.0, 1.0]])
        rpe1 = np.array([[0.0, 1.0], [1.0, 0.0]])
        di = exp12.compute_same_gene_di(k562, rpe1)
        assert di[0] == pytest.approx(1.0, abs=1e-10)
        assert di[1] == pytest.approx(1.0, abs=1e-10)


class TestCohensD:
    def test_identical_groups_give_zero(self):
        g = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        assert exp12.cohens_d(g, g) == pytest.approx(0.0, abs=1e-10)

    def test_known_separation(self):
        np.random.seed(None)
        g1 = np.random.randn(10000)
        g2 = np.random.randn(10000) + 1.0
        d = exp12.cohens_d(g1, g2)
        assert d == pytest.approx(-1.0, abs=0.05)

    def test_sign_convention(self):
        g1 = np.array([10.0, 11.0, 9.0, 10.5])
        g2 = np.array([0.0, 1.0, -1.0, 0.5])
        d = exp12.cohens_d(g1, g2)
        assert d > 0

    def test_matches_across_scripts(self):
        g1 = np.random.randn(50)
        g2 = np.random.randn(50) + 1.0
        d12 = exp12.cohens_d(g1, g2)
        d14 = exp14.cohens_d(g1, g2)
        assert d12 == pytest.approx(d14, abs=1e-10)


class TestAurocForPredictor:
    def test_perfect_negative_predictor(self):
        predictor = np.array([10.0, 9.0, 8.0, 1.0, 0.5, 0.1])
        outcome = np.array([0.0, 0.1, 0.2, 0.8, 0.9, 1.0])
        auc = exp13.auroc_for_predictor(predictor, outcome, 0.5, higher_is_worse=True)
        assert auc == pytest.approx(1.0, abs=1e-10)

    def test_perfect_positive_predictor(self):
        predictor = np.array([0.1, 0.2, 0.3, 0.8, 0.9, 1.0])
        outcome = np.array([0.0, 0.1, 0.2, 0.8, 0.9, 1.0])
        auc = exp13.auroc_for_predictor(predictor, outcome, 0.5, higher_is_worse=False)
        assert auc == pytest.approx(1.0, abs=1e-10)

    def test_higher_is_worse_flips_direction(self):
        predictor = np.array([0.1, 0.5, 0.9, 0.2, 0.6, 1.0])
        outcome = np.array([0.8, 0.7, 0.1, 0.9, 0.3, 0.05])
        auc_worse = exp13.auroc_for_predictor(predictor, outcome, 0.5, higher_is_worse=True)
        auc_better = exp13.auroc_for_predictor(predictor, outcome, 0.5, higher_is_worse=False)
        assert auc_worse == pytest.approx(1.0 - auc_better, abs=1e-10)

    def test_all_same_label_gives_nan(self):
        predictor = np.array([1.0, 2.0, 3.0])
        outcome = np.array([1.0, 1.0, 1.0])
        assert np.isnan(exp13.auroc_for_predictor(predictor, outcome, 0.5))


class TestResidualize:
    def test_residuals_uncorrelated_with_covariate(self):
        np.random.seed(None)
        n = 500
        cov = np.random.randn(n)
        x = 2.0 * cov + np.random.randn(n) * 0.5
        resid = exp13.residualize(x, cov)
        rho, _ = __import__("scipy").stats.spearmanr(resid, cov)
        assert abs(rho) < 0.15

    def test_residualize_constant_covariate_preserves_mean(self):
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        cov = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        resid = exp13.residualize(x, cov)
        assert np.mean(resid) == pytest.approx(0.0, abs=1e-10)


class TestClassifyMoa:
    def test_hdac_is_machinery(self):
        assert exp14.classify_moa("HDAC inhibitor") == "machinery"

    def test_dopamine_agonist_is_receptor(self):
        assert exp14.classify_moa("dopamine receptor agonist") == "receptor"

    def test_egfr_is_kinase(self):
        assert exp14.classify_moa("EGFR inhibitor") == "kinase"

    def test_nan_returns_none(self):
        assert exp14.classify_moa(float("nan")) is None

    def test_empty_returns_none(self):
        assert exp14.classify_moa("") is None

    def test_unknown_moa_returns_other(self):
        assert exp14.classify_moa("bacterial cell wall synthesis inhibitor") == "other"

    def test_case_insensitive(self):
        assert exp14.classify_moa("hdac inhibitor") == "machinery"
        assert exp14.classify_moa("DOPAMINE RECEPTOR AGONIST") == "receptor"

    def test_pipe_separated_moa_first_match_wins(self):
        assert exp14.classify_moa("HDAC inhibitor|serotonin receptor antagonist") == "machinery"

    def test_machinery_priority_over_receptor(self):
        assert exp14.classify_moa("topoisomerase inhibitor") == "machinery"
        assert exp14.classify_moa("proteasome inhibitor") == "machinery"
        assert exp14.classify_moa("tubulin inhibitor") == "machinery"
