"""Tests for Experiment 15: Positive control validation."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "experiments"))

from positive_control_validation import (
    is_nuclear_receptor,
    is_constitutive_machinery,
    is_narrow_machinery,
    compute_di,
    generate_spike_in,
    bootstrap_mean_diff,
)


class TestClassification:

    def test_nhr_glucocorticoid(self):
        assert is_nuclear_receptor("glucocorticoid receptor agonist")

    def test_nhr_estrogen(self):
        assert is_nuclear_receptor("estrogen receptor antagonist")

    def test_nhr_ppar(self):
        assert is_nuclear_receptor("PPAR receptor agonist")

    def test_not_nhr_dopamine(self):
        assert not is_nuclear_receptor("dopamine receptor antagonist")

    def test_not_nhr_proteasome(self):
        assert not is_nuclear_receptor("proteasome inhibitor")

    def test_not_nhr_empty(self):
        assert not is_nuclear_receptor("")

    def test_not_nhr_nan(self):
        assert not is_nuclear_receptor(float("nan"))

    def test_machinery_proteasome(self):
        assert is_constitutive_machinery("proteasome inhibitor")

    def test_machinery_tubulin(self):
        assert is_constitutive_machinery("tubulin inhibitor")

    def test_machinery_hdac(self):
        assert is_constitutive_machinery("HDAC inhibitor")

    def test_machinery_topoisomerase(self):
        assert is_constitutive_machinery("topoisomerase inhibitor")

    def test_machinery_excludes_nhr(self):
        assert not is_constitutive_machinery("glucocorticoid receptor agonist")

    def test_machinery_not_kinase(self):
        assert not is_constitutive_machinery("EGFR inhibitor")

    def test_narrow_proteasome(self):
        assert is_narrow_machinery("proteasome inhibitor")

    def test_narrow_tubulin(self):
        assert is_narrow_machinery("tubulin inhibitor")

    def test_narrow_not_hdac(self):
        assert not is_narrow_machinery("HDAC inhibitor")

    def test_narrow_not_topoisomerase(self):
        assert not is_narrow_machinery("topoisomerase inhibitor")

    def test_nhr_compound_moa_excluded_from_machinery(self):
        assert not is_constitutive_machinery("PPAR receptor agonist|HDAC inhibitor")


class TestComputeDI:

    def test_identical_signatures_zero_di(self):
        v = np.random.default_rng(0).standard_normal(100)
        sigs = np.tile(v, (10, 1))
        assert compute_di(sigs) == pytest.approx(0.0, abs=1e-10)

    def test_orthogonal_signatures_high_di(self):
        sigs = np.eye(10, 100)
        di = compute_di(sigs)
        assert di == pytest.approx(1.0, abs=1e-10)

    def test_di_increases_with_spread(self):
        rng = np.random.default_rng(42)
        v = rng.standard_normal(200)
        v /= np.linalg.norm(v)
        dis = []
        for sigma in [0.01, 0.1, 0.5, 1.0]:
            sigs = []
            for _ in range(20):
                n = rng.standard_normal(200)
                n -= np.dot(n, v) * v
                n /= np.linalg.norm(n)
                theta = abs(rng.normal(0, sigma))
                sigs.append(np.cos(theta) * v + np.sin(theta) * n)
            dis.append(compute_di(np.array(sigs)))
        for i in range(len(dis) - 1):
            assert dis[i] < dis[i + 1], f"DI should increase: {dis}"


class TestSpikeIn:

    def test_zero_sigma_zero_tau_gives_zero_di(self):
        rng = np.random.default_rng(99)
        di = generate_spike_in(0.0, 20, 0.0, rng)
        assert di == pytest.approx(0.0, abs=1e-10)

    def test_zero_sigma_nonzero_tau_gives_nonzero_di(self):
        rng = np.random.default_rng(99)
        dis = [generate_spike_in(0.0, 20, 0.1, np.random.default_rng(i)) for i in range(100)]
        assert np.mean(dis) > 0.001

    def test_large_sigma_gives_high_di(self):
        rng = np.random.default_rng(99)
        di = generate_spike_in(1.5, 20, 0.01, rng)
        assert di > 0.3

    def test_monotonicity_across_sigma(self):
        sigmas = [0.0, 0.1, 0.3, 0.6, 1.0]
        mean_dis = []
        for sigma in sigmas:
            dis = [generate_spike_in(sigma, 15, 0.05, np.random.default_rng(i + int(sigma * 1000)))
                   for i in range(200)]
            mean_dis.append(np.mean(dis))
        for i in range(len(mean_dis) - 1):
            assert mean_dis[i] < mean_dis[i + 1], f"Should be monotonic: {list(zip(sigmas, mean_dis))}"


class TestBootstrapMeanDiff:

    def test_separated_groups_significant(self):
        rng = np.random.default_rng(42)
        a = rng.normal(10, 1, 100)
        b = rng.normal(5, 1, 100)
        result = bootstrap_mean_diff(a, b, n_bootstrap=5000, seed=42)
        assert result["ci_lo"] > 0
        assert result["point_diff"] > 0

    def test_identical_groups_not_significant(self):
        rng = np.random.default_rng(42)
        a = rng.normal(5, 1, 100)
        b = rng.normal(5, 1, 100)
        result = bootstrap_mean_diff(a, b, n_bootstrap=5000, seed=42)
        assert result["ci_lo"] < 0 < result["ci_hi"]


class TestOnRealData:

    @pytest.fixture
    def drug_df(self):
        csv_path = (Path(__file__).resolve().parent.parent.parent /
                    "drug-perturbation-geometry" / "zenodo_v1" / "drug_instability_8949.csv")
        if not csv_path.exists():
            pytest.skip("Drug CSV not available")
        return pd.read_csv(csv_path)

    def test_nhr_count(self, drug_df):
        n = drug_df["moa"].apply(is_nuclear_receptor).sum()
        assert n == 125

    def test_constitutive_machinery_count(self, drug_df):
        n = drug_df["moa"].apply(is_constitutive_machinery).sum()
        assert n == 59

    def test_narrow_machinery_count(self, drug_df):
        n = drug_df["moa"].apply(is_narrow_machinery).sum()
        assert n == 22

    def test_no_overlap(self, drug_df):
        nhr = drug_df["moa"].apply(is_nuclear_receptor)
        mach = drug_df["moa"].apply(is_constitutive_machinery)
        assert (nhr & mach).sum() == 0

    def test_narrow_subset_of_broad(self, drug_df):
        broad = set(drug_df[drug_df["moa"].apply(is_constitutive_machinery)].index)
        narrow = set(drug_df[drug_df["moa"].apply(is_narrow_machinery)].index)
        assert narrow.issubset(broad)
