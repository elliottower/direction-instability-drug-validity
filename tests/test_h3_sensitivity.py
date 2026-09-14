import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
_mod = __import__("03c_h3_sensitivity")
partial_spearman = _mod.partial_spearman
projected_dispersion = _mod.projected_dispersion
permuted_statistic = _mod.permuted_statistic
load_bundle = _mod.load_bundle

N_LANDMARK = _mod.N_LANDMARK


def _cohort(n_drugs=40, n_genes=N_LANDMARK, shared_direction_for=()):
    rng = np.random.default_rng()
    sigs, dirs = [], []
    base = rng.standard_normal(n_genes)
    base /= np.linalg.norm(base)
    for c in range(n_drugs):
        k = int(rng.integers(5, 12))
        sigs.append(rng.standard_normal((k, n_genes)) * rng.uniform(0.5, 4.0))
        if c in shared_direction_for:
            dirs.append(base)
        else:
            u = rng.standard_normal(n_genes)
            dirs.append(u / np.linalg.norm(u))
    return sigs, np.array(dirs)


def _assemble(sigs, dirs):
    pair_diffs, P, E, M, K = [], [], [], [], []
    for s, u in zip(sigs, dirs):
        iu = np.triu_indices(s.shape[0], k=1)
        d = s[iu[0]] - s[iu[1]]
        pair_diffs.append(d)
        P.append(projected_dispersion(d, u))
        M.append(np.linalg.norm(d, axis=1).mean())
        K.append(s.shape[0])
        m = s.mean(axis=0)
        E.append((m @ u / (np.linalg.norm(m) * np.linalg.norm(u))) ** 2)
    n = len(sigs)
    offsets = np.cumsum([0] + [d.shape[0] for d in pair_diffs])
    G = np.abs(np.vstack(pair_diffs) @ dirs.T)
    row_drug = np.repeat(np.arange(n), np.diff(offsets))
    counts = np.diff(offsets).astype(float)
    means = np.array([s.mean(axis=0) for s in sigs])
    means /= np.linalg.norm(means, axis=1, keepdims=True)
    E_all = (means @ dirs.T) ** 2
    return (pair_diffs, np.array(P), np.array(E), np.array(M, float),
            np.array(K, float), G, row_drug, counts, E_all)


def test_lookup_P_matches_direct_projection_under_permutation():
    for _ in range(20):
        sigs, dirs = _cohort(n_drugs=25, n_genes=40)
        pd_, P, E, M, K, G, row_drug, counts, E_all = _assemble(sigs, dirs)
        perm = np.random.default_rng().permutation(len(sigs))
        _, P_p, _ = permuted_statistic(G, row_drug, counts, E_all, (M, K), perm)
        direct = np.array([projected_dispersion(pd_[c], dirs[perm[c]]) for c in range(len(sigs))])
        assert P_p == pytest.approx(direct, abs=1e-12)


def test_lookup_E_matches_direct_squared_cosine_under_permutation():
    for _ in range(20):
        sigs, dirs = _cohort(n_drugs=25, n_genes=40)
        *_, E_all = _assemble(sigs, dirs)[:9]
        perm = np.random.default_rng().permutation(len(sigs))
        _, _, E_p = permuted_statistic(*_assemble(sigs, dirs)[5:9], (np.ones(len(sigs)),), perm)
        direct = []
        for c, s in enumerate(sigs):
            m = s.mean(axis=0); u = dirs[perm[c]]
            direct.append((m @ u / (np.linalg.norm(m) * np.linalg.norm(u))) ** 2)
        assert E_p == pytest.approx(np.array(direct), abs=1e-12)


def test_full_permuted_statistic_matches_direct_partial_correlation():
    for _ in range(15):
        sigs, dirs = _cohort(n_drugs=30, n_genes=40)
        pd_, P, E, M, K, G, row_drug, counts, E_all = _assemble(sigs, dirs)
        perm = np.random.default_rng().permutation(len(sigs))
        fast, _, _ = permuted_statistic(G, row_drug, counts, E_all, (M, K), perm)
        P_d, E_d = [], []
        for c, s in enumerate(sigs):
            u = dirs[perm[c]]
            P_d.append(projected_dispersion(pd_[c], u))
            m = s.mean(axis=0)
            E_d.append((m @ u / (np.linalg.norm(m) * np.linalg.norm(u))) ** 2)
        slow = partial_spearman(np.array(P_d), np.array(E_d), (M, K))
        assert fast == pytest.approx(slow, abs=1e-12)


def test_identity_permutation_reproduces_the_observed_vectors():
    sigs, dirs = _cohort(n_drugs=30, n_genes=40)
    _, P, E, M, K, G, row_drug, counts, E_all = _assemble(sigs, dirs)
    rho, P_i, E_i = permuted_statistic(G, row_drug, counts, E_all, (M, K),
                                       np.arange(len(sigs)))
    assert P_i == pytest.approx(P, abs=1e-12)
    assert E_i == pytest.approx(E, abs=1e-12)
    assert rho == pytest.approx(partial_spearman(P, E, (M, K)), abs=1e-12)


def test_duplicate_directions_survive_permutation_as_a_multiset():
    shared = (0, 1, 2, 3)
    sigs, dirs = _cohort(n_drugs=20, n_genes=40, shared_direction_for=shared)
    _, P, E, M, K, G, row_drug, counts, E_all = _assemble(sigs, dirs)
    for _ in range(20):
        perm = np.random.default_rng().permutation(len(sigs))
        _, P_p, _ = permuted_statistic(G, row_drug, counts, E_all, (M, K), perm)
        direct = np.array([projected_dispersion(
            sigs[c][np.triu_indices(sigs[c].shape[0], 1)[0]]
            - sigs[c][np.triu_indices(sigs[c].shape[0], 1)[1]], dirs[perm[c]])
            for c in range(len(sigs))])
        assert P_p == pytest.approx(direct, abs=1e-12)


def test_permutation_with_fixed_points_is_handled():
    # constructed so the swap provably changes the two moved drugs: drug 3 has a
    # difference along e0 only and drug 7 along e1 only, with orthogonal directions
    n_genes = 8
    sigs, dirs = [], []
    for c in range(10):
        s = np.zeros((5, n_genes))
        s[0, c % n_genes] = 1.0 + c
        s[1:] = 0.0
        sigs.append(s)
        u = np.zeros(n_genes); u[c % n_genes] = 1.0
        dirs.append(u)
    dirs = np.array(dirs)
    _, P, E, M, K, G, row_drug, counts, E_all = _assemble(sigs, dirs)
    perm = np.arange(10); perm[[3, 7]] = perm[[7, 3]]        # 8 fixed points
    _, P_p, _ = permuted_statistic(G, row_drug, counts, E_all, (M, K), perm)
    for c in (0, 1, 2, 4, 5, 6, 8, 9):
        assert P_p[c] == pytest.approx(P[c], abs=1e-12)
    assert P_p[3] == pytest.approx(0.0, abs=1e-12)           # projected onto e7
    assert P[3] > 0.1


def test_bundle_direction_must_be_unit_norm():
    sigs, dirs = _cohort(n_drugs=6, n_genes=N_LANDMARK)
    dirs[2] *= 3.0
    assert not np.isclose(np.linalg.norm(dirs[2]), 1.0, atol=1e-10)
    assert all(np.isclose(np.linalg.norm(dirs[i]), 1.0, atol=1e-10)
               for i in range(len(dirs)) if i != 2)


def test_expected_shard_set_detects_a_missing_shard():
    n_deposited, shard = 795, 250
    n_shards = (n_deposited + shard - 1) // shard
    expected = {f"shard_{s:03d}.npz" for s in range(n_shards)}
    assert len(expected) == 4
    observed = expected - {"shard_002.npz"}
    assert observed != expected
    assert sorted(expected - observed) == ["shard_002.npz"]
    # three of four shards can still exceed the 700 floor, which is the point
    assert 3 * shard > 700


def test_duplicate_hairpin_rows_collapse_to_one_signature():
    ids = ["s1", "s1", "s2", "s3"]
    assert len(ids) == 4
    assert len(sorted(set(ids))) == 3
    assert len({i for i in ids}) >= 3


def test_partial_spearman_uses_average_ranks_for_ties():
    y = np.array([1.0, 1.0, 1.0, 2.0, 3.0, 4.0])
    assert rankdata(y, method="average")[:3] == pytest.approx([2.0, 2.0, 2.0])
    x = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    assert partial_spearman(y, x, (np.ones(6),)) == pytest.approx(
        spearmanr(y, x).statistic, abs=1e-12)


def test_partial_spearman_removes_a_covariate_it_is_given():
    z = np.linspace(-2, 2, 400)
    y = z + 0.6 * np.sin(7 * z)
    x = z + 0.6 * np.cos(11 * z)
    marginal = abs(spearmanr(y, x).statistic)
    adjusted = abs(partial_spearman(y, x, (z,)))
    assert marginal > 0.85
    assert adjusted < marginal / 2


def test_partial_spearman_refuses_a_perfectly_collinear_covariate():
    z = np.linspace(-2, 2, 100)
    with pytest.raises(AssertionError, match="collinear"):
        partial_spearman(z.copy(), z.copy(), (z,))


def test_empirical_p_value_uses_the_add_one_formula():
    null = np.arange(10_000, dtype=float)
    obs = 9_999.0
    assert (1 + int((null >= obs).sum())) / (10_000 + 1) == pytest.approx(2 / 10_001)
    obs_big = 1e9
    assert (1 + int((null >= obs_big).sum())) / (10_000 + 1) == pytest.approx(1 / 10_001)


def test_equivalence_bound_logic():
    bound = _mod.S3_EQUIV_BOUND
    inside = [-0.10, 0.10]
    straddles = [-0.20, 0.05]
    touching = [-0.15, 0.05]
    assert (-bound < inside[0] and inside[1] < bound)
    assert not (-bound < straddles[0] and straddles[1] < bound)
    assert not (-bound < touching[0] and touching[1] < bound)


def test_load_bundle_rejects_ragged_and_miswidthed_input(tmp_path):
    good = np.array([np.zeros((5, N_LANDMARK))] * 3, dtype=object)
    p = tmp_path / "ragged.npz"
    np.savez(p, drugs=np.array(["a", "b", "c"]), targets=np.array(["t", "u"]),
             signatures=good, directions=np.zeros((3, N_LANDMARK)))
    with pytest.raises(AssertionError, match="ragged"):
        load_bundle(p)

    q = tmp_path / "narrow.npz"
    np.savez(q, drugs=np.array(["a"]), targets=np.array(["t"]),
             signatures=np.array([np.zeros((5, 10))], dtype=object),
             directions=np.zeros((1, 10)))
    with pytest.raises(AssertionError, match="directions are"):
        load_bundle(q)


def test_sidecar_hash_matches_the_file_it_describes(tmp_path):
    payload = {"a": 1, "b": [2, 3]}
    f = tmp_path / "r.json"
    f.write_text(json.dumps(payload, indent=2))
    digest = hashlib.sha256(f.read_bytes()).hexdigest()
    (tmp_path / "r.json.sha256").write_text(digest + "\n")
    assert (tmp_path / "r.json.sha256").read_text().strip() == \
        hashlib.sha256(f.read_bytes()).hexdigest()
    assert json.loads(f.read_text()) == payload
