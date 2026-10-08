"""Power of S2-T: the unique-target permutation must fire on a planted signal.

The existing tests establish that the permutation is built correctly, that it keeps
drugs of one target together and that it is not a record permutation. None of them
establishes that it can refuse anything. A failing S2-T is evidence only if the same
machinery returns a small p-value when the association really is target-specific, so
both directions are exercised here against the code the run uses.

The two cohorts differ in one respect. In the first, a drug's dispersion along its
own target axis tracks its alignment to that axis, and reading any other axis finds
isotropic noise instead; reassigning axes destroys the association. In the second,
every target shares one axis, so reassigning axes changes nothing and the observed
statistic must sit inside its own null.
"""
import site
from pathlib import Path

import numpy as np

site.addsitedir(str(Path(__file__).resolve().parents[1] / "experiments"))
_mod = __import__("03c_h3_sensitivity")
partial_spearman = _mod.partial_spearman
permuted_statistic = _mod.permuted_statistic

from geometry.inference import unique_target_permutations

N_GENES = 200
N_PERM = 2000
NOISE = 1.0


def _cohort(target_specific):
    """Signatures, directions and target labels for one of the two cohorts.

    Each drug carries a quality q. Its mean sits q along the axis it is scored
    against and its contexts jitter by q along the same axis, so projected dispersion
    and alignment both rise with q and the two move together. The isotropic noise is
    large enough that the mean pairwise difference norm is effectively constant,
    which keeps the magnitude covariate from absorbing the planted signal.
    """
    rng = np.random.default_rng()
    n_targets = 30
    axes = np.linalg.qr(rng.standard_normal((N_GENES, n_targets)))[0].T
    shared = axes[0]

    sigs, dirs, targets = [], [], []
    for t in range(n_targets):
        axis = axes[t] if target_specific else shared
        for _ in range(int(rng.integers(2, 6))):
            q = float(rng.uniform(0.3, 3.0))
            k = int(rng.integers(6, 12))
            jitter = rng.standard_normal(k)[:, None] * q * axis
            noise = rng.standard_normal((k, N_GENES)) * NOISE
            sigs.append(q * axis + jitter + noise)
            dirs.append(axis)
            targets.append(f"T{t:02d}")
    return sigs, np.asarray(dirs, dtype=np.float64), np.asarray(targets)


def _permutation_test(sigs, dirs, targets):
    """The observed partial correlation and its unique-target permutation p-value.

    Built the way the script builds it, so the test exercises the same lookup tables
    and the same permutation generator rather than a restatement of them.
    """
    n = len(sigs)
    pair_diffs, P, E, M_delta, K = [], np.empty(n), np.empty(n), np.empty(n), np.empty(n)
    for i, (S, u) in enumerate(zip(sigs, dirs)):
        iu = np.triu_indices(S.shape[0], k=1)
        diffs = S[iu[0]] - S[iu[1]]
        pair_diffs.append(diffs)
        P[i] = np.abs(diffs @ u).mean()
        M_delta[i] = np.linalg.norm(diffs, axis=1).mean()
        K[i] = S.shape[0]
        mean_sig = S.mean(axis=0)
        E[i] = (mean_sig @ u / (np.linalg.norm(mean_sig) * np.linalg.norm(u))) ** 2

    offsets = np.cumsum([0] + [d.shape[0] for d in pair_diffs])
    G = np.abs(np.vstack(pair_diffs) @ dirs.T)
    row_drug = np.repeat(np.arange(n), np.diff(offsets))
    counts = np.diff(offsets).astype(np.float64)
    unit_means = np.array([S.mean(axis=0) for S in sigs])
    unit_means /= np.linalg.norm(unit_means, axis=1, keepdims=True)
    E_all = (unit_means @ dirs.T) ** 2
    covars = (M_delta, K)

    observed = partial_spearman(P, E, covars)
    _, first_row = np.unique(targets, return_index=True)
    null = np.array([
        permuted_statistic(G, row_drug, counts, E_all, covars, first_row[assigned])[0]
        for assigned in unique_target_permutations(targets, N_PERM, int(np.random.default_rng().integers(1 << 30)))])
    p_perm = (1 + int((null >= observed).sum())) / (N_PERM + 1)
    return observed, null, p_perm


def test_s2t_fires_when_the_association_is_target_specific():
    observed, null, p_perm = _permutation_test(*_cohort(target_specific=True))
    assert p_perm < 0.01, (
        f"a planted target-specific association went undetected: observed {observed:.4f}, "
        f"null median {np.median(null):.4f}, p_perm {p_perm:.4f}")
    assert np.median(null) < observed / 2, (
        f"reassigning the axes left most of the association standing, so the cohort "
        f"is not target-specific enough to establish power: observed {observed:.4f}, "
        f"null median {np.median(null):.4f}")


def test_s2t_stays_quiet_when_every_target_shares_one_axis():
    observed, null, p_perm = _permutation_test(*_cohort(target_specific=False))
    assert p_perm > 0.05, (
        f"the permutation refused a cohort carrying no target-specific signal: "
        f"observed {observed:.4f}, null median {np.median(null):.4f}, p_perm {p_perm:.4f}")


def test_the_two_cohorts_separate():
    """Both arms in one run, so a change that breaks the contrast fails here."""
    specific, _, p_specific = _permutation_test(*_cohort(target_specific=True))
    _, _, p_shared = _permutation_test(*_cohort(target_specific=False))
    assert p_specific < 0.01 < p_shared, (
        f"the permutation does not separate the two cohorts: target-specific "
        f"p_perm {p_specific:.4f}, shared-axis p_perm {p_shared:.4f}")
