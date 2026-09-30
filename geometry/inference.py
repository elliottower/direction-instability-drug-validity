"""Inference machinery for the H3 registrations frozen at 7f57136.

Every rule here is fixed by `experiments/PREREG_H3_S1S3_CORRECTED_BASELINE.md` and
`experiments/PREREG_H3_REFERENCE_DISCORDANCE.md`. Nothing in this module reads a
result to decide what to do with it.

The unit of resampling and of permutation is the target, not the drug: drugs that
share an annotated target inherit the same reference direction, so treating them
as independent overstates the information they carry.
"""
import numpy as np
from scipy.stats import rankdata, spearmanr

RATIO = 0.5          # what counts as a substantive drop in the comparison rule
MIN_REFERENCE = 0.15  # below this the comparison rule does not apply
PRACTICAL_NULL = 0.15  # the equivalence bound for "practically null"


def rank_partial_correlation(y: np.ndarray, x: np.ndarray, covars) -> float:
    """Partial Spearman by rank residualization, the procedure frozen at f288507.

    Midrank each variable, regress ranked y and ranked x separately on the ranked
    covariates with an explicit intercept, and correlate the residuals.
    """
    ry = rankdata(y, method="average")
    rx = rankdata(x, method="average")
    rc = np.column_stack([rankdata(c, method="average") for c in covars])
    design = np.column_stack([np.ones(len(ry)), rc])

    def resid(v):
        beta, *_ = np.linalg.lstsq(design, v, rcond=None)
        return v - design @ beta

    a, b = resid(ry), resid(rx)
    # Under near-perfect collinearity with the covariates both residual vectors are
    # rounding error, and their correlation is a confident number computed from
    # noise. Refuse rather than report it.
    for name, resid_vec, orig in (("y", a, ry), ("x", b, rx)):
        scale = np.linalg.norm(orig - orig.mean())
        assert np.linalg.norm(resid_vec) > 1e-8 * max(scale, 1.0), (
            f"{name} is collinear with the covariates; the partial correlation "
            "would be computed from residual noise")
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


def target_members(targets: np.ndarray):
    """Row indices of each unique target, in first-appearance order."""
    uniq, inverse = np.unique(targets, return_inverse=True)
    return uniq, inverse, [np.flatnonzero(inverse == t) for t in range(len(uniq))]


def target_balanced_spearman(x: np.ndarray, y: np.ndarray, targets: np.ndarray) -> float:
    """Spearman across targets of the within-target medians of each variable."""
    _, _, members = target_members(targets)
    mx = np.array([np.median(x[m]) for m in members])
    my = np.array([np.median(y[m]) for m in members])
    return float(spearmanr(mx, my).statistic)


def target_balanced_mean(values: np.ndarray, targets: np.ndarray) -> float:
    """Mean within each target, then mean across targets with equal weight."""
    _, _, members = target_members(targets)
    return float(np.mean([values[m].mean() for m in members]))


def own_target_percentile(alignment: np.ndarray, own: np.ndarray) -> np.ndarray:
    """Rank of a drug's own target among all candidate targets, in [0, 1].

    Args:
        alignment: (n_drugs, n_targets); entry (c, t) is drug c's alignment with
            target t's direction.
        own: (n_drugs,) column index of each drug's annotated target.

    Returns:
        (n_drugs,) fraction of the other targets the drug aligns with less well,
        counting ties as half. 0.5 under no target identity.
    """
    n_drugs, n_targets = alignment.shape
    assert n_targets > 1, "an own-target percentile needs at least two candidate targets"
    own_value = alignment[np.arange(n_drugs), own][:, None]
    below = (alignment < own_value).sum(axis=1)
    ties = (alignment == own_value).sum(axis=1) - 1      # the own column ties itself
    return (below + 0.5 * ties) / (n_targets - 1)


def percentile_table(alignment: np.ndarray) -> np.ndarray:
    """(n_drugs, n_targets) own-target percentile for every candidate assignment.

    Entry (c, j) is what `own_target_percentile` returns for drug c if target j
    were its target. A permutation then becomes a lookup instead of a recount,
    which is what makes 10,000 permutations affordable.
    """
    n_drugs, n_targets = alignment.shape
    assert n_targets > 1, "an own-target percentile needs at least two candidate targets"
    order = np.argsort(alignment, axis=1, kind="stable")
    ranks = np.empty_like(order)
    rows = np.arange(n_drugs)[:, None]
    ranks[rows, order] = np.arange(n_targets)[None, :]
    # ties: count strictly below and half the ties, as own_target_percentile does
    table = np.empty_like(alignment, dtype=float)
    for c in range(n_drugs):
        values = alignment[c]
        sorted_values = values[order[c]]
        below = np.searchsorted(sorted_values, values, side="left")
        at_or_below = np.searchsorted(sorted_values, values, side="right")
        ties = at_or_below - below - 1
        table[c] = (below + 0.5 * ties) / (n_targets - 1)
    return table


def cluster_bootstrap(targets: np.ndarray, statistic, n_replicates: int, seed: int) -> np.ndarray:
    """Resample targets with replacement, carrying every drug of a drawn target.

    `statistic` takes an index array of rows and returns a float or a 1-D array;
    paired quantities are returned together so that they share a replicate.
    """
    _, _, members = target_members(targets)
    n_targets = len(members)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_replicates):
        drawn = rng.integers(0, n_targets, n_targets)
        idx = np.concatenate([members[t] for t in drawn])
        draws.append(np.atleast_1d(np.asarray(statistic(idx), dtype=float)))
    out = np.vstack(draws)
    assert np.isfinite(out).all(), "nonfinite bootstrap replicate"
    return out


def unique_target_permutations(targets: np.ndarray, n_permutations: int, seed: int):
    """Permute direction assignment among unique targets, not among drug records.

    Yields (n_drugs,) arrays of target column indices. Every drug of one target
    receives the same reassigned direction, so target group sizes and the
    within-target dependence of the observed mapping are preserved. Fixed points
    are allowed.
    """
    _, inverse, members = target_members(targets)
    n_targets = len(members)
    rng = np.random.default_rng(seed)
    for _ in range(n_permutations):
        yield rng.permutation(n_targets)[inverse]


def percentile_interval(draws: np.ndarray, level: float = 95.0):
    """Two-sided percentile interval of bootstrap draws."""
    tail = (100.0 - level) / 2.0
    return float(np.percentile(draws, tail)), float(np.percentile(draws, 100.0 - tail))


def excludes_zero(interval) -> bool:
    lo, hi = interval
    return lo > 0.0 or hi < 0.0


def is_practically_null(draws: np.ndarray, bound: float = PRACTICAL_NULL) -> bool:
    """True when the 90% interval lies inside (-bound, bound)."""
    lo, hi = percentile_interval(draws, level=90.0)
    return -bound < lo and hi < bound


def comparison_reading(rho_reference: float, rho_alternative: float,
                       difference_draws: np.ndarray, alternative_draws: np.ndarray) -> str:
    """The frozen comparison rule for an alternative estimate of one association.

    Every reading requires a paired interval, so a large apparent drop whose
    interval spans zero is inconclusive rather than attenuated.
    """
    if abs(rho_reference) < MIN_REFERENCE:
        return "reference too small to attenuate"
    difference_excludes_zero = excludes_zero(percentile_interval(difference_draws, 95.0))
    alternative_excludes_zero = excludes_zero(percentile_interval(alternative_draws, 95.0))
    same_sign = rho_reference * rho_alternative > 0

    if not same_sign and difference_excludes_zero:
        return "reversed"
    if same_sign and abs(rho_alternative) <= RATIO * abs(rho_reference) and difference_excludes_zero:
        return "attenuated"
    if same_sign and abs(rho_alternative) >= RATIO * abs(rho_reference) and alternative_excludes_zero:
        return "retained"
    return "inconclusive"


def permutation_reading(observed: float, null_draws: np.ndarray) -> str:
    """Assignment-dependent, permutation-resistant, or indeterminate.

    The observed statistic is compared with a null that reassigns whole target
    directions. Resistance is never described as a generic or shared component
    here; that reading needs the other modules.
    """
    median = float(np.median(null_draws))
    if observed >= 0:
        beyond = observed > np.percentile(null_draws, 99.0)
    else:
        beyond = observed < np.percentile(null_draws, 1.0)
    inside_central = (np.percentile(null_draws, 2.5) <= observed <= np.percentile(null_draws, 97.5))
    retains_magnitude = median * observed > 0 and abs(median) >= RATIO * abs(observed)

    if beyond and abs(median) <= RATIO * abs(observed):
        return "assignment-dependent"
    if inside_central and retains_magnitude:
        return "permutation-resistant"
    return "indeterminate"


def two_sided_permutation_p(observed: float, null_draws: np.ndarray) -> float:
    """(1 + #{|null| >= |observed|}) / (n + 1)."""
    n = len(null_draws)
    return float((1 + int((np.abs(null_draws) >= abs(observed)).sum())) / (n + 1))
