"""Does rank residualization bias a partial Spearman correlation, or did a test omit an intercept?

Two papers in this family disagree about one choice. The atlas analysis
residualizes with Pearson and then ranks, and its methods state that
rank-then-residualize "produce[s] spurious partial correlations of rho = 0.05 to
0.25 under the null at large n", citing its own regression test. The H3
registration frozen at f288507 registered rank-then-residualize with an explicit
intercept and ran it.

The disagreement is settled by reading the test rather than by arguing about the
methods. Its comparator regresses rank vectors on ranked covariates with no
intercept column, and ranks are strictly positive with mean (n+1)/2, so the fit
is forced through the origin and leaves a large common offset in both residual
vectors. Its own reference arm omits the intercept too, and is unharmed only
because the simulated raw variables are already centered. So the test compares a
no-intercept fit on uncentered data against a no-intercept fit on centered data,
which is not a comparison between the two procedures.

Four estimators run on one set of draws from the test's own generating model, so
the arms differ in the estimator and in nothing else.

    uv run --no-project --with numpy --with scipy \
        python experiments/03k_partial_spearman_estimators.py
"""
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy import stats

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
OUT = REPO / "results" / "03k_partial_spearman_estimators"

# the atlas test's own settings, reused so the arms are comparable to its claim
N = 5000
N_REPLICATES = 200
COVARIATE_SCALE = 0.3
SEED = 20261008
# the band the atlas test asserts against: correct under 0.02, biased over 0.05
ATLAS_UNBIASED_BAND = 0.02
ATLAS_BIASED_FLOOR = 0.05


def _residualize(design, v):
    beta, *_ = np.linalg.lstsq(design, v, rcond=None)
    return v - design @ beta


def rank_residualize_with_intercept(x, y, Z):
    """The estimator registered at f288507 and used by the H3 analyses."""
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    rZ = np.column_stack([stats.rankdata(Z[:, i]) for i in range(Z.shape[1])])
    design = np.column_stack([np.ones(len(rx)), rZ])
    return float(np.corrcoef(_residualize(design, rx), _residualize(design, ry))[0, 1])


def rank_residualize_no_intercept(x, y, Z):
    """The atlas test's comparator, reproduced exactly, as a software diagnostic."""
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    rZ = np.column_stack([stats.rankdata(Z[:, i]) for i in range(Z.shape[1])])
    return float(np.corrcoef(_residualize(rZ, rx), _residualize(rZ, ry))[0, 1])


def pearson_residualize_no_intercept(x, y, Z):
    """The atlas estimator as its power_analysis.partial_spearman implements it."""
    return float(stats.spearmanr(_residualize(Z, x), _residualize(Z, y))[0])


def pearson_residualize_with_intercept(x, y, Z):
    """The atlas procedure with the intercept its implementation leaves out."""
    design = np.column_stack([np.ones(len(x)), Z])
    return float(stats.spearmanr(_residualize(design, x), _residualize(design, y))[0])


ESTIMATORS = {
    "rank_residualize_with_intercept": rank_residualize_with_intercept,
    "rank_residualize_no_intercept": rank_residualize_no_intercept,
    "pearson_residualize_no_intercept": pearson_residualize_no_intercept,
    "pearson_residualize_with_intercept": pearson_residualize_with_intercept,
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    rng = np.random.default_rng(SEED)
    draws = {name: np.empty(N_REPLICATES) for name in ESTIMATORS}
    # the offset the missing intercept leaves behind, measured rather than argued
    residual_means = np.empty((N_REPLICATES, 2))

    for r in range(N_REPLICATES):
        Z = rng.standard_normal((N, 2))
        x = Z @ rng.standard_normal(2) * COVARIATE_SCALE + rng.standard_normal(N)
        y = Z @ rng.standard_normal(2) * COVARIATE_SCALE + rng.standard_normal(N)
        for name, estimator in ESTIMATORS.items():
            draws[name][r] = estimator(x, y, Z)

        rx, ry = stats.rankdata(x), stats.rankdata(y)
        rZ = np.column_stack([stats.rankdata(Z[:, i]) for i in range(2)])
        residual_means[r] = (_residualize(rZ, rx).mean(), _residualize(rZ, ry).mean())

    summary = {
        "question": ("whether the registered rank-then-residualize estimator is biased "
                     "under conditional independence, or whether the test that reported "
                     "the bias omitted an intercept"),
        "status": ("software diagnostic of a test's comparator; no data are analyzed and "
                   "no registered statistic is recomputed"),
        "generating_model": {
            "source": "direction-instability-atlas tests/test_partial_spearman.py, "
                      "test_rank_then_residualize_is_biased",
            "n": N, "n_replicates": N_REPLICATES, "n_covariates": 2,
            "covariate_scale": COVARIATE_SCALE,
            "description": ("Z standard normal; x and y each a linear function of Z plus "
                            "independent standard normal noise, so x and y are "
                            "conditionally independent given Z and the true partial "
                            "correlation is zero"),
            "seed": SEED,
        },
        "atlas_assertions_under_test": {
            "correct_method_abs_mean_below": ATLAS_UNBIASED_BAND,
            "buggy_method_mean_above": ATLAS_BIASED_FLOOR,
        },
        "estimators": {
            name: {
                "mean": float(d.mean()),
                "median": float(np.median(d)),
                "sd": float(d.std(ddof=1)),
                "min": float(d.min()),
                "max": float(d.max()),
                "passes_atlas_unbiased_band": bool(abs(d.mean()) < ATLAS_UNBIASED_BAND),
                "clears_atlas_biased_floor": bool(d.mean() > ATLAS_BIASED_FLOOR),
            }
            for name, d in draws.items()
        },
        "rank_residual_mean_without_intercept": {
            "description": ("the mean of each residual vector under the no-intercept rank "
                            "fit; a fit through the origin on strictly positive ranks of "
                            f"mean {(N + 1) / 2} cannot remove the level, and the common "
                            "offset that survives is what correlates"),
            "x_residual_mean": float(residual_means[:, 0].mean()),
            "y_residual_mean": float(residual_means[:, 1].mean()),
        },
        "environment": {"numpy": np.__version__, "scipy": scipy.__version__,
                        "python": platform.python_version()},
    }

    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "estimator_comparison.json"
    if target.exists():
        raise AssertionError(f"{target} exists; a result file is never overwritten")
    target.write_text(json.dumps(summary, indent=2))
    (OUT / "estimator_comparison.json.sha256").write_text(sha256_file(target) + "\n")

    for name, stat in summary["estimators"].items():
        print(f"{name:38s} mean {stat['mean']:+.4f}  sd {stat['sd']:.4f}")
    print(f"\nno-intercept rank residual means: "
          f"x {summary['rank_residual_mean_without_intercept']['x_residual_mean']:.1f}, "
          f"y {summary['rank_residual_mean_without_intercept']['y_residual_mean']:.1f}")
    print(f"written to {target}")


if __name__ == "__main__":
    main()
