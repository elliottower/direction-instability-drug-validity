"""How distinguishable are the shRNA reference directions from each other?

Descriptive, and deliberately so: no hypothesis, no criterion, no null. S2-T
reassigns target directions among targets and asks whether the P-E association
survives. The test can only be sensitive to target specificity if reassignment
actually changes the axis a drug is scored against, so the mutual geometry of the
258 unique directions sets a ceiling on what the test could have detected,
whatever the biology is. That geometry is a property of the reference, measurable
without reference to any result, and it is measured here.

The second half reads the detection threshold off the nulls the registered run
recorded: the smallest attenuation that would have returned p < 0.01.

    PYTHONPATH=. uv run --no-project --with numpy --with scipy \
        python experiments/03j_reference_axis_geometry.py
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
SENSITIVITY = REPO / "results" / "03c_h3_sensitivity"
OUT = REPO / "results" / "03j_reference_axis_geometry"

# The bundle the registered S1-S3 run read, named by the hash its result file
# records under provenance.bundle_sha256. The local results/03c_h3_sensitivity
# copy is a different, superseded build (Deviation 17), so the path is not enough
# to identify the right file and the hash is checked rather than assumed.
EXPECTED_BUNDLE_SHA256 = "bfa02ff0b991c99efc40310dca28e03978a93f1913664982f04fc5c3d657e8f6"
DEFAULT_BUNDLE = SENSITIVITY / "gctx_rebuild" / "cohort_bundle.npz"
N_LANDMARK = 978

# an artifact that carries its own expected hash attests to itself, so the
# expected values for the recorded run live here and not in the files they pin
EXPECTED_DRAWS_SHA256 = "e3b26398e2f423e96539f7e365ebe44f36b6e28a4035a2a232b0c856e638febf"
EXPECTED_RESULTS_SHA256 = "6f60f1ad82de2c24a020eac8d30f807420590591ada2c959e1fa5e584ec4538f"

# the criterion the registration set for S2 and S2-T
P_CRITERION = 0.01


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def unique_directions(targets, dirs):
    """One direction per target, and a refusal if a target carries two of them.

    Built the way the analysis builds it: the first row of each target in sorted
    unique order. Every drug of a target is supposed to share one direction, so
    that is checked rather than relied on.
    """
    unique_targets, representative = np.unique(targets, return_index=True)
    for name, row in zip(unique_targets, representative):
        members = dirs[targets == name]
        if not (np.abs(members - dirs[row]).max() == 0.0):
            raise AssertionError(f"target {name} carries more than one direction")
    return unique_targets, dirs[representative]


def cosine_structure(U):
    """Pairwise |cos| among unit rows, and how many directions the set spans."""
    norms = np.linalg.norm(U, axis=1)
    if not (np.abs(norms - 1.0).max() < 1e-9):
        raise AssertionError(f"directions are not unit vectors; worst norm {norms.max()}")
    gram = U @ U.T
    off = np.abs(gram[np.triu_indices(len(U), k=1)])

    # participation ratio of the Gram spectrum: the number of directions the set
    # behaves as, which equals len(U) for an orthonormal set and 1 for a set that
    # is one direction repeated
    eigenvalues = np.linalg.eigvalsh(gram)
    eigenvalues = np.clip(eigenvalues, 0.0, None)
    total = float(eigenvalues.sum())
    participation = total ** 2 / float((eigenvalues ** 2).sum())
    descending = np.sort(eigenvalues)[::-1]
    cumulative = np.cumsum(descending) / total

    return {
        "n_directions": int(len(U)),
        "n_pairs": int(off.size),
        "abs_cosine": {
            "mean": float(off.mean()),
            "median": float(np.median(off)),
            "p25": float(np.quantile(off, 0.25)),
            "p75": float(np.quantile(off, 0.75)),
            "p90": float(np.quantile(off, 0.90)),
            "p99": float(np.quantile(off, 0.99)),
            "max": float(off.max()),
            "fraction_above_0.3": float((off > 0.3).mean()),
            "fraction_above_0.5": float((off > 0.5).mean()),
            "fraction_above_0.7": float((off > 0.7).mean()),
        },
        "expected_abs_cosine_for_isotropic_directions": float(
            np.sqrt(2.0 / (np.pi * N_LANDMARK))),
        "effective_rank_participation_ratio": float(participation),
        "n_components_for_90pct_of_trace": int(np.searchsorted(cumulative, 0.90) + 1),
        "n_components_for_50pct_of_trace": int(np.searchsorted(cumulative, 0.50) + 1),
        "leading_eigenvalue_share": float(descending[0] / total),
    }


def detection_threshold(draws_path, observed):
    """The smallest observed statistic each recorded null would have refused.

    S2 and S2-T refuse when the one-sided permutation p falls below 0.01, so the
    threshold is the null's 99th percentile and the detectable effect is the
    distance from the null's centre to that percentile.
    """
    draws = np.load(draws_path)
    out = {}
    for key, label in (("permutation_null_target", "S2T_unique_target"),
                       ("permutation_null", "S2_drug_record")):
        null = draws[key]
        threshold = float(np.quantile(null, 1.0 - P_CRITERION))
        median = float(np.median(null))
        out[label] = {
            "n_permutations": int(null.size),
            "null_median": median,
            "null_sd": float(null.std(ddof=1)),
            "threshold_at_p_lt_0.01": threshold,
            "minimum_detectable_contribution": threshold - median,
            "observed": observed,
            "observed_minus_threshold": observed - threshold,
            "observed_percentile_in_null": float((null < observed).mean()),
        }
    return out


def main(bundle_path):
    bundle_path = Path(bundle_path)
    digest = sha256_file(bundle_path)
    if not (digest == EXPECTED_BUNDLE_SHA256):
        raise AssertionError(
            f"{bundle_path} hashes to {digest}; the registered run read "
            f"{EXPECTED_BUNDLE_SHA256}")

    results_path = SENSITIVITY / "h3_sensitivity_results.json"
    if not (sha256_file(results_path) == EXPECTED_RESULTS_SHA256):
        raise AssertionError("the recorded S1-S3 result file is not the one this pins")
    draws_path = SENSITIVITY / "h3_sensitivity_draws.npz"
    if not (sha256_file(draws_path) == EXPECTED_DRAWS_SHA256):
        raise AssertionError("the recorded draws file is not the one this pins")

    recorded = json.loads(results_path.read_text())
    observed = float(recorded["S2T_unique_target_permutation"]["rho_observed"])
    if not (recorded["provenance"]["bundle_sha256"] == EXPECTED_BUNDLE_SHA256):
        raise AssertionError("the recorded run did not read the bundle this pins")

    z = np.load(bundle_path, allow_pickle=True)
    targets = np.asarray([str(t) for t in z["targets"]])
    dirs = np.asarray(z["directions"], dtype=np.float64)
    if not (dirs.shape == (len(targets), N_LANDMARK)):
        raise AssertionError(f"directions are {dirs.shape}")

    unique_targets, U = unique_directions(targets, dirs)
    if not (len(unique_targets) == recorded["S2T_unique_target_permutation"]["n_unique_targets"]):
        raise AssertionError(
            f"{len(unique_targets)} unique targets here, "
            f"{recorded['S2T_unique_target_permutation']['n_unique_targets']} in the run")

    summary = {
        "question": ("whether reassigning a target direction to another target changes "
                     "the axis a drug is scored against, which bounds what S2-T could "
                     "have detected"),
        "status": ("descriptive diagnostic of the reference, computed after the "
                   "registered verdicts were recorded; no criterion and no null"),
        "reference": "shRNA consensus knockdown directions, one per unique target",
        "pins": {
            "bundle_sha256": digest,
            "results_sha256": EXPECTED_RESULTS_SHA256,
            "draws_sha256": EXPECTED_DRAWS_SHA256,
        },
        "cohort": {
            "n_drugs": int(len(targets)),
            "n_unique_targets": int(len(unique_targets)),
            "drugs_per_target": {
                "min": int(np.unique(targets, return_counts=True)[1].min()),
                "median": float(np.median(np.unique(targets, return_counts=True)[1])),
                "max": int(np.unique(targets, return_counts=True)[1].max()),
            },
        },
        "direction_geometry": cosine_structure(U),
        "detection_threshold": detection_threshold(draws_path, observed),
        "environment": {"numpy": np.__version__, "scipy": scipy.__version__,
                        "python": platform.python_version()},
    }

    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "reference_axis_geometry.json"
    if target.exists():
        raise AssertionError(f"{target} exists; a result file is never overwritten")
    target.write_text(json.dumps(summary, indent=2))
    (OUT / "reference_axis_geometry.json.sha256").write_text(sha256_file(target) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"\nwritten to {target}")


if __name__ == "__main__":
    a = argparse.ArgumentParser()
    a.add_argument("--bundle", default=str(DEFAULT_BUNDLE))
    main(a.parse_args().bundle)
