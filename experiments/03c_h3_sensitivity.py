"""H3 sensitivity analyses S1-S3, registered in PREREG_H3_MAGNITUDE_AND_SHARED_AXIS.md
(frozen at f288507, before the LINCS matrices were rebuilt).

S1  the association survives adjustment for pairwise-difference magnitude and coverage
S2  it exceeds a null that permutes target assignment while keeping the shared axis
S3  the raw score stays practically equivalent to zero under the same adjustment

The amendment frozen at 7f57136 adds S1-T, S2-T and S3-T: the same three
statistics with the target rather than the drug as the unit, resampled as clusters
and permuted among unique targets. All three target-level versions are required
and none compensates for another; the drug-level versions are reported as
registered, and where the two disagree the target-level verdict governs.

Input is the bundle written by the rebuild stage: per-drug signature matrices, the
unit target direction, and the identifiers, for drugs that reproduced the deposited
artifact. Nothing here re-derives the cohort; it fails closed if the bundle and the
deposited artifact disagree.

    PYTHONPATH=. uv run --no-project --with numpy --with scipy \
        python experiments/03c_h3_sensitivity.py \
        --bundle results/03c_h3_sensitivity/cohort_bundle.npz \
        --manifest results/03c_h3_sensitivity/rebuild_manifest.json

PYTHONPATH is what puts `geometry` on the path under --no-project.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import rankdata, spearmanr

from geometry.inference import (cluster_bootstrap, percentile_interval,
                                unique_target_permutations)

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
DEPOSITED = REPO / "results" / "03_phenotype_projection" / "phenotype_projection_results.json"
OUT = REPO / "results" / "03c_h3_sensitivity"

# frozen in the registration
MIN_COMMON_RECORDS = 700
DEPOSITED_RECORDS = 795
RECON_TOL = 1e-6
# the corrected artifact, which replaces the superseded 65e5d10e2720... of f288507
EXPECTED_REFERENCE_SHA256 = "fd69e26fc9a3917323065b631688baeab8b283f735c8bf5b16210ba67bd21425"
# The amendment's input table pins lincs_subset.npz and lincs_shrna.npz, which the
# extraction retired under Deviation 12 produced. The rebuild that builds this
# bundle reads GSE92742 directly and never opens either file, so pinning them here
# would assert a dependency the run does not have. The sources it does read are
# pinned instead. Deviation 17.
EXPECTED_SIG_INFO_SHA256 = "19da29c0ee12ddf27f9698cd0da40beaff58657dcde9d382aae068737e831299"
EXPECTED_SHRNA_SIG_INFO_SHA256 = "bd396fa0e1a2f00c1b5f2c8d2b35f9a056f5e5353382475655869038037ec014"
# Stamped on first retrieval rather than registered in advance, as
# results/03c_h3_sensitivity/input_pin_check.json records. It is trust on first
# use and is named as such rather than presented as a registered pin.
EXPECTED_GCTX_SHA256 = "b293f3fb7c2298a60526de727e5400d8400af4b77a23c4ed2116f86199fb45e8"
S1_MIN_EFFECT = 0.20
S3_EQUIV_BOUND = 0.15
N_LANDMARK = 978
N_BOOT = 10_000
N_PERM = 10_000
SEED_BOOT = 20260913          # drug-level, as registered at f288507
SEED_PERM = 20260914          # drug-record permutation, as registered at f288507
SEED_BOOT_TARGET = 20260926   # target-level, added by the amendment
SEED_PERM_TARGET = 20260927   # unique-target permutation, added by the amendment


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def partial_spearman(y, x, covars):
    """Rank residualization, exactly as registered.

    Midrank each variable, regress ranked y and ranked x separately on the ranked
    covariates with an explicit intercept, and correlate the residuals. Not a
    library's partial_corr, whose tie handling and intercept are unstated.
    """
    ry = rankdata(y, method="average")
    rx = rankdata(x, method="average")
    rc = np.column_stack([rankdata(c, method="average") for c in covars])
    design = np.column_stack([np.ones(len(ry)), rc])

    def resid(v):
        beta, *_ = np.linalg.lstsq(design, v, rcond=None)
        return v - design @ beta

    a, b = resid(ry), resid(rx)
    # Under near-perfect collinearity with the covariates both residual vectors
    # are rounding error, and their correlation is a confident number computed
    # from noise. Refuse rather than report it.
    for name, resid_vec, orig in (("y", a, ry), ("x", b, rx)):
        scale = np.linalg.norm(orig - orig.mean())
        if not (np.linalg.norm(resid_vec) > 1e-8 * max(scale, 1.0)):
            raise AssertionError(
                f"{name} is collinear with the covariates; the partial correlation "
                "would be computed from residual noise")
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / denom)


def projected_dispersion(pair_diffs, direction):
    """P for one drug: mean |(s_i - s_j) . u| over distinct context pairs."""
    return float(np.abs(pair_diffs @ direction).mean())


def permuted_statistic(G, row_drug, counts, E_all, covars, perm):
    """The adjusted statistic under one target-assignment permutation.

    A permutation only changes which frozen direction each drug is assigned, so
    projecting every pair difference onto every direction once (G) makes this a
    column lookup rather than a recomputation. Equivalence with direct
    recomputation is asserted at runtime and covered by the tests.
    """
    n = len(counts)
    P_p = np.bincount(row_drug, weights=G[np.arange(len(row_drug)), perm[row_drug]],
                      minlength=n) / counts
    E_p = E_all[np.arange(n), perm]
    return partial_spearman(P_p, E_p, covars), P_p, E_p


def target_representatives(targets, dirs):
    """One drug index per unique target, whose direction stands for that target.

    The permutation reassigns directions by target, while the lookup tables are
    indexed by drug, so each target needs a drug that carries its direction. Every
    drug of a target must already share that direction, which is asserted here
    rather than assumed: the bundle builds a drug's direction from its target.
    """
    unique, first_row = np.unique(targets, return_index=True)
    for target, row in zip(unique, first_row):
        rows = np.flatnonzero(targets == target)
        if not (np.allclose(dirs[rows], dirs[row], atol=1e-12)):
            raise AssertionError(f"drugs annotated to {target} do not share one direction")
    return unique, first_row


def load_bundle(path):
    """The bundle is an untrusted input; every assertion here is repeated from
    the stage that wrote it, because a short or ragged bundle must not reach the
    statistics."""
    z = np.load(path, allow_pickle=True)
    drugs = [str(d) for d in z["drugs"]]
    targets = [str(t) for t in z["targets"]]
    sigs = list(z["signatures"])          # one (K_c, 978) array per drug
    dirs = np.asarray(z["directions"], dtype=np.float64)
    if not (len(drugs) == len(targets) == len(sigs) == len(dirs)):
        raise AssertionError(
            f"ragged bundle: {len(drugs)} drugs, {len(targets)} targets, "
            f"{len(sigs)} matrices, {len(dirs)} directions")
    if not (dirs.shape == (len(drugs), N_LANDMARK)):
        raise AssertionError(f"directions are {dirs.shape}")
    return drugs, targets, sigs, dirs


REQUIRED_MANIFEST_KEYS = ("stage_fingerprint", "fingerprint_parts", "n_bundled",
                          "n_deposited", "cohort_identifier_sha256", "bundle_sha256")
REQUIRED_FINGERPRINT_PARTS = ("gctx", "sig_info", "shrna_sig_info", "deposited",
                              "landmark_order", "loader", "extract_code")


def main(bundle_path, manifest_path):
    bundle_path = Path(bundle_path)
    manifest = json.loads(Path(manifest_path).read_text())
    missing = [k for k in REQUIRED_MANIFEST_KEYS if k not in manifest]
    if not (not missing):
        raise AssertionError(f"rebuild manifest is missing {missing}")
    if not (manifest["bundle_sha256"] == sha256_file(bundle_path)):
        raise AssertionError("the manifest does not describe this bundle")
    if not (manifest["n_deposited"] == DEPOSITED_RECORDS):
        raise AssertionError(f"manifest records {manifest['n_deposited']} deposited, expected {DEPOSITED_RECORDS}")
    # the fingerprint must be the hash of the parts it claims to summarize,
    # and those parts must be the registered ones
    parts = manifest["fingerprint_parts"]
    recomputed_fp = hashlib.sha256(json.dumps(parts, sort_keys=True).encode()).hexdigest()
    if not (manifest["stage_fingerprint"] == recomputed_fp):
        raise AssertionError("the manifest fingerprint is not the hash of its own component hashes")
    absent = [k for k in REQUIRED_FINGERPRINT_PARTS if k not in parts]
    if not (not absent):
        raise AssertionError(f"the manifest fingerprint does not name {absent}")
    if not (parts["deposited"] == EXPECTED_REFERENCE_SHA256):
        raise AssertionError("the manifest pins an artifact other than the corrected one")
    if not (parts["sig_info"] == EXPECTED_SIG_INFO_SHA256):
        raise AssertionError("the manifest pins different compound signature metadata")
    if not (parts["shrna_sig_info"] == EXPECTED_SHRNA_SIG_INFO_SHA256):
        raise AssertionError("the manifest pins different shRNA signature metadata")
    if not (parts["gctx"] == EXPECTED_GCTX_SHA256):
        raise AssertionError("the manifest pins a different GSE92742 release")
    if not (sha256_file(DEPOSITED) == EXPECTED_REFERENCE_SHA256):
        raise AssertionError("the reference artifact is not the corrected one pinned in the amendment")
    drugs, targets, sigs, dirs = load_bundle(bundle_path)
    n = len(drugs)

    # --- the cohort is fixed here, once, before any statistic is computed
    if not (n == len(set(drugs))):
        raise AssertionError("drug identifiers are not unique; the sampling unit is the drug")
    _dep_records = json.loads(DEPOSITED.read_text())
    if not (len(_dep_records) == DEPOSITED_RECORDS):
        raise AssertionError(f"deposited holds {len(_dep_records)} records")
    deposited = {r["drug"]: r for r in _dep_records}
    if not (len(deposited) == DEPOSITED_RECORDS):
        raise AssertionError("deposited drug identifiers are not unique")
    if not (not (set(drugs) - set(deposited))):
        raise AssertionError("rebuilt identifiers absent from the deposited artifact")
    if not (n >= MIN_COMMON_RECORDS):
        raise AssertionError(f"{n} common records, registration requires {MIN_COMMON_RECORDS}")
    if not (manifest["n_bundled"] == n):
        raise AssertionError(f"manifest says {manifest['n_bundled']} bundled, the bundle holds {n}")
    cohort_hash = hashlib.sha256("\n".join(sorted(drugs)).encode()).hexdigest()
    if not (manifest["cohort_identifier_sha256"] == cohort_hash):
        raise AssertionError("the manifest cohort hash does not match the bundle")

    # --- per-drug quantities, and the per-drug reproduction check
    pair_diffs = []
    P, E, D = np.empty(n), np.empty(n), np.empty(n)
    M_delta, M_mean, K = np.empty(n), np.empty(n), np.empty(n)
    worst = {"D": 0.0, "P": 0.0, "E": 0.0}
    for i, (drug, S, u) in enumerate(zip(drugs, sigs, dirs)):
        S = np.asarray(S, dtype=np.float64)
        if not (S.ndim == 2 and S.shape[1] == N_LANDMARK):
            raise AssertionError(f"{drug}: signatures are {S.shape}")
        if not (S.shape[0] >= 5):
            raise AssertionError(f"{drug}: {S.shape[0]} contexts, cohort requires >= 5")
        if not (np.isfinite(S).all() and np.isfinite(u).all()):
            raise AssertionError(f"{drug}: nonfinite input")
        if not (np.isclose(np.linalg.norm(u), 1.0, atol=1e-10)):
            raise AssertionError(f"{drug}: direction is not unit")
        if not (np.all(np.linalg.norm(S, axis=1) > 0)):
            raise AssertionError(f"{drug}: zero-norm signature")
        iu = np.triu_indices(S.shape[0], k=1)
        diffs = S[iu[0]] - S[iu[1]]
        pair_diffs.append(diffs)
        P[i] = projected_dispersion(diffs, u)
        M_delta[i] = float(np.linalg.norm(diffs, axis=1).mean())
        M_mean[i] = float(np.linalg.norm(S, axis=1).mean())
        K[i] = S.shape[0]
        unit = S / np.linalg.norm(S, axis=1, keepdims=True)
        cos = unit @ unit.T
        D[i] = 1.0 - float(cos[iu].mean())
        mean_sig = S.mean(axis=0)
        E[i] = float((mean_sig @ u / (np.linalg.norm(mean_sig) * np.linalg.norm(u))) ** 2)
        if not (np.isfinite([P[i], E[i], D[i], M_delta[i]]).all()):
            raise AssertionError(f"{drug}: nonfinite statistic")
        dep = deposited[drug]
        if not (np.isfinite([dep["raw_instability"], dep["projected_instability"],
                             dep["on_target_enrichment"]]).all()):
            raise AssertionError(f"{drug}: nonfinite deposited value")
        if not (dep["target"] == targets[i]):
            raise AssertionError(f"target assignment differs for {drug}")
        if not (int(dep["n_celllines"]) == S.shape[0]):
            raise AssertionError(f"n_celllines differs for {drug}")
        for key, got, want in (("D", D[i], dep["raw_instability"]),
                               ("P", P[i], dep["projected_instability"]),
                               ("E", E[i], dep["on_target_enrichment"])):
            worst[key] = max(worst[key], abs(got - want))
    for key, w in worst.items():
        if not (w < RECON_TOL):
            raise AssertionError(f"{key} reproduces to {w:.3g}, tolerance {RECON_TOL}")

    covars = (M_delta, K)
    rho_obs = partial_spearman(P, E, covars)
    rho_raw = partial_spearman(D, E, covars)
    rho_meannorm = partial_spearman(P, E, (M_mean, K))     # declared secondary

    # --- S2. Each permutation reassigns which frozen direction a drug gets, so
    # projecting every drug's pair differences onto every direction once turns a
    # permutation into a lookup instead of a recomputation.
    offsets = np.cumsum([0] + [d.shape[0] for d in pair_diffs])
    all_diffs = np.vstack(pair_diffs)
    G = np.abs(all_diffs @ dirs.T)                    # (total_pairs, n_drugs)
    row_drug = np.repeat(np.arange(n), np.diff(offsets))
    counts = np.diff(offsets).astype(np.float64)
    unit_means = np.array([np.asarray(S).mean(axis=0) for S in sigs])
    unit_means /= np.linalg.norm(unit_means, axis=1, keepdims=True)
    E_all = (unit_means @ dirs.T) ** 2                # (n_drugs, n_drugs)

    def stat_under(perm):
        return permuted_statistic(G, row_drug, counts, E_all, covars, perm)[0]

    ident = np.arange(n)
    _, P_ident, E_ident = permuted_statistic(G, row_drug, counts, E_all, covars, ident)
    dP = float(np.abs(P_ident - P).max())
    dE = float(np.abs(E_ident - E).max())
    if not (dP < 1e-10):
        raise AssertionError(f"lookup P departs from direct P by {dP:.3g}")
    if not (dE < 1e-12):
        raise AssertionError(f"lookup E departs from direct E by {dE:.3g}")
    identity = stat_under(ident)
    if not (abs(identity - rho_obs) < 1e-9):
        raise AssertionError(f"identity permutation gives {identity:.12f}, observed {rho_obs:.12f}")

    rng_perm = np.random.default_rng(SEED_PERM)
    null = np.array([stat_under(rng_perm.permutation(n)) for _ in range(N_PERM)])
    if not (np.isfinite(null).all()):
        raise AssertionError("nonfinite permutation replicate")
    p_perm = (1 + int((null >= rho_obs).sum())) / (N_PERM + 1)

    # --- S2-T. The same statistic under a permutation of unique targets: every
    # drug of one target receives the same reassigned direction, so target group
    # sizes and the within-target dependence of the observed mapping survive.
    targets_array = np.asarray(targets)
    unique_targets, representative = target_representatives(targets_array, dirs)
    covariate_fingerprint = (M_delta.tobytes(), K.tobytes())
    null_target = np.array([
        stat_under(representative[assigned])
        for assigned in unique_target_permutations(targets_array, N_PERM, SEED_PERM_TARGET)])
    if not (np.isfinite(null_target).all()):
        raise AssertionError("nonfinite unique-target permutation replicate")
    p_perm_target = (1 + int((null_target >= rho_obs).sum())) / (N_PERM + 1)
    # M_delta is the mean pairwise difference norm and K the number of cell lines;
    # neither uses the target direction, which is why the permutation holds them
    # fixed. Assert it rather than trust it: a covariate that did depend on u would
    # have to be recomputed inside the loop.
    if not ((M_delta.tobytes(), K.tobytes()) == covariate_fingerprint):
        raise AssertionError(
            "a covariate changed while directions were reassigned; it depends on the "
            "target direction and cannot be held fixed under permutation")

    # --- S1 and S3. Percentile bootstrap; ranks and both regressions refit inside.
    rng_boot = np.random.default_rng(SEED_BOOT)
    boot_p = np.empty(N_BOOT)
    boot_raw = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx = rng_boot.integers(0, n, n)
        boot_p[b] = partial_spearman(P[idx], E[idx], (M_delta[idx], K[idx]))
        boot_raw[b] = partial_spearman(D[idx], E[idx], (M_delta[idx], K[idx]))
    if not (np.isfinite(boot_p).all() and np.isfinite(boot_raw).all()):
        raise AssertionError("nonfinite bootstrap replicate")

    # --- S1-T and S3-T. The same two statistics, resampling targets rather than
    # drugs, so that drugs sharing a reference direction travel together.
    def _paired_partials(idx):
        return (partial_spearman(P[idx], E[idx], (M_delta[idx], K[idx])),
                partial_spearman(D[idx], E[idx], (M_delta[idx], K[idx])))

    boot_target = cluster_bootstrap(targets_array, _paired_partials, N_BOOT, SEED_BOOT_TARGET)
    boot_p_target, boot_raw_target = boot_target[:, 0], boot_target[:, 1]

    # a zero denominator makes that drug's ratio missing, not the whole description
    valid = M_delta > 0
    ratio = P[valid] / M_delta[valid]
    ratio_stats = {
        "n": int(valid.sum()), "n_missing_M_delta_zero": int((~valid).sum()),
        "mean": float(ratio.mean()), "median": float(np.median(ratio)),
        "min": float(ratio.min()), "max": float(ratio.max()),
        "spearman_with_E": float(spearmanr(ratio, E[valid]).statistic)}

    ci_s1 = [float(np.percentile(boot_p, 2.5)), float(np.percentile(boot_p, 97.5))]
    ci_s3 = [float(np.percentile(boot_raw, 5.0)), float(np.percentile(boot_raw, 95.0))]
    s1 = bool(rho_obs >= S1_MIN_EFFECT and ci_s1[0] > 0)
    s2 = bool(p_perm < 0.01)
    s3 = bool(-S3_EQUIV_BOUND < ci_s3[0] and ci_s3[1] < S3_EQUIV_BOUND)

    ci_s1_target = list(percentile_interval(boot_p_target, 95.0))
    ci_s3_target = list(percentile_interval(boot_raw_target, 90.0))
    s1_target = bool(rho_obs >= S1_MIN_EFFECT and ci_s1_target[0] > 0)
    s2_target = bool(p_perm_target < 0.01)
    s3_target = bool(-S3_EQUIV_BOUND < ci_s3_target[0] and ci_s3_target[1] < S3_EQUIV_BOUND)

    out = {
        "registration": {"file": "experiments/PREREG_H3_MAGNITUDE_AND_SHARED_AXIS.md",
                         "frozen_at": "f288507"},
        "cohort": {"n": n, "deposited": DEPOSITED_RECORDS,
                   "restricted_to_rebuilt_subset": n < DEPOSITED_RECORDS,
                   "missing_from_rebuild": sorted(set(deposited) - set(drugs)),
                   "identifier_hash": hashlib.sha256("\n".join(sorted(drugs)).encode()).hexdigest(),
                   "max_abs_reproduction_error": {k: float(v) for k, v in worst.items()}},
        "S1_magnitude_and_coverage": {
            "partial_rho": rho_obs, "ci95_percentile": ci_s1,
            "display": f"{rho_obs:.4f} [{ci_s1[0]:.4f}, {ci_s1[1]:.4f}]",
            "min_effect": S1_MIN_EFFECT, "holds": s1,
            "unit": "drug, as registered at f288507"},
        "S2_shared_axis_permutation": {
            "rho_observed": rho_obs, "p_perm": p_perm,
            "null_pct": {"p2.5": float(np.percentile(null, 2.5)),
                         "p50": float(np.percentile(null, 50)),
                         "p97.5": float(np.percentile(null, 97.5))},
            "n_permutations": N_PERM, "holds": s2,
            "unit": "drug record, as registered at f288507"},
        "S3_raw_equivalence": {
            "partial_rho": rho_raw, "ci90_percentile": ci_s3,
            "display": f"{rho_raw:.4f} [{ci_s3[0]:.4f}, {ci_s3[1]:.4f}]",
            "equivalence_bound": S3_EQUIV_BOUND, "holds": s3,
            "unit": "drug, as registered at f288507"},
        "S1T_magnitude_and_coverage": {
            "partial_rho": rho_obs, "ci95_percentile": ci_s1_target,
            "display": f"{rho_obs:.4f} [{ci_s1_target[0]:.4f}, {ci_s1_target[1]:.4f}]",
            "min_effect": S1_MIN_EFFECT, "holds": s1_target,
            "unit": "target cluster"},
        "S2T_unique_target_permutation": {
            "rho_observed": rho_obs, "p_perm": p_perm_target,
            "null_pct": {"p2.5": float(np.percentile(null_target, 2.5)),
                         "p50": float(np.percentile(null_target, 50)),
                         "p97.5": float(np.percentile(null_target, 97.5))},
            "attenuation_vs_null_median": float(rho_obs - np.median(null_target)),
            "n_permutations": N_PERM, "n_unique_targets": int(len(unique_targets)),
            "holds": s2_target, "unit": "unique target"},
        "S3T_raw_equivalence": {
            "partial_rho": rho_raw, "ci90_percentile": ci_s3_target,
            "display": f"{rho_raw:.4f} [{ci_s3_target[0]:.4f}, {ci_s3_target[1]:.4f}]",
            "equivalence_bound": S3_EQUIV_BOUND, "holds": s3_target,
            "unit": "target cluster"},
        "identity_check": {"max_abs_P": dP, "max_abs_E": dE,
                           "abs_partial_rho": abs(identity - rho_obs)},
        # The amendment frozen at 7f57136 puts the gate on the target-level
        # versions: where a drug-level verdict and its target-level counterpart
        # disagree, the target-level verdict governs.
        "gate": {"all_three_required": True,
                 "governing_unit": "target",
                 "h3_interpretation_survives": bool(s1_target and s2_target and s3_target),
                 "drug_level_verdict": bool(s1 and s2 and s3),
                 "units_agree": bool((s1, s2, s3) == (s1_target, s2_target, s3_target))},
        "secondary_no_criterion": {
            "partial_rho_given_mean_signature_norm": rho_meannorm,
            "P_over_M_delta": ratio_stats},
        "provenance": {
            "rebuild_manifest": manifest,
            "bundle_sha256": sha256_file(bundle_path),
            "deposited_artifact_sha256": sha256_file(DEPOSITED),
            "script_sha256": sha256_file(Path(__file__)),
            "seeds": {"bootstrap": SEED_BOOT, "permutation": SEED_PERM,
                      "bootstrap_target": SEED_BOOT_TARGET,
                      "permutation_target": SEED_PERM_TARGET},
            "numpy": np.__version__, "scipy": scipy.__version__,
            "python": platform.python_version()},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    # the replicate arrays travel with the result, so the intervals and the
    # p-value can be recomputed without rerunning
    draws = OUT / "h3_sensitivity_draws.npz"
    np.savez_compressed(draws, permutation_null=null, bootstrap_P=boot_p, bootstrap_raw=boot_raw,
                        permutation_null_target=null_target,
                        bootstrap_P_target=boot_p_target, bootstrap_raw_target=boot_raw_target)
    out["provenance"]["draws_sha256"] = sha256_file(draws)

    path = OUT / "h3_sensitivity_results.json"
    path.write_text(json.dumps(out, indent=2))
    # a file cannot contain its own digest; the hash goes in a sidecar
    (OUT / "h3_sensitivity_results.json.sha256").write_text(sha256_file(path) + "\n")
    print(json.dumps({k: out[k] for k in
                      ("cohort", "S1_magnitude_and_coverage", "S2_shared_axis_permutation",
                       "S3_raw_equivalence", "S1T_magnitude_and_coverage",
                       "S2T_unique_target_permutation", "S3T_raw_equivalence", "gate")}, indent=2))
    print(f"\nwritten to {path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--manifest", required=True,
                    help="manifest written by the rebuild stage")
    a = ap.parse_args()
    main(a.bundle, a.manifest)
