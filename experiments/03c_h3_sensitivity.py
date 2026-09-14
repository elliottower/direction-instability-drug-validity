"""H3 sensitivity analyses S1-S3, registered in PREREG_H3_MAGNITUDE_AND_SHARED_AXIS.md
(frozen at f288507, before the LINCS matrices were rebuilt).

S1  the association survives adjustment for pairwise-difference magnitude and coverage
S2  it exceeds a null that permutes target assignment while keeping the shared axis
S3  the raw score stays practically equivalent to zero under the same adjustment

All three are required; none compensates for another.

Input is the bundle written by the rebuild stage: per-drug signature matrices, the
unit target direction, and the identifiers, for drugs that reproduced the deposited
artifact. Nothing here re-derives the cohort; it fails closed if the bundle and the
deposited artifact disagree.

    uv run --no-project --with numpy --with scipy python experiments/03c_h3_sensitivity.py \
        --bundle results/03c_h3_sensitivity/cohort_bundle.npz \
        --manifest results/03c_h3_sensitivity/rebuild_manifest.json
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import rankdata, spearmanr

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
DEPOSITED = REPO / "results" / "03_phenotype_projection" / "phenotype_projection_results.json"
OUT = REPO / "results" / "03c_h3_sensitivity"

# frozen in the registration
MIN_COMMON_RECORDS = 700
DEPOSITED_RECORDS = 795
RECON_TOL = 1e-6
EXPECTED_DEPOSITED_SHA256 = "65e5d10e272037987384f89e6208de478fa17f6b8fb946add893e5a24c2d80c4"
EXPECTED_LOADER_SHA256 = "b7ec2cd46be4a800dc546594e82491f785fdf504d4bd8bc783be293b343a1611"
EXPECTED_SHARD_SIZE = 250
EXPECTED_MIN_HAIRPINS = 3
S1_MIN_EFFECT = 0.20
S3_EQUIV_BOUND = 0.15
N_LANDMARK = 978
N_BOOT = 10_000
N_PERM = 10_000
SEED_BOOT = 20260913
SEED_PERM = 20260914


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
        assert np.linalg.norm(resid_vec) > 1e-8 * max(scale, 1.0), (
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


def load_bundle(path):
    """The bundle is an untrusted input; every assertion here is repeated from
    the stage that wrote it, because a short or ragged bundle must not reach the
    statistics."""
    z = np.load(path, allow_pickle=True)
    drugs = [str(d) for d in z["drugs"]]
    targets = [str(t) for t in z["targets"]]
    sigs = list(z["signatures"])          # one (K_c, 978) array per drug
    dirs = np.asarray(z["directions"], dtype=np.float64)
    assert len(drugs) == len(targets) == len(sigs) == len(dirs), (
        f"ragged bundle: {len(drugs)} drugs, {len(targets)} targets, "
        f"{len(sigs)} matrices, {len(dirs)} directions")
    assert dirs.shape == (len(drugs), N_LANDMARK), f"directions are {dirs.shape}"
    return drugs, targets, sigs, dirs


REQUIRED_MANIFEST_KEYS = ("stage_fingerprint", "fingerprint_parts", "n_bundled",
                          "n_deposited", "cohort_identifier_sha256", "bundle_sha256",
                          "shrna_consensus_sha256", "shrna_hairpins_sha256")


def main(bundle_path, manifest_path):
    bundle_path = Path(bundle_path)
    manifest = json.loads(Path(manifest_path).read_text())
    missing = [k for k in REQUIRED_MANIFEST_KEYS if k not in manifest]
    assert not missing, f"rebuild manifest is missing {missing}"
    assert manifest["bundle_sha256"] == sha256_file(bundle_path), (
        "the manifest does not describe this bundle")
    assert manifest["n_deposited"] == DEPOSITED_RECORDS, (
        f"manifest records {manifest['n_deposited']} deposited, expected {DEPOSITED_RECORDS}")
    # the fingerprint must be the hash of the parts it claims to summarize,
    # and those parts must be the registered ones
    parts = manifest["fingerprint_parts"]
    recomputed_fp = hashlib.sha256(json.dumps(parts, sort_keys=True).encode()).hexdigest()
    assert manifest["stage_fingerprint"] == recomputed_fp, (
        "the manifest fingerprint is not the hash of its own component hashes")
    assert parts["deposited"] == EXPECTED_DEPOSITED_SHA256, "manifest pins a different deposited artifact"
    assert parts["loader"] == EXPECTED_LOADER_SHA256, "manifest pins a loader other than the registered one"
    assert parts["shard_size"] == EXPECTED_SHARD_SIZE, f"manifest shard size {parts['shard_size']}"
    assert parts["min_hairpins"] == EXPECTED_MIN_HAIRPINS, f"manifest hairpin minimum {parts['min_hairpins']}"
    assert sha256_file(DEPOSITED) == EXPECTED_DEPOSITED_SHA256, (
        "the deposited artifact is not the one pinned in the registration")
    drugs, targets, sigs, dirs = load_bundle(bundle_path)
    n = len(drugs)

    # --- the cohort is fixed here, once, before any statistic is computed
    assert n == len(set(drugs)), "drug identifiers are not unique; the sampling unit is the drug"
    _dep_records = json.loads(DEPOSITED.read_text())
    assert len(_dep_records) == DEPOSITED_RECORDS, f"deposited holds {len(_dep_records)} records"
    deposited = {r["drug"]: r for r in _dep_records}
    assert len(deposited) == DEPOSITED_RECORDS, "deposited drug identifiers are not unique"
    assert not (set(drugs) - set(deposited)), "rebuilt identifiers absent from the deposited artifact"
    assert n >= MIN_COMMON_RECORDS, f"{n} common records, registration requires {MIN_COMMON_RECORDS}"
    assert manifest["n_bundled"] == n, (
        f"manifest says {manifest['n_bundled']} bundled, the bundle holds {n}")
    cohort_hash = hashlib.sha256("\n".join(sorted(drugs)).encode()).hexdigest()
    assert manifest["cohort_identifier_sha256"] == cohort_hash, (
        "the manifest cohort hash does not match the bundle")

    # --- per-drug quantities, and the per-drug reproduction check
    pair_diffs = []
    P, E, D = np.empty(n), np.empty(n), np.empty(n)
    M_delta, M_mean, K = np.empty(n), np.empty(n), np.empty(n)
    worst = {"D": 0.0, "P": 0.0, "E": 0.0}
    for i, (drug, S, u) in enumerate(zip(drugs, sigs, dirs)):
        S = np.asarray(S, dtype=np.float64)
        assert S.ndim == 2 and S.shape[1] == N_LANDMARK, f"{drug}: signatures are {S.shape}"
        assert S.shape[0] >= 5, f"{drug}: {S.shape[0]} contexts, cohort requires >= 5"
        assert np.isfinite(S).all() and np.isfinite(u).all(), f"{drug}: nonfinite input"
        assert np.isclose(np.linalg.norm(u), 1.0, atol=1e-10), f"{drug}: direction is not unit"
        assert np.all(np.linalg.norm(S, axis=1) > 0), f"{drug}: zero-norm signature"
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
        assert np.isfinite([P[i], E[i], D[i], M_delta[i]]).all(), f"{drug}: nonfinite statistic"
        dep = deposited[drug]
        assert np.isfinite([dep["raw_instability"], dep["projected_instability"],
                            dep["on_target_enrichment"]]).all(), f"{drug}: nonfinite deposited value"
        assert dep["target"] == targets[i], f"target assignment differs for {drug}"
        assert int(dep["n_celllines"]) == S.shape[0], f"n_celllines differs for {drug}"
        for key, got, want in (("D", D[i], dep["raw_instability"]),
                               ("P", P[i], dep["projected_instability"]),
                               ("E", E[i], dep["on_target_enrichment"])):
            worst[key] = max(worst[key], abs(got - want))
    for key, w in worst.items():
        assert w < RECON_TOL, f"{key} reproduces to {w:.3g}, tolerance {RECON_TOL}"

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
    assert dP < 1e-10, f"lookup P departs from direct P by {dP:.3g}"
    assert dE < 1e-12, f"lookup E departs from direct E by {dE:.3g}"
    identity = stat_under(ident)
    assert abs(identity - rho_obs) < 1e-9, (
        f"identity permutation gives {identity:.12f}, observed {rho_obs:.12f}")

    rng_perm = np.random.default_rng(SEED_PERM)
    null = np.array([stat_under(rng_perm.permutation(n)) for _ in range(N_PERM)])
    assert np.isfinite(null).all(), "nonfinite permutation replicate"
    p_perm = (1 + int((null >= rho_obs).sum())) / (N_PERM + 1)

    # --- S1 and S3. Percentile bootstrap; ranks and both regressions refit inside.
    rng_boot = np.random.default_rng(SEED_BOOT)
    boot_p = np.empty(N_BOOT)
    boot_raw = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx = rng_boot.integers(0, n, n)
        boot_p[b] = partial_spearman(P[idx], E[idx], (M_delta[idx], K[idx]))
        boot_raw[b] = partial_spearman(D[idx], E[idx], (M_delta[idx], K[idx]))
    assert np.isfinite(boot_p).all() and np.isfinite(boot_raw).all(), "nonfinite bootstrap replicate"

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
            "min_effect": S1_MIN_EFFECT, "holds": s1},
        "S2_shared_axis_permutation": {
            "rho_observed": rho_obs, "p_perm": p_perm,
            "null_pct": {"p2.5": float(np.percentile(null, 2.5)),
                         "p50": float(np.percentile(null, 50)),
                         "p97.5": float(np.percentile(null, 97.5))},
            "n_permutations": N_PERM, "holds": s2},
        "S3_raw_equivalence": {
            "partial_rho": rho_raw, "ci90_percentile": ci_s3,
            "display": f"{rho_raw:.4f} [{ci_s3[0]:.4f}, {ci_s3[1]:.4f}]",
            "equivalence_bound": S3_EQUIV_BOUND, "holds": s3},
        "identity_check": {"max_abs_P": dP, "max_abs_E": dE,
                           "abs_partial_rho": abs(identity - rho_obs)},
        "gate": {"all_three_required": True,
                 "h3_interpretation_survives": bool(s1 and s2 and s3)},
        "secondary_no_criterion": {
            "partial_rho_given_mean_signature_norm": rho_meannorm,
            "P_over_M_delta": ratio_stats},
        "provenance": {
            "rebuild_manifest": manifest,
            "bundle_sha256": sha256_file(bundle_path),
            "deposited_artifact_sha256": sha256_file(DEPOSITED),
            "script_sha256": sha256_file(Path(__file__)),
            "seeds": {"bootstrap": SEED_BOOT, "permutation": SEED_PERM},
            "numpy": np.__version__, "scipy": scipy.__version__,
            "python": platform.python_version()},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    # the replicate arrays travel with the result, so the intervals and the
    # p-value can be recomputed without rerunning
    draws = OUT / "h3_sensitivity_draws.npz"
    np.savez_compressed(draws, permutation_null=null, bootstrap_P=boot_p, bootstrap_raw=boot_raw)
    out["provenance"]["draws_sha256"] = sha256_file(draws)

    path = OUT / "h3_sensitivity_results.json"
    path.write_text(json.dumps(out, indent=2))
    # a file cannot contain its own digest; the hash goes in a sidecar
    (OUT / "h3_sensitivity_results.json.sha256").write_text(sha256_file(path) + "\n")
    print(json.dumps({k: out[k] for k in
                      ("cohort", "S1_magnitude_and_coverage", "S2_shared_axis_permutation",
                       "S3_raw_equivalence", "gate")}, indent=2))
    print(f"\nwritten to {path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--manifest", required=True,
                    help="manifest written by the rebuild stage")
    a = ap.parse_args()
    main(a.bundle, a.manifest)
