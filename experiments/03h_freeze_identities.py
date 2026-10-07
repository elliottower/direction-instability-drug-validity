"""Regenerate the frozen identities and the implementation manifest from the files.

Hashes that are patched by hand go stale in places nobody looks. The implementation
manifest had `03e` at the superseded gate's hash in one block and the live hash in
another, and its commands still named an interface the code no longer accepts. Both
artifacts are therefore generated in one deterministic step from whatever is actually
on disk, and an artifact that does not exist yet is written as a refusal to pretend,
not as an absent key.

Two outputs, deliberately apart:

  registry/frozen/expected_identities.json   what the gate compares against
  results/03d_.../implementation_manifest_amendment_3.json   what the freeze records

The first exists because an artifact carrying its own expected hash attests to itself.
The gate reads its expectations from here, and Amendment 3's B8 pins this file.
"""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

SCHEMA_VERSION = 1
PENDING = "PENDING"

REBUILD_ARTIFACTS = ("landmark_gene_ids.json", "shrna_consensus.npz",
                     "shrna_signatures.npz", "rebuild_manifest.json")
MANIFESTS = ("cohort_795_shrna_paired.json", "cohort_812_compound.json",
             "shrna_eligible_targets.json", "shrna_excluded_records.json")
CODE = ("experiments/03e_reconstruction_gate.py",
        "experiments/03f_shrna_coverage_identity.py",
        "experiments/03g_build_cohort_manifests.py",
        "experiments/03h_freeze_identities.py",
        "experiments/03i_freeze_analysis_bases.py",
        "experiments/03d_h3_reference_discordance.py",
        "experiments/modal_h3_execute.py",
        "experiments/modal_h3_rebuild.py",
        "geometry/inference.py", "geometry/references.py",
        "geometry/single_cell.py", "geometry/direction_instability.py")


class FreezeError(AssertionError):
    """A refusal, raised so `python -O` cannot strip it."""


def require(condition, message) -> None:
    if not condition:
        raise FreezeError(message)


def _relative(path, repo) -> str:
    """The path as written where it sits under the repo, and as given where it does not."""
    try:
        return str(Path(path).resolve().relative_to(Path(repo).resolve()))
    except ValueError:
        return str(path)


def canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def sha256_file(path) -> str:
    path = Path(path)
    if not path.exists():
        return PENDING
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            sha.update(block)
    return sha.hexdigest()


def git_head() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                             check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return PENDING
    dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True,
                           text=True).stdout.strip()
    return out.stdout.strip() + ("+dirty" if dirty else "")


# the four constructions Amendment 4 freezes a basis for, and the fields of each
# record that are pinned here rather than inside the artifact
BASIS_CONSTRUCTIONS = ("C1-K562", "C1-RPE1", "C1-GW", "C1-GW-phenotype-positive")
BASIS_PINNED_FIELDS = ("file", "file_sha256", "basis_sha256", "n_landmarks_matched",
                       "n_landmarks_in_basis")


def analysis_bases_pin(repo: Path, bases: Path) -> dict:
    """Amendment 4's bases, pinned from outside the file that holds them.

    The artifact already hashes each ordered basis, but a list edited together
    with its own recomputed digest passes that check, so these values live apart
    from it and Amendment 4's text pins this file's own sha256.

    They are not added to `expected_identities.json`: Amendment 3 pins that file
    by sha256 at `3de491d6234c9dc9f56072227b50e7a513fdc25720a9bc113ef56c0939e3b9b2`
    and a frozen pin is not rewritten to carry a later amendment's values.
    """
    generator = repo / "experiments/03i_freeze_analysis_bases.py"
    if not Path(bases).exists():
        return {"file_sha256": PENDING, "generator_sha256": sha256_file(generator),
                "schema_version": PENDING, "landmark_axis_sha256": PENDING,
                "constructions": {name: PENDING for name in BASIS_CONSTRUCTIONS}}
    document = json.loads(Path(bases).read_text())
    require(sorted(document["constructions"]) == sorted(BASIS_CONSTRUCTIONS),
            f"the bases artifact holds {sorted(document['constructions'])}, not "
            f"{sorted(BASIS_CONSTRUCTIONS)}")
    return {
        "name": "analysis_bases_pin", "schema_version": SCHEMA_VERSION,
        "what_this_is": ("the expected identity of Amendment 4's analysis bases, held "
                         "outside the artifact it describes, because an artifact carrying "
                         "its own expected hash attests to itself"),
        "registered_by": "experiments/PREREG_H3_S1S3_AMENDMENT_4.md",
        "file": _relative(bases, repo),
        "file_sha256": sha256_file(bases),
        "generator_sha256": sha256_file(generator),
        "schema_version": document["schema_version"],
        "landmark_axis_sha256": document["gene_info"]["landmark_axis_sha256"],
        "landmark_axis_file_sha256": sha256_file(Path(bases).parent / "landmark_gene_ids.json"),
        "constructions": {
            name: {field: document["constructions"][name][field]
                   for field in BASIS_PINNED_FIELDS}
            for name in BASIS_CONSTRUCTIONS},
        "pooled_construction": {
            field: document["pooled_construction"][field]
            for field in ("file", "file_sha256", "n_landmarks_declared",
                          "n_landmarks_measured", "declared_sha256", "measured_sha256")},
        "registered_comparisons": {key: row["n_shared"] for key, row
                                   in document["registered_comparisons"].items()},
    }


def frozen_identities(repo: Path, rebuilt: Path, coverage: Path, axis_map: Path,
                      deposited: Path, legacy: Path, registry: Path) -> dict:
    mapping = {}
    if Path(axis_map).exists():
        payload = json.loads(Path(axis_map).read_text())
        order = payload["declared_index_of_each_actual_column"]
        mapping = {"axis_map_file_sha256": sha256_file(axis_map),
                   "axis_map_content_sha256": hashlib.sha256(
                       canonical(order).encode()).hexdigest()}
        require(mapping["axis_map_content_sha256"] == payload["map_sha256"],
                "the map's permutation does not hash to its own map_sha256")
    else:
        mapping = {"axis_map_file_sha256": PENDING, "axis_map_content_sha256": PENDING}

    return {
        "name": "frozen_identities", "schema_version": SCHEMA_VERSION,
        "what_this_is": ("the expected hashes the gate compares against, held outside "
                         "every artifact they describe, because an artifact carrying its "
                         "own expected hash attests to itself"),
        "coverage_report_sha256": sha256_file(coverage),
        "coverage_code_sha256": sha256_file(repo / "experiments/03f_shrna_coverage_identity.py"),
        "deposited_sha256": sha256_file(deposited),
        "legacy_extraction_sha256": sha256_file(legacy),
        **mapping,
        "manifest_sha256": {name: sha256_file(registry / name) for name in MANIFESTS},
        "rebuild_sha256": {
            **{name: sha256_file(Path(rebuilt) / name) for name in REBUILD_ARTIFACTS},
            **{shard.name: sha256_file(shard)
               for shard in sorted(Path(rebuilt).glob("shard_*.npz"))}},
    }


def passing_reports(repo: Path) -> dict:
    """Both gate v2 reports, with the predicate each one held and its own hash.

    A report is recorded here only when its release condition held. A failed report
    stays on disk and in the ledger as a failed report, and is not listed as passing.
    """
    found = {}
    for directory in sorted((repo / "results/03c_h3_sensitivity").glob("gate_v2_*")):
        report = directory / "reconstruction_gate.json"
        if not report.exists():
            continue
        payload = json.loads(report.read_text())
        if not payload.get("release_authorized"):
            continue
        found[payload["cohort_scope"]] = {
            "report": _relative(report, repo),
            "sha256": sha256_file(report),
            "n_drugs": payload["cohort"]["n_drugs"],
            "cohort_pair_sha256": payload["cohort"]["pair_sha256"],
            "gate_code_sha256": payload["gate_code_sha256"],
            "predicate_held": payload["predicate_held"],
            "tolerances": {"rtol": payload["compound_source_reconstruction"]["rtol"],
                           "atol": payload["compound_source_reconstruction"]["atol"]},
            "compound_max_abs_difference":
                payload["compound_source_reconstruction"]["max_abs_difference"],
            "shrna": {k: payload["shrna_reconstruction_on_the_eligible_set"][k]
                      for k in ("n_signatures", "n_targets",
                                "max_abs_signature_difference",
                                "max_abs_consensus_difference",
                                "max_abs_stored_direction_difference")},
            "legacy": {
                "raw_exact_over_n_signatures":
                    payload["legacy_artifact_integrity"]["raw_signature_reproduction"],
                "under_the_map":
                    payload["legacy_artifact_integrity"]["under_the_map"],
                "under_the_declared_labels":
                    payload["legacy_artifact_integrity"]["under_the_declared_labels"]},
        }
    return found


def implementation_manifest(repo: Path, frozen: dict, registry: Path,
                            coverage: Path) -> dict:
    index = json.loads((registry / "INDEX.json").read_text())
    gate_v1 = repo / "experiments/superseded/03e_reconstruction_gate_v1.py"
    amendment = repo / "experiments/PREREG_H3_S1S3_AMENDMENT_3_DRAFT.md"
    report = json.loads(Path(coverage).read_text()) if Path(coverage).exists() else {}

    return {
        "what_this_binds": (
            "the production identity Amendment 3's B8 requires. Generated from the files "
            "by experiments/03h_freeze_identities.py, never patched by hand."),
        "generated_from_git": git_head(),
        "supersedes_in_part": {
            "file": "results/03d_h3_reference_discordance/implementation_manifest.json",
            "sha256": sha256_file(repo / "results/03d_h3_reference_discordance"
                                       / "implementation_manifest.json"),
            "what": ("its `inputs` entries naming lincs_subset.npz and lincs_shrna.npz as "
                     "production inputs, and its `commands_in_order` steps 3 and 6, whose "
                     "commands name an interface gate v2 does not accept")},
        "amendment": {
            "file": "experiments/PREREG_H3_S1S3_AMENDMENT_3_DRAFT.md",
            "draft_sha256": sha256_file(amendment),
            "freeze_commit": PENDING,
            "amends": {"file": "experiments/PREREG_H3_S1S3_AMENDMENT_2.md",
                       "commit": "fd1ae8d", "clause_replaced": "A2"}},
        "frozen_identities": {
            "file": "registry/frozen/expected_identities.json",
            "sha256": PENDING,       # written after the file itself, below
            "why": ("the gate's expectations live here and not inside the artifacts they "
                    "describe")},
        "cohort_manifests": {
            "directory": "registry/cohorts",
            "file_sha256": index["manifest_sha256"],
            "identity": index["identity"],
            "built_by": "experiments/03g_build_cohort_manifests.py",
            "builder_sha256": sha256_file(repo / "experiments/03g_build_cohort_manifests.py"),
            "why": ("the cohorts were results files; phenotype_projection_results.json "
                    "carries projected instability and on-target enrichment, so pinning it "
                    "coupled eligibility to derived values")},
        "shrna_coverage": {
            "measured_in": _relative(coverage, repo),
            "measurement_sha256": frozen["coverage_report_sha256"],
            "code_sha256": frozen["coverage_code_sha256"],
            "eligible_targets": report.get("counts", {}).get("eligible_targets"),
            "paired_records": report.get("counts", {}).get("valued"),
            "excluded_records": report.get("counts", {}).get("excluded_records"),
            "excluded_targets": report.get("excluded_targets"),
            "consensus_used": report.get("inputs", {}).get("consensus"),
            "how_the_quantities_are_missing": (
                "the two shRNA keys are absent from the seventeen records rather than "
                "present and null; the file holds no explicit nulls, and eligibility means "
                "both quantities exist and are numeric, whichever way missingness is "
                "represented"),
            "superseded": {
                "file": ("results/superseded/"
                         "shrna_coverage_identity_v2_wrong_consensus_2026-10-06.json"),
                "why": ("computed against results/03c_h3_sensitivity/shrna_consensus.npz, a "
                        "partial copy without gene_ids, rather than the authoritative "
                        "rebuild consensus B1 names. The 258-target set, its hash and every "
                        "set equality are identical; only the provenance was wrong.")}},
        "gate": {
            "version": 2,
            "code": "experiments/03e_reconstruction_gate.py",
            "code_sha256": sha256_file(repo / "experiments/03e_reconstruction_gate.py"),
            "v1_preserved": "experiments/superseded/03e_reconstruction_gate_v1.py",
            "v1_sha256": sha256_file(gate_v1),
            "release_predicates": ["cohort_identity", "compound_source_reconstruction",
                                   "shrna_reconstruction_on_the_eligible_set",
                                   "shrna_eligibility_partition",
                                   "legacy_artifact_integrity"],
            "axis_map": {"file_sha256": frozen["axis_map_file_sha256"],
                         "content_sha256": frozen["axis_map_content_sha256"]},
            "passing_reports": passing_reports(repo) or PENDING,
            "output_layout": ("one immutable directory per run, "
                              "gate_v2_<scope>_<timestamp>/reconstruction_gate.json, "
                              "because the report filename is fixed and two scopes sharing "
                              "a directory would overwrite one another")},
        "code_sha256": {name: sha256_file(repo / name) for name in CODE},
        "commands_in_order": [
            {"step": 1, "what": "the cohort manifests",
             "command": ("python experiments/03g_build_cohort_manifests.py --shrna-paired "
                         "results/03_phenotype_projection/phenotype_projection_results.json "
                         "--compound results/03b_h3_crispri/crispri_cohort.json --out-dir "
                         "registry/cohorts --index registry/cohorts/INDEX.json")},
            {"step": 2, "what": "the shRNA eligibility measurement",
             "command": ("python experiments/03f_shrna_coverage_identity.py --deposited "
                         "results/03b_h3_crispri/h3_crispri_results.json --consensus "
                         "results/03c_h3_sensitivity/gctx_rebuild/shrna_consensus.npz "
                         "--eligible registry/cohorts/shrna_eligible_targets.json "
                         "--excluded registry/cohorts/shrna_excluded_records.json "
                         "--compound-cohort registry/cohorts/cohort_812_compound.json --out "
                         "results/03c_h3_sensitivity/shrna_coverage_identity_v2.json")},
            {"step": 3, "what": "the legacy coordinate map, where the pinned GCTX is",
             "command": "modal run experiments/modal_h3_execute.py::stage_pin_the_legacy_axis_map"},
            {"step": 4, "what": "the frozen identities and this manifest",
             "command": "python experiments/03h_freeze_identities.py (see --help)"},
            {"step": 5, "what": "gate v2, once per scope, into separate directories",
             "command": ("python experiments/03e_reconstruction_gate.py --rebuilt "
                         "results/03c_h3_sensitivity/gctx_rebuild --extraction "
                         "data/lincs_extraction --cohort-manifest "
                         "registry/cohorts/cohort_812_compound.json --paired-manifest "
                         "registry/cohorts/cohort_795_shrna_paired.json --eligible-targets "
                         "registry/cohorts/shrna_eligible_targets.json --excluded-records "
                         "registry/cohorts/shrna_excluded_records.json --axis-map "
                         "registry/frozen/legacy_axis_map.json --coverage-report "
                         "results/03c_h3_sensitivity/shrna_coverage_identity_v2.json "
                         "--deposited-records results/03b_h3_crispri/h3_crispri_results.json "
                         "--frozen-identities registry/frozen/expected_identities.json "
                         "--cohort-scope {compound|shrna_paired} --gctx <gctx> "
                         "--shrna-siginfo <s> --compound-siginfo <c> --gene-info <g> "
                         "--output results/03c_h3_sensitivity"),
             "note": "run on Modal, via experiments/modal_h3_execute.py"},
            {"step": 6, "what": "R0-R7, only after both scopes pass",
             "command": ("python experiments/03d_h3_reference_discordance.py --r5-stage "
                         "response --rebuilt results/03c_h3_sensitivity/gctx_rebuild "
                         "--expected-mapping-sha256 152361cb3174a5fb7aae0229c3e3d049dc00d49"
                         "d9e442925156a9fe0564b89d3 ...")}],
        "state": (
            "Gate v1 failed on both registered cohort invocations. Its production compound "
            "comparison passed on both cohorts; its shRNA comparison passed on the 795 "
            "cohort and was not evaluated on the 812 cohort, because the gate treated the "
            "compound-complete target universe as shRNA-complete; the retained legacy "
            "extraction failed the v1 equality predicate on both cohorts. Before freeze and "
            "before any R0-R7 statistic, Amendment 3 specifies corrected eligibility and "
            "legacy-custody predicates. No v1 failure is erased or relabeled."),
        "retained_not_consumed": {
            "cohort_bundle.npz": ("present in the rebuild and read by no gate predicate; "
                                  "recorded here as provenance rather than pinned as a "
                                  "gate input"),
            "results/03c_h3_sensitivity/shrna_consensus.npz": (
                "a partial copy without gene_ids; the authoritative consensus is the one "
                "under gctx_rebuild/ that B1 names")},
        "open_before_the_freeze": (
            ["this file's own sha256 written into the freeze commit message"]
            if len(passing_reports(repo)) == 2 else
            ["a passing gate v2 report on each of the two cohort scopes"]),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--rebuilt", type=Path)
    parser.add_argument("--registry", type=Path)
    parser.add_argument("--coverage", type=Path)
    parser.add_argument("--axis-map", type=Path)
    parser.add_argument("--deposited", type=Path)
    parser.add_argument("--legacy", type=Path)
    parser.add_argument("--analysis-bases", type=Path,
                        help="registry/frozen/analysis_bases.json, pinned from outside itself")
    parser.add_argument("--bases-pin-out", type=Path,
                        help="where to write that pin. It is a separate file because "
                             "Amendment 3 pins expected_identities.json by sha256.")
    parser.add_argument("--frozen-out", type=Path)
    parser.add_argument("--manifest-out", type=Path)
    args = parser.parse_args(argv)

    if args.frozen_out is None and args.bases_pin_out is None:
        raise FreezeError("nothing to write: name --frozen-out, --bases-pin-out, or both")

    # Amendment 4's pin is written on its own, so regenerating it cannot disturb
    # expected_identities.json, which Amendment 3 pins by sha256
    pin = None
    if args.bases_pin_out is not None:
        require(args.analysis_bases is not None, "--bases-pin-out needs --analysis-bases")
        pin = analysis_bases_pin(args.repo, args.analysis_bases)
        args.bases_pin_out.parent.mkdir(parents=True, exist_ok=True)
        args.bases_pin_out.write_text(json.dumps(pin, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"analysis_bases_pin": str(args.bases_pin_out),
                          "analysis_bases_pin_sha256": sha256_file(args.bases_pin_out),
                          "pins_bases_file_sha256": pin["file_sha256"]}, indent=2))
    if args.frozen_out is None:
        return pin, None

    for name in ("rebuilt", "registry", "coverage", "axis_map", "deposited", "legacy",
                 "manifest_out"):
        require(getattr(args, name) is not None,
                f"--frozen-out regenerates the whole identity and needs --{name.replace('_', '-')}")
    frozen = frozen_identities(args.repo, args.rebuilt, args.coverage, args.axis_map,
                               args.deposited, args.legacy, args.registry)
    args.frozen_out.parent.mkdir(parents=True, exist_ok=True)
    args.frozen_out.write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")

    manifest = implementation_manifest(args.repo, frozen, args.registry, args.coverage)
    manifest["frozen_identities"]["sha256"] = sha256_file(args.frozen_out)
    args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    pending = sorted(k for k, v in frozen.items() if v == PENDING) + sorted(
        f"rebuild/{k}" for k, v in frozen["rebuild_sha256"].items() if v == PENDING) + \
        sorted(f"manifest/{k}" for k, v in frozen["manifest_sha256"].items() if v == PENDING)
    print(json.dumps({"frozen_identities": str(args.frozen_out),
                      "frozen_identities_sha256": manifest["frozen_identities"]["sha256"],
                      "implementation_manifest": str(args.manifest_out),
                      "implementation_manifest_sha256": sha256_file(args.manifest_out),
                      "still_pending": pending}, indent=2))
    return frozen, manifest


if __name__ == "__main__":
    main()
