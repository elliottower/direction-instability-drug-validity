"""Whether the shRNA eligibility boundary is the boundary the deposited analysis had.

The gate on the 812-drug cohort refused with "8 registered targets are absent from
the rebuild". That reads as a rebuild defect. It is a gate defect: the gate took the
compound cohort's 266 target symbols as the shRNA universe, and the shRNA-paired
universe is the 258 targets pinned in `registry/cohorts/shrna_eligible_targets.json`.

What makes that a coverage boundary rather than a rebuild that silently lost targets
is a conjunction, and every clause is asserted here against a pinned manifest rather
than inferred from whatever the files happen to hold:

  the rebuild's shRNA targets equal the pinned eligible set
  the deposited records carrying both shRNA quantities have exactly that target set
  the records carrying neither are exactly the pinned excluded records
  no record carries one shRNA quantity and not the other

A quantity is missing either because its key is absent or because the key holds null,
and the two are not the same fact about the file. Both count as missing, the mixed case
is refused, and which mechanism the file actually uses is recorded rather than assumed:
the 17 excluded records omit the keys and hold no explicit null.
  the deposited records and the compound cohort correspond one to one on (drug, target)
  the eligible and excluded manifests partition the compound cohort exactly

Cardinality is not identity: "258" is satisfied by a set with one target dropped and
another substituted, so every comparison is on a set or a canonical hash.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

SHRNA_FIELDS = ("proj_shrna", "enrich_shrna")


class CoverageError(AssertionError):
    """A refusal, raised so `python -O` cannot strip it."""


def require(condition, message) -> None:
    if not condition:
        raise CoverageError(message)


def listing(items) -> str:
    items = sorted(map(str, items))
    return f"{len(items)}: " + ", ".join(items)


def canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def sha256_file(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pair_key(record) -> tuple:
    return (str(record["drug"]), str(record["target"]))


def load_manifest(path: Path, expected_name: str) -> dict:
    manifest = json.loads(Path(path).read_text())
    require(manifest.get("name") == expected_name,
            f"{path} names {manifest.get('name')!r}, not {expected_name!r}")
    require(manifest.get("schema_version") == 1,
            f"{path} is schema version {manifest.get('schema_version')}, not 1")
    return manifest


def pairs_of(manifest: dict) -> set:
    records = manifest["records"]
    pairs = {pair_key(r) for r in records}
    require(len(pairs) == len(records),
            f"{manifest['name']}: a (drug, target) pair repeats in the manifest")
    require(sha256_text(canonical(sorted(pairs))) == manifest["pair_sha256"],
            f"{manifest['name']}: the records do not hash to the manifest's pair_sha256")
    return pairs


def measure(deposited_path: Path, consensus_path: Path, eligible_path: Path,
            excluded_path: Path, compound_path: Path) -> dict:
    deposited = json.loads(Path(deposited_path).read_text())
    require(isinstance(deposited, list) and deposited,
            f"{deposited_path} does not hold a list of records")
    for i, record in enumerate(deposited):
        require(isinstance(record, dict) and "drug" in record and "target" in record,
                f"deposited record {i} has no drug and target")
        for field in SHRNA_FIELDS:
            value = record.get(field)
            require(value is None or isinstance(value, (int, float)),
                    f"deposited record {i}: {field} is {type(value).__name__}")

    eligible = load_manifest(eligible_path, "shrna_eligible_targets")
    excluded = load_manifest(excluded_path, "shrna_excluded_records")
    compound = load_manifest(compound_path, "cohort_812_compound")

    eligible_targets = set(map(str, eligible["targets"]))
    require(sha256_text(canonical(sorted(eligible_targets))) == eligible["target_sha256"],
            "the eligible manifest's targets do not hash to its target_sha256")
    excluded_pairs = pairs_of(excluded)
    compound_pairs = pairs_of(compound)

    deposited_pairs = [pair_key(r) for r in deposited]
    require(len(set(deposited_pairs)) == len(deposited_pairs),
            "a (drug, target) pair repeats in the deposited records, so the record key "
            "is not unique and a one-to-one correspondence is not defined")
    deposited_pairs = set(deposited_pairs)

    # one to one on the record key, because equal counts are satisfied by a cohort that
    # drops one record and adds another
    require(deposited_pairs == compound_pairs,
            "the deposited records and the compound cohort do not correspond: "
            f"only deposited {listing(deposited_pairs - compound_pairs)}; "
            f"only cohort {listing(compound_pairs - deposited_pairs)}")

    # a record with one quantity and not the other is neither valued nor excluded, so it
    # would fall out of the partition silently
    def missing(record, field) -> bool:
        return record.get(field) is None

    half = [pair_key(r) for r in deposited
            if missing(r, "proj_shrna") != missing(r, "enrich_shrna")]
    require(not half,
            "a deposited record carries one shRNA quantity and not the other, so the gap "
            f"is per-field and not a coverage boundary: {listing(half)}")

    valued = [r for r in deposited if not missing(r, "proj_shrna")]
    unvalued = [r for r in deposited if missing(r, "proj_shrna")]
    valued_pairs = {pair_key(r) for r in valued}
    unvalued_pairs = {pair_key(r) for r in unvalued}
    valued_targets = {t for _, t in valued_pairs}
    unvalued_targets = {t for _, t in unvalued_pairs}

    consensus = np.load(consensus_path, allow_pickle=True)
    rebuild_targets = {str(gene) for gene in consensus["genes"].tolist()}
    require(len(rebuild_targets) == len(consensus["genes"]),
            "the rebuild consensus repeats a target identifier")

    result = {
        "inputs": {
            "deposited": str(deposited_path), "deposited_sha256": sha256_file(deposited_path),
            "consensus": str(consensus_path), "consensus_sha256": sha256_file(consensus_path),
            "eligible_manifest_sha256": sha256_file(eligible_path),
            "excluded_manifest_sha256": sha256_file(excluded_path),
            "compound_manifest_sha256": sha256_file(compound_path),
        },
        "counts": {"deposited_records": len(deposited), "valued": len(valued),
                   "unvalued": len(unvalued), "compound_cohort": len(compound_pairs),
                   "eligible_targets": len(eligible_targets),
                   "excluded_records": len(excluded_pairs),
                   "rebuild_shrna_targets": len(rebuild_targets),
                   "compound_cohort_targets": len({t for _, t in compound_pairs})},
        "set_equalities": {
            "rebuild_targets_equal_the_eligible_set": rebuild_targets == eligible_targets,
            "valued_targets_equal_the_eligible_set": valued_targets == eligible_targets,
            "unvalued_records_equal_the_excluded_manifest": unvalued_pairs == excluded_pairs,
            "deposited_corresponds_to_the_compound_cohort": deposited_pairs == compound_pairs,
            "no_record_half_valued": not half,
            "eligible_and_excluded_targets_disjoint": eligible_targets.isdisjoint(
                unvalued_targets),
            "partition_is_exact": valued_pairs | excluded_pairs == compound_pairs
                                  and not (valued_pairs & excluded_pairs),
        },
        "differences": {
            "rebuild_minus_eligible": sorted(rebuild_targets - eligible_targets),
            "eligible_minus_rebuild": sorted(eligible_targets - rebuild_targets),
            "valued_targets_symmetric_difference": sorted(valued_targets ^ eligible_targets),
            "unvalued_symmetric_difference": sorted(
                f"{d}|{t}" for d, t in unvalued_pairs ^ excluded_pairs),
        },
        "excluded_targets": sorted(unvalued_targets),
        "how_the_quantities_are_missing": {
            field: {"key_absent": sum(1 for r in deposited if field not in r),
                    "key_present_and_null": sum(1 for r in deposited
                                                if field in r and r[field] is None)}
            for field in SHRNA_FIELDS},
        "record_shapes": sorted(
            ({"keys": list(keys), "n_records": n} for keys, n in
             {tuple(sorted(r)): sum(1 for s in deposited if tuple(sorted(s)) == tuple(sorted(r)))
              for r in deposited}.items()),
            key=lambda entry: -entry["n_records"]),
        "hashes": {
            "rebuild_target_sha256": sha256_text(canonical(sorted(rebuild_targets))),
            "eligible_target_sha256": eligible["target_sha256"],
            "unvalued_pair_sha256": sha256_text(canonical(sorted(unvalued_pairs))),
            "excluded_manifest_pair_sha256": excluded["pair_sha256"],
        },
    }
    result["holds"] = all(result["set_equalities"].values())
    result["reading"] = (
        f"the shRNA-eligible universe is the {len(eligible_targets)} targets and "
        f"{len(valued)} drug-target records pinned in the manifests; the complement in "
        f"the compound cohort is exactly the {len(excluded_pairs)} pinned records, whose "
        f"targets are {', '.join(sorted(unvalued_targets))}; the boundary is coverage in "
        "the source, carried identically by the deposited analysis, and not an artifact "
        "of the rebuild"
        if result["holds"] else
        "the boundary does not hold as registered; see set_equalities and differences")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--deposited", type=Path, required=True)
    parser.add_argument("--consensus", type=Path, required=True)
    parser.add_argument("--eligible", type=Path, required=True)
    parser.add_argument("--excluded", type=Path, required=True)
    parser.add_argument("--compound-cohort", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    result = measure(args.deposited, args.consensus, args.eligible, args.excluded,
                     args.compound_cohort)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")   # written before any refusal

    for name, held in result["set_equalities"].items():
        require(held, f"{name} is false; the shRNA eligibility boundary is not as pinned")
    print(json.dumps({k: v for k, v in result.items() if k != "inputs"}, indent=2))
    return result


if __name__ == "__main__":
    main()
