"""Minimal cohort manifests: stable identifiers only, each with its own identity fields.

The registered cohorts were results files. `phenotype_projection_results.json` carries
projected instability and on-target enrichment alongside the drug and target, so
eligibility was coupled to derived values, and a gate pinning it pinned numbers it had
no business pinning. These manifests carry `(drug, target)` and nothing else.

Four objects, because the gate needs all four distinguished:

  cohort_795          the shRNA-paired cohort
  cohort_812          the compound cohort, which is the 795 plus the extension
  shrna_eligible      the exact target-identifier set shRNA quantities exist for
  shrna_excluded      the exact drug-target records the 812 cohort holds and 795 does not

Cardinality is not identity. "258 targets" is satisfied by a set with one target
dropped and another substituted, so every manifest carries a canonical hash of its
sorted contents and the gate compares sets, never counts.

Written canonically — sorted, fixed separators, newline-terminated — so the hash of a
manifest is a function of its contents and not of how it was serialized.
"""
import argparse
import hashlib
import json
from pathlib import Path

SCHEMA_VERSION = 1


class ManifestError(AssertionError):
    """A refusal, raised so `python -O` cannot strip it."""


def require(condition, message) -> None:
    if not condition:
        raise ManifestError(message)


def listing(items) -> str:
    items = sorted(map(str, items))
    return f"{len(items)}: " + ", ".join(items)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def pair_key(record) -> tuple:
    return (str(record["drug"]), str(record["target"]))


def identity(records, name: str) -> dict:
    """Every field the gate compares a cohort on, so no comparison rests on a count."""
    pairs = sorted(pair_key(r) for r in records)
    require(len(set(pairs)) == len(pairs),
            f"{name}: a drug-target pair repeats, so the record key is not unique")
    drugs = sorted({d for d, _ in pairs})
    targets = sorted({t for _, t in pairs})
    require(len(drugs) == len(pairs),
            f"{name}: {len(pairs)} records hold {len(drugs)} distinct drugs; the unit is the drug")
    return {
        "schema_version": SCHEMA_VERSION,
        "name": name,
        "n_records": len(pairs),
        "n_unique_drugs": len(drugs),
        "n_unique_targets": len(targets),
        "pair_sha256": sha256_text(canonical(pairs)),
        "drug_sha256": sha256_text(canonical(drugs)),
        "target_sha256": sha256_text(canonical(targets)),
        "uniqueness": "the drug is unique across records; the (drug, target) pair is unique",
        "records": [{"drug": d, "target": t} for d, t in pairs],
    }


def target_set_identity(targets, name: str) -> dict:
    ordered = sorted(map(str, targets))
    require(len(set(ordered)) == len(ordered),
            f"{name}: a target repeats, so the set is not well defined")
    return {"schema_version": SCHEMA_VERSION, "name": name, "n_targets": len(ordered),
            "target_sha256": sha256_text(canonical(ordered)), "targets": ordered}


def write(path: Path, payload: dict) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(paired_path: Path, compound_path: Path, out_dir: Path) -> dict:
    paired = json.loads(Path(paired_path).read_text())
    compound = json.loads(Path(compound_path).read_text())
    require(isinstance(paired, list) and isinstance(compound, list),
            "both cohort sources must be lists of records")

    paired_pairs = {pair_key(r) for r in paired}
    compound_pairs = {pair_key(r) for r in compound}
    require(paired_pairs <= compound_pairs,
            "the shRNA-paired cohort is not a subset of the compound cohort: "
            f"{listing(paired_pairs - compound_pairs)}")
    excluded = sorted(compound_pairs - paired_pairs)

    m795 = identity(paired, "cohort_795_shrna_paired")
    m812 = identity(compound, "cohort_812_compound")
    eligible = target_set_identity({t for _, t in paired_pairs}, "shrna_eligible_targets")
    excluded_m = identity([{"drug": d, "target": t} for d, t in excluded],
                          "shrna_excluded_records")
    excluded_m["targets"] = sorted({t for _, t in excluded})

    # the partition is asserted, because 812 = 795 + 17 by count is also satisfied by a
    # cohort that drops one record and adds another
    require(m795["n_records"] + excluded_m["n_records"] == m812["n_records"],
            f"{m795['n_records']} + {excluded_m['n_records']} does not equal "
            f"{m812['n_records']}")
    require(set(eligible["targets"]).isdisjoint(excluded_m["targets"]),
            "a target is both shRNA-eligible and excluded: "
            f"{listing(set(eligible['targets']) & set(excluded_m['targets']))}")

    m812["extension_over_795"] = {
        "n_records": excluded_m["n_records"],
        "pair_sha256": excluded_m["pair_sha256"],
        "manifest": "shrna_excluded_records.json",
    }

    written = {}
    for stem, payload in (("cohort_795_shrna_paired", m795),
                          ("cohort_812_compound", m812),
                          ("shrna_eligible_targets", eligible),
                          ("shrna_excluded_records", excluded_m)):
        written[f"{stem}.json"] = write(out_dir / f"{stem}.json", payload)

    return {
        "built_from": {"shrna_paired_source": str(paired_path),
                       "compound_source": str(compound_path),
                       "shrna_paired_source_sha256": hashlib.sha256(
                           Path(paired_path).read_bytes()).hexdigest(),
                       "compound_source_sha256": hashlib.sha256(
                           Path(compound_path).read_bytes()).hexdigest()},
        "manifest_sha256": written,
        "identity": {"cohort_795": {k: v for k, v in m795.items() if k != "records"},
                     "cohort_812": {k: v for k, v in m812.items() if k != "records"},
                     "shrna_eligible_targets": {k: v for k, v in eligible.items()
                                                if k != "targets"},
                     "shrna_excluded_records": {k: v for k, v in excluded_m.items()
                                                if k != "records"}},
        "excluded_records": [{"drug": d, "target": t} for d, t in excluded],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--shrna-paired", type=Path, required=True)
    parser.add_argument("--compound", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    args = parser.parse_args(argv)

    summary = build(args.shrna_paired, args.compound, args.out_dir)
    args.index.parent.mkdir(parents=True, exist_ok=True)
    args.index.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


if __name__ == "__main__":
    main()
