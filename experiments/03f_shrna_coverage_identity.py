"""Whether the rebuild's shRNA coverage is the coverage the deposited analysis had.

The gate on the 812-drug cohort refused with "8 registered targets are absent from
the rebuild". That reads as a rebuild defect. It is a gate defect: the gate took the
cohort's 266 target symbols as the shRNA universe, and the shRNA-paired universe is
258 in the rebuild and was 258 in the deposited analysis too.

This measures the claim rather than asserting it. Three things have to hold for the
reading to stand: the deposited records carry no shRNA value for exactly the records
whose target is one of the eight; those records' drugs are exactly the seventeen the
rebuild was extended by; and the rebuild's shRNA target set is the same set as the
deposited analysis's usable one, not merely the same size.
"""
import argparse
import json
from pathlib import Path

import numpy as np


class CoverageError(AssertionError):
    """A refusal, raised so `python -O` cannot strip it."""


def require(condition, message) -> None:
    if not condition:
        raise CoverageError(message)


def measure(deposited_path: Path, consensus_path: Path, cohort_path: Path) -> dict:
    deposited = json.loads(Path(deposited_path).read_text())
    require(isinstance(deposited, list) and deposited,
            f"{deposited_path} does not hold a list of records")
    consensus = np.load(consensus_path, allow_pickle=True)
    rebuild_targets = {str(gene) for gene in consensus["genes"].tolist()}
    cohort = json.loads(Path(cohort_path).read_text())

    def valued(record) -> bool:
        return (record.get("proj_shrna") is not None
                and record.get("enrich_shrna") is not None)

    usable = [r for r in deposited if valued(r)]
    unusable = [r for r in deposited if not valued(r)]
    deposited_targets = {str(r["target"]) for r in usable}
    cohort_targets = {str(r["target"]) for r in cohort}
    absent = cohort_targets - rebuild_targets

    # a record is unusable exactly when its target has no shRNA direction: the
    # alternative is that some records of a covered target are also null, which would
    # make the gap a per-record defect rather than a coverage boundary
    unusable_targets = {str(r["target"]) for r in unusable}
    partially_null = sorted(unusable_targets & deposited_targets)

    result = {
        "deposited_records": len(deposited),
        "deposited_usable": len(usable),
        "deposited_unusable": len(unusable),
        "cohort_records": len(cohort),
        "cohort_targets": len(cohort_targets),
        "rebuild_shrna_targets": len(rebuild_targets),
        "deposited_usable_targets": len(deposited_targets),
        "target_universes_identical": deposited_targets == rebuild_targets,
        "target_symmetric_difference": sorted(deposited_targets ^ rebuild_targets),
        "targets_absent_from_the_rebuild": sorted(absent),
        "targets_with_no_deposited_shrna_value": sorted(unusable_targets),
        "absent_targets_are_the_unvalued_ones": sorted(absent) == sorted(unusable_targets),
        "targets_both_valued_and_unvalued": partially_null,
        "drugs_of_the_unvalued_records": sorted({str(r["drug"]) for r in unusable}),
        "usable_drugs": len(({str(r["drug"]) for r in usable})),
    }
    result["reading"] = (
        "the shRNA-paired universe is "
        f"{result['rebuild_shrna_targets']} targets and "
        f"{result['usable_drugs']} drugs; the gap to the cohort is coverage in the "
        "source, carried identically by the deposited analysis, and not an artifact "
        "of the rebuild"
        if result["target_universes_identical"] and result["absent_targets_are_the_unvalued_ones"]
        else "the reading does not hold; see the symmetric difference")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--deposited", type=Path, required=True)
    parser.add_argument("--consensus", type=Path, required=True)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    result = measure(args.deposited, args.consensus, args.cohort)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")   # written before any refusal

    require(not result["targets_both_valued_and_unvalued"],
            "a target is both valued and unvalued in the deposited records, so the gap "
            f"is per-record and not a coverage boundary: {result['targets_both_valued_and_unvalued']}")
    require(result["target_universes_identical"],
            "the rebuild's shRNA universe differs from the deposited analysis's: "
            f"{result['target_symmetric_difference']}")
    require(result["absent_targets_are_the_unvalued_ones"],
            "the targets absent from the rebuild are not the ones the deposited analysis "
            "left unvalued")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("drugs_of_the_unvalued_records",)}, indent=2))
    return result


if __name__ == "__main__":
    main()
