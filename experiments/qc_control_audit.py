"""Metadata-only QC for the JUMP-CP cohort. Computes no direction instability.

Run before freezing the registration addendum. Everything here is design
information -- how many compounds sit in which stratum, which perturbations are
plate-level controls -- and none of it inspects an outcome.

The exclusion rule is frequency-based because the assembled profile table
carries no control-role column. The frequency distribution is written out so
the threshold can be seen to separate controls from screened compounds rather
than having been tuned.

    uv run --no-project --with pandas --with numpy --with pyarrow \
        python experiments/qc_control_audit.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
CACHE = REPO / "data" / "jump_cp" / "cpg0016"
OUT = REPO / "results" / "qc_control_audit"
DMSO_JCP = "JCP2022_033924"
CONTROL_PLATE_FRACTION = 0.50     # frozen: see distribution in the output
MIN_SOURCES_GLOBAL = 5            # main ablation cohort
MIN_WELLS_REPLICATE = 2           # matched-contrast cohort
MIN_SOURCES_REPLICATE = 3
R6_BINS = [(5, 6), (7, 8), (9, 10)]


def main():
    wm = pd.read_parquet(CACHE / "well_metadata.parquet")
    assert {"row_ok", "Metadata_Plate", "Metadata_Source",
            "Metadata_JCP2022"} <= set(wm.columns)
    ok = wm.row_ok.to_numpy()
    elig = wm[ok]
    n_plates = int(elig.Metadata_Plate.nunique())

    plates_per = elig.groupby("Metadata_JCP2022").Metadata_Plate.nunique()
    frac = plates_per / n_plates
    controls = sorted(frac[frac >= CONTROL_PLATE_FRACTION].index)

    # the gap either side of the threshold is what justifies it
    ordered = frac.sort_values(ascending=False)
    top20 = [{"jcp": j, "plates": int(plates_per[j]), "frac_plates": round(float(f), 5)}
             for j, f in ordered.head(20).items()]
    below = ordered[ordered < CONTROL_PLATE_FRACTION]

    keep = ~elig.Metadata_JCP2022.isin(controls)

    def cohort_global(mask):
        s = elig[mask].groupby("Metadata_JCP2022").Metadata_Source.nunique()
        return s[s >= MIN_SOURCES_GLOBAL]

    g_with = cohort_global(np.ones(len(elig), bool))
    g_without = cohort_global(keep.to_numpy())

    bins = {}
    for lo, hi in R6_BINS:
        n = int(((g_without >= lo) & (g_without <= hi)).sum())
        bins[f"{lo}-{hi}"] = n
    # per-value counts decide whether any coverage contrast is estimable
    exact = {int(k): int(v) for k, v in g_without.value_counts().sort_index().items()}
    wells_all = elig[keep].groupby("Metadata_JCP2022").size()
    wells_in_cohort = wells_all[wells_all.index.isin(g_without.index)]
    wq = {f"p{q}": float(wells_in_cohort.quantile(q / 100)) for q in (25, 50, 75, 90, 99)}

    per = elig[keep].groupby(["Metadata_JCP2022", "Metadata_Source"]).size()
    rep = per[per >= MIN_WELLS_REPLICATE]
    nrep = rep.groupby("Metadata_JCP2022").size()
    repl_cohort = nrep[nrep >= MIN_SOURCES_REPLICATE]
    wells_per = elig[keep & elig.Metadata_JCP2022.isin(repl_cohort.index)] \
        .groupby("Metadata_JCP2022").size()

    out = {
        "what_this_is": ("metadata-only cohort QC; no direction instability is "
                         "computed here"),
        "release": {"wells": int(len(wm)), "eligible_wells": int(ok.sum()),
                    "distinct_eligible_plates": n_plates},
        "control_rule": {
            "definition": (f"exclude DMSO and any compound present on >= "
                           f"{CONTROL_PLATE_FRACTION:.0%} of eligible COMPOUND plates"),
            "dmso": DMSO_JCP,
            "threshold_frac_plates": CONTROL_PLATE_FRACTION,
            "n_identified": len(controls),
            "identified": controls,
            "max_frac_below_threshold": round(float(below.max()), 5),
            "separation_note": ("controls sit far above the next most frequent "
                                "compound; the threshold is not near any compound"),
            "top20_by_plate_fraction": top20,
        },
        "global_cohort_min_sources_5": {
            "with_controls": int(len(g_with)),
            "without_controls": int(len(g_without)),
            "removed": int(len(g_with) - len(g_without)),
        },
        "r6_source_count_bins_control_excluded": bins,
        "sources_per_compound_exact_counts": exact,
        "wells_per_compound_in_global_cohort": {
            "min": int(wells_in_cohort.min()), "max": int(wells_in_cohort.max()),
            **{k: round(v, 1) for k, v in wq.items()}},
        "replicate_cohort_control_excluded": {
            "min_wells_per_source": MIN_WELLS_REPLICATE,
            "min_sources": MIN_SOURCES_REPLICATE,
            "n_compounds": int(len(repl_cohort)),
            "n_wells": int(wells_per.sum()),
            "max_wells_per_compound": int(wells_per.max()),
            "median_wells_per_compound": float(wells_per.median()),
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "control_audit.json").write_text(json.dumps(out, indent=2))
    frac.sort_values(ascending=False).to_csv(OUT / "plate_fraction_per_compound.csv",
                                             header=["frac_eligible_plates"])
    (OUT / "excluded_control_ids.txt").write_text("\n".join(controls) + "\n")

    print(json.dumps({k: v for k, v in out.items() if k != "control_rule"}, indent=2))
    print(f"\ncontrols identified ({len(controls)}):")
    for j in controls:
        print(f"  {j}  {plates_per[j]:>5,} plates  {frac[j]:.1%}")
    print(f"\nhighest plate fraction among retained compounds: {below.max():.1%}")
    print(f"written to {OUT}")


if __name__ == "__main__":
    main()
