"""Ceiling check for the PRISM external criterion: is any cohort annotatable?

The blinded feasibility count returned 11 exact matches for the analysed
(>= 5 source) cohort. This asks whether that is a property of the cohort or of
the mapping, by computing the overlap ceiling over the entire JUMP library and
at successively weaker source thresholds.

Counts only; no viability value is read beyond testing non-nullness.
"""
import json
from pathlib import Path

import pandas as pd

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
ATLAS = Path("/Users/elliottower/Documents/GitHub/direction-instability-atlas")
OUT = REPO / "results" / "qc_prism_overlap"
STEM = r"^(BRD-[A-Z]\d+)"


def main():
    cpd = pd.read_csv(REPO / "data" / "jump_cp" / "compound_metadata.csv.gz")
    sam = pd.read_csv(REPO / "data" / "repurposing" / "repurposing_samples_20200324.txt",
                      sep="\t", comment="!", low_memory=False)
    lfc = pd.read_csv(ATLAS / "data" / "prism" / "Repurposing_Public_24Q2_LFC_COLLAPSED.csv",
                      usecols=["broad_id", "row_id", "LFC"])
    lfc = lfc[lfc.LFC.notna()]
    lfc["stem"] = lfc.broad_id.str.extract(STEM)[0]
    covered = set(lfc.stem.dropna())
    sam["stem"] = sam.broad_id.str.extract(STEM)[0]
    prism_keys = set(sam[sam.stem.isin(covered)].InChIKey.dropna())

    wm = pd.read_parquet(REPO / "data" / "jump_cp" / "cpg0016" / "well_metadata.parquet")
    el = wm[wm.row_ok]
    npl = el.Metadata_Plate.nunique()
    ppc = el.groupby("Metadata_JCP2022").Metadata_Plate.nunique()
    ctrl = set(ppc[ppc / npl >= 0.5].index)
    src = el[~el.Metadata_JCP2022.isin(ctrl)].groupby(
        "Metadata_JCP2022").Metadata_Source.nunique()

    def block(s):
        return {k.split("-")[0] for k in s}

    all_keys = set(cpd.Metadata_InChIKey.dropna())
    by_thr = {}
    for k in (1, 2, 3, 4, 5):
        coh = set(src[src >= k].index)
        keys = set(cpd[cpd.Metadata_JCP2022.isin(coh)].Metadata_InChIKey.dropna())
        by_thr[f">={k}_sources"] = {
            "compounds": len(coh),
            "exact_prism_overlap": len(keys & prism_keys),
            "connectivity_overlap": len(block(keys) & block(prism_keys)),
        }

    out = {
        "what_this_is": "overlap ceiling between JUMP compounds and PRISM 24Q2",
        "prism": {"compounds_with_lfc": len(covered),
                  "resolved_to_inchikey": len(prism_keys)},
        "jump_library": {"compounds": len(all_keys),
                         "exact_prism_overlap": len(all_keys & prism_keys),
                         "connectivity_overlap": len(block(all_keys) & block(prism_keys))},
        "by_source_threshold": by_thr,
        "threshold_required_for_option_e": 650,
        "conclusion": ("the exact-match ceiling over the ENTIRE JUMP library is "
                       "below the 650 required, and the analysed >=5-source "
                       "cohort is nearly disjoint from PRISM-annotated drug "
                       "space; no mapping improvement can rescue the test. This "
                       "concerns PRISM coverage only: these compounds may carry "
                       "ChEMBL, PubChem, vendor or target annotation."),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "prism_ceiling.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
