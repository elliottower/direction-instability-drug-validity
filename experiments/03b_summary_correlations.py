"""Summary correlations of the corrected CRISPRi convergent-validity records.

A deterministic extraction from `results/03b_h3_crispri/h3_crispri_results.json`,
not a new analysis: every value here was printed by the corrected run of
`03b_h3_crispri_ground_truth.py` (commit f822fb1). The raw-instability
correlations on the matched cohort are deliberately not computed; they are
registered in `experiments/PREREG_H3_REFERENCE_DISCORDANCE.md`.

    uv run --no-project --with scipy python experiments/03b_summary_correlations.py
"""
import hashlib
import json
from pathlib import Path

import scipy
from scipy.stats import spearmanr

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
RECORDS = REPO / "results" / "03b_h3_crispri" / "h3_crispri_results.json"
OUT = REPO / "results" / "03b_h3_crispri" / "summary_correlations.json"
EXPECTED_RECORDS_SHA256 = "f00f8428071178bd2317139c5a1f534e458931684fb6da65c750f9f1e5456eb3"


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rho(records, x, y):
    r = spearmanr([rec[x] for rec in records], [rec[y] for rec in records])
    return [float(r.statistic), float(r.pvalue)]


def main():
    assert sha256_file(RECORDS) == EXPECTED_RECORDS_SHA256, "records are not the corrected artifact"
    records = json.loads(RECORDS.read_text())
    shrna = [r for r in records if "proj_shrna" in r]
    crispri = [r for r in records if "proj_crispri" in r]
    matched = [r for r in crispri if "proj_shrna" in r]

    out = {
        "note": ("deterministic extraction from the corrected per-drug records; every value "
                 "was printed by the corrected run. Matched raw correlations are not computed."),
        "shrna_all": {"n": len(shrna),
                      "projected_vs_enrichment": rho(shrna, "proj_shrna", "enrich_shrna"),
                      "raw_vs_enrichment": rho(shrna, "raw_instability", "enrich_shrna")},
        "crispri_all": {"n": len(crispri),
                        "projected_vs_enrichment": rho(crispri, "proj_crispri", "enrich_crispri"),
                        "raw_vs_enrichment": rho(crispri, "raw_instability", "enrich_crispri")},
        "matched": {"n": len(matched),
                    "shrna_projected": rho(matched, "proj_shrna", "enrich_shrna"),
                    "crispri_projected": rho(matched, "proj_crispri", "enrich_crispri")},
        "crispri_targets": {"n_distinct": len({r["target"] for r in crispri})},
        "provenance": {"records": str(RECORDS.relative_to(REPO)),
                       "records_sha256": EXPECTED_RECORDS_SHA256,
                       "script_sha256": sha256_file(__file__),
                       "scipy": scipy.__version__,
                       "statistic": "scipy.stats.spearmanr, two-sided p"},
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
