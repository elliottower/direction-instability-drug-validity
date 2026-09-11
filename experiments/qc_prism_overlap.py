"""Blinded feasibility count for the PRISM external-criterion test (Option E).

Answers one question: how many analysed JUMP-CP compounds carry a usable PRISM
viability annotation? It runs BEFORE the hypothesis is frozen, so it must not
reveal anything the hypothesis could be tuned to.

It reports counts only. It never emits a compound identity, a viability value or
any summary of one, a rank shift, or any association. PRISM `LFC` is used only
to test non-nullness and to count distinct cell lines: `n_lines` is computed,
`median(LFC)` is not. That boundary is what keeps this blind.

Mapping is exact-InChIKey and one-to-one in BOTH directions. The connectivity
fallback is counted separately and can never by itself flip the decision to
proceed, because block 1 of an InChIKey drops stereochemistry and protonation
and does not desalt a disconnected counterion.

    # one-time, see SAMPLES_URL below
    uv run --no-project --with pandas --with numpy --with pyarrow \
        python experiments/qc_prism_overlap.py
"""
import hashlib
import json
from pathlib import Path

import pandas as pd

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
ATLAS = Path("/Users/elliottower/Documents/GitHub/direction-instability-atlas")
CACHE = REPO / "data" / "jump_cp" / "cpg0016"
JUMP_CPD = REPO / "data" / "jump_cp" / "compound_metadata.csv.gz"
PRISM_LFC = ATLAS / "data" / "prism" / "Repurposing_Public_24Q2_LFC_COLLAPSED.csv"
# PRISM keys on broad_id and ships no structure key; the Repurposing Hub
# *samples* table is the public bridge. The *drugs* table already in the tree is
# keyed on pert_iname and cannot be used. This 2020 table predates PRISM 24Q2,
# so compounds added since are unmapped: that loss is counted, not assumed away.
SAMPLES = REPO / "data" / "repurposing" / "repurposing_samples_20200324.txt"
# the repo-hub path serves an HTML page, not the file; this LINCS mirror is
# the working source and carries broad_id, InChIKey and smiles
SAMPLES_URL = ("https://raw.githubusercontent.com/broadinstitute/"
               "lincs-cell-painting/master/metadata/moa/clue/"
               "repurposing_samples_20200324.txt")
OUT = REPO / "results" / "qc_prism_overlap"

CONTROL_PLATE_FRACTION = 0.50
MIN_SOURCES_GLOBAL = 5
MIN_PRISM_LINES = 100          # frozen: minimum cell lines per compound
MIN_OVERLAP_TO_PROCEED = 650   # frozen before running: below this, abandon E
BRD_STEM = r"^(BRD-[A-Z]\d+)"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def connectivity(key):
    """First InChIKey block.

    Drops the stereochemical and protonation layers. It is NOT salt stripping:
    a disconnected counterion changes the connectivity layer, so this is not a
    substitute for parent standardization.
    """
    return key.split("-")[0] if isinstance(key, str) else None


def one_to_one(df, a, b):
    """Rows of df whose values are unique in both directions between a and b."""
    pairs = df[[a, b]].dropna().drop_duplicates()
    ab = pairs.groupby(a)[b].nunique()
    ba = pairs.groupby(b)[a].nunique()
    return pairs[pairs[a].isin(ab[ab == 1].index) & pairs[b].isin(ba[ba == 1].index)]


def jump_cohort():
    wm = pd.read_parquet(CACHE / "well_metadata.parquet")
    assert "row_ok" in wm.columns, "well metadata lacks row_ok"
    elig = wm[wm.row_ok]
    n_plates = elig.Metadata_Plate.nunique()
    ppc = elig.groupby("Metadata_JCP2022").Metadata_Plate.nunique()
    controls = set(ppc[ppc / n_plates >= CONTROL_PLATE_FRACTION].index)
    src = elig[~elig.Metadata_JCP2022.isin(controls)] \
        .groupby("Metadata_JCP2022").Metadata_Source.nunique()
    return set(src[src >= MIN_SOURCES_GLOBAL].index), len(controls)


def main():
    for p in (JUMP_CPD, PRISM_LFC, CACHE / "well_metadata.parquet"):
        assert p.exists(), f"missing input: {p}"
    assert SAMPLES.exists(), (
        f"missing {SAMPLES}\nDownload once from {SAMPLES_URL}. It is the only "
        "public broad_id <-> InChIKey bridge; PRISM ships no structure key.")

    cohort, n_controls = jump_cohort()

    # ---- JUMP side: one JCP <-> one InChIKey, both directions
    cpd = pd.read_csv(JUMP_CPD)
    assert {"Metadata_JCP2022", "Metadata_InChIKey"} <= set(cpd.columns)
    jp = cpd[cpd.Metadata_JCP2022.isin(cohort)][
        ["Metadata_JCP2022", "Metadata_InChIKey"]].dropna().drop_duplicates()
    keys_per_jcp = jp.groupby("Metadata_JCP2022").Metadata_InChIKey.nunique()
    jcps_per_key = jp.groupby("Metadata_InChIKey").Metadata_JCP2022.nunique()
    n_jcp_multi_key = int((keys_per_jcp > 1).sum())
    n_key_multi_jcp = int((jcps_per_key > 1).sum())
    jump = one_to_one(jp, "Metadata_JCP2022", "Metadata_InChIKey")

    # ---- bridge side: one InChIKey <-> one BRD stem, both directions
    samples = pd.read_csv(SAMPLES, sep="\t", comment="!", low_memory=False)
    # the shipped header is "InChIKey"; normalize by stripping separators so a
    # future release spelling it "inchi_key" still resolves
    sc = {c.lower().replace("_", "").replace("-", ""): c for c in samples.columns}
    assert "broadid" in sc and "inchikey" in sc, \
        f"samples file lacks broad_id/InChIKey; has {list(samples.columns)}"
    bridge = samples[[sc["broadid"], sc["inchikey"]]].dropna()
    bridge.columns = ["broad_id", "inchi_key"]
    bridge["stem"] = bridge.broad_id.str.extract(BRD_STEM)[0]
    n_bridge_unparsed = int(bridge.stem.isna().sum())
    bridge = bridge.dropna(subset=["stem"])

    # ---- PRISM side: coverage only. LFC is counted, never summarized.
    lfc = pd.read_csv(PRISM_LFC, usecols=["broad_id", "row_id", "LFC"])
    lfc = lfc[lfc.LFC.notna()]
    lfc["stem"] = lfc.broad_id.str.extract(BRD_STEM)[0]
    n_lfc_unparsed = int(lfc.stem.isna().sum())
    lfc = lfc.dropna(subset=["stem"])
    lines = lfc.groupby("stem").row_id.nunique()
    covered = set(lines[lines >= MIN_PRISM_LINES].index)
    n_covered_absent_from_bridge = len(covered - set(bridge.stem))

    b11 = one_to_one(bridge[bridge.stem.isin(covered)], "inchi_key", "stem")

    # ---- exact join, the only route that can license proceeding
    exact = jump.merge(b11, left_on="Metadata_InChIKey", right_on="inchi_key",
                       how="inner")
    n_exact = int(exact.Metadata_JCP2022.nunique())

    # ---- connectivity fallback, reported as sensitivity only
    rest = jump[~jump.Metadata_JCP2022.isin(exact.Metadata_JCP2022)].copy()
    rest["conn"] = rest.Metadata_InChIKey.map(connectivity)
    bconn = b11.assign(conn=b11.inchi_key.map(connectivity))
    conn_ok = bconn.groupby("conn").stem.nunique()
    jconn = rest.groupby("conn").Metadata_JCP2022.nunique()
    n_conn = int(rest[rest.conn.isin(conn_ok[conn_ok == 1].index)
                      & rest.conn.isin(jconn[jconn == 1].index)].Metadata_JCP2022.nunique())
    n_conn_ambiguous = int(rest[rest.conn.isin(
        conn_ok[conn_ok > 1].index)].Metadata_JCP2022.nunique())

    out = {
        "what_this_is": ("blinded feasibility count; no viability value, rank "
                         "shift, compound identity or association is computed"),
        "frozen_before_running": {
            "min_prism_cell_lines": MIN_PRISM_LINES,
            "min_overlap_to_proceed": MIN_OVERLAP_TO_PROCEED,
            "control_plate_fraction": CONTROL_PLATE_FRACTION,
            "min_sources": MIN_SOURCES_GLOBAL,
            "primary_mapping": "exact InChIKey, one-to-one in both directions",
            "connectivity_role": "sensitivity only; cannot license proceeding",
        },
        "inputs": {
            "script_sha256": sha256(Path(__file__)),
            "prism_lfc_sha256": sha256(PRISM_LFC),
            "samples_sha256": sha256(SAMPLES),
            "jump_compound_metadata_sha256": sha256(JUMP_CPD),
        },
        "jump": {
            "controls_excluded": n_controls,
            "cohort_compounds": len(cohort),
            "with_inchikey": int(jp.Metadata_JCP2022.nunique()),
            "jcp_with_multiple_inchikeys_excluded": n_jcp_multi_key,
            "inchikeys_shared_by_multiple_jcp_excluded": n_key_multi_jcp,
            "one_to_one_compounds": int(jump.Metadata_JCP2022.nunique()),
        },
        "prism": {
            "compounds_with_lfc": int(lines.size),
            f"compounds_with_ge_{MIN_PRISM_LINES}_lines": len(covered),
            "covered_but_absent_from_bridge": n_covered_absent_from_bridge,
            "bridge_ids_unparsed": n_bridge_unparsed,
            "lfc_ids_unparsed": n_lfc_unparsed,
            "bridge_one_to_one_keys": int(b11.inchi_key.nunique()),
        },
        "mapping": {
            "exact_one_to_one": n_exact,
            "connectivity_only_unique_both_sides": n_conn,
            "connectivity_ambiguous_excluded": n_conn_ambiguous,
        },
        "decision": {
            "threshold": MIN_OVERLAP_TO_PROCEED,
            "exact_overlap": n_exact,
            "exact_proceed": bool(n_exact >= MIN_OVERLAP_TO_PROCEED),
            "expanded_overlap": n_exact + n_conn,
            "expanded_proceed_requires_parent_standardization": bool(
                n_exact < MIN_OVERLAP_TO_PROCEED
                and n_exact + n_conn >= MIN_OVERLAP_TO_PROCEED),
        },
    }
    blob = json.dumps(out)
    for token in ("JCP2022_", "BRD-"):
        assert token not in blob, f"identity leaked into output: {token}"
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "prism_overlap.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    d = out["decision"]
    if d["exact_proceed"]:
        print("\nPROCEED: exact one-to-one overlap clears the frozen threshold.")
    elif d["expanded_proceed_requires_parent_standardization"]:
        print("\nDO NOT PROCEED on this count. Exact mapping is short of the "
              "threshold and only the connectivity fallback would cross it. "
              "Attempt RDKit parent standardization from structures, or abandon "
              "Option E as confirmatory.")
    else:
        print("\nABANDON Option E as confirmatory: overlap is below threshold "
              "by every route.")


if __name__ == "__main__":
    main()
