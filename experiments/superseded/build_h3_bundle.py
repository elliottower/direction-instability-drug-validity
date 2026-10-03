"""Build the H3 cohort bundle from the original LINCS extraction.

The deposited artifact was produced by `03_phenotype_projection.py` from
`lincs_subset.npz` and `lincs_shrna.npz`, which are held on the Modal volume
`drug-perturbation-vol` from that run. Reconstructing those matrices from the
GEO GCTX reproduced their values exactly but not their coordinate order, and the
squared cosine between a drug mean and a target direction is only invariant to a
coordinate permutation applied to both. Using the original extraction removes
the ordering question rather than answering it, and is faithful to the artifact
by construction. Deviation 9 records this.

The aggregation below is the one in `03_phenotype_projection.py`, not the
loader's consensus function: per drug and cell line, the mean over every
signature; per target gene, the mean over every hairpin, at least three.

    uv run --no-project --with numpy --with pandas python \
        experiments/build_h3_bundle.py --data <dir holding the two npz files>
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
DPG = Path("/Users/elliottower/Documents/GitHub/drug-perturbation-geometry/data")
DEPOSITED = REPO / "results" / "03_phenotype_projection" / "phenotype_projection_results.json"
OUT = REPO / "results" / "03c_h3_sensitivity"
MIN_HAIRPINS = 3
MIN_CONTEXTS = 5
RECON_TOL = 1e-6


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(data_dir):
    data_dir = Path(data_dir)
    comp = np.load(data_dir / "lincs_subset.npz", allow_pickle=True)
    shr = np.load(data_dir / "lincs_shrna.npz", allow_pickle=True)
    assert [str(g) for g in comp["gene_ids"]] == [str(g) for g in shr["gene_ids"]], \
        "the two extractions do not share a gene axis"
    n_genes = comp["signatures"].shape[1]

    comp_sigs = comp["signatures"]
    comp_idx = {str(s): i for i, s in enumerate(comp["sig_ids"])}
    shr_sigs = shr["signatures"]
    shr_idx = {str(s): i for i, s in enumerate(shr["sig_ids"])}

    siginfo = pd.read_csv(DPG / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t", low_memory=False)
    filtered = siginfo[siginfo.sig_id.astype(str).isin(comp_idx) & siginfo.pert_iname.notna()].copy()
    filtered["_idx"] = filtered.sig_id.astype(str).map(comp_idx)

    drug_cell = {}
    for (drug, cell), grp in filtered.groupby(["pert_iname", "cell_id"]):
        drug_cell.setdefault(drug, {})[cell] = comp_sigs[grp._idx.values].mean(axis=0)

    shr_info = pd.read_csv(DPG / "lincs_shrna_siginfo.csv.gz")
    shr_f = shr_info[shr_info.sig_id.astype(str).isin(shr_idx)].copy()
    shr_f["_idx"] = shr_f.sig_id.astype(str).map(shr_idx)
    consensus, hairpins = {}, {}
    for gene, grp in shr_f.groupby("pert_iname"):
        if len(grp) < MIN_HAIRPINS:
            continue
        consensus[gene] = shr_sigs[grp._idx.values].mean(axis=0)
        hairpins[gene] = {"n": int(len(grp)), "sig_ids": sorted(grp.sig_id.astype(str))}

    deposited = json.loads(DEPOSITED.read_text())
    drugs, targets, mats, dirs = [], [], [], []
    for rec in deposited:
        d, t = rec["drug"], rec["target"]
        cells = drug_cell.get(d, {})
        if t not in consensus or len(cells) < MIN_CONTEXTS:
            continue
        u = consensus[t]
        drugs.append(d); targets.append(t)
        mats.append(np.array(list(cells.values()), dtype=np.float64))
        dirs.append(u / np.linalg.norm(u))

    worst = {"D": 0.0, "P": 0.0, "E": 0.0, "K": 0}
    dep = {r["drug"]: r for r in deposited}
    for d, S, u in zip(drugs, mats, dirs):
        r = dep[d]
        iu = np.triu_indices(S.shape[0], 1)
        unit = S / np.linalg.norm(S, axis=1, keepdims=True)
        D = 1.0 - float((unit @ unit.T)[iu].mean())
        P = float(np.abs((S[iu[0]] - S[iu[1]]) @ u).mean())
        m = S.mean(axis=0)
        E = float((m @ u / (np.linalg.norm(m) * np.linalg.norm(u))) ** 2)
        worst["D"] = max(worst["D"], abs(D - r["raw_instability"]))
        worst["P"] = max(worst["P"], abs(P - r["projected_instability"]))
        worst["E"] = max(worst["E"], abs(E - r["on_target_enrichment"]))
        worst["K"] = max(worst["K"], abs(S.shape[0] - r["n_celllines"]))

    print(f"rebuilt {len(drugs)} of {len(deposited)} deposited drugs, {n_genes} genes")
    for k in ("D", "P", "E"):
        print(f"  max |{k}_rebuilt - {k}_deposited| = {worst[k]:.3g}")
    print(f"  n_celllines mismatches = {worst['K']}")

    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "cohort_bundle.npz", drugs=np.array(drugs),
                        targets=np.array(targets),
                        signatures=np.array(mats, dtype=object), directions=np.array(dirs))
    manifest = {
        "source": "original extraction from Modal volume drug-perturbation-vol",
        "deviation": "Deviation 9",
        "stage_fingerprint": "", "fingerprint_parts": {},
        "n_bundled": len(drugs), "n_deposited": len(deposited),
        "reproduction_max_abs_error": {k: float(v) for k, v in worst.items()},
        "cohort_identifier_sha256": hashlib.sha256("\n".join(sorted(drugs)).encode()).hexdigest(),
        "bundle_sha256": sha256_file(OUT / "cohort_bundle.npz"),
        "lincs_subset_sha256": sha256_file(data_dir / "lincs_subset.npz"),
        "lincs_shrna_sha256": sha256_file(data_dir / "lincs_shrna.npz"),
        "deposited_sha256": sha256_file(DEPOSITED),
        "n_targets_with_consensus": len(consensus),
    }
    parts = {k: manifest[k] for k in ("lincs_subset_sha256", "lincs_shrna_sha256",
                                      "deposited_sha256", "n_targets_with_consensus")}
    manifest["fingerprint_parts"] = parts
    manifest["stage_fingerprint"] = hashlib.sha256(
        json.dumps(parts, sort_keys=True).encode()).hexdigest()
    (OUT / "rebuild_manifest.json").write_text(json.dumps(manifest, indent=2))
    (OUT / "shrna_hairpins.json").write_text(json.dumps(hairpins, indent=2, sort_keys=True))
    print(f"\nwritten to {OUT}")
    for k in ("D", "P", "E"):
        assert worst[k] < RECON_TOL, f"{k} reproduces to {worst[k]:.3g}, tolerance {RECON_TOL}"
    assert worst["K"] == 0
    print("reproduction within tolerance on every drug")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", required=True)
    main(ap.parse_args().data)
