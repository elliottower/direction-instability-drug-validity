"""R0-R7 of the reference-discordance registration frozen at 7f57136.

Why H3's pattern reverses between the shRNA and CRISPRi references. Everything
here is fixed by `experiments/PREREG_H3_REFERENCE_DISCORDANCE.md`: the cohorts,
the constructions, the inference machinery, every threshold, and which outcome
counts as which reading.

The single-cell modules R0.3, R0.4 and R0.6 run on Modal
(`experiments/modal_03d_single_cell.py`) and arrive here as a JSON file.

    PYTHONPATH=. uv run --no-project --with numpy --with scipy --with pandas \\
        --with anndata python experiments/03d_h3_reference_discordance.py \\
        --data ../drug-perturbation-geometry/data \\
        --perturbseq data/scperturb/ReplogleWeissman2022_K562_essential.h5ad \\
        --replogle data/replogle2022 \\
        [--prism ../direction-instability-atlas/data/prism_19q4] \\
        [--single-cell results/03d_h3_reference_discordance/single_cell_audit.json]
"""
import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from geometry.direction_instability import direction_instability
from geometry.inference import (cluster_bootstrap, comparison_reading, excludes_zero,
                                percentile_table,
                                is_practically_null, own_target_percentile, percentile_interval,
                                permutation_reading, rank_partial_correlation,
                                target_balanced_mean, target_balanced_spearman,
                                two_sided_permutation_p, unique_target_permutations)
from geometry.single_cell import split_half_reliability_by_unit
from geometry.references import (MIN_LANDMARKS, N_LANDMARK, Reference, alignment_matrix,
                                 landmark_symbols, load_replogle_bulk, pooled_crispri_reference,
                                 projected_dispersion, shared_space, unit)

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
CRISPRI_RECORDS = REPO / "results" / "03b_h3_crispri" / "h3_crispri_results.json"
OUT = REPO / "results" / "03d_h3_reference_discordance"

# pinned in the registration
EXPECTED_RECORDS_SHA256 = "f00f8428071178bd2317139c5a1f534e458931684fb6da65c750f9f1e5456eb3"
EXPECTED_LINCS_SUBSET_SHA256 = "2ad0f5d30ab826f9ec0cfe37f6b53b2829bfef1d73b7920c3441de6407adb4ec"
EXPECTED_LINCS_SHRNA_SHA256 = "4a990e5072a43f59fdda60a2ff040f355f06be947c14bcfa87056c66332b58e7"
EXPECTED_PERTURBSEQ_SHA256 = "412fd0df8c4ccea9f4db91cd88033c49200838b29d40945e48574be588b48789"
EXPECTED_BULK_SHA256 = {
    "K562_essential_normalized_bulk_01.h5ad":
        "c1ca6456c9c9f1aa2b02c496eb64d1dc3e6a852edbd744d682b8d2c95fd36829",
    "rpe1_normalized_bulk_01.h5ad":
        "a3c5bfd0f15d63938bc80c9b8874b9cd761e3a23caf5ffe7966bae4e887ec89d",
    "K562_gwps_normalized_bulk_01.h5ad":
        "37e48c474d8b5dead4151f96ea8f5fe7bbe6beb10eeea48685b740c3f74490a2",
}

MIN_CELL_LINES = 5
MIN_SIGNATURES = 3        # distinct signature ids, as the corrected artifact counts them
RECON_TOL = 1e-6
N_BOOT = 10_000
N_PERM = 10_000
SEED_BOOT = 20260922
SEED_PERM = 20260923
SEED_SPLIT = 20260924   # split-half reliability, as on the single-cell side
ENERGY_P = 0.05                  # R4's phenotype-positive rule
R4_MIN_DRUGS, R4_MIN_TARGETS = 60, 30
R5_MIN_DRUGS, R5_MIN_TARGETS = 60, 20
R5_MIN_LINES = 100
R5_RATIO_DENOMINATOR = 0.15
R7C_MIN_RELIABILITY = 0.5
R7C_MIN_CRISPRI_TARGETS, R7C_MIN_SHRNA_TARGETS = 20, 100
Q_EQUIVALENCE = (0.40, 0.60)
HEMATOPOIETIC = ["HL60", "JURKAT", "NOMO1", "PL21", "SKM1", "THP1", "U266", "U937", "WSUDLCL2"]
MYELOID = ["HL60", "THP1", "U937", "NOMO1", "PL21", "SKM1"]
BRD_STEM = re.compile(r"(BRD-[A-Z]\d{8})")


class Draws(dict):
    """Every resampling distribution the run produces, keyed once.

    The registrations require the replicate arrays to travel with the result, so
    the intervals and p-values can be recomputed without rerunning. A key is
    written once: a silent overwrite would lose a distribution.
    """

    def keep(self, key, array):
        assert key not in self, f"two distributions claim the key {key}"
        self[key] = np.asarray(array)
        return array


def log(message):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


# --------------------------------------------------------------------------- cohorts


@dataclass
class Arm:
    """One set of drugs, their per-cell-line signatures, and their targets."""

    name: str
    drugs: list
    targets: np.ndarray
    signatures: list           # (K_c, 978) per drug
    cell_lines: list           # cell-line identifiers per drug, same order as rows

    def __len__(self):
        return len(self.drugs)

    def subset(self, keep) -> "Arm":
        keep = np.asarray(keep)
        return Arm(self.name, [self.drugs[i] for i in keep], self.targets[keep],
                   [self.signatures[i] for i in keep], [self.cell_lines[i] for i in keep])

    def means(self, positions=None) -> np.ndarray:
        return np.vstack([S.mean(axis=0) if positions is None else S[:, positions].mean(axis=0)
                          for S in self.signatures])


def build_drug_signatures(data_dir: Path):
    """Per drug and cell line, the mean over every signature, as the artifact does."""
    compounds = np.load(data_dir / "lincs_subset.npz", allow_pickle=True)
    signatures = compounds["signatures"]
    position = {str(sig_id): i for i, sig_id in enumerate(compounds["sig_ids"])}
    assert len(position) == len(compounds["sig_ids"]), (
        "duplicate signature ids in the compound extraction")
    siginfo = pd.read_csv(data_dir / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(position) & siginfo.pert_iname.notna()].copy()
    siginfo["_row"] = siginfo.sig_id.astype(str).map(position)

    per_drug = {}
    for (drug, cell), group in siginfo.groupby(["pert_iname", "cell_id"]):
        per_drug.setdefault(drug, {})[cell] = signatures[group._row.values].mean(axis=0)
    gene_ids = [str(g) for g in compounds["gene_ids"]]
    return per_drug, gene_ids


def build_shrna_reference(data_dir: Path, name="shRNA") -> Reference:
    """Consensus shRNA direction per target, indexed by signature id.

    The defect of Deviation 9 was to index the signature matrix by the metadata's
    row order. The map below is built from the matrix's own identifiers.
    """
    shrna = np.load(data_dir / "lincs_shrna.npz", allow_pickle=True)
    signatures = shrna["signatures"]
    position = {str(sig_id): i for i, sig_id in enumerate(shrna["sig_ids"])}
    assert len(position) == len(shrna["sig_ids"]), "duplicate signature ids in the shRNA extraction"
    siginfo = pd.read_csv(data_dir / "lincs_shrna_siginfo.csv.gz")
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(position)].copy()
    siginfo["_row"] = siginfo.sig_id.astype(str).map(position)

    directions, signature_counts = {}, {}
    for gene, group in siginfo.groupby("pert_iname"):
        rows = sorted(set(group.sig_id.astype(str)))
        if len(rows) < MIN_SIGNATURES:
            continue
        consensus = signatures[[position[r] for r in rows]].mean(axis=0)
        if np.linalg.norm(consensus) == 0:
            continue
        directions[gene] = unit(consensus)
        signature_counts[gene] = len(rows)
    # a signature id identifies one experiment, not one reagent: a hairpin measured
    # in several cell lines contributes several signatures. The corrected artifact
    # counts distinct signatures, so the eligibility rule is kept and the field is
    # named for what it counts.
    audit = {"n_targets": len(directions), "min_distinct_signatures": MIN_SIGNATURES,
             "n_distinct_signatures": signature_counts}
    return Reference(name, directions, np.arange(N_LANDMARK), audit)


def drug_targets(data_dir: Path):
    """First listed target per drug, and whether the annotation lists more."""
    labels = json.loads((data_dir / "frozen_drug_labels.json").read_text())
    first, multiple, moa = {}, {}, {}
    for entry in labels["drugs"]:
        target = entry.get("target")
        if not target:
            continue
        first[entry["pert_iname"]] = target.split("|")[0].strip()
        multiple[entry["pert_iname"]] = "|" in target
        moa[entry["pert_iname"]] = str(entry.get("moa", ""))
    return first, multiple, moa


def assemble_arm(name, per_drug, targets, reference: Reference, drugs=None) -> Arm:
    """Drugs with a direction under `reference` and enough cell lines."""
    chosen, chosen_targets, signatures, cells = [], [], [], []
    for drug, per_cell in sorted(per_drug.items()):
        if drugs is not None and drug not in drugs:
            continue
        target = targets.get(drug)
        if target is None or target not in reference.directions or len(per_cell) < MIN_CELL_LINES:
            continue
        chosen.append(drug)
        chosen_targets.append(target)
        signatures.append(np.array([per_cell[c] for c in sorted(per_cell)], dtype=np.float64))
        cells.append(sorted(per_cell))
    return Arm(name, chosen, np.array(chosen_targets), signatures, cells)


# ------------------------------------------------------------------- per-drug quantities


def quantities(arm: Arm, reference: Reference) -> dict:
    """D, P and E for every drug of an arm under one reference, on its gene space."""
    positions = reference.positions
    whole_space = len(positions) == N_LANDMARK and positions[0] == 0 and positions[-1] == N_LANDMARK - 1
    restricted = arm.signatures if whole_space else [S[:, positions] for S in arm.signatures]
    P, E, D_h = np.empty(len(arm)), np.empty(len(arm)), np.empty(len(arm))
    for i, (S, target) in enumerate(zip(restricted, arm.targets)):
        u = reference.directions[target]
        upper = np.triu_indices(S.shape[0], k=1)
        differences = S[upper[0]] - S[upper[1]]
        P[i] = projected_dispersion(differences, u)
        mean = S.mean(axis=0)
        E[i] = float((mean @ u / (np.linalg.norm(mean) * np.linalg.norm(u))) ** 2)
        D_h[i] = direction_instability(S)
    D = np.array([direction_instability(S) for S in arm.signatures])
    return {"P": P, "E": E, "D": D, "D_h": D_h,
            "K": np.array([S.shape[0] for S in arm.signatures], dtype=float)}


def own_target_q(arm: Arm, reference: Reference):
    """Each drug's own-target percentile over the arm's distinct targets."""
    genes = sorted(set(arm.targets))
    column = {gene: i for i, gene in enumerate(genes)}
    matrix = alignment_matrix(arm.means(reference.positions), reference, genes)
    own = np.array([column[t] for t in arm.targets])
    return own_target_percentile(matrix, own), matrix, genes


# ------------------------------------------------------------------------- readings


def association(x, y, targets, label):
    """A Spearman correlation with its target-cluster interval and readings."""
    rho = float(spearmanr(x, y).statistic)
    draws = cluster_bootstrap(targets, lambda idx: spearmanr(x[idx], y[idx]).statistic,
                              N_BOOT, SEED_BOOT)[:, 0]
    ci95 = list(percentile_interval(draws, 95.0))
    return {"label": label, "n": int(len(x)), "n_targets": int(len(set(targets))),
            "rho": rho, "ci95": ci95, "ci90": list(percentile_interval(draws, 90.0)),
            "excludes_zero": excludes_zero(ci95), "practically_null": is_practically_null(draws),
            "target_balanced_rho": target_balanced_spearman(x, y, targets)}, draws


def paired_comparison(reference_rho, alternative_rho, targets, statistic_pair, label):
    """Compare two estimates of one association on the same drugs, as the rule fixes."""
    draws = cluster_bootstrap(targets, statistic_pair, N_BOOT, SEED_BOOT)
    difference = draws[:, 1] - draws[:, 0]
    reading = comparison_reading(reference_rho, alternative_rho, difference, draws[:, 1])
    return {"label": label, "rho_reference": reference_rho, "rho_alternative": alternative_rho,
            "difference": alternative_rho - reference_rho,
            "difference_ci95": list(percentile_interval(difference, 95.0)),
            "alternative_ci95": list(percentile_interval(draws[:, 1], 95.0)),
            "reading": reading}, draws


def pair_projection_table(arm: Arm, reference: Reference, genes) -> np.ndarray:
    """(n_drugs, n_genes) projected dispersion of each drug onto every direction.

    Projecting every pair difference onto every direction once turns a permutation
    into a lookup. Ten thousand permutations of 795 drugs over 258 targets is not
    affordable any other way.
    """
    directions = np.vstack([reference.directions[gene] for gene in genes])
    rows, counts = [], []
    for S in arm.signatures:
        restricted = S[:, reference.positions]
        upper = np.triu_indices(restricted.shape[0], k=1)
        rows.append(restricted[upper[0]] - restricted[upper[1]])
        counts.append(len(upper[0]))
    projections = np.abs(np.vstack(rows) @ directions.T)
    owner = np.repeat(np.arange(len(arm)), counts)
    totals = np.zeros((len(arm), len(genes)))
    np.add.at(totals, owner, projections)
    return totals / np.array(counts)[:, None]


def permutation_null(arm: Arm, values, genes, projection_table, percentiles_or_matrix,
                     statistic, label):
    """A unique-target permutation null for a statistic built from the directions.

    `statistic(P_perm, E_perm)` returns the quantity under one reassignment.
    """
    column = {gene: i for i, gene in enumerate(genes)}
    own = np.array([column[t] for t in arm.targets])
    observed = statistic(values["P"], values["E"])
    rows = np.arange(len(arm))
    null = np.array([statistic(projection_table[rows, assigned],
                               percentiles_or_matrix[rows, assigned])
                     for assigned in unique_target_permutations(arm.targets, N_PERM, SEED_PERM)])
    identity = statistic(projection_table[rows, own], percentiles_or_matrix[rows, own])
    assert abs(identity - observed) < 1e-9, (
        f"{label}: the lookup gives {identity:.12f} where the direct statistic gives {observed:.12f}")
    return {"label": label, "observed": observed, "null_median": float(np.median(null)),
            "null_pct": {"p1": float(np.percentile(null, 1)), "p2.5": float(np.percentile(null, 2.5)),
                         "p50": float(np.percentile(null, 50)),
                         "p97.5": float(np.percentile(null, 97.5)),
                         "p99": float(np.percentile(null, 99))},
            "attenuation": float(observed - np.median(null)),
            "p_two_sided": two_sided_permutation_p(observed, null),
            "reading": permutation_reading(observed, null)}, null


def q_reading(arm: Arm, reference: Reference, label):
    """Own-target percentile, its permutation p-value, and its interval."""
    q, matrix, genes = own_target_q(arm, reference)
    balanced = target_balanced_mean(q, arm.targets)
    draws = cluster_bootstrap(arm.targets, lambda idx: target_balanced_mean(q[idx], arm.targets[idx]),
                              N_BOOT, SEED_BOOT)[:, 0]
    ci95, ci90 = percentile_interval(draws, 95.0), percentile_interval(draws, 90.0)

    table = percentile_table(matrix)
    rows = np.arange(len(arm))
    null = np.array([target_balanced_mean(table[rows, assigned], arm.targets)
                     for assigned in unique_target_permutations(arm.targets, N_PERM, SEED_PERM)])
    p_value = float((1 + int((null >= balanced).sum())) / (N_PERM + 1))

    if p_value < 0.01 and ci95[0] > 0.5:
        reading = "carries target identity"
    elif Q_EQUIVALENCE[0] < ci90[0] and ci90[1] < Q_EQUIVALENCE[1]:
        reading = "lacks target identity"
    else:
        reading = "inconclusive"

    # a target whose direction carries a component common to many drugs scores well
    # for all of them; z-scoring within target removes that advantage
    z_matrix = (matrix - matrix.mean(axis=0)) / np.maximum(matrix.std(axis=0), 1e-12)
    column = {gene: i for i, gene in enumerate(genes)}
    own = np.array([column[t] for t in arm.targets])
    q_z = own_target_percentile(z_matrix, own)

    return {"label": label, "q_bar": balanced, "ci95": list(ci95), "ci90": list(ci90),
            "p_permutation": p_value, "null_median": float(np.median(null)),
            "equivalence_region": list(Q_EQUIVALENCE), "reading": reading,
            "q_bar_within_target_z": target_balanced_mean(q_z, arm.targets),
            "per_target_q": {t: float(q[arm.targets == t].mean())
                             for t in sorted(set(arm.targets))}}, q, draws, null


def leave_one_target_out(x, y, targets, label):
    """Every primary correlation re-estimated with each target omitted."""
    estimates = {}
    for target in sorted(set(targets)):
        keep = targets != target
        estimates[target] = float(spearmanr(x[keep], y[keep]).statistic)
    values = np.array(list(estimates.values()))
    full = float(spearmanr(x, y).statistic)
    worst = max(estimates, key=lambda t: abs(estimates[t] - full))
    return {"label": label, "full": full, "min": float(values.min()), "max": float(values.max()),
            "largest_change_target": worst, "largest_change": float(estimates[worst] - full),
            "sign_stable": bool(np.all(np.sign(values) == np.sign(full)))}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--perturbseq", type=Path, required=True)
    parser.add_argument("--replogle", type=Path, required=True)
    parser.add_argument("--prism", type=Path)
    parser.add_argument("--single-cell", type=Path)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    run(args)


def verify_inputs(args):
    """Every pinned input, checked before anything is computed."""
    pins = {"crispri_records": (CRISPRI_RECORDS, EXPECTED_RECORDS_SHA256),
            "lincs_subset": (args.data / "lincs_subset.npz", EXPECTED_LINCS_SUBSET_SHA256),
            "lincs_shrna": (args.data / "lincs_shrna.npz", EXPECTED_LINCS_SHRNA_SHA256),
            "perturbseq": (args.perturbseq, EXPECTED_PERTURBSEQ_SHA256)}
    for name, path in EXPECTED_BULK_SHA256.items():
        pins[name] = (args.replogle / name, EXPECTED_BULK_SHA256[name])
    recorded = {}
    for name, (path, expected) in pins.items():
        digest = sha256_file(path)
        assert digest == expected, f"{name}: sha256 {digest}, pinned {expected}"
        recorded[name] = digest
    return recorded


def reproduction_gate(arm: Arm, shrna_values, crispri_values, records):
    """The arms must reproduce the corrected per-drug records before anything runs."""
    by_drug = {r["drug"]: r for r in records}
    worst = {"D": 0.0, "P_shrna": 0.0, "E_shrna": 0.0, "P_crispri": 0.0, "E_crispri": 0.0}
    for i, drug in enumerate(arm.drugs):
        record = by_drug[drug]
        assert record["target"] == arm.targets[i], f"{drug}: target differs"
        assert int(record["n_celllines"]) == len(arm.cell_lines[i]), f"{drug}: n_celllines differs"
        reference_D = (shrna_values or crispri_values)["D"]
        worst["D"] = max(worst["D"], abs(reference_D[i] - record["raw_instability"]))
        if "proj_shrna" in record and shrna_values is not None:
            worst["P_shrna"] = max(worst["P_shrna"], abs(shrna_values["P"][i] - record["proj_shrna"]))
            worst["E_shrna"] = max(worst["E_shrna"], abs(shrna_values["E"][i] - record["enrich_shrna"]))
        if "proj_crispri" in record and crispri_values is not None:
            worst["P_crispri"] = max(worst["P_crispri"],
                                     abs(crispri_values["P"][i] - record["proj_crispri"]))
            worst["E_crispri"] = max(worst["E_crispri"],
                                     abs(crispri_values["E"][i] - record["enrich_crispri"]))
    for quantity, error in worst.items():
        assert error < RECON_TOL, f"{quantity} reproduces to {error:.3g}, tolerance {RECON_TOL}"
    return {k: float(v) for k, v in worst.items()}


DOSE_WINDOW = (2.0, 3.0)         # Amendment 1, A1: the registered 2.5 uM with a tolerance


def prism_mapping(prism_dir: Path, drugs, pert_ids) -> dict:
    """Stage one of R5: which treatment stands for which drug, and nothing else.

    Reads the treatment table and the missingness of the response matrix. No
    response value is read, summarized or returned, so the mapping can be frozen
    and hashed before any association exists. Amendment 1, A4.
    """
    treatment_path = prism_dir / "primary-screen-replicate-collapsed-treatment-info.csv"
    lfc_path = prism_dir / "primary-screen-replicate-collapsed-logfold-change.csv"
    treatments = pd.read_csv(treatment_path, low_memory=False)
    responses = pd.read_csv(lfc_path, index_col=0)
    lines_measured = responses.notna().sum(axis=0)       # missingness, not response
    del responses

    treatments["stem"] = treatments.broad_id.str.extract(BRD_STEM)[0]
    treatments = treatments[treatments.column_name.isin(lines_measured.index)].copy()
    treatments["n_lines"] = treatments.column_name.map(lines_measured)
    eligible = treatments[(treatments.n_lines >= R5_MIN_LINES)
                          & (treatments.dose >= DOSE_WINDOW[0])
                          & (treatments.dose <= DOSE_WINDOW[1])]

    def normalize(text):
        return " ".join(str(text).lower().split())

    by_stem, by_name = {}, {}
    for _, row in eligible.iterrows():
        by_stem.setdefault(row.stem, []).append(row.column_name)
        by_name.setdefault(normalize(row["name"]), []).append(row.column_name)

    accepted, rejected = {}, {}
    for drug in drugs:
        stems = {stem for pid in pert_ids.get(drug, []) for stem in BRD_STEM.findall(str(pid))}
        candidates = [c for stem in sorted(stems) for c in by_stem.get(stem, [])]
        route = "broad identifier"
        if not candidates:
            candidates = by_name.get(normalize(drug), [])
            route = "exact name"
        if not candidates:
            rejected[drug] = (f"no treatment inside {DOSE_WINDOW[0]}-{DOSE_WINDOW[1]} uM with "
                              f"log-fold change in at least {R5_MIN_LINES} cell lines, by "
                              "identifier or by exact name in the release")
            continue
        best = max(sorted(candidates), key=lambda column: (lines_measured[column], column))
        row = eligible[eligible.column_name == best].iloc[0]
        accepted[drug] = {"column": best, "stem": row.stem, "route": route,
                          "n_lines": int(lines_measured[best]), "dose": float(row.dose),
                          "screen": row.screen_id, "n_candidates": len(set(candidates))}

    return {"amendment": "experiments/PREREG_H3_REFERENCE_DISCORDANCE_AMENDMENT_1.md",
            "dose_window": list(DOSE_WINDOW), "min_lines": R5_MIN_LINES,
            "accepted": accepted, "rejected": rejected,
            "n_accepted": len(accepted), "n_rejected": len(rejected),
            "provenance": {"lfc_file": str(lfc_path), "treatment_file": str(treatment_path),
                           "lfc_sha256": sha256_file(lfc_path),
                           "treatment_sha256": sha256_file(treatment_path),
                           "no_response_value_was_read": True}}


def prism_response(prism_dir: Path, mapping_path: Path):
    """Stage two of R5: the toxicity measure, from a mapping it did not choose.

    The stage accepts only a mapping written to disk by stage one, and refuses to
    build one of its own, so no mapping can be revised after a response value has
    been summarized.
    """
    assert mapping_path.exists(), (
        f"{mapping_path} is missing: the mapping must be frozen before any response "
        "value is summarized")
    frozen = json.loads(mapping_path.read_text())
    assert frozen.get("accepted"), "the frozen mapping accepts no drug"
    lfc_path = prism_dir / "primary-screen-replicate-collapsed-logfold-change.csv"
    assert sha256_file(lfc_path) == frozen["provenance"]["lfc_sha256"], (
        "the response matrix is not the one the mapping was frozen against")

    columns = sorted({entry["column"] for entry in frozen["accepted"].values()})
    responses = pd.read_csv(lfc_path, index_col=0, usecols=["Unnamed: 0"] + columns)
    toxicity = {column: -float(responses[column].median(skipna=True)) for column in columns}
    killed = {column: float((responses[column] < -1).sum() / responses[column].notna().sum())
              for column in columns}
    return toxicity, killed, {"mapping_file": str(mapping_path),
                              "mapping_sha256": sha256_file(mapping_path),
                              "n_columns_read": len(columns)}


def run(args):
    OUT.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    draws = Draws()
    result = {"registration": {"file": "experiments/PREREG_H3_REFERENCE_DISCORDANCE.md",
                               "frozen_at": "7f57136"},
              "inputs": verify_inputs(args), "modules": {}}

    log("building drug signatures and references")
    per_drug, gene_ids = build_drug_signatures(args.data)
    symbols = landmark_symbols(args.data / "GSE92742_Broad_LINCS_gene_info.txt.gz", gene_ids)
    targets, multi_target, moa = drug_targets(args.data)
    records = json.loads(CRISPRI_RECORDS.read_text())

    shrna = build_shrna_reference(args.data)
    c0 = pooled_crispri_reference(args.perturbseq, symbols, name="C0")
    c1 = load_replogle_bulk(args.replogle / "K562_essential_normalized_bulk_01.h5ad",
                            symbols, "C1-K562")
    c1_rpe1 = load_replogle_bulk(args.replogle / "rpe1_normalized_bulk_01.h5ad", symbols, "C1-RPE1")
    c1_gw = load_replogle_bulk(args.replogle / "K562_gwps_normalized_bulk_01.h5ad", symbols, "C1-GW")
    c1_gw_positive = load_replogle_bulk(
        args.replogle / "K562_gwps_normalized_bulk_01.h5ad", symbols, "C1-GW-phenotype-positive",
        qualifying_rows=lambda obs: obs["energy_test_p_value"] < ENERGY_P)

    result["modules"]["R0.2_gene_space"] = {
        reference.name: {k: v for k, v in reference.audit.items()
                         if k not in ("rows_per_gene", "rows_used_per_gene", "cells_per_gene")}
        for reference in (c0, c1, c1_rpe1, c1_gw, c1_gw_positive)}
    assert len(c1.positions) >= MIN_LANDMARKS, (
        f"C1 covers {len(c1.positions)} landmarks, below the registered floor of {MIN_LANDMARKS}")
    result["modules"]["R0.2_gene_space"]["landmark_floor"] = {
        "floor": MIN_LANDMARKS, "C1_landmarks": int(len(c1.positions)),
        "note": "a usability floor, not evidence that C1 covers the landmark space"}

    log("assembling cohorts")
    crispri_drugs = {r["drug"] for r in records if "proj_crispri" in r}
    shrna_arm = assemble_arm("shRNA arm", per_drug, targets, shrna,
                             drugs={r["drug"] for r in records if "proj_shrna" in r})
    crispri_arm = assemble_arm("CRISPRi arm", per_drug, targets, c0, drugs=crispri_drugs)
    matched_drugs = set(shrna_arm.drugs) & set(crispri_arm.drugs)
    gw_arm = assemble_arm("genome-wide arm", per_drug, targets, c1_gw)

    shrna_values = quantities(shrna_arm, shrna)
    c0_values = quantities(crispri_arm, c0)
    # the CRISPRi arm reproduces its CRISPRi fields, and its shRNA fields wherever
    # the drug also carries a shRNA direction
    crispri_with_shrna = crispri_arm.subset(
        [i for i, t in enumerate(crispri_arm.targets) if t in shrna.directions])
    result["gate"] = {
        "shrna": reproduction_gate(shrna_arm, shrna_values, None, records),
        "crispri": reproduction_gate(crispri_arm, None, c0_values, records),
        "crispri_shrna_fields": reproduction_gate(
            crispri_with_shrna, quantities(crispri_with_shrna, shrna), None, records)}
    log(f"reproduction holds; {len(shrna_arm)} shRNA drugs, {len(crispri_arm)} CRISPRi drugs")

    # the harmonized comparators: the shRNA quantities and raw instability on C1's
    # gene space, so that a comparison between references is not also a comparison
    # between gene spaces
    space = shared_space(c1, c1_rpe1)
    c1_shared, c1_rpe1_shared = c1.restricted_to(space), c1_rpe1.restricted_to(space)
    shrna_h = shrna.restricted_to(c1.positions)
    c0h = c0.restricted_to(c1.positions)

    c1_values = quantities(crispri_arm, c1)
    c0h_values = quantities(crispri_arm, c0h)

    # ---- R0.5 agreement between constructions, R0.8 concentration
    common = sorted(set(crispri_arm.targets) & set(shrna.directions))
    agreement = {}
    for name, reference in (("C0_vs_C1", (c0h, c1)), ("C1_vs_RPE1", (c1_shared, c1_rpe1_shared))):
        left, right = reference
        cosines = {gene: float(abs(left.directions[gene] @ right.directions[gene]))
                   for gene in sorted(set(left.directions) & set(right.directions) & set(crispri_arm.targets))}
        agreement[name] = {"median": float(np.median(list(cosines.values()))),
                           "min": float(min(cosines.values())), "max": float(max(cosines.values())),
                           "per_target": cosines}
    result["modules"]["R0.5_agreement"] = agreement

    def leading_share(reference, genes):
        matrix = np.vstack([reference.directions[g] for g in genes])
        similarity = matrix @ matrix.T
        eigenvalues = np.linalg.eigvalsh(similarity)
        return float(eigenvalues[-1] / eigenvalues.sum())

    concentration_draws = cluster_bootstrap(
        np.array(common),
        lambda idx: leading_share(c1, [common[i] for i in idx]) -
                    leading_share(shrna, [common[i] for i in idx]),
        1000, SEED_BOOT)[:, 0]
    draws.keep("R0.8_concentration_difference", concentration_draws)
    result["modules"]["R0.8_concentration"] = {
        "n_targets": len(common),
        "C1_leading_share": leading_share(c1, common),
        "shRNA_leading_share": leading_share(shrna, common),
        "difference_ci95": list(percentile_interval(concentration_draws, 95.0)),
        "reading": ("CRISPRi directions more concentrated"
                    if percentile_interval(concentration_draws, 95.0)[0] > 0 else "inconclusive")}

    # ---- R0.9 the headline correlations under each construction
    headline = {}
    for name, values in (("C0", c0_values), ("C0h", c0h_values), ("C1", c1_values)):
        for quantity in ("P", "D"):
            summary, drawn = association(values[quantity], values["E"], crispri_arm.targets,
                                         f"{name}: rho({quantity}, E)")
            headline[f"{name}_{quantity}_E"] = summary
            draws.keep(f"R0.9_{name}_{quantity}_E", drawn)
    raw_comparison, raw_comparison_draws = paired_comparison(
        headline["C0_D_E"]["rho"], headline["C1_D_E"]["rho"], crispri_arm.targets,
        lambda idx: (spearmanr(c0_values["D"][idx], c0_values["E"][idx]).statistic,
                     spearmanr(c1_values["D"][idx], c1_values["E"][idx]).statistic),
        "C1 against C0, rho(D, E)")
    projected_comparison, projected_comparison_draws = paired_comparison(
        headline["C0_P_E"]["rho"], headline["C1_P_E"]["rho"], crispri_arm.targets,
        lambda idx: (spearmanr(c0_values["P"][idx], c0_values["E"][idx]).statistic,
                     spearmanr(c1_values["P"][idx], c1_values["E"][idx]).statistic),
        "C1 against C0, rho(P, E)")
    draws.keep("R0.9_C1_vs_C0_raw_paired", raw_comparison_draws)
    draws.keep("R0.9_C1_vs_C0_projected_paired", projected_comparison_draws)
    reversal = {"attenuated": "construction-dependent", "reversed": "construction-dependent",
                "retained": "construction-robust"}.get(raw_comparison["reading"], "inconclusive")
    result["modules"]["R0.9_constructions"] = {"headline": headline, "raw": raw_comparison,
                                               "projected": projected_comparison,
                                               "raw_reversal_reading": reversal}

    # ---- R1 dependence, weighting, composition
    r1 = {"dependence": {}, "weighting": {}, "composition": {}}
    for arm, values, name in ((shrna_arm, shrna_values, "shRNA"),
                              (crispri_arm, c1_values, "C1"), (crispri_arm, c0_values, "C0")):
        for quantity in ("P", "D"):
            summary, drawn = association(values[quantity], values["E"], arm.targets,
                                         f"{name}: rho({quantity}, E)")
            r1["dependence"][f"{name}_{quantity}_E"] = summary
            draws.keep(f"R1_{name}_{quantity}_E", drawn)
            weighting, weighting_draws = paired_comparison(
                summary["rho"], summary["target_balanced_rho"], arm.targets,
                lambda idx, v=values, q=quantity: (
                    spearmanr(v[q][idx], v["E"][idx]).statistic,
                    target_balanced_spearman(v[q][idx], v["E"][idx], arm.targets[idx])),
                f"{name}: target-balanced against drug-weighted, rho({quantity}, E)")
            r1["weighting"][f"{name}_{quantity}_E"] = weighting
            draws.keep(f"R1_weighting_{name}_{quantity}_E_paired", weighting_draws)

    matched_index = [i for i, d in enumerate(crispri_arm.drugs) if d in matched_drugs]
    matched_arm = crispri_arm.subset(matched_index)
    matched_c1 = quantities(matched_arm, c1)
    matched_shrna_h = quantities(matched_arm, shrna_h)
    matched_shrna = quantities(matched_arm, shrna)
    for quantity in ("P", "D"):
        comparison, composition_draws = paired_comparison(
            float(spearmanr(matched_shrna_h[quantity], matched_shrna_h["E"]).statistic),
            float(spearmanr(matched_c1[quantity], matched_c1["E"]).statistic),
            matched_arm.targets,
            lambda idx, q=quantity: (
                spearmanr(matched_shrna_h[q][idx], matched_shrna_h["E"][idx]).statistic,
                spearmanr(matched_c1[q][idx], matched_c1["E"][idx]).statistic),
            f"matched: C1 against shRNA-h, rho({quantity}, E)")
        r1["composition"][f"{quantity}_E_harmonized"] = comparison
        draws.keep(f"R1_composition_{quantity}_E_paired", composition_draws)
        r1["composition"][f"{quantity}_E_unharmonized"] = {
            "rho_shRNA_full_space": float(spearmanr(matched_shrna[quantity],
                                                    matched_shrna["E"]).statistic),
            "note": "reported beside the harmonized comparison; the gap is the gene-space part"}
    r1["composition"]["n_matched"] = len(matched_arm)
    result["modules"]["R1_dependence_weighting_composition"] = r1

    # ---- R2 target identity
    r2 = {}
    for arm, reference, name in ((shrna_arm, shrna, "shRNA"), (crispri_arm, c1, "C1"),
                                 (crispri_arm, c1_rpe1, "C1-RPE1"), (crispri_arm, c0, "C0")):
        r2[name], _, q_draws, q_null = q_reading(arm, reference, f"{name}: own-target percentile")
        draws.keep(f"R2_{name}_qbar_bootstrap", q_draws)
        draws.keep(f"R2_{name}_qbar_null", q_null)
    result["modules"]["R2_target_identity"] = r2

    # ---- R3 permutation resistance
    r3 = {}
    for arm, reference, values, name in ((crispri_arm, c1, c1_values, "C1"),
                                         (crispri_arm, c0, c0_values, "C0"),
                                         (shrna_arm, shrna, shrna_values, "shRNA")):
        _, matrix, genes = own_target_q(arm, reference)
        projections = pair_projection_table(arm, reference, genes)
        for quantity in ("D", "P"):
            def statistic(P_perm, E_perm, v=values, q=quantity):
                return float(spearmanr(v["D"] if q == "D" else P_perm, E_perm).statistic)
            summary, null = permutation_null(arm, values, genes, projections, matrix, statistic,
                                             f"{name}: rho({quantity}, E) under unique-target permutation")
            r3[f"{name}_{quantity}_E"] = summary
            draws.keep(f"R3_{name}_{quantity}_E_null", null)
    result["modules"]["R3_permutation_resistance"] = r3
    result["modules"]["R3_permutation_resistance"]["note"] = (
        "permutation resistance is not attributed to any program here")

    # ---- R4 essential-gene composition
    essential_targets = set(c1.directions)
    r4 = {}
    for label, reference in (("phenotype_positive", c1_gw_positive), ("unrestricted", c1_gw)):
        arm = assemble_arm(f"genome-wide {label}", per_drug, targets, reference)
        values = quantities(arm, reference)
        is_essential = np.array([t in essential_targets for t in arm.targets])
        subsets = {}
        for subset_name, mask in (("essential", is_essential), ("other", ~is_essential)):
            if mask.sum() < 10:
                subsets[subset_name] = {"n": int(mask.sum()), "reading": "too few drugs"}
                continue
            for quantity in ("P", "D"):
                summary, _ = association(values[quantity][mask], values["E"][mask],
                                         arm.targets[mask], f"{label}/{subset_name}: rho({quantity}, E)")
                subsets[f"{subset_name}_{quantity}_E"] = summary
        other_mask = ~is_essential
        gate_met = bool(other_mask.sum() >= R4_MIN_DRUGS and
                        len(set(arm.targets[other_mask])) >= R4_MIN_TARGETS)
        r4[label] = {"n_drugs": len(arm), "n_essential": int(is_essential.sum()),
                     "n_other": int(other_mask.sum()),
                     "n_other_targets": int(len(set(arm.targets[other_mask]))),
                     "gate_met": gate_met, "subsets": subsets,
                     "note": ("" if gate_met else
                              "the registered minimum cohort was not reached, so the estimates are "
                              "descriptive and neither support nor count against an "
                              "essential-gene-composition explanation")}
    result["modules"]["R4_essential_composition"] = r4

    # ---- R5 toxicity
    if args.prism is None:
        result["modules"]["R5_toxicity"] = {"runnable": False, "reason": "no PRISM release given"}
    else:
        siginfo = pd.read_csv(args.data / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                              low_memory=False, usecols=["pert_iname", "pert_id"])
        pert_ids = {drug: sorted(set(group.pert_id.astype(str)))
                    for drug, group in siginfo.groupby("pert_iname")}

        mapping_path = args.output / "r5_prism_mapping.json"
        if not mapping_path.exists():
            mapping = prism_mapping(args.prism, crispri_arm.drugs, pert_ids)
            mapping_path.write_text(json.dumps(mapping, indent=2))
            log(f"mapping frozen at {mapping_path}; commit it before the response stage")
        frozen = json.loads(mapping_path.read_text())
        accepted = frozen["accepted"]

        mapped = [i for i, drug in enumerate(crispri_arm.drugs) if drug in accepted]
        n_targets = len(set(crispri_arm.targets[mapped]))
        gate_met = len(mapped) >= R5_MIN_DRUGS and n_targets >= R5_MIN_TARGETS
        r5 = {"runnable": gate_met, "n_mapped": len(mapped), "n_targets": n_targets,
              "gate": {"min_drugs": R5_MIN_DRUGS, "min_targets": R5_MIN_TARGETS},
              "dose_window": frozen["dose_window"], "mapping_table": mapping_path.name,
              "mapping_sha256": sha256_file(mapping_path), "provenance": frozen["provenance"],
              "reading": None}
        if not gate_met:
            r5["reading"] = ("the registered minimum cohort was not reached; R5 is reported "
                             "descriptively and labeled exploratory, and it neither supports "
                             "nor counts against the toxicity explanation")
        else:
            toxicity, killed, response_provenance = prism_response(args.prism, mapping_path)
            r5["response_provenance"] = response_provenance
            idx = np.array(mapped)
            T = np.array([toxicity[accepted[crispri_arm.drugs[i]]["column"]] for i in idx])
            fraction_killed = np.array([killed[accepted[crispri_arm.drugs[i]]["column"]]
                                        for i in idx])
            sub_targets = crispri_arm.targets[idx]
            e_t, e_t_draws = association(c1_values["E"][idx], T, sub_targets, "rho(E, toxicity)")
            d_e, d_e_draws = association(c1_values["D"][idx], c1_values["E"][idx], sub_targets,
                                         "rho(D, E) on the complete-case cohort")
            draws.keep("R5_rho_E_toxicity", e_t_draws)
            draws.keep("R5_rho_D_E_complete_case", d_e_draws)
            partial = rank_partial_correlation(c1_values["D"][idx], c1_values["E"][idx], (T,))
            partial_draws = cluster_bootstrap(
                sub_targets,
                lambda rows: (spearmanr(c1_values["D"][idx][rows], c1_values["E"][idx][rows]).statistic,
                              rank_partial_correlation(c1_values["D"][idx][rows],
                                                       c1_values["E"][idx][rows], (T[rows],))),
                N_BOOT, SEED_BOOT)
            draws.keep("R5_partial_paired", partial_draws)
            attenuation = partial_draws[:, 0] - partial_draws[:, 1]
            reading_rule = comparison_reading(d_e["rho"], partial,
                                              partial_draws[:, 1] - partial_draws[:, 0],
                                              partial_draws[:, 1])
            if e_t["ci95"][0] > 0 and (is_practically_null(partial_draws[:, 1])
                                       or reading_rule in ("attenuated", "reversed")):
                reading = "toxicity accounts for the raw association"
            elif e_t["practically_null"] or reading_rule == "retained":
                reading = "counts against"
            else:
                reading = "inconclusive"
            r5.update({"rho_E_toxicity": e_t, "rho_D_E": d_e, "partial_rho_D_E_given_T": partial,
                       "partial_ci95": list(percentile_interval(partial_draws[:, 1], 95.0)),
                       "attenuation": float(d_e["rho"] - partial),
                       "attenuation_ci95": list(percentile_interval(attenuation, 95.0)),
                       "ratio": (float(partial / d_e["rho"]) if d_e["rho"] >= R5_RATIO_DENOMINATOR
                                 else None),
                       "fraction_of_lines_below_minus_one": {
                           "median": float(np.median(fraction_killed))},
                       "comparison_reading": reading_rule, "reading": reading,
                       "measure": ("viability response at a near-common exposure of "
                                   f"{DOSE_WINDOW[0]}-{DOSE_WINDOW[1]} uM; not a dose-response "
                                   "summary and not a potency estimate")})
        result["modules"]["R5_toxicity"] = r5

    # ---- R6 cell context
    r6 = {}
    for quantity in ("P", "D"):
        k562 = quantities(crispri_arm, c1_shared)
        rpe1 = quantities(crispri_arm, c1_rpe1_shared)
        comparison, context_draws = paired_comparison(
            float(spearmanr(k562[quantity], k562["E"]).statistic),
            float(spearmanr(rpe1[quantity], rpe1["E"]).statistic),
            crispri_arm.targets,
            lambda idx, q=quantity: (spearmanr(k562[q][idx], k562["E"][idx]).statistic,
                                     spearmanr(rpe1[q][idx], rpe1["E"][idx]).statistic),
            f"RPE1 against K562, rho({quantity}, E)")
        r6[f"R6a_{quantity}_E"] = comparison
        draws.keep(f"R6a_{quantity}_E_paired", context_draws)
    r6["R6a_reading"] = {"retained": "raw reversal replicates in RPE1",
                         "attenuated": "raw reversal is context-dependent",
                         "reversed": "raw reversal is context-dependent"}.get(
        r6["R6a_D_E"]["reading"], "inconclusive")

    cell_info = pd.read_csv(args.data / "GSE92742_Broad_LINCS_cell_info.txt.gz", sep="\t",
                            low_memory=False)
    heme_lines = set(HEMATOPOIETIC) & set(cell_info.cell_id)
    r6["R6b_exploratory"] = hematopoietic_probe(crispri_arm, c1, heme_lines, MYELOID)
    result["modules"]["R6_cell_context"] = r6

    # ---- R7 sensitivities
    r7 = {"R7a_single_target": {}, "R7b_signed": {}, "R7d_leave_one_target_out": {},
          "R7e_exclusions": {}, "R7f_target_gene": {}}
    for arm, values, reference, name in ((shrna_arm, shrna_values, shrna, "shRNA"),
                                         (crispri_arm, c1_values, c1, "C1")):
        single = np.array([not multi_target.get(d, False) for d in arm.drugs])
        r7["R7a_single_target"][name] = {
            "n_single": int(single.sum()), "n_multi": int((~single).sum())}
        for quantity in ("P", "D"):
            if single.sum() >= 10:
                summary, _ = association(values[quantity][single], values["E"][single],
                                         arm.targets[single], f"{name}: single-target rho({quantity}, E)")
                r7["R7a_single_target"][f"{name}_{quantity}_E"] = summary
            r7["R7d_leave_one_target_out"][f"{name}_{quantity}_E"] = leave_one_target_out(
                values[quantity], values["E"], arm.targets, f"{name}: rho({quantity}, E)")

        inhibitors = np.array(["inhibitor" in moa.get(d, "").lower() for d in arm.drugs])
        if inhibitors.sum() >= 10:
            means = arm.means(reference.positions)
            signed = np.array([float(means[i] @ reference.directions[arm.targets[i]] /
                                     np.linalg.norm(means[i]))
                               for i in range(len(arm))])
            r7["R7b_signed"][name] = {
                "n_inhibitors": int(inhibitors.sum()),
                "mean_signed_cosine": target_balanced_mean(signed[inhibitors],
                                                           arm.targets[inhibitors]),
                "rho_D_signed": float(spearmanr(values["D"][inhibitors], signed[inhibitors]).statistic),
                "rho_P_signed": float(spearmanr(values["P"][inhibitors], signed[inhibitors]).statistic)}

        without_gene = {}
        # one cut reference per target, not one per drug
        cuts = {target: reference.without_gene(target, symbols) for target in set(arm.targets)}
        changed_targets = {t for t, cut in cuts.items() if cut is not reference}
        changed = sum(t in changed_targets for t in arm.targets)
        P_no, E_no = np.empty(len(arm)), np.empty(len(arm))
        for i, target in enumerate(arm.targets):
            cut = cuts[target]
            S = arm.signatures[i][:, cut.positions]
            u = cut.directions[target]
            upper = np.triu_indices(S.shape[0], k=1)
            P_no[i] = projected_dispersion(S[upper[0]] - S[upper[1]], u)
            mean = S.mean(axis=0)
            E_no[i] = float((mean @ u / (np.linalg.norm(mean) * np.linalg.norm(u))) ** 2)
        for quantity, series in (("P", P_no), ("D", values["D"])):
            summary, _ = association(series, E_no, arm.targets,
                                     f"{name}: rho({quantity}, E) without the target's own gene")
            without_gene[f"{quantity}_E"] = summary
        without_gene["n_drugs_changed"] = int(changed)
        without_gene["n_targets_changed"] = int(len(changed_targets))
        r7["R7f_target_gene"][name] = without_gene

    eligible = {d for d, t in targets.items() if d in per_drug and len(per_drug[d]) >= MIN_CELL_LINES}
    r7["R7e_exclusions"] = {
        "n_eligible_drugs": len(eligible),
        "n_without_shrna_reference": len(eligible - set(shrna_arm.drugs)),
        "n_without_crispri_reference": len(eligible - set(crispri_arm.drugs)),
        "note": "no drug is filtered on alignment"}

    shrna_reliabilities = shrna_reliability(args.data, set(shrna_arm.targets))
    result["modules"]["R0.6_shrna_reliability"] = {
        "n_targets": len(shrna_reliabilities),
        "median": float(np.median(list(shrna_reliabilities.values()))) if shrna_reliabilities else None,
        "per_target": shrna_reliabilities}
    if args.single_cell is not None:
        audit = json.loads(args.single_cell.read_text())
        audit.setdefault("R0.6_reliability", {})["shRNA"] = shrna_reliabilities
        result["modules"]["R0_single_cell"] = audit
        r7["R7c_reliability"] = reliability_restricted(audit, crispri_arm, c1, shrna_arm, shrna)
    else:
        r7["R7c_reliability"] = {"runnable": False,
                                 "reason": "no single-cell audit given, so CRISPRi reliability is unknown"}
    result["modules"]["R7_sensitivities"] = r7

    draws_path = args.output / "reference_discordance_draws.npz"
    np.savez_compressed(draws_path, **draws)
    result["draws"] = {"file": draws_path.name, "keys": sorted(draws)}
    result["provenance"] = {"script_sha256": sha256_file(Path(__file__)),
                            "draws_sha256": sha256_file(draws_path),
                            "seeds": {"bootstrap": SEED_BOOT, "permutation": SEED_PERM},
                            "n_bootstrap": N_BOOT, "n_permutations": N_PERM,
                            "numpy": np.__version__, "pandas": pd.__version__}
    path = args.output / "reference_discordance_results.json"
    path.write_text(json.dumps(result, indent=2, default=float))
    (args.output / "reference_discordance_results.json.sha256").write_text(sha256_file(path) + "\n")
    log(f"written to {path}")


def hematopoietic_probe(arm: Arm, reference: Reference, heme_lines, myeloid_lines):
    """R6b, exploratory: own-target percentile in blood lines against the others."""
    out = {"note": "exploratory; carries no reading in either direction"}
    for label, lines in (("hematopoietic", set(heme_lines)), ("myeloid", set(myeloid_lines))):
        rows, targets_kept = [], []
        for i in range(len(arm)):
            in_lines = np.array([c in lines for c in arm.cell_lines[i]])
            if in_lines.any() and (~in_lines).any():
                rows.append((i, in_lines))
                targets_kept.append(arm.targets[i])
        if len(rows) < 10:
            out[label] = {"n": len(rows), "reading": "too few drugs"}
            continue
        genes = sorted(set(arm.targets))
        column = {gene: i for i, gene in enumerate(genes)}
        inside, outside = [], []
        for i, in_lines in rows:
            S = arm.signatures[i][:, reference.positions]
            for group, store in ((in_lines, inside), (~in_lines, outside)):
                mean = S[group].mean(axis=0)[None, :]
                matrix = alignment_matrix(mean, reference, genes)
                store.append(own_target_percentile(matrix, np.array([column[arm.targets[i]]]))[0])
        difference = np.array(inside) - np.array(outside)
        targets_kept = np.array(targets_kept)
        draws = cluster_bootstrap(targets_kept,
                                 lambda idx: target_balanced_mean(difference[idx], targets_kept[idx]),
                                 N_BOOT, SEED_BOOT)[:, 0]
        out[label] = {"n_drugs": len(rows), "n_targets": int(len(set(targets_kept))),
                      "mean_difference": target_balanced_mean(difference, targets_kept),
                      "ci95": list(percentile_interval(draws, 95.0))}
    return out


def shrna_reliability(data_dir: Path, targets_wanted) -> dict:
    """Split-half reliability of each shRNA consensus, split by hairpin.

    A target with one hairpin measured in several cell lines has no hairpin split
    to make, so its signatures are halved at random instead.
    """
    shrna = np.load(data_dir / "lincs_shrna.npz", allow_pickle=True)
    signatures = shrna["signatures"]
    position = {str(sig_id): i for i, sig_id in enumerate(shrna["sig_ids"])}
    siginfo = pd.read_csv(data_dir / "lincs_shrna_siginfo.csv.gz")
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(position)
                      & siginfo.pert_iname.isin(set(targets_wanted))].copy()
    siginfo["_row"] = siginfo.sig_id.astype(str).map(position)

    out = {}
    for gene, group in siginfo.groupby("pert_iname"):
        rows = group.drop_duplicates("sig_id")
        if len(rows) < MIN_SIGNATURES:
            continue
        matrix = signatures[rows._row.values]
        hairpins = rows.pert_id.astype(str).to_numpy() if "pert_id" in rows else \
            np.array([str(i) for i in range(len(rows))])
        out[str(gene)] = split_half_reliability_by_unit(matrix, hairpins, SEED_SPLIT)
    return out


def reliability_restricted(single_cell, crispri_arm, c1, shrna_arm, shrna):
    """R7c: the primary correlations on targets whose split-half reliability holds up."""
    reliability = single_cell.get("R0.6_reliability", {})
    crispri_reliable = {t for t, v in reliability.get("C1-K562", {}).items()
                        if v >= R7C_MIN_RELIABILITY}
    shrna_reliable = {t for t, v in reliability.get("shRNA", {}).items()
                      if v >= R7C_MIN_RELIABILITY}
    out = {"threshold": R7C_MIN_RELIABILITY,
           "n_crispri_targets": len(crispri_reliable), "n_shrna_targets": len(shrna_reliable)}
    if len(crispri_reliable) < R7C_MIN_CRISPRI_TARGETS or len(shrna_reliable) < R7C_MIN_SHRNA_TARGETS:
        out["runnable"] = False
        out["reason"] = "too few targets survive the reliability threshold"
        return out
    out["runnable"] = True
    for arm, reference, keep, name in ((crispri_arm, c1, crispri_reliable, "C1"),
                                       (shrna_arm, shrna, shrna_reliable, "shRNA")):
        mask = np.array([t in keep for t in arm.targets])
        values = quantities(arm.subset(np.flatnonzero(mask)), reference)
        for quantity in ("P", "D"):
            summary, _ = association(values[quantity], values["E"], arm.targets[mask],
                                     f"{name}: reliable targets, rho({quantity}, E)")
            out[f"{name}_{quantity}_E"] = summary
    return out


if __name__ == "__main__":
    main()
