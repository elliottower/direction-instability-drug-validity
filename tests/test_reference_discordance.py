import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from geometry.direction_instability import direction_instability
from geometry.references import Reference, unit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
_mod = __import__("03d_h3_reference_discordance")
Arm = _mod.Arm


def _arm(n_drugs=24, n_targets=6, n_genes=40, contexts=(5, 9)):
    rng = np.random.default_rng()
    targets = np.array([f"T{i % n_targets}" for i in range(n_drugs)])
    signatures, cells = [], []
    for _ in range(n_drugs):
        k = int(rng.integers(*contexts))
        signatures.append(rng.standard_normal((k, n_genes)) * rng.uniform(0.5, 3.0))
        cells.append([f"LINE{j}" for j in range(k)])
    drugs = [f"drug{i}" for i in range(n_drugs)]
    return Arm("synthetic", drugs, targets, signatures, cells)


def _reference(arm, n_genes=40, positions=None, aligned_with=None):
    rng = np.random.default_rng()
    positions = np.arange(n_genes) if positions is None else positions
    directions = {}
    for target in sorted(set(arm.targets)):
        if aligned_with is None:
            directions[target] = unit(rng.standard_normal(len(positions)))
        else:
            member = np.flatnonzero(arm.targets == target)[0]
            directions[target] = unit(arm.signatures[member][:, positions].mean(axis=0)
                                      + aligned_with * rng.standard_normal(len(positions)))
    return Reference("synthetic", directions, positions)


def test_quantities_are_the_definitions_not_an_approximation():
    arm = _arm()
    reference = _reference(arm)
    values = _mod.quantities(arm, reference)

    for i in range(len(arm)):
        S = arm.signatures[i]
        u = reference.directions[arm.targets[i]]
        upper = np.triu_indices(S.shape[0], k=1)
        expected_P = np.abs((S[upper[0]] - S[upper[1]]) @ u).mean()
        mean = S.mean(axis=0)
        expected_E = (mean @ u / (np.linalg.norm(mean) * np.linalg.norm(u))) ** 2
        assert values["P"][i] == pytest.approx(expected_P)
        assert values["E"][i] == pytest.approx(expected_E)
        assert values["D"][i] == pytest.approx(direction_instability(S))
        assert values["K"][i] == S.shape[0]


def test_raw_instability_on_a_narrower_gene_space_is_reported_separately():
    arm = _arm(n_genes=40)
    narrow = _reference(arm, positions=np.arange(0, 40, 2))
    values = _mod.quantities(arm, narrow)

    assert values["D"] == pytest.approx([direction_instability(S) for S in arm.signatures])
    assert values["D_h"] == pytest.approx(
        [direction_instability(S[:, narrow.positions]) for S in arm.signatures])
    assert not np.allclose(values["D"], values["D_h"])


def test_pair_projection_table_equals_the_direct_projection_for_every_direction():
    arm = _arm(n_drugs=10, n_targets=4)
    reference = _reference(arm)
    genes = sorted(set(arm.targets))
    table = _mod.pair_projection_table(arm, reference, genes)

    for i in range(len(arm)):
        S = arm.signatures[i]
        upper = np.triu_indices(S.shape[0], k=1)
        differences = S[upper[0]] - S[upper[1]]
        for j, gene in enumerate(genes):
            assert table[i, j] == pytest.approx(
                np.abs(differences @ reference.directions[gene]).mean())


def test_permutation_null_identity_reproduces_the_observed_statistic():
    arm = _arm(n_drugs=18, n_targets=6)
    reference = _reference(arm, aligned_with=0.5)
    values = _mod.quantities(arm, reference)
    q, matrix, genes = _mod.own_target_q(arm, reference)
    table = _mod.pair_projection_table(arm, reference, genes)

    def statistic(P_perm, E_perm):
        return float(spearmanr(P_perm, E_perm).statistic)

    summary, null = _mod.permutation_null(arm, values, genes, table, matrix, statistic, "test")
    assert summary["observed"] == pytest.approx(
        spearmanr(values["P"], values["E"]).statistic)
    assert len(null) == _mod.N_PERM
    assert summary["p_two_sided"] > 0


def test_own_target_percentile_is_one_when_directions_are_the_drug_means():
    arm = _arm(n_drugs=12, n_targets=12)          # one drug per target
    reference = _reference(arm, aligned_with=0.0)  # each direction is its drug's mean
    q, matrix, genes = _mod.own_target_q(arm, reference)
    assert q == pytest.approx(np.ones(len(arm)))


def test_reproduction_gate_accepts_the_values_it_was_built_from_and_rejects_a_drift():
    arm = _arm()
    reference = _reference(arm)
    values = _mod.quantities(arm, reference)
    records = [{"drug": arm.drugs[i], "target": arm.targets[i],
                "n_celllines": len(arm.cell_lines[i]),
                "raw_instability": values["D"][i],
                "proj_shrna": values["P"][i], "enrich_shrna": values["E"][i]}
               for i in range(len(arm))]

    worst = _mod.reproduction_gate(arm, values, None, records)
    assert max(worst.values()) < 1e-12

    records[3]["proj_shrna"] += 1e-4
    with pytest.raises(AssertionError):
        _mod.reproduction_gate(arm, values, None, records)

    records[3]["proj_shrna"] -= 1e-4
    records[5]["n_celllines"] += 1
    with pytest.raises(AssertionError):
        _mod.reproduction_gate(arm, values, None, records)


def test_association_reports_the_target_balanced_estimate_beside_the_drug_weighted_one():
    arm = _arm(n_drugs=30, n_targets=5)
    rng = np.random.default_rng()
    x, y = rng.standard_normal(30), rng.standard_normal(30)

    summary, draws = _mod.association(x, y, arm.targets, "test")
    assert summary["rho"] == pytest.approx(spearmanr(x, y).statistic)
    assert summary["n_targets"] == 5
    assert len(draws) == _mod.N_BOOT
    assert summary["ci95"][0] <= summary["rho"] <= summary["ci95"][1]


def test_leave_one_target_out_finds_the_target_that_carries_the_estimate():
    rng = np.random.default_rng()
    targets = np.array(["A"] * 10 + ["B"] * 10 + ["C"] * 10)
    x = rng.standard_normal(30)
    y = rng.standard_normal(30)
    x[:10] = np.arange(10)          # A alone carries a strong association
    y[:10] = np.arange(10)

    out = _mod.leave_one_target_out(x, y, targets, "test")
    assert out["largest_change_target"] == "A"
    assert abs(out["largest_change"]) > 0.1


def _prism_files(tmp_path, treatments, columns_with_values):
    """A miniature PRISM 19Q4 primary screen: a wide matrix plus its treatment table."""
    lines = [f"ACH-{i:06d}" for i in range(120)]
    matrix = pd.DataFrame({column: values(lines) for column, values in columns_with_values.items()},
                          index=lines)
    matrix.to_csv(tmp_path / "primary-screen-replicate-collapsed-logfold-change.csv")
    pd.DataFrame(treatments).to_csv(
        tmp_path / "primary-screen-replicate-collapsed-treatment-info.csv", index=False)
    return tmp_path


def test_prism_mapping_uses_only_the_pinned_files(tmp_path):
    columns = {
        "BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-2.0] * len(lines),
        "BRD-K22222222-001-01-9::2.5::HTS": lambda lines: [-0.1] * len(lines),
        # measured in too few cell lines
        "BRD-K33333333-001-01-9::2.5::HTS": lambda lines: [-3.0] * 20 + [np.nan] * (len(lines) - 20),
    }
    treatments = [
        {"column_name": "BRD-K11111111-001-01-9::2.5::HTS", "broad_id": "BRD-K11111111-001-01-9",
         "name": "Alpha", "dose": 2.5, "screen_id": "HTS"},
        {"column_name": "BRD-K22222222-001-01-9::2.5::HTS", "broad_id": "BRD-K22222222-001-01-9",
         "name": "Beta", "dose": 2.5, "screen_id": "HTS"},
        {"column_name": "BRD-K33333333-001-01-9::2.5::HTS", "broad_id": "BRD-K33333333-001-01-9",
         "name": "Gamma", "dose": 2.5, "screen_id": "HTS"},
    ]
    prism = _prism_files(tmp_path, treatments, columns)

    drugs = ["alpha", "beta", "gamma", "delta", "byid"]
    pert_ids = {"byid": ["BRD-K11111111"], "delta": ["BRD-K99999999"]}
    mapping, rejected, toxicity, killed, provenance = _mod.prism_toxicity(prism, drugs, pert_ids)

    assert mapping["alpha"]["route"] == "exact name"
    assert mapping["byid"]["route"] == "broad identifier"
    assert "gamma" in rejected and "delta" in rejected      # too few lines; no entry at all
    assert toxicity[mapping["alpha"]["column"]] == pytest.approx(2.0)
    assert toxicity[mapping["beta"]["column"]] == pytest.approx(0.1)
    assert killed[mapping["alpha"]["column"]] == pytest.approx(1.0)
    assert killed[mapping["beta"]["column"]] == pytest.approx(0.0)
    assert mapping["alpha"]["dose"] == 2.5
    assert provenance["lfc_sha256"] and provenance["treatment_sha256"]
    assert provenance["departures_from_the_registered_description"]


def test_prism_mapping_refuses_a_name_that_is_not_in_the_release(tmp_path):
    columns = {"BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-1.5] * len(lines)}
    treatments = [{"column_name": "BRD-K11111111-001-01-9::2.5::HTS",
                   "broad_id": "BRD-K11111111-001-01-9", "name": "Alpha", "dose": 2.5,
                   "screen_id": "HTS"}]
    prism = _prism_files(tmp_path, treatments, columns)

    mapping, rejected, _, _, _ = _mod.prism_toxicity(prism, ["alpha-two"], {})
    assert mapping == {}
    assert "alpha-two" in rejected


def test_prism_mapping_prefers_the_treatment_measured_in_more_lines(tmp_path):
    columns = {
        "BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-1.0] * 110 + [np.nan] * 10,
        "BRD-K44444444-001-01-9::2.5::MTS004": lambda lines: [-4.0] * len(lines),
    }
    treatments = [
        {"column_name": "BRD-K11111111-001-01-9::2.5::HTS", "broad_id": "BRD-K11111111-001-01-9",
         "name": "Shared", "dose": 2.5, "screen_id": "HTS"},
        {"column_name": "BRD-K44444444-001-01-9::2.5::MTS004", "broad_id": "BRD-K44444444-001-01-9",
         "name": "Shared", "dose": 2.5, "screen_id": "MTS004"},
    ]
    prism = _prism_files(tmp_path, treatments, columns)

    mapping, _, toxicity, _, _ = _mod.prism_toxicity(prism, ["shared"], {})
    assert mapping["shared"]["n_lines"] == 120
    assert mapping["shared"]["n_candidates"] == 2
    assert toxicity[mapping["shared"]["column"]] == pytest.approx(4.0)


def test_hematopoietic_probe_says_when_there_are_too_few_drugs():
    arm = _arm(n_drugs=6, n_targets=3)
    reference = _reference(arm)
    out = _mod.hematopoietic_probe(arm, reference, {"LINE0"}, ["LINE0"])
    assert out["hematopoietic"]["reading"] == "too few drugs"
