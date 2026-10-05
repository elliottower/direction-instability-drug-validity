import json
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
    assert max(v for v in worst.values() if isinstance(v, float)) < 1e-12
    assert worst["crispri_fields_are_a_correction_not_a_reproduction"] is True

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
    tmp_path.mkdir(parents=True, exist_ok=True)
    lines = [f"ACH-{i:06d}" for i in range(120)]
    matrix = pd.DataFrame({column: values(lines) for column, values in columns_with_values.items()},
                          index=lines)
    matrix.to_csv(tmp_path / "primary-screen-replicate-collapsed-logfold-change.csv")
    pd.DataFrame(treatments).to_csv(
        tmp_path / "primary-screen-replicate-collapsed-treatment-info.csv", index=False)
    return tmp_path


def test_prism_mapping_uses_only_the_pinned_files_and_the_dose_window(tmp_path):
    columns = {
        "BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-2.0] * len(lines),
        "BRD-K22222222-001-01-9::2.5::HTS": lambda lines: [-0.1] * len(lines),
        # outside the 2.0-3.0 window fixed by Amendment 1
        "BRD-K55555555-001-01-9::0.08::HTS": lambda lines: [-5.0] * len(lines),
        # measured in too few cell lines
        "BRD-K33333333-001-01-9::2.5::HTS": lambda lines: [-3.0] * 20 + [np.nan] * (len(lines) - 20),
    }
    treatments = [
        {"column_name": "BRD-K11111111-001-01-9::2.5::HTS", "broad_id": "BRD-K11111111-001-01-9",
         "name": "Alpha", "dose": 2.5, "screen_id": "HTS"},
        {"column_name": "BRD-K22222222-001-01-9::2.5::HTS", "broad_id": "BRD-K22222222-001-01-9",
         "name": "Beta", "dose": 2.5, "screen_id": "HTS"},
        {"column_name": "BRD-K55555555-001-01-9::0.08::HTS", "broad_id": "BRD-K55555555-001-01-9",
         "name": "LowDose", "dose": 0.08, "screen_id": "HTS"},
        {"column_name": "BRD-K33333333-001-01-9::2.5::HTS", "broad_id": "BRD-K33333333-001-01-9",
         "name": "Gamma", "dose": 2.5, "screen_id": "HTS"},
    ]
    prism = _prism_files(tmp_path, treatments, columns)
    pert_ids = {"byid": ["BRD-K11111111"], "delta": ["BRD-K99999999"]}

    mapping = _mod.prism_mapping(prism, ["alpha", "beta", "gamma", "delta", "byid", "lowdose"],
                                 pert_ids)
    accepted, rejected = mapping["accepted"], mapping["rejected"]
    assert accepted["alpha"]["route"] == "exact name"
    assert accepted["byid"]["route"] == "broad identifier"
    assert accepted["alpha"]["dose"] == 2.5
    assert "lowdose" in rejected          # outside the frozen dose window
    assert "gamma" in rejected            # too few cell lines
    assert "delta" in rejected            # absent from the release
    assert mapping["dose_window"] == [2.0, 3.0]
    assert mapping["provenance"]["response_values_used_for_mapping"] is False
    assert mapping["provenance"]["response_missingness_used_for_mapping"] is True


def test_prism_response_refuses_a_mapping_it_was_not_given(tmp_path):
    columns = {"BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-2.0] * len(lines)}
    treatments = [{"column_name": "BRD-K11111111-001-01-9::2.5::HTS",
                   "broad_id": "BRD-K11111111-001-01-9", "name": "Alpha", "dose": 2.5,
                   "screen_id": "HTS"}]
    prism = _prism_files(tmp_path, treatments, columns)

    with pytest.raises(AssertionError, match="frozen before any response value"):
        _mod.prism_response(prism, tmp_path / "does_not_exist.json")

    mapping = _mod.prism_mapping(prism, ["alpha"], {})
    path = tmp_path / "r5_prism_mapping.json"
    path.write_text(json.dumps(mapping))
    toxicity, killed, provenance = _mod.prism_response(prism, path)
    assert toxicity[mapping["accepted"]["alpha"]["column"]] == pytest.approx(2.0)
    assert killed[mapping["accepted"]["alpha"]["column"]] == pytest.approx(1.0)
    assert provenance["mapping_sha256"]


def test_prism_response_refuses_a_response_matrix_that_moved(tmp_path):
    columns = {"BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-2.0] * len(lines)}
    treatments = [{"column_name": "BRD-K11111111-001-01-9::2.5::HTS",
                   "broad_id": "BRD-K11111111-001-01-9", "name": "Alpha", "dose": 2.5,
                   "screen_id": "HTS"}]
    prism = _prism_files(tmp_path, treatments, columns)
    mapping = _mod.prism_mapping(prism, ["alpha"], {})
    path = tmp_path / "r5_prism_mapping.json"
    path.write_text(json.dumps(mapping))

    _prism_files(tmp_path, treatments,
                 {"BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-1.0] * len(lines)})
    with pytest.raises(AssertionError, match="not the one the mapping was frozen against"):
        _mod.prism_response(prism, path)


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

    mapping = _mod.prism_mapping(prism, ["shared"], {})
    assert mapping["accepted"]["shared"]["n_lines"] == 120
    assert mapping["accepted"]["shared"]["n_candidates"] == 2


def test_hematopoietic_probe_says_when_there_are_too_few_drugs():
    arm = _arm(n_drugs=6, n_targets=3)
    reference = _reference(arm)
    out = _mod.hematopoietic_probe(arm, reference, {"LINE0"}, ["LINE0"])
    assert out["hematopoietic"]["reading"] == "too few drugs"


def _minimal_prism(tmp_path):
    columns = {"BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-2.0] * len(lines),
               "BRD-K22222222-001-01-9::2.5::HTS": lambda lines: [-0.4] * len(lines)}
    treatments = [
        {"column_name": "BRD-K11111111-001-01-9::2.5::HTS", "broad_id": "BRD-K11111111-001-01-9",
         "name": "Alpha", "dose": 2.5, "screen_id": "HTS"},
        {"column_name": "BRD-K22222222-001-01-9::2.5::HTS", "broad_id": "BRD-K22222222-001-01-9",
         "name": "Beta", "dose": 2.5, "screen_id": "HTS"}]
    return _prism_files(tmp_path, treatments, columns)


def test_the_mapping_is_written_canonically_and_its_hash_is_stable(tmp_path):
    prism = _minimal_prism(tmp_path)
    mapping = _mod.prism_mapping(prism, ["alpha", "beta"], {})

    first = _mod.write_mapping(mapping, tmp_path / "one.json")
    shuffled = {k: mapping[k] for k in reversed(list(mapping))}
    second = _mod.write_mapping(shuffled, tmp_path / "two.json")
    assert first == second
    assert (tmp_path / "one.json").read_text() == (tmp_path / "two.json").read_text()


def test_mapping_depends_on_missingness_and_not_on_the_response_values(tmp_path):
    treatments = [{"column_name": "BRD-K11111111-001-01-9::2.5::HTS",
                   "broad_id": "BRD-K11111111-001-01-9", "name": "Alpha", "dose": 2.5,
                   "screen_id": "HTS"}]
    quiet = _prism_files(tmp_path / "quiet", treatments,
                         {"BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-0.01] * len(lines)})
    loud = _prism_files(tmp_path / "loud", treatments,
                        {"BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-9.0] * len(lines)})

    a = _mod.prism_mapping(quiet, ["alpha"], {})
    b = _mod.prism_mapping(loud, ["alpha"], {})
    for mapping in (a, b):
        mapping["provenance"].pop("lfc_file")
        mapping["provenance"].pop("lfc_sha256")
        mapping["provenance"].pop("treatment_file")
        mapping["provenance"].pop("treatment_sha256")
    assert a == b
    assert a["provenance"]["response_values_used_for_mapping"] is False


def test_the_response_stage_refuses_a_mapping_whose_hash_is_not_the_expected_one(tmp_path):
    prism = _minimal_prism(tmp_path)
    mapping = _mod.prism_mapping(prism, ["alpha", "beta"], {})
    path = tmp_path / "r5_prism_mapping.json"
    expected = _mod.write_mapping(mapping, path)

    _mod.validate_mapping(json.loads(path.read_text()), ["alpha", "beta"], expected, path)

    edited = json.loads(path.read_text())
    edited["rejected"]["beta"] = "removed by hand"
    edited["accepted"].pop("beta")
    edited["n_accepted"], edited["n_rejected"] = 1, 1
    _mod.write_mapping(edited, path)
    with pytest.raises(AssertionError, match="not a frozen mapping"):
        _mod.validate_mapping(json.loads(path.read_text()), ["alpha", "beta"], expected, path)


def test_validate_mapping_refuses_a_mapping_built_for_another_cohort(tmp_path):
    prism = _minimal_prism(tmp_path)
    path = tmp_path / "r5_prism_mapping.json"
    expected = _mod.write_mapping(_mod.prism_mapping(prism, ["alpha"], {}), path)
    frozen = json.loads(path.read_text())

    with pytest.raises(AssertionError, match="built for a different cohort"):
        _mod.validate_mapping(frozen, ["alpha", "beta"], expected, path)


def test_validate_mapping_refuses_an_unregistered_route_or_a_dose_outside_the_window(tmp_path):
    prism = _minimal_prism(tmp_path)
    path = tmp_path / "r5_prism_mapping.json"
    _mod.write_mapping(_mod.prism_mapping(prism, ["alpha"], {}), path)

    frozen = json.loads(path.read_text())
    frozen["accepted"]["alpha"]["route"] = "fuzzy name"
    expected = _mod.write_mapping(frozen, path)
    with pytest.raises(AssertionError, match="unregistered route"):
        _mod.validate_mapping(json.loads(path.read_text()), ["alpha"], expected, path)

    frozen = json.loads(path.read_text())
    frozen["accepted"]["alpha"]["route"] = "exact name"
    frozen["accepted"]["alpha"]["dose"] = 0.05
    expected = _mod.write_mapping(frozen, path)
    with pytest.raises(AssertionError, match="dose outside the window"):
        _mod.validate_mapping(json.loads(path.read_text()), ["alpha"], expected, path)


def test_the_response_stage_reads_only_the_columns_the_mapping_names(tmp_path):
    prism = _minimal_prism(tmp_path)
    path = tmp_path / "r5_prism_mapping.json"
    mapping = _mod.prism_mapping(prism, ["alpha"], {})      # beta is never offered
    _mod.write_mapping(mapping, path)

    toxicity, killed, provenance = _mod.prism_response(prism, path)
    assert provenance["n_columns_read"] == 1
    assert set(toxicity) == {mapping["accepted"]["alpha"]["column"]}


def test_the_identifier_route_wins_when_both_routes_have_a_candidate(tmp_path):
    columns = {"BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-2.0] * len(lines),
               "BRD-K77777777-001-01-9::2.5::HTS": lambda lines: [-8.0] * len(lines)}
    treatments = [
        # the drug's own identifier, under a different name in the release
        {"column_name": "BRD-K11111111-001-01-9::2.5::HTS", "broad_id": "BRD-K11111111-001-01-9",
         "name": "ReleaseName", "dose": 2.5, "screen_id": "HTS"},
        # a different compound that happens to carry the drug's name
        {"column_name": "BRD-K77777777-001-01-9::2.5::HTS", "broad_id": "BRD-K77777777-001-01-9",
         "name": "Alpha", "dose": 2.5, "screen_id": "HTS"}]
    prism = _prism_files(tmp_path, treatments, columns)

    mapping = _mod.prism_mapping(prism, ["alpha"], {"alpha": ["BRD-K11111111-001-01-9"]})
    entry = mapping["accepted"]["alpha"]
    assert entry["route"] == "broad identifier"
    assert entry["stem"] == "BRD-K11111111"       # not the name match, which had more effect


def test_the_frozen_regex_matches_stems_and_nothing_broader(tmp_path):
    columns = {"BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-2.0] * len(lines)}
    treatments = [{"column_name": "BRD-K11111111-001-01-9::2.5::HTS",
                   "broad_id": "BRD-K11111111-001-01-9", "name": "Alpha", "dose": 2.5,
                   "screen_id": "HTS"}]
    prism = _prism_files(tmp_path, treatments, columns)

    # a suffixed identifier still yields the stem; a malformed one yields nothing
    by_suffix = _mod.prism_mapping(prism, ["x"], {"x": ["BRD-K11111111-999-99-9"]})
    assert by_suffix["accepted"]["x"]["route"] == "broad identifier"
    malformed = _mod.prism_mapping(prism, ["y"], {"y": ["BRD-K1111111", "BRDK11111111"]})
    assert "y" in malformed["rejected"]


def test_validate_mapping_refuses_two_drugs_claiming_one_treatment(tmp_path):
    prism = _minimal_prism(tmp_path)
    path = tmp_path / "r5_prism_mapping.json"
    _mod.write_mapping(_mod.prism_mapping(prism, ["alpha", "beta"], {}), path)

    frozen = json.loads(path.read_text())
    frozen["accepted"]["beta"]["column"] = frozen["accepted"]["alpha"]["column"]
    expected = _mod.write_mapping(frozen, path)
    with pytest.raises(AssertionError, match="same treatment column"):
        _mod.validate_mapping(json.loads(path.read_text()), ["alpha", "beta"], expected, path)


def test_validate_mapping_refuses_a_partition_that_does_not_account_for_every_drug(tmp_path):
    prism = _minimal_prism(tmp_path)
    path = tmp_path / "r5_prism_mapping.json"
    _mod.write_mapping(_mod.prism_mapping(prism, ["alpha", "beta"], {}), path)

    frozen = json.loads(path.read_text())
    frozen["accepted"].pop("beta")
    frozen["n_accepted"] = 1
    expected = _mod.write_mapping(frozen, path)
    with pytest.raises(AssertionError, match="account for every drug"):
        _mod.validate_mapping(json.loads(path.read_text()), ["alpha", "beta"], expected, path)


def test_an_exact_tie_goes_to_the_lexically_first_column(tmp_path):
    lines = [f"ACH-{i:06d}" for i in range(120)]
    columns = {"BRD-K99999999-001-01-9::2.5::HTS": lambda lines: [-1.0] * len(lines),
               "BRD-K11111111-001-01-9::2.5::HTS": lambda lines: [-7.0] * len(lines)}
    treatments = [
        {"column_name": "BRD-K99999999-001-01-9::2.5::HTS", "broad_id": "BRD-K99999999-001-01-9",
         "name": "Tied", "dose": 2.5, "screen_id": "HTS"},
        {"column_name": "BRD-K11111111-001-01-9::2.5::HTS", "broad_id": "BRD-K11111111-001-01-9",
         "name": "Tied", "dose": 2.5, "screen_id": "HTS"}]
    prism = _prism_files(tmp_path, treatments, columns)

    mapping = _mod.prism_mapping(prism, ["tied"], {})
    entry = mapping["accepted"]["tied"]
    assert entry["n_candidates"] == 2
    assert entry["column"] == "BRD-K11111111-001-01-9::2.5::HTS"   # lexically first, both at 120



def test_the_discordance_script_holds_no_assert_statements():
    # `assert` disappears under -O. The gate was converted first and this file was
    # left, which reads as finished when it is not.
    import ast

    source = (Path(__file__).resolve().parents[1] / "experiments"
              / "03d_h3_reference_discordance.py").read_text()
    remaining = [node.lineno for node in ast.walk(ast.parse(source))
                 if isinstance(node, ast.Assert)]
    assert remaining == [], f"assert statements at lines {remaining}"


def test_the_rebuild_loader_refuses_what_it_used_to_tolerate(tmp_path):
    import json

    import numpy as np

    rng = np.random.default_rng()
    genes = [str(i) for i in range(_mod.N_LANDMARK)]
    (tmp_path / "landmark_gene_ids.json").write_text(json.dumps(genes))

    with pytest.raises(AssertionError, match="no shards"):
        _mod.build_drug_signatures_from_rebuild(tmp_path)

    def shard(path, matrices, cells):
        np.savez_compressed(path, fingerprint=np.array("fp"),
                            __cells__=np.array(json.dumps(cells)), **matrices)

    # a matrix whose rows contradict its declared cell lines: zip() would truncate
    shard(tmp_path / "shard_000.npz",
          {"drugA": rng.standard_normal((1, _mod.N_LANDMARK))},
          {"drugA": ["MCF7", "PC3"]})
    with pytest.raises(AssertionError, match=r"the matrix is \(1, 978\) against 2"):
        _mod.build_drug_signatures_from_rebuild(tmp_path)

    # the same drug from two shards
    shard(tmp_path / "shard_000.npz", {"drugA": rng.standard_normal((2, _mod.N_LANDMARK))},
          {"drugA": ["MCF7", "PC3"]})
    shard(tmp_path / "shard_001.npz", {"drugA": rng.standard_normal((2, _mod.N_LANDMARK))},
          {"drugA": ["MCF7", "PC3"]})
    with pytest.raises(AssertionError, match="more than one shard"):
        _mod.build_drug_signatures_from_rebuild(tmp_path)

    # a cell line declared twice
    (tmp_path / "shard_001.npz").unlink()
    shard(tmp_path / "shard_000.npz", {"drugA": rng.standard_normal((2, _mod.N_LANDMARK))},
          {"drugA": ["MCF7", "MCF7"]})
    with pytest.raises(AssertionError, match="declared twice"):
        _mod.build_drug_signatures_from_rebuild(tmp_path)

    # consistent input loads
    shard(tmp_path / "shard_000.npz", {"drugA": rng.standard_normal((2, _mod.N_LANDMARK))},
          {"drugA": ["MCF7", "PC3"]})
    per_drug, axis = _mod.build_drug_signatures_from_rebuild(tmp_path)
    assert sorted(per_drug["drugA"]) == ["MCF7", "PC3"]
    assert axis == genes



def test_the_gate_measures_the_crispri_correction_rather_than_requiring_it_to_vanish():
    # Amendment 3: the deposited C0 fields placed the reference by the extraction's
    # declared labels while the signatures followed another axis, so the corrected
    # run is deliberately different from them and must not be asserted equal
    arm = _arm()
    reference = _reference(arm)
    values = _mod.quantities(arm, reference)
    records = [{"drug": arm.drugs[i], "target": arm.targets[i],
                "n_celllines": len(arm.cell_lines[i]),
                "raw_instability": values["D"][i],
                "proj_crispri": values["P"][i] + 0.5,
                "enrich_crispri": values["E"][i] + 0.5}
               for i in range(len(arm))]

    worst = _mod.reproduction_gate(arm, None, values, records)
    assert worst["P_crispri"] == pytest.approx(0.5, abs=1e-9)
    assert worst["D"] < 1e-12
    assert worst["crispri_fields_are_a_correction_not_a_reproduction"] is True

    # the legacy rerun of the deposited configuration still demands equality
    with pytest.raises(AssertionError, match="P_crispri reproduces"):
        _mod.reproduction_gate(arm, None, values, records, crispri_reproduces=True)


def test_the_invariant_quantities_are_still_required_to_reproduce():
    arm = _arm()
    reference = _reference(arm)
    values = _mod.quantities(arm, reference)
    records = [{"drug": arm.drugs[i], "target": arm.targets[i],
                "n_celllines": len(arm.cell_lines[i]),
                "raw_instability": values["D"][i] + (0.5 if i == 2 else 0.0),
                "proj_shrna": values["P"][i], "enrich_shrna": values["E"][i]}
               for i in range(len(arm))]

    with pytest.raises(AssertionError, match="invariant under the shared permutation"):
        _mod.reproduction_gate(arm, values, None, records)
