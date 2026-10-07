"""Amendment 5's two-layer reproduction gate.

The frozen flat absolute tolerance stays on the comparison it was registered for,
the retired route against the deposited records. The elementwise float32 rule
governs the comparison Amendment 3 created, the authoritative route against those
same records. These check that each rule is applied where it belongs, that neither
can be satisfied in aggregate, and that the gate still refuses a wrong value.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
MEASUREMENT = REPO / "results/03d_h3_reference_discordance/two_layer_gate_measurement.json"


def load():
    spec = importlib.util.spec_from_file_location(
        "discordance", REPO / "experiments/03d_h3_reference_discordance.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def mod():
    return load()


def an_arm(mod, drugs, targets, n_contexts=3, n_genes=12):
    signatures, cells = [], []
    for index, _ in enumerate(drugs):
        base = np.arange(1, n_contexts * n_genes + 1, dtype=np.float64)
        signatures.append(base.reshape(n_contexts, n_genes) + index)
        cells.append([f"CELL{j}" for j in range(n_contexts)])
    return mod.Arm(name="test", drugs=np.array(drugs), targets=np.array(targets),
                   signatures=signatures, cell_lines=cells)


def a_reference(mod, targets, n_genes=12):
    directions = {t: mod.unit(np.arange(1.0, n_genes + 1) + i)
                  for i, t in enumerate(sorted(set(targets)))}
    return mod.Reference("shRNA", directions, np.arange(n_genes))


def records_for(mod, arm, reference, nudge=None):
    values = mod.quantities(arm, reference)
    out = []
    for i, drug in enumerate(arm.drugs):
        record = {"drug": str(drug), "target": str(arm.targets[i]),
                  "n_celllines": len(arm.cell_lines[i]),
                  "raw_instability": float(values["D"][i]),
                  "proj_shrna": float(values["P"][i]),
                  "enrich_shrna": float(values["E"][i])}
        if nudge and str(drug) == nudge[0]:
            record[nudge[1]] = record[nudge[1]] + nudge[2]
        out.append(record)
    return out, values


PROVENANCE = {"deposited_records_sha256": "d" * 64,
              "retired_extraction_sha256": "e" * 64,
              "analysis_code_sha256": "c" * 64}


def a_legacy_audit(tmp_path, mod, holds=True, n_drugs=795, provenance=None):
    """A sealed audit and the external pin that authenticates it.

    The pin is a separate object on purpose: the gate must not be able to read the
    expected identity out of the artifact it is checking.
    """
    table = tmp_path / "legacy_reproduction_per_drug.json"
    table.write_text(json.dumps({"rows": [], "n_rows": n_drugs * 3}))
    table_sha = mod.sha256_file(table)
    provenance = dict(PROVENANCE if provenance is None else provenance)
    layer = {quantity: {"max_absolute_error": 9.7e-07 if holds else 2e-06,
                        "max_normalized_residual": 0.122,
                        "flat_1e-6_holds": holds,
                        "n_over_flat_1e-6": 0 if holds else 3,
                        "drugs_over_flat_1e-6": [] if holds else ["a", "b", "c"]}
             for quantity in ("P_shrna", "E_shrna", "D")}
    audit = {"schema": "legacy_reproduction_audit/1",
             "flat_rule": "max |a-b| < 1e-6",
             "flat_tolerance": 1e-6,
             "quantities": ["P_shrna", "E_shrna", "D"],
             "n_drugs": n_drugs,
             "layers": {"legacy_reproduction": layer},
             "per_drug_table": {"file": table.name, "sha256": table_sha,
                                "n_rows": n_drugs * 3},
             "provenance": provenance}
    path = tmp_path / "legacy_reproduction_audit.json"
    path.write_text(json.dumps(audit, indent=2))
    pin = {"audit_sha256": mod.sha256_file(path),
           "flat_rule": "max |a-b| < 1e-6",
           "n_drugs": n_drugs,
           "per_drug_sha256": table_sha,
           "provenance": dict(provenance)}
    return path, pin


def rewritten(path, mod, pin, **changes):
    """Change the audit's content and re-pin it, so a refusal is not the hash check."""
    audit = json.loads(path.read_text())
    audit.update(changes)
    path.write_text(json.dumps(audit, indent=2))
    return {**pin, "audit_sha256": mod.sha256_file(path)}


def test_the_gate_passes_when_both_layers_hold(tmp_path, mod):
    drugs, targets = ["a", "b", "c"], ["T1", "T2", "T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)

    report = mod.reproduction_gate(arm, values, None, records,
                                   legacy_record=legacy, legacy_expected=pin)
    assert report["n_failures"] == 0
    assert report["legacy_reproduction"]["flat_tolerance_holds_on_every_quantity"]
    assert sorted(report["production_equivalence"]) == ["D", "E_shrna", "P_shrna"]
    assert len(report["per_drug"]) == 9


def test_a_difference_inside_the_elementwise_allowance_passes(tmp_path, mod):
    drugs, targets = ["a", "b"], ["T1", "T2"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    deposited = records[0]["proj_shrna"]
    # half the allowance at this magnitude, which a flat 1e-6 would reject when
    # the magnitude is large enough
    records[0]["proj_shrna"] = deposited + 0.5 * mod.elementwise_allowance(deposited)

    report = mod.reproduction_gate(arm, values, None, records,
                                   legacy_record=legacy, legacy_expected=pin)
    assert report["n_failures"] == 0
    assert report["production_equivalence"]["P_shrna"]["max_normalized_residual"] \
        == pytest.approx(0.5, rel=1e-6)


def test_a_difference_over_the_elementwise_allowance_refuses(tmp_path, mod):
    drugs, targets = ["a", "b"], ["T1", "T2"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    deposited = records[1]["proj_shrna"]
    records[1]["proj_shrna"] = deposited + 1.5 * mod.elementwise_allowance(deposited)

    with pytest.raises(mod.DiscordanceError, match="production equivalence fails"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin)


def test_one_drug_cannot_be_excused_by_the_others(tmp_path, mod):
    drugs = [f"d{i}" for i in range(40)]
    targets = ["T1"] * 40
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    deposited = records[7]["proj_shrna"]
    records[7]["proj_shrna"] = deposited + 2.0 * mod.elementwise_allowance(deposited)

    # 39 of 40 drugs reproduce exactly, and the rule is elementwise, so the mean
    # residual being far below 1 does not carry the one that is over
    with pytest.raises(mod.DiscordanceError, match="d7 P_shrna"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin)


def test_a_failing_legacy_layer_refuses_before_the_production_layer(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod, holds=False)

    with pytest.raises(mod.DiscordanceError, match="is a finding and is not amended away"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin)


def test_the_gate_refuses_without_the_legacy_record(mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)

    with pytest.raises(mod.DiscordanceError, match="sealed audit"):
        mod.reproduction_gate(arm, values, None, records, legacy_record=None)


def test_a_record_missing_a_shrna_field_is_a_cohort_defect(tmp_path, mod):
    drugs, targets = ["a", "b"], ["T1", "T2"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    del records[1]["enrich_shrna"]

    # the old gate checked these fields only where present, so a record missing
    # one was silently not compared
    with pytest.raises(mod.DiscordanceError, match="carries no enrich_shrna"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin)


def test_a_non_finite_value_never_reaches_the_comparison(tmp_path, mod):
    drugs, targets = ["a", "b"], ["T1", "T2"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    records[0]["proj_shrna"] = float("nan")

    with pytest.raises(mod.DiscordanceError, match="deposited P_shrna is not finite"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin)

    records[0]["proj_shrna"] = 1.0
    broken = {**values, "P": np.array([np.inf] + list(values["P"][1:]))}
    with pytest.raises(mod.DiscordanceError, match="computed P_shrna is not finite"):
        mod.reproduction_gate(arm, broken, None, records,
                              legacy_record=legacy, legacy_expected=pin)


def test_every_drug_and_quantity_is_reported_with_its_allowance(tmp_path, mod):
    drugs, targets = ["a", "b", "c"], ["T1", "T2", "T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)

    report = mod.reproduction_gate(arm, values, None, records,
                                   legacy_record=legacy, legacy_expected=pin)
    for row in report["per_drug"]:
        assert set(row) >= {"drug", "quantity", "computed", "deposited", "absolute_error",
                            "allowed_elementwise", "normalized_residual", "within_flat",
                            "within_elementwise", "enforced"}
        assert row["allowed_elementwise"] == pytest.approx(
            mod.RECON_ATOL + mod.RECON_RTOL * abs(row["deposited"]))
    for quantity in ("P_shrna", "E_shrna", "D"):
        summary = report["production_equivalence"][quantity]
        assert summary["n_drugs"] == 3
        assert "worst_drug" in summary and "drugs_over_flat" in summary


def test_the_measured_legacy_layer_holds_on_the_real_cohort():
    layers = json.loads(MEASUREMENT.read_text())["layers"]
    legacy, production = layers["legacy_reproduction"], layers["production_equivalence"]

    assert legacy["n_drugs"] == 795 and production["n_drugs"] == 795
    for quantity in ("P_shrna", "E_shrna", "D"):
        # the retired route reproduces the deposited analysis under the frozen
        # flat tolerance, so that tolerance is kept rather than loosened
        assert legacy[quantity]["flat_1e-6_holds"], quantity
        assert legacy[quantity]["n_over_flat_1e-6"] == 0, quantity
        assert production[quantity]["elementwise_float32_holds"], quantity
        assert production[quantity]["max_normalized_residual"] < 1.0, quantity
    # and the production route is the one the flat rule rejects, for one drug
    assert production["P_shrna"]["n_over_flat_1e-6"] == 1
    assert production["P_shrna"]["drugs_over_flat_1e-6"] == ["levofloxacin"]


# ---- the legacy audit is a trust boundary, so it is authenticated, not believed


def test_a_legacy_audit_whose_hash_is_not_the_pinned_one_refuses(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)

    with pytest.raises(mod.DiscordanceError, match="Amendment 5 pins"):
        mod.reproduction_gate(arm, values, None, records, legacy_record=legacy,
                              legacy_expected={**pin, "audit_sha256": "0" * 64})


def test_one_altered_byte_in_the_legacy_audit_refuses(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    legacy.write_text(legacy.read_text().replace("legacy_reproduction_audit/1",
                                                 "legacy_reproduction_audit/l"))

    with pytest.raises(mod.DiscordanceError, match="Amendment 5 pins"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin)


def test_a_fabricated_audit_that_says_everything_holds_refuses(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    # a replacement that would satisfy every predicate the gate reads, written by
    # someone who never ran the retired route
    forged = tmp_path / "forged.json"
    forged.write_text(legacy.read_text())

    with pytest.raises(mod.DiscordanceError, match="Amendment 5 pins"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=forged,
                              legacy_expected={**pin, "audit_sha256": "f" * 64})


def test_the_gate_refuses_when_no_pin_is_supplied(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, _ = a_legacy_audit(tmp_path, mod)

    with pytest.raises(mod.DiscordanceError, match="attests to itself"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=None)


def test_a_legacy_audit_covering_the_wrong_number_of_drugs_refuses(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    # re-pinned, so the refusal is the count and not the hash
    pin = rewritten(legacy, mod, pin, n_drugs=794)

    with pytest.raises(mod.DiscordanceError, match="covers 794 drugs"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy,
                              legacy_expected={**pin, "n_drugs": 795})


def test_a_legacy_audit_that_contradicts_itself_refuses(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)

    # the boolean, the count and the list are three statements of one fact
    for field, value, message in (
            ("n_over_flat_1e-6", 2, "while counting 2 drugs over"),
            ("drugs_over_flat_1e-6", ["levofloxacin"], "while naming"),
            ("max_absolute_error", 4e-06, "is not below")):
        directory = tmp_path / field
        directory.mkdir()
        legacy, pin = a_legacy_audit(directory, mod)
        audit = json.loads(legacy.read_text())
        audit["layers"]["legacy_reproduction"]["P_shrna"][field] = value
        legacy.write_text(json.dumps(audit, indent=2))
        pin = {**pin, "audit_sha256": mod.sha256_file(legacy)}

        with pytest.raises(mod.DiscordanceError, match=message):
            mod.reproduction_gate(arm, values, None, records,
                                  legacy_record=legacy, legacy_expected=pin)


def test_a_legacy_audit_measured_on_other_inputs_refuses(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    # the audit holds, and it was measured against a different deposited cohort
    pin = {**pin, "provenance": {**pin["provenance"],
                                 "deposited_records_sha256": "9" * 64}}

    with pytest.raises(mod.DiscordanceError, match="not evidence about this cohort"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin)


def test_a_legacy_audit_without_its_complete_table_refuses(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)

    with pytest.raises(mod.DiscordanceError, match="not accepted without the complete"):
        mod.reproduction_gate(arm, values, None, records, legacy_record=legacy,
                              legacy_expected={**pin, "per_drug_sha256": "a" * 64})


def test_a_legacy_audit_declaring_its_own_rule_refuses(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    pin = rewritten(legacy, mod, pin, flat_rule="max |a-b| < 1e-4")

    with pytest.raises(mod.DiscordanceError, match="restate its own rule"):
        mod.reproduction_gate(arm, values, None, records, legacy_record=legacy,
                              legacy_expected={**pin, "flat_rule": "max |a-b| < 1e-6"})


# ---- cohort identity: a record nobody compared is not a reproduction


def test_a_duplicated_deposited_drug_refuses(tmp_path, mod):
    drugs, targets = ["a", "b"], ["T1", "T2"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    # keying by drug would keep the last of these two and compare against it
    records.append({**records[0], "proj_shrna": records[0]["proj_shrna"] + 1.0})

    with pytest.raises(mod.DiscordanceError, match="duplicated drug identifiers"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin)


def test_a_deposited_record_the_arm_never_compares_refuses(tmp_path, mod):
    drugs, targets = ["a", "b"], ["T1", "T2"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)

    with pytest.raises(mod.DiscordanceError, match="would go uncompared"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin,
                              expected_drugs={"a", "b", "c"})


def test_a_drug_with_no_deposited_record_refuses(tmp_path, mod):
    drugs, targets = ["a", "b"], ["T1", "T2"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    records = [r for r in records if r["drug"] != "b"]

    # previously an incidental KeyError, which is not a stated refusal
    with pytest.raises(mod.DiscordanceError, match="no deposited record"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin)


def test_a_crispri_field_without_its_pair_refuses(tmp_path, mod):
    drugs, targets = ["a", "b"], ["T1", "T2"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    for record in records:
        record["proj_crispri"] = record["proj_shrna"]
        record["enrich_crispri"] = record["enrich_shrna"]
    # enrichment serialized without projection was silently skipped, and the
    # opposite case raised an incidental key error
    del records[1]["proj_crispri"]

    with pytest.raises(mod.DiscordanceError, match="carries no proj_crispri"):
        mod.reproduction_gate(arm, None, values, records,
                              legacy_record=legacy, legacy_expected=pin)

    records[1]["proj_crispri"] = records[1]["proj_shrna"]
    del records[1]["enrich_crispri"]
    with pytest.raises(mod.DiscordanceError, match="carries no enrich_crispri"):
        mod.reproduction_gate(arm, None, values, records,
                              legacy_record=legacy, legacy_expected=pin)


# ---- the flat rule is strict, at the boundary and in the audit


def test_the_flat_rule_is_strict_at_exactly_the_tolerance(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    deposited = records[0]["proj_shrna"]
    records[0]["proj_shrna"] = deposited + mod.RECON_TOL

    report = mod.reproduction_gate(arm, values, None, records,
                                   legacy_record=legacy, legacy_expected=pin)
    row = next(r for r in report["per_drug"] if r["quantity"] == "P_shrna")
    assert row["absolute_error"] == pytest.approx(mod.RECON_TOL, rel=1e-9)
    assert not row["within_flat"]


def test_a_legacy_audit_at_exactly_the_tolerance_refuses(tmp_path, mod):
    drugs, targets = ["a"], ["T1"]
    arm = an_arm(mod, drugs, targets)
    reference = a_reference(mod, targets)
    records, values = records_for(mod, arm, reference)
    legacy, pin = a_legacy_audit(tmp_path, mod)
    audit = json.loads(legacy.read_text())
    audit["layers"]["legacy_reproduction"]["D"]["max_absolute_error"] = mod.RECON_TOL
    legacy.write_text(json.dumps(audit, indent=2))
    pin = {**pin, "audit_sha256": mod.sha256_file(legacy)}

    with pytest.raises(mod.DiscordanceError, match="is not below"):
        mod.reproduction_gate(arm, values, None, records,
                              legacy_record=legacy, legacy_expected=pin)
