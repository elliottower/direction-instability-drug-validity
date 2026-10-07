"""The real frozen bases, against the pin that lives outside them.

The synthetic tests in `test_analysis_basis.py` check the rule. These check the
artifact the registered run will actually consume: its identities, its counts, the
coordinates each construction excludes, the intersections the registered
comparisons run on, and that an edit to any of it is refused.
"""
import hashlib
import json
from pathlib import Path

import pytest

from geometry.references import MIN_LANDMARKS, N_LANDMARK, basis_sha256, validated_analysis_bases

REPO = Path(__file__).resolve().parent.parent
BASES = REPO / "registry/frozen/analysis_bases.json"
PIN = REPO / "registry/frozen/analysis_bases_pin.json"
AXIS = REPO / "registry/frozen/landmark_gene_ids.json"
GENERATOR = REPO / "experiments/03i_freeze_analysis_bases.py"

# what Amendment 4 states, restated here so a silent change to the artifact fails
EXPECTED = {
    "C1-K562": {"matched": 728, "basis": 728, "excluded": []},
    "C1-RPE1": {"matched": 813, "basis": 812, "excluded": ["CCL2"]},
    "C1-GW": {"matched": 721, "basis": 715,
              "excluded": ["BAMBI", "ICAM1", "MEST", "PXN", "SLC25A14", "TCTN1"]},
    "C1-GW-phenotype-positive": {"matched": 721, "basis": 715,
                                 "excluded": ["BAMBI", "ICAM1", "MEST", "PXN", "SLC25A14",
                                              "TCTN1"]},
}
EXPECTED_SHARED = {
    "R0.5: C0 against C1-K562": 728,
    "R0.5: C1-K562 against C1-RPE1": 688,
    "R0.5: C1-K562 against C1-GW": 712,
    "R6a: C1-K562 against C1-RPE1": 688,
}


@pytest.fixture
def gene_ids():
    return json.loads(AXIS.read_text())


@pytest.fixture
def pin():
    return json.loads(PIN.read_text())


def test_the_frozen_bases_pass_their_own_external_pin(gene_ids, pin):
    bases = validated_analysis_bases(BASES, gene_ids, pin, generator=GENERATOR)
    assert sorted(bases) == sorted(EXPECTED)


def test_each_construction_holds_the_counts_amendment_4_states(gene_ids, pin):
    bases = validated_analysis_bases(BASES, gene_ids, pin, generator=GENERATOR)
    for name, expected in EXPECTED.items():
        record = bases[name]
        assert record["n_landmarks_matched"] == expected["matched"], name
        assert record["n_landmarks_in_basis"] == expected["basis"], name
        assert len(record["landmark_gene_ids"]) == expected["basis"], name
        assert sorted(gene["symbol"] for gene in record["excluded_as_undefined"]) \
            == expected["excluded"], name
        assert record["n_landmarks_in_basis"] >= MIN_LANDMARKS, name


def test_the_excluded_coordinates_leave_only_the_arms_reading_their_source(gene_ids, pin):
    bases = validated_analysis_bases(BASES, gene_ids, pin, generator=GENERATOR)
    held = {name: set(record["landmark_positions"]) for name, record in bases.items()}
    excluded = {name: {gene["landmark_position"]: gene["symbol"]
                       for gene in record["excluded_as_undefined"]}
                for name, record in bases.items()}

    ccl2 = next(position for position, symbol in excluded["C1-RPE1"].items() if symbol == "CCL2")
    assert ccl2 not in held["C1-RPE1"]
    # the K562-essential release does not carry CCL2 among its matched landmarks,
    # so its absence there is the deposit and not this rule
    assert ccl2 not in held["C1-K562"] and ccl2 not in excluded["C1-K562"]

    for position, symbol in excluded["C1-GW"].items():
        assert position not in held["C1-GW"], symbol
        assert position not in held["C1-GW-phenotype-positive"], symbol
        assert position in held["C1-K562"] or position in held["C1-RPE1"], symbol
        assert position not in excluded["C1-K562"] and position not in excluded["C1-RPE1"]


def test_the_two_genome_wide_arms_share_one_basis(gene_ids, pin):
    bases = validated_analysis_bases(BASES, gene_ids, pin, generator=GENERATOR)
    left, right = bases["C1-GW"], bases["C1-GW-phenotype-positive"]

    assert left["file_sha256"] == right["file_sha256"]
    assert left["basis_sha256"] == right["basis_sha256"]
    assert left["landmark_positions"] == right["landmark_positions"]


def test_the_registered_comparisons_run_on_the_intersections_amendment_4_states():
    comparisons = json.loads(BASES.read_text())["registered_comparisons"]
    assert sorted(comparisons) == sorted(EXPECTED_SHARED)
    for key, n_shared in EXPECTED_SHARED.items():
        assert comparisons[key]["n_shared"] == n_shared, key
    assert comparisons["R0.5: C1-K562 against C1-GW"]["coordinates_lost_to_the_rule"] == 6
    assert comparisons["R6a: C1-K562 against C1-RPE1"]["coordinates_lost_to_the_rule"] == 0


def test_the_pooled_construction_declares_the_whole_axis_and_measures_728():
    document = json.loads(BASES.read_text())
    pooled = document["pooled_construction"]

    assert pooled["n_landmarks_declared"] == N_LANDMARK
    assert pooled["n_landmarks_measured"] == 728
    assert len(pooled["landmark_positions"]) == N_LANDMARK
    assert len(pooled["measured_positions"]) == 728
    assert pooled["declared_sha256"] != pooled["measured_sha256"]
    # declared and measured are different spaces, so the record names them apart
    assert "n_landmarks_matched" not in pooled


def test_the_coordinates_c0_measures_are_c1_k562s_basis():
    document = json.loads(BASES.read_text())
    pooled = document["pooled_construction"]
    c1 = document["constructions"]["C1-K562"]

    # both come from the K562-essential gene space, so R0.5's C0-against-C1
    # intersection is C1-K562's basis and not a third set of coordinates
    assert pooled["measured_sha256"] == c1["basis_sha256"]
    assert pooled["measured_gene_ids"] == c1["landmark_gene_ids"]
    assert pooled["measured_positions"] == c1["landmark_positions"]
    assert document["registered_comparisons"]["R0.5: C0 against C1-K562"]["n_shared"] == \
        c1["n_landmarks_in_basis"]


def test_the_pooled_source_is_pinned_from_outside_the_artifact(pin):
    pooled = json.loads(BASES.read_text())["pooled_construction"]
    pinned = pin["pooled_construction"]

    assert pooled["file"] == pinned["file"] == "ReplogleWeissman2022_K562_essential.h5ad"
    assert pooled["file_sha256"] == pinned["file_sha256"]
    assert len(pooled["file_sha256"]) == 64
    assert pooled["measured_sha256"] == pinned["measured_sha256"]


def test_r0_5_registers_three_pairs_and_none_can_be_dropped():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "discordance", REPO / "experiments/03d_h3_reference_discordance.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    comparisons = json.loads(BASES.read_text())["registered_comparisons"]

    # the registration names C0 and C1; C1 and C1-RPE1; C1 and C1-GW
    assert sorted(module.REGISTERED_R05_PAIRS) == ["C0_vs_C1", "C1_vs_GW", "C1_vs_RPE1"]
    for key, comparison in module.REGISTERED_R05_PAIRS.items():
        assert comparison in comparisons, key
    assert comparisons[module.REGISTERED_R05_PAIRS["C1_vs_GW"]]["n_shared"] == 712
    assert comparisons[module.REGISTERED_R05_PAIRS["C1_vs_RPE1"]]["n_shared"] == 688
    assert comparisons[module.REGISTERED_R05_PAIRS["C0_vs_C1"]]["n_shared"] == 728


def test_r0_5_takes_the_signed_cosine():
    source = (REPO / "experiments/03d_h3_reference_discordance.py").read_text()
    block = source[source.index("R0.5 agreement between constructions"):
                   source.index('result["modules"]["R0.5_agreement"]')]

    # the registration calls E = cos^2 signless where it wants signlessness and
    # asks for "the cosine" here, so an absolute value would read two
    # constructions that assign a target opposite directions as agreeing
    assert "abs(" not in block
    assert "left.directions[gene] @ right.directions[gene]" in block
    assert '"cosine_is_signed": True' in block
    assert '"n_negative"' in block


def test_the_axis_file_is_the_axis_the_bases_were_frozen_against(gene_ids, pin):
    document = json.loads(BASES.read_text())
    assert len(gene_ids) == N_LANDMARK
    assert basis_sha256(gene_ids) == document["gene_info"]["landmark_axis_sha256"]
    assert basis_sha256(gene_ids) == pin["landmark_axis_sha256"]
    assert hashlib.sha256(AXIS.read_bytes()).hexdigest() == pin["landmark_axis_file_sha256"]


def test_amendment_3s_pinned_identity_file_is_unchanged():
    frozen = REPO / "registry/frozen/expected_identities.json"
    assert hashlib.sha256(frozen.read_bytes()).hexdigest() == \
        "3de491d6234c9dc9f56072227b50e7a513fdc25720a9bc113ef56c0939e3b9b2"
    # Amendment 4's pins are a separate artifact, so freezing them leaves that one alone
    assert "analysis_bases" not in json.loads(frozen.read_text())


def test_an_edit_to_the_bases_is_refused_however_consistent_it_is(tmp_path, gene_ids, pin):
    document = json.loads(BASES.read_text())
    record = document["constructions"]["C1-RPE1"]
    dropped = record["landmark_gene_ids"][:-1]
    # the edit recomputes the artifact's own digest, which is why the pin is external
    record["landmark_gene_ids"] = dropped
    record["landmark_positions"] = record["landmark_positions"][:-1]
    record["landmark_symbols"] = record["landmark_symbols"][:-1]
    record["n_landmarks_in_basis"] = len(dropped)
    record["basis_sha256"] = basis_sha256(dropped)
    edited = tmp_path / "analysis_bases.json"
    edited.write_text(json.dumps(document, indent=2, sort_keys=True))

    with pytest.raises(AssertionError, match="has changed since it was pinned"):
        validated_analysis_bases(edited, gene_ids, pin, generator=GENERATOR)


def test_a_pin_naming_a_different_source_digest_is_refused(tmp_path, gene_ids, pin):
    tampered = json.loads(json.dumps(pin))
    tampered["constructions"]["C1-RPE1"]["file_sha256"] = "0" * 64
    tampered["file_sha256"] = hashlib.sha256(BASES.read_bytes()).hexdigest()

    with pytest.raises(AssertionError, match="file_sha256"):
        validated_analysis_bases(BASES, gene_ids, tampered, generator=GENERATOR)


def test_a_basis_frozen_against_another_axis_is_refused(gene_ids, pin):
    with pytest.raises(AssertionError, match="landmark axis"):
        validated_analysis_bases(BASES, list(reversed(gene_ids)), pin, generator=GENERATOR)


def test_a_changed_generator_is_refused(tmp_path, gene_ids, pin):
    other = tmp_path / "03i_freeze_analysis_bases.py"
    other.write_text(GENERATOR.read_text() + "\n# a change to the code that measures them\n")

    with pytest.raises(AssertionError, match="has changed since they were pinned"):
        validated_analysis_bases(BASES, gene_ids, pin, generator=other)


def test_no_analysis_code_states_a_contract_as_a_removable_assertion():
    import ast

    offending = {}
    for directory in ("geometry", "experiments"):
        for path in sorted((REPO / directory).rglob("*.py")):
            # a superseded script is the record of what ran, not code that will
            # run again, and editing it would falsify that record
            if "superseded" in path.parts:
                continue
            lines = [node.lineno for node in ast.walk(ast.parse(path.read_text()))
                     if isinstance(node, ast.Assert)]
            if lines:
                offending[str(path.relative_to(REPO))] = lines

    # python -O strips `assert`, so a contract written as one is not a contract.
    # Every refusal in the analysis code raises instead.
    assert offending == {}, offending
