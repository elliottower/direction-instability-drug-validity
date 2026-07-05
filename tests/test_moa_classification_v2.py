"""Tests for the broad MOA classification used in experiment 13f v2."""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "experiments"))
from importlib import import_module

mod = import_module("13f_moa_classification_auroc_v2")


class TestClassifyBroad:
    def test_hdac_is_machinery(self):
        assert mod.classify_broad("HDAC inhibitor") == "machinery"

    def test_parp_is_machinery(self):
        assert mod.classify_broad("PARP inhibitor") == "machinery"

    def test_bet_is_machinery(self):
        assert mod.classify_broad("BET inhibitor") == "machinery"

    def test_tubulin_polymerization_is_machinery(self):
        assert mod.classify_broad("tubulin polymerization inhibitor") == "machinery"

    def test_dopamine_receptor_is_receptor(self):
        assert mod.classify_broad("dopamine receptor antagonist") == "receptor"

    def test_calcium_channel_is_receptor(self):
        assert mod.classify_broad("calcium channel blocker") == "receptor"

    def test_adenosine_receptor_is_receptor(self):
        assert mod.classify_broad("adenosine receptor agonist") == "receptor"

    def test_egfr_is_kinase(self):
        assert mod.classify_broad("EGFR inhibitor") == "kinase"

    def test_aurora_kinase_is_kinase(self):
        assert mod.classify_broad("Aurora kinase inhibitor") == "kinase"

    def test_akt_is_kinase(self):
        assert mod.classify_broad("AKT inhibitor") == "kinase"

    def test_cox_is_excluded(self):
        assert mod.classify_broad("cyclooxygenase inhibitor") == "excluded"

    def test_pde_is_excluded(self):
        assert mod.classify_broad("phosphodiesterase inhibitor") == "excluded"

    def test_mao_is_excluded(self):
        assert mod.classify_broad("monoamine oxidase inhibitor") == "excluded"

    def test_glucocorticoid_is_excluded(self):
        assert mod.classify_broad("glucocorticoid receptor agonist") == "excluded"

    def test_estrogen_is_excluded(self):
        assert mod.classify_broad("estrogen receptor agonist") == "excluded"

    def test_ppar_is_excluded(self):
        assert mod.classify_broad("PPAR receptor agonist") == "excluded"

    def test_bacterial_is_excluded(self):
        assert mod.classify_broad("bacterial cell wall synthesis inhibitor") == "excluded"

    def test_ssri_is_excluded(self):
        assert mod.classify_broad("selective serotonin reuptake inhibitor (SSRI)") == "excluded"

    def test_nan_is_no_annotation(self):
        assert mod.classify_broad(float("nan")) == "no_annotation"

    def test_empty_is_no_annotation(self):
        assert mod.classify_broad("") == "no_annotation"

    def test_machinery_priority_over_kinase(self):
        assert mod.classify_broad("CDK inhibitor|cell cycle inhibitor") == "machinery"

    def test_machinery_priority_over_receptor(self):
        assert mod.classify_broad("proteasome inhibitor|dopamine receptor antagonist") == "machinery"

    def test_kinase_priority_over_receptor(self):
        assert mod.classify_broad("EGFR inhibitor|calcium channel blocker") == "kinase"

    def test_nhr_excluded_even_with_compound_moa(self):
        moa = "cannabinoid receptor agonist|PPAR receptor agonist"
        assert mod.classify_broad(moa) == "excluded"

    def test_nhr_excluded_before_machinery_check(self):
        moa = "estrogen receptor agonist|proteasome inhibitor"
        assert mod.classify_broad(moa) == "excluded"


class TestIsNuclearReceptor:
    def test_glucocorticoid(self):
        assert mod.is_nuclear_receptor("glucocorticoid receptor agonist")

    def test_ppar(self):
        assert mod.is_nuclear_receptor("PPAR receptor agonist")

    def test_dopamine_is_not(self):
        assert not mod.is_nuclear_receptor("dopamine receptor antagonist")

    def test_hdac_is_not(self):
        assert not mod.is_nuclear_receptor("HDAC inhibitor")

    def test_compound_with_nhr(self):
        assert mod.is_nuclear_receptor("cannabinoid receptor agonist|PPAR receptor agonist")


DRUG_REPO = Path(__file__).resolve().parent.parent.parent / "drug-perturbation-geometry"
INSTABILITY_CSV = DRUG_REPO / "zenodo_v1" / "drug_instability_8949.csv"


@pytest.mark.skipif(not INSTABILITY_CSV.exists(), reason="LINCS CSV not available")
class TestOnRealData:
    @pytest.fixture(scope="class")
    def classified_df(self):
        df = pd.read_csv(INSTABILITY_CSV)
        df["broad_class"] = df["moa"].apply(mod.classify_broad)
        df["is_nuclear_receptor"] = df["moa"].apply(mod.is_nuclear_receptor)
        return df

    def test_counts_sum(self, classified_df):
        annotated = classified_df[classified_df["broad_class"] != "no_annotation"]
        n_mach = (annotated["broad_class"] == "machinery").sum()
        n_kin = (annotated["broad_class"] == "kinase").sum()
        n_recv = (annotated["broad_class"] == "receptor").sum()
        n_excl = (annotated["broad_class"] == "excluded").sum()
        assert n_mach + n_kin + n_recv + n_excl == len(annotated)
        assert len(annotated) == classified_df["moa"].notna().sum()

    def test_expected_counts(self, classified_df):
        counts = classified_df["broad_class"].value_counts()
        assert counts["machinery"] == 180
        assert counts["receptor"] == 585
        assert counts["kinase"] == 153
        assert counts["excluded"] == 864

    def test_classes_disjoint(self, classified_df):
        annotated = classified_df[classified_df["broad_class"] != "no_annotation"]
        for _, row in annotated.iterrows():
            moa = row["moa"]
            cls = row["broad_class"]
            assert cls in {"machinery", "kinase", "receptor", "excluded"}, \
                f"Unexpected class {cls} for {moa}"

    def test_zero_nhr_in_binary_set(self, classified_df):
        binary = classified_df[classified_df["broad_class"].isin(["machinery", "receptor"])]
        nhr_in_binary = binary[binary["is_nuclear_receptor"]]
        assert len(nhr_in_binary) == 0, \
            f"NHR drugs in binary set: {nhr_in_binary['drug_name'].tolist()}"

    def test_all_nhr_in_excluded(self, classified_df):
        nhr = classified_df[classified_df["is_nuclear_receptor"]]
        assert (nhr["broad_class"] == "excluded").all(), \
            f"NHR drugs not all excluded: {nhr[nhr['broad_class'] != 'excluded'][['drug_name', 'broad_class']].to_dict()}"

    def test_nhr_total_count(self, classified_df):
        assert classified_df["is_nuclear_receptor"].sum() == 125

    def test_sensitivity_moves_all_nhr(self, classified_df):
        binary_df = classified_df[classified_df["broad_class"].isin(["machinery", "receptor"])].copy()
        annotated = classified_df[classified_df["broad_class"] != "no_annotation"]
        nhr_excluded = annotated[
            (annotated["broad_class"] == "excluded") & annotated["is_nuclear_receptor"]
        ].copy()
        assert len(nhr_excluded) == 125, \
            f"Sensitivity should move all 125 NHR drugs, got {len(nhr_excluded)}"
        nhr_as_mach = nhr_excluded.copy()
        nhr_as_mach["broad_class"] = "machinery"
        sens_df = pd.concat([binary_df, nhr_as_mach])
        assert (sens_df["broad_class"] == "machinery").sum() == 180 + 125
