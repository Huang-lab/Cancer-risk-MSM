"""Unit tests for phenotype roster building against real BioMe schemas.

Fixture layout (synthetic, TEST_### IDs):
    TEST_001  CRC case (C18.7, encounter + problem_list); no famhx.       WES: yes
    TEST_002  Breast case (C50.911);                       no famhx.       WES: yes
    TEST_003  no cancer dx;   Mother w/ Breast Cancer  (1deg breast)       WES: yes
    TEST_004  no cancer dx;   Father w/ Colon Cancer   (1deg colorectal)   WES: no
    TEST_005  CRC case (C20); Brother w/ Diabetes    (no cancer famhx)     WES: no
    TEST_006  Prostate case (C61); Brother w/ Prostate Cancer              WES: yes
    TEST_007  no cancer dx;   Mother w/ "Cancer" (unclassified any)        WES: yes

Cancer types (config): colorectal, breast, prostate, lung, ovarian, pancreatic.
"""
from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.phenotype.cases import build_cases, read_dx_ehr
from src.phenotype.famhx import (
    classify_cancer_type,
    famhx_flags,
    is_cancer_free_text,
    read_family_history,
    relationship_degree,
)
from src.phenotype.icd_mapping import IcdToPhecodeX, infer_icd_version
from src.phenotype.roster import assign_group, build
from src.phenotype.vcf_samples import load_wes_samples


FIXTURES = Path(__file__).parent / "synthetic"
CANCER_PHECODES = {
    "colorectal":  ["CA_101.1", "CA_101.11"],
    "breast":      ["CA_105.1"],
    "prostate":    ["CA_105.2"],
    "lung":        ["CA_101.6"],
    "ovarian":     ["CA_105.3"],
    "pancreatic":  ["CA_101.2"],
}
TYPE_PATTERNS = {
    "colorectal":  ["colon", "colorectal", "rectal"],
    "breast":      ["breast"],
    "prostate":    ["prostate"],
    "lung":        ["lung", "bronchus"],
    "ovarian":     ["ovary", "ovarian"],
    "pancreatic":  ["pancreas", "pancreatic"],
}
FREE_TEXT_PATTERNS = ["cancer", "malignant", "carcinoma", "lymphoma", "melanoma"]


# ---------- ICD mapping ----------

def test_infer_icd_version_9_vs_10():
    assert infer_icd_version("C18.7") == "10"
    assert infer_icd_version("Z00.00") == "10"
    assert infer_icd_version("153.0") == "9"
    assert infer_icd_version("174") == "9"


def test_icd_to_phecodex_expanded_map():
    m = IcdToPhecodeX.load(None)
    assert m.map_code("C18.7") == "CA_101.1"
    assert m.map_code("C20") == "CA_101.11"
    assert m.map_code("C50.911") == "CA_105.1"
    assert m.map_code("C61") == "CA_105.2"
    assert m.map_code("C34.11") == "CA_101.6"
    assert m.map_code("C56.9") == "CA_105.3"
    assert m.map_code("C25.0") == "CA_101.2"
    assert m.map_code("Z00.00") is None
    # ICD-9 fallback
    assert m.map_code("185", "9") == "CA_105.2"


# ---------- Relationship / famhx text ----------

def test_relationship_degree_common_variants():
    assert relationship_degree("Mother") == 1
    assert relationship_degree("Father") == 1
    assert relationship_degree("Brother") == 1
    assert relationship_degree("Maternal Grandfather") == 2
    assert relationship_degree("MGF") == 2
    assert relationship_degree("Uncle") == 2
    assert relationship_degree("Half Brother") == 2
    assert relationship_degree("Cousin") is None


def test_cancer_free_text_patterns():
    assert is_cancer_free_text("Breast cancer", FREE_TEXT_PATTERNS) is True
    assert is_cancer_free_text("Non-Hodgkin lymphoma", FREE_TEXT_PATTERNS) is True
    assert is_cancer_free_text("Hypertension", FREE_TEXT_PATTERNS) is False


def test_classify_cancer_type_truth_table():
    assert classify_cancer_type("Prostate Cancer", TYPE_PATTERNS) == "prostate"
    assert classify_cancer_type("Breast Cancer", TYPE_PATTERNS) == "breast"
    assert classify_cancer_type("Colon cancer", TYPE_PATTERNS) == "colorectal"
    assert classify_cancer_type("colorectal adenocarcinoma", TYPE_PATTERNS) == "colorectal"
    assert classify_cancer_type("Lung Cancer", TYPE_PATTERNS) == "lung"
    assert classify_cancer_type("Cancer", TYPE_PATTERNS) is None      # unclassified
    assert classify_cancer_type("Diabetes", TYPE_PATTERNS) is None
    assert classify_cancer_type("", TYPE_PATTERNS) is None


# ---------- Case building ----------

def _encounter_schema() -> dict:
    return {
        "person_id_col": "sample_name", "icd_col": "dx_code",
        "icd_version_col": "dx_code_type", "date_col": "encounter_date",
        "date_format": "%m/%d/%Y", "sep": "|",
    }


def test_read_dx_ehr_and_build_cases_from_fixtures():
    m = IcdToPhecodeX.load(None)
    rows = read_dx_ehr(FIXTURES / "encounter_diagnosis.tsv", _encounter_schema(),
                       m, source="encounter_diagnosis")
    # Cancer-mapping rows: TEST_001 (C18.7), TEST_002 (C50.911), TEST_005 (C20),
    # TEST_006 (C61). Non-cancer rows filtered out.
    pids = [r.person_id for r in rows]
    assert set(pids) == {"TEST_001", "TEST_002", "TEST_005", "TEST_006"}

    cases = build_cases(rows, CANCER_PHECODES,
                        enrollment_by_pid={f"TEST_{i:03d}": "2015-01-01"
                                           for i in range(1, 8)})
    by_key = {(c.person_id, c.cancer): c for c in cases}
    assert by_key[("TEST_001", "colorectal")].is_case
    assert by_key[("TEST_001", "colorectal")].first_dx_date == "2020-05-10"
    assert by_key[("TEST_001", "colorectal")].incident_flag is True
    assert by_key[("TEST_006", "prostate")].is_case
    assert ("TEST_003", "colorectal") not in by_key


# ---------- Family history parsing ----------

def test_read_family_history_pipe_delim_and_classification():
    schema = {
        "person_id_col": "sample_name", "relationship_col": "relation",
        "condition_col": "problem_description", "condition_code_col": None,
        "degree_col": None, "sep": "|", "date_format": "%m/%d/%Y",
    }
    m = IcdToPhecodeX.load(None)
    fam = read_family_history(
        FIXTURES / "family_history.tsv", schema,
        cancer_free_text_patterns=FREE_TEXT_PATTERNS,
        cancer_type_patterns=TYPE_PATTERNS,
        cancer_phecodes=CANCER_PHECODES, mapper=m,
        relative_degrees=[1, 2],
    )
    # TEST_003 -> mother w/ Breast Cancer
    flags3 = famhx_flags(fam["TEST_003"], CANCER_PHECODES)
    assert flags3["has_famhx_1deg_any_cancer"] is True
    assert flags3["has_famhx_1deg_breast"] is True
    assert flags3["has_famhx_1deg_colorectal"] is False
    # TEST_004 -> father w/ Colon Cancer
    flags4 = famhx_flags(fam["TEST_004"], CANCER_PHECODES)
    assert flags4["has_famhx_1deg_any_cancer"] is True
    assert flags4["has_famhx_1deg_colorectal"] is True
    # TEST_006 -> brother w/ Prostate Cancer
    flags6 = famhx_flags(fam["TEST_006"], CANCER_PHECODES)
    assert flags6["has_famhx_1deg_prostate"] is True
    # TEST_007 -> mother w/ generic "Cancer" (unclassified)
    flags7 = famhx_flags(fam["TEST_007"], CANCER_PHECODES)
    assert flags7["has_famhx_1deg_any_cancer"] is True
    assert flags7["has_famhx_1deg_breast"] is False
    assert flags7["has_famhx_1deg_colorectal"] is False
    # TEST_005 -> brother w/ diabetes (not a cancer entry)
    assert "TEST_005" not in fam


# ---------- Group assignment ----------

def test_assign_group_truth_table():
    assert assign_group(is_case=True,  has_famhx_any=True)  == "case"
    assert assign_group(is_case=False, has_famhx_any=True)  == "famhx"
    assert assign_group(is_case=False, has_famhx_any=False) == "control"


# ---------- WES sample loader ----------

def test_load_wes_samples_reads_ids():
    cfg = {"inputs": {"wes": {"sample_list_file": str(FIXTURES / "wes_samples.txt")}}}
    ws = load_wes_samples(cfg)
    assert ws == {"TEST_001", "TEST_002", "TEST_003", "TEST_006", "TEST_007"}

    cfg_missing = {"inputs": {"wes": {"sample_list_file": "/no/such/file"}}}
    assert load_wes_samples(cfg_missing) == set()

    cfg_null = {"inputs": {"wes": {"sample_list_file": None}}}
    assert load_wes_samples(cfg_null) == set()


# ---------- End-to-end + cancer_counts.tsv ----------

def _make_cfg(root: Path, wes_list: str | None) -> dict:
    return {
        "inputs": {
            "msm_data_root": str(root),
            "wes": {"sample_list_file": wes_list},
            "phenotypes": {
                "date": "synthetic",
                "files": {
                    "encounter_diagnosis": "encounter_diagnosis.tsv",
                    "problem_list": "problem_list.tsv",
                    "family_history": "family_history.tsv",
                    "demographics": "demographics.tsv",
                    "phecodex": "phecodex_missing.tsv",
                },
            },
        },
        "phenotype": {
            "icd_to_phecodex_map": None,
            "cancer_phecodes": CANCER_PHECODES,
            "ehr_schema": {
                "encounter_diagnosis": _encounter_schema(),
                "problem_list": {
                    "person_id_col": "sample_name", "icd_col": "icd9_code",
                    "icd_version_col": "code_type", "noted_date_col": "noted_date",
                    "resolved_col": "resolved_date", "status_col": "status",
                    "date_format": "%m/%d/%Y", "sep": "|",
                },
                "family_history": {
                    "person_id_col": "sample_name", "relationship_col": "relation",
                    "condition_col": "problem_description", "condition_code_col": None,
                    "degree_col": None, "date_format": "%m/%d/%Y", "sep": "|",
                },
                "demographics": {
                    "person_id_col": "sample_name", "birth_year_col": "birth_year",
                    "sex_col": "sex", "enrollment_date_col": "enrollment_date",
                    "date_format": "%m/%d/%Y", "sep": "|",
                },
                "phecodex_roster": {
                    "format": "auto", "person_id_col": "sample_name",
                    "phecodex_col": "phecodeX", "value_col": "case_control", "sep": "|",
                },
            },
            "roster": {
                "famhx_scope": "any_cancer", "control_exclusion": "any_cancer",
                "emit_prevalent_cases": True, "min_control_encounters": 0,
            },
            "famhx": {
                "relative_degrees": [1, 2],
                "cancer_free_text_patterns": FREE_TEXT_PATTERNS,
                "cancer_type_patterns": TYPE_PATTERNS,
            },
        },
    }


def _stage_fixtures(tmp_path: Path) -> Path:
    ph = tmp_path / "phenotypes" / "synthetic"
    ph.mkdir(parents=True)
    for f in ["encounter_diagnosis.tsv", "problem_list.tsv",
              "family_history.tsv", "demographics.tsv"]:
        shutil.copy(FIXTURES / f, ph / f)
    (ph / "phecodex_missing.tsv").write_text("sample_name\n")
    return ph


def _read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_end_to_end_roster_and_cancer_counts(tmp_path: Path):
    _stage_fixtures(tmp_path)
    cfg = _make_cfg(tmp_path, wes_list=str(FIXTURES / "wes_samples.txt"))
    out = tmp_path / "phenotype_out"
    build(cfg, out)

    # 1. Roster truth table (spot-check)
    roster = _read_tsv(out / "roster.tsv")
    by_key = {(r["person_id"], r["cancer"]): r for r in roster}
    assert by_key[("TEST_001", "colorectal")]["group"] == "case"
    assert by_key[("TEST_001", "colorectal")]["has_wes"] == "1"
    assert by_key[("TEST_006", "prostate")]["group"] == "case"
    assert by_key[("TEST_003", "colorectal")]["group"] == "famhx"    # has any-cancer famhx
    assert by_key[("TEST_003", "breast")]["group"] == "famhx"
    assert by_key[("TEST_004", "colorectal")]["group"] == "famhx"
    assert by_key[("TEST_007", "lung")]["group"] == "famhx"           # unclassified "Cancer" counts as any
    # People with a different cancer are excluded from other-cancer controls
    assert by_key[("TEST_001", "breast")]["group"] == "excluded"
    assert by_key[("TEST_006", "colorectal")]["group"] == "excluded"

    # 2. cancer_counts.tsv — the deliverable table
    counts = {r["cancer"]: r for r in _read_tsv(out / "cancer_counts.tsv")}

    # Colorectal: cases TEST_001 (WES) + TEST_005 -> 2 total, 1 with WES.
    assert counts["colorectal"]["n_ehr_cases"] == "2"
    assert counts["colorectal"]["n_ehr_cases_with_wes"] == "1"
    # 4 people with any 1deg cancer famhx (TEST_003/4/6/7); 3 in WES (TEST_003/6/7).
    assert counts["colorectal"]["n_famhx_1deg_any_cancer"] == "4"
    assert counts["colorectal"]["n_famhx_1deg_any_cancer_with_wes"] == "3"
    # Only TEST_004 has 1deg colorectal famhx (father-colon); not in WES.
    assert counts["colorectal"]["n_famhx_1deg_this_cancer"] == "1"
    assert counts["colorectal"]["n_famhx_1deg_this_cancer_with_wes"] == "0"

    # Breast: TEST_002 case (WES); TEST_003 breast-specific famhx (WES).
    assert counts["breast"]["n_ehr_cases"] == "1"
    assert counts["breast"]["n_ehr_cases_with_wes"] == "1"
    assert counts["breast"]["n_famhx_1deg_this_cancer"] == "1"
    assert counts["breast"]["n_famhx_1deg_this_cancer_with_wes"] == "1"

    # Prostate: TEST_006 case (WES) + brother-prostate famhx (WES).
    assert counts["prostate"]["n_ehr_cases"] == "1"
    assert counts["prostate"]["n_ehr_cases_with_wes"] == "1"
    assert counts["prostate"]["n_famhx_1deg_this_cancer"] == "1"
    assert counts["prostate"]["n_famhx_1deg_this_cancer_with_wes"] == "1"

    # Lung / ovarian / pancreatic: no cancer-specific cases or famhx in fixture.
    for cancer in ["lung", "ovarian", "pancreatic"]:
        assert counts[cancer]["n_ehr_cases"] == "0"
        assert counts[cancer]["n_famhx_1deg_this_cancer"] == "0"
        # But everyone with any 1deg famhx still counts as famhx_any:
        assert counts[cancer]["n_famhx_1deg_any_cancer"] == "4"


def test_end_to_end_without_wes_list_zeros_wes_columns(tmp_path: Path):
    _stage_fixtures(tmp_path)
    cfg = _make_cfg(tmp_path, wes_list=None)
    out = tmp_path / "phenotype_out"
    build(cfg, out)
    counts = {r["cancer"]: r for r in _read_tsv(out / "cancer_counts.tsv")}
    for cancer, row in counts.items():
        assert row["n_ehr_cases_with_wes"] == "0"
        assert row["n_famhx_1deg_any_cancer_with_wes"] == "0"
        assert row["n_famhx_1deg_this_cancer_with_wes"] == "0"
