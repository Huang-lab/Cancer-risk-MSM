"""Unit tests for phenotype roster building.

Truth table (from synthetic fixtures):

    TEST_001: CRC case (C18.7); no famhx     -> CRC:case, breast:control
    TEST_002: Breast case (C50.911); no famhx -> CRC:control, breast:case
    TEST_003: no cancer dx; mother w/ breast cancer -> CRC:famhx, breast:famhx
    TEST_004: no cancer dx; maternal grandmother w/ colon cancer -> CRC:famhx, breast:famhx
    TEST_005: CRC case (C20); brother w/ diabetes -> CRC:case, breast:excluded
              (excluded because control_exclusion=any_cancer)
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.phenotype.cases import (
    build_cases,
    read_dx_ehr,
)
from src.phenotype.famhx import (
    famhx_flags,
    is_cancer_free_text,
    read_family_history,
    relationship_degree,
)
from src.phenotype.icd_mapping import IcdToPhecodeX, infer_icd_version
from src.phenotype.roster import assign_group, build


FIXTURES = Path(__file__).parent / "synthetic"
CANCER_PHECODES = {"colorectal": ["CA_101.1", "CA_101.11"], "breast": ["CA_105.1"]}


# ---------- ICD mapping ----------

def test_infer_icd_version_9_vs_10():
    assert infer_icd_version("C18.7") == "10"
    assert infer_icd_version("Z00.00") == "10"
    assert infer_icd_version("153.0") == "9"
    assert infer_icd_version("174") == "9"


def test_icd_to_phecodex_builtin_crc_and_breast():
    m = IcdToPhecodeX.load(None)
    assert m.map_code("C18.7") == "CA_101.1"
    assert m.map_code("C19") == "CA_101.11"
    assert m.map_code("C50.911") == "CA_105.1"
    assert m.map_code("D05.11") == "CA_105.1"
    assert m.map_code("Z00.00") is None
    assert m.map_code("C18.7", "10") == "CA_101.1"
    assert m.map_code("174", "9") == "CA_105.1"


# ---------- Relationship / famhx text detection ----------

def test_relationship_degree_common_variants():
    assert relationship_degree("MOTHER") == 1
    assert relationship_degree("Father") == 1
    assert relationship_degree("sibling") == 1
    assert relationship_degree("MATERNAL_GRANDMOTHER") == 2
    assert relationship_degree("maternal grandmother") == 2
    assert relationship_degree("MGF") == 2
    assert relationship_degree("uncle") == 2
    assert relationship_degree("HALF_BROTHER") == 2
    assert relationship_degree("cousin") is None
    assert relationship_degree("") is None


def test_cancer_free_text_patterns():
    pats = ["cancer", "malignant", "carcinoma", "lymphoma", "melanoma"]
    assert is_cancer_free_text("Breast cancer", pats) is True
    assert is_cancer_free_text("Malignant neoplasm of colon", pats) is True
    assert is_cancer_free_text("Non-Hodgkin lymphoma", pats) is True
    assert is_cancer_free_text("Hypertension", pats) is False
    assert is_cancer_free_text("", pats) is False


# ---------- Case building ----------

def test_read_dx_ehr_and_build_cases_from_fixtures():
    schema = {
        "person_id_col": "MRN", "icd_col": "ICD_CODE",
        "icd_version_col": "ICD_VERSION", "date_col": "ENCOUNTER_DATE",
        "date_format": "%Y-%m-%d", "sep": "\t",
    }
    m = IcdToPhecodeX.load(None)
    rows = read_dx_ehr(FIXTURES / "encounter_diagnosis.tsv", schema, m,
                      source="encounter_diagnosis")
    # 4 cancer-mapping rows: TEST_001 (C18.7), TEST_002 (C50.911), TEST_005 (C20)
    assert len(rows) == 3
    cases = build_cases(rows, CANCER_PHECODES,
                        enrollment_by_pid={"TEST_001": "2015-01-01",
                                           "TEST_002": "2015-01-01",
                                           "TEST_005": "2015-01-01"})
    by_key = {(c.person_id, c.cancer): c for c in cases}
    assert by_key[("TEST_001", "colorectal")].is_case
    assert by_key[("TEST_001", "colorectal")].first_dx_date == "2020-05-10"
    assert by_key[("TEST_001", "colorectal")].incident_flag is True
    assert by_key[("TEST_002", "breast")].is_case
    assert by_key[("TEST_005", "colorectal")].is_case
    assert ("TEST_003", "colorectal") not in by_key


# ---------- Family history ----------

def test_read_family_history_and_flags():
    schema = {
        "person_id_col": "MRN", "relationship_col": "RELATIONSHIP",
        "condition_col": "CONDITION", "condition_code_col": None,
        "degree_col": None, "onset_age_col": "AGE_AT_DIAGNOSIS", "sep": "\t",
    }
    m = IcdToPhecodeX.load(None)
    fam = read_family_history(
        FIXTURES / "family_history.tsv", schema,
        cancer_free_text_patterns=["cancer", "malignant", "carcinoma"],
        cancer_phecodes=CANCER_PHECODES, mapper=m,
        relative_degrees=[1, 2],
    )
    # TEST_003: mother w/ breast cancer -> 1st deg, any-cancer positive
    assert "TEST_003" in fam
    flags3 = famhx_flags(fam["TEST_003"], CANCER_PHECODES)
    assert flags3["has_famhx_1deg_any_cancer"] is True
    assert flags3["has_famhx_2deg_any_cancer"] is False
    # TEST_004: MGM w/ colon cancer -> 2nd deg, any-cancer positive
    flags4 = famhx_flags(fam["TEST_004"], CANCER_PHECODES)
    assert flags4["has_famhx_1deg_any_cancer"] is False
    assert flags4["has_famhx_2deg_any_cancer"] is True
    # TEST_005: brother w/ diabetes -> not in famhx map
    assert "TEST_005" not in fam


# ---------- Group assignment ----------

def test_assign_group_truth_table():
    assert assign_group(is_case=True,  has_famhx_any=True)  == "case"
    assert assign_group(is_case=True,  has_famhx_any=False) == "case"
    assert assign_group(is_case=False, has_famhx_any=True)  == "famhx"
    assert assign_group(is_case=False, has_famhx_any=False) == "control"


# ---------- End-to-end synthetic build ----------

def _make_cfg(root: Path) -> dict:
    """Build a config pointing msm_data_root at `root`. Expects the phenotypes
    directory laid out as <root>/phenotypes/synthetic/*.tsv."""
    return {
        "inputs": {
            "msm_data_root": str(root),
            "phenotypes": {
                "date": "synthetic",
                "files": {
                    "encounter_diagnosis": "encounter_diagnosis.tsv",
                    "problem_list": "problem_list.tsv",
                    "family_history": "family_history.tsv",
                    "demographics": "demographics.tsv",
                    "phecodex": "phecodex_missing.tsv",   # optional; not present
                },
            },
        },
        "phenotype": {
            "icd_to_phecodex_map": None,
            "cancer_phecodes": CANCER_PHECODES,
            "ehr_schema": {
                "encounter_diagnosis": {
                    "person_id_col": "MRN", "icd_col": "ICD_CODE",
                    "icd_version_col": "ICD_VERSION", "date_col": "ENCOUNTER_DATE",
                    "date_format": "%Y-%m-%d", "sep": "\t",
                },
                "problem_list": {
                    "person_id_col": "MRN", "icd_col": "ICD_CODE",
                    "icd_version_col": "ICD_VERSION", "noted_date_col": "NOTED_DATE",
                    "resolved_col": "RESOLVED_DATE", "status_col": "STATUS",
                    "date_format": "%Y-%m-%d", "sep": "\t",
                },
                "family_history": {
                    "person_id_col": "MRN", "relationship_col": "RELATIONSHIP",
                    "condition_col": "CONDITION", "condition_code_col": None,
                    "degree_col": None, "onset_age_col": "AGE_AT_DIAGNOSIS",
                    "sep": "\t",
                },
                "demographics": {
                    "person_id_col": "MRN", "birth_year_col": "BIRTH_YEAR",
                    "sex_col": "SEX", "enrollment_date_col": "ENROLLMENT_DATE",
                    "date_format": "%Y-%m-%d", "sep": "\t",
                },
                "phecodex_roster": {
                    "format": "auto", "person_id_col": "MRN",
                    "phecodex_col": "phecodeX", "value_col": "case_control",
                    "sep": "\t",
                },
            },
            "roster": {
                "famhx_scope": "any_cancer",
                "control_exclusion": "any_cancer",
                "emit_prevalent_cases": True,
                "min_control_encounters": 0,
            },
            "famhx": {
                "relative_degrees": [1, 2],
                "cancer_free_text_patterns": ["cancer", "malignant", "carcinoma"],
            },
        },
    }


def test_end_to_end_roster_build(tmp_path: Path):
    # Lay out <tmp>/phenotypes/synthetic/*.tsv so `_phenotypes_dir` resolves.
    import shutil
    ph = tmp_path / "phenotypes" / "synthetic"
    ph.mkdir(parents=True)
    for f in ["encounter_diagnosis.tsv", "problem_list.tsv",
              "family_history.tsv", "demographics.tsv"]:
        shutil.copy(FIXTURES / f, ph / f)
    (ph / "phecodex_missing.tsv").write_text("MRN\n")
    cfg = _make_cfg(tmp_path)
    out = tmp_path / "phenotype_out"
    build(cfg, out)

    # Roster: read back and assert truth-table assignments.
    rows = list(csv_read(out / "roster.tsv"))
    by_key = {(r["person_id"], r["cancer"]): r for r in rows}

    assert by_key[("TEST_001", "colorectal")]["group"] == "case"
    assert by_key[("TEST_002", "breast")]["group"] == "case"
    assert by_key[("TEST_005", "colorectal")]["group"] == "case"

    # famhx (no cancer dx, has famhx of any cancer)
    assert by_key[("TEST_003", "colorectal")]["group"] == "famhx"
    assert by_key[("TEST_003", "breast")]["group"] == "famhx"
    assert by_key[("TEST_004", "colorectal")]["group"] == "famhx"

    # control_exclusion=any_cancer -> TEST_005 has CRC, so their "breast" row is excluded
    assert by_key[("TEST_005", "breast")]["group"] == "excluded"
    assert by_key[("TEST_005", "breast")]["excluded_from_control_reason"] == "has_other_cancer"

    # Same for TEST_001/002 non-case cancers.
    assert by_key[("TEST_001", "breast")]["group"] == "excluded"
    assert by_key[("TEST_002", "colorectal")]["group"] == "excluded"

    # Summary file exists and mentions each cancer
    summary = (out / "roster_summary.md").read_text()
    assert "colorectal" in summary and "breast" in summary


def csv_read(path: Path):
    import csv
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        yield from rdr
