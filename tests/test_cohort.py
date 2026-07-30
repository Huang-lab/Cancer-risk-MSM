"""Tests for the EHR cohort build: cases, controls, famhx, lifestyle.

Fixture truth (all synthetic TEST_### IDs):

  cancers_combined.tsv     TEST_001 colon, TEST_002 breast (x2 rows),
                           TEST_005 rectal, TEST_006 prostate,
                           TEST_009 unclassifiable type (still a case)
  demographics_cohort.tsv  roster = TEST_001..TEST_008  (note: 009 absent)
  encounter_diagnosis.tsv  TEST_001 C18.7, TEST_002 C50.911, TEST_005 C20,
                           TEST_006 C61  (cancer ICD codes)
  family_history.tsv       TEST_003 mother-breast, TEST_004 father-colon,
                           TEST_006 brother-prostate, TEST_007 mother-"Cancer",
                           TEST_005 brother-diabetes (not cancer)
"""
from __future__ import annotations

import csv
import shutil
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.phenotype import covariates as cov
from src.phenotype.cancer_cases import (
    describe_schema,
    load_cancer_cases,
    normalize_cancer_type,
)
from src.phenotype.cohort import build
from src.phenotype.controls import assign_groups, is_cancer_phecode

FIXTURES = Path(__file__).parent / "synthetic"


# ---------- cancer case loading ----------

def test_normalize_cancer_type():
    assert normalize_cancer_type("Colon adenocarcinoma") == "colorectal"
    assert normalize_cancer_type("Rectal cancer") == "colorectal"
    assert normalize_cancer_type("Invasive ductal carcinoma of breast") == "breast"
    assert normalize_cancer_type("Prostate cancer") == "prostate"
    assert normalize_cancer_type("Some Rare Thing") is None
    assert normalize_cancer_type("") is None
    assert normalize_cancer_type("unknown") is None


def test_describe_schema_reports_columns_only():
    info = describe_schema(FIXTURES / "cancers_combined.tsv")
    assert info["n_rows"] == 6
    assert "sample_name" in info["columns"]
    assert info["detected_cancer_type_col"] == "cancer_type"
    assert info["detected_genetics_col"] == "has_genetics"


def test_load_cancer_cases_aggregates_multirow_people():
    cases, rep = load_cancer_cases(FIXTURES / "cancers_combined.tsv")
    # 6 rows -> 5 distinct people (TEST_002 appears twice)
    assert rep["n_rows"] == 6
    assert rep["n_people"] == 5
    assert set(cases) == {"TEST_001", "TEST_002", "TEST_005", "TEST_006", "TEST_009"}
    # Both TEST_002 rows normalize to the same type -> one entry, two raw strings
    assert cases["TEST_002"].cancer_types == {"breast"}
    assert cases["TEST_002"].n_rows == 2
    assert len(cases["TEST_002"].raw_types) == 2
    assert cases["TEST_005"].cancer_types == {"colorectal"}
    assert cases["TEST_006"].has_genetics is False or cases["TEST_006"].has_genetics is True


def test_unclassifiable_type_still_counts_as_a_case():
    """A type string we can't map must not silently drop the person."""
    cases, rep = load_cancer_cases(FIXTURES / "cancers_combined.tsv")
    assert "TEST_009" in cases
    assert cases["TEST_009"].cancer_types == set()      # no type assigned
    assert cases["TEST_009"].n_cancers == 1             # but still one cancer
    assert rep["n_rows_unclassified_type"] == 1
    assert "Some Rare Thing" in rep["unclassified_type_examples"]


def test_missing_id_column_reports_rather_than_crashing(tmp_path: Path):
    p = tmp_path / "bad.tsv"
    p.write_text("foo\tbar\n1\t2\n")
    cases, rep = load_cancer_cases(p, person_id_col="sample_name")
    assert cases == {}
    assert any("not found" in w for w in rep["warnings"])


# ---------- control assignment ----------

def test_is_cancer_phecode():
    assert is_cancer_phecode("CA_101.1") is True
    assert is_cancer_phecode("153") is True
    assert is_cancer_phecode("174.1") is True
    assert is_cancer_phecode("401.1") is False
    assert is_cancer_phecode("") is False


def test_assign_groups_basic_set_difference():
    roster = {f"TEST_{i:03d}" for i in range(1, 9)}
    cases = {"TEST_001", "TEST_002"}
    groups, rep = assign_groups(roster, cases)
    assert groups["TEST_001"] == "case"
    assert groups["TEST_003"] == "control"
    assert rep.n_cases == 2
    assert rep.n_controls == 6


def test_icd_screen_keeps_undocumented_cancers_out_of_controls():
    """Someone with a cancer ICD code missing from the case file must not
    become a 'non-cancer' control."""
    roster = {"TEST_001", "TEST_002", "TEST_003"}
    cases = {"TEST_001"}
    icd_cancer = {"TEST_002"}          # has cancer ICD but absent from case file
    groups, rep = assign_groups(roster, cases, icd_cancer)
    assert groups["TEST_001"] == "case"
    assert groups["TEST_002"] == "excluded_icd_cancer"     # NOT a control
    assert groups["TEST_003"] == "control"
    assert rep.n_controls == 1
    assert rep.n_excluded_icd_cancer == 1


def test_case_absent_from_roster_is_kept_and_reported():
    roster = {"TEST_001", "TEST_002"}
    cases = {"TEST_001", "TEST_009"}      # 009 not in roster
    groups, rep = assign_groups(roster, cases)
    assert groups["TEST_009"] == "case"   # kept, not dropped
    assert rep.n_cases_not_in_roster == 1
    assert any("absent from the roster" in w for w in rep.warnings)


def test_min_encounters_filter():
    roster = {"TEST_001", "TEST_002"}
    groups, rep = assign_groups(roster, set(), None,
                                {"TEST_001": 5, "TEST_002": 1}, min_encounters=2)
    assert groups["TEST_001"] == "control"
    assert groups["TEST_002"] == "excluded_few_encounters"


# ---------- lifestyle covariates ----------

def test_normalize_category_exact_beats_substring():
    smoking = {
        "never": ["never", "never smoker"],
        "former": ["former", "former smoker"],
        "current": ["current", "current smoker", "smoker"],
    }
    # 'never smoker' must map to never, not current via the 'smoker' pattern
    assert cov.normalize_category("Never smoker", smoking) == "never"
    assert cov.normalize_category("Former smoker", smoking) == "former"
    assert cov.normalize_category("Current Smoker", smoking) == "current"
    assert cov.normalize_category("", smoking) == "unknown"
    assert cov.normalize_category("gibberish", smoking) == "unknown"


def test_derive_bmi_units_and_plausibility():
    assert round(cov.derive_bmi(170, 82), 1) == 28.4        # cm + kg
    assert round(cov.derive_bmi(68, 150), 1) == 22.8        # inches + lb
    assert cov.derive_bmi(170, 5) is None                   # implausible -> None
    assert cov.derive_bmi(None, 80) is None
    assert cov.derive_bmi(0, 80) is None


def test_pick_measurement_respects_index_date():
    rows = [(date(2019, 1, 5), 28.4), (date(2021, 9, 20), 31.2)]
    # No index date -> most recent
    v, d = cov.pick_measurement(rows, None)
    assert v == 31.2 and d == date(2021, 9, 20)
    # Index date in 2020 -> only the pre-index measurement is admissible
    v, d = cov.pick_measurement(rows, date(2020, 6, 1))
    assert v == 28.4 and d == date(2019, 1, 5)
    # Index date before everything -> nothing usable, not a silent fallback
    v, d = cov.pick_measurement(rows, date(2018, 1, 1))
    assert v is None and d is None


def test_read_bmi_drops_implausible_and_derives_missing():
    schema = {"person_id_col": "sample_name", "bmi_col": "bmi",
              "height_col": "height", "weight_col": "weight",
              "date_col": "encounter_date", "date_format": "%m/%d/%Y", "sep": "|"}
    bmis = cov.read_bmi(FIXTURES / "vitals.tsv", schema)
    assert round(bmis["TEST_001"][0], 1) == 31.2       # most recent
    assert round(bmis["TEST_003"][0], 1) == 22.0       # derived from h/w
    assert "TEST_004" not in bmis                      # bmi=999 dropped


def test_read_social_history_latest_value():
    schema = {"person_id_col": "sample_name", "smoking_status_col": "smoking_status",
              "alcohol_col": "alcohol_use", "date_col": "encounter_date",
              "date_format": "%m/%d/%Y", "sep": "|"}
    smap = {"never": ["never", "never smoker"], "former": ["former smoker"],
            "current": ["current smoker"]}
    amap = {"never": ["never", "no"], "current": ["occasional", "social", "heavy"]}
    social = cov.read_social_history(FIXTURES / "social_history.tsv", schema, smap, amap)
    # TEST_001 has 2018 former + 2021 current -> latest wins
    assert social["TEST_001"] == ("current", "current")
    assert social["TEST_003"] == ("never", "current")


def test_read_ob_history_max_parity_min_afb():
    schema = {"person_id_col": "sample_name", "parity_col": "parity",
              "age_at_first_birth_col": "age_at_first_birth", "sep": "|"}
    obs = cov.read_ob_history(FIXTURES / "ob_history.tsv", schema)
    assert obs["TEST_002"] == (3, 27.0)     # parity cumulative -> max
    assert obs["TEST_003"] == (0, None)


# ---------- end-to-end ----------

def _cfg(root: Path) -> dict:
    common = {"person_id_col": "sample_name", "sep": "|"}
    return {
        "project": {"output_root": str(root / "out")},
        "inputs": {
            "msm_data_root": str(root),
            "cancer_cases": {
                "path": str(FIXTURES / "cancers_combined.tsv"),
                "person_id_col": "sample_name", "cancer_type_col": None, "sep": "\t",
            },
            "phenotypes": {
                "date": "synthetic",
                "files": {
                    "demographics": "demographics_cohort.tsv",
                    "family_history": "family_history.tsv",
                    "social_history": "social_history.tsv",
                    "vitals": "vitals.tsv",
                    "ob_history": "ob_history.tsv",
                    "encounter_diagnosis": "encounter_diagnosis.tsv",
                    "problem_list": "problem_list.tsv",
                },
            },
        },
        "phenotype": {
            "icd_to_phecodex_map": None,
            "cancer_phecodes": {"colorectal": ["CA_101.1"], "breast": ["CA_105.1"]},
            "ehr_schema": {
                "demographics": {**common, "birth_year_col": "birth_year", "sex_col": "sex"},
                "family_history": {**common, "relationship_col": "relation",
                                   "condition_col": "problem_description",
                                   "condition_code_col": None, "degree_col": None},
                "social_history": {**common, "smoking_status_col": "smoking_status",
                                   "alcohol_col": "alcohol_use",
                                   "date_col": "encounter_date", "date_format": "%m/%d/%Y"},
                "vitals": {**common, "bmi_col": "bmi", "height_col": "height",
                           "weight_col": "weight", "date_col": "encounter_date",
                           "date_format": "%m/%d/%Y"},
                "ob_history": {**common, "parity_col": "parity",
                               "age_at_first_birth_col": "age_at_first_birth"},
                "encounter_diagnosis": {**common, "icd_col": "dx_code",
                                        "icd_version_col": "dx_code_type",
                                        "date_col": "encounter_date",
                                        "date_format": "%m/%d/%Y"},
                "problem_list": {**common, "icd_col": "icd9_code",
                                 "icd_version_col": "code_type",
                                 "noted_date_col": "noted_date",
                                 "date_format": "%m/%d/%Y"},
            },
            "famhx": {
                "relative_degrees": [1],
                "cancer_free_text_patterns": ["cancer", "malignant", "carcinoma"],
                "cancer_type_patterns": {
                    "colorectal": ["colon", "rectal"], "breast": ["breast"],
                    "prostate": ["prostate"],
                },
            },
        },
        "cohort": {
            "roster_source": "demographics",
            "screen_icd_for_cancer": True,
            "min_control_encounters": 0,
            "smoking_map": {"never": ["never smoker"], "former": ["former smoker"],
                            "current": ["current smoker"]},
            "alcohol_map": {"never": ["never", "no"],
                            "current": ["occasional", "social", "heavy"]},
        },
    }


def _stage(tmp_path: Path) -> None:
    ph = tmp_path / "phenotypes" / "synthetic"
    ph.mkdir(parents=True)
    for f in ["demographics_cohort.tsv", "family_history.tsv", "social_history.tsv",
              "vitals.tsv", "ob_history.tsv", "encounter_diagnosis.tsv",
              "problem_list.tsv"]:
        shutil.copy(FIXTURES / f, ph / f)


def _read(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_end_to_end_cohort(tmp_path: Path, capsys):
    _stage(tmp_path)
    out = tmp_path / "cohort_out"
    assert build(_cfg(tmp_path), out) == 0

    rows = {r["sample_name"]: r for r in _read(out / "cohort.tsv")}

    # Cases from the case file
    assert rows["TEST_001"]["group"] == "case"
    assert rows["TEST_001"]["cancer_types"] == "colorectal"
    assert rows["TEST_002"]["cancer_types"] == "breast"
    assert rows["TEST_006"]["group"] == "case"
    # Case absent from the roster is still emitted as a case
    assert rows["TEST_009"]["group"] == "case"

    # Clean controls: no cancer in the case file AND no cancer ICD code
    assert rows["TEST_003"]["group"] == "control"
    assert rows["TEST_004"]["group"] == "control"
    assert rows["TEST_007"]["group"] == "control"
    assert rows["TEST_008"]["group"] == "control"

    # First-degree family history, everyone (cases and controls alike)
    assert rows["TEST_003"]["has_famhx_1deg_any_cancer"] == "1"
    assert rows["TEST_003"]["has_famhx_1deg_breast"] == "1"
    assert rows["TEST_004"]["has_famhx_1deg_colorectal"] == "1"
    assert rows["TEST_006"]["has_famhx_1deg_prostate"] == "1"
    # Unclassified "Cancer" counts for any-cancer but no specific type
    assert rows["TEST_007"]["has_famhx_1deg_any_cancer"] == "1"
    assert rows["TEST_007"]["has_famhx_1deg_breast"] == "0"
    # Brother with diabetes is not cancer famhx
    assert rows["TEST_005"]["has_famhx_1deg_any_cancer"] == "0"

    # Lifestyle
    assert rows["TEST_001"]["smoking_status"] == "current"
    assert rows["TEST_003"]["smoking_status"] == "never"
    assert rows["TEST_001"]["bmi"] == "31.2"
    assert rows["TEST_002"]["parity"] == "3"
    assert rows["TEST_002"]["age_at_first_birth"] == "27"
    assert rows["TEST_001"]["sex"] == "M"

    summary = (out / "cohort_summary.md").read_text()
    assert "case" in summary and "control" in summary
    assert "Covariate missingness" in summary


def test_missingness_reflects_merged_lifestyle_not_demographics_only(tmp_path: Path):
    """Regression: lifestyle values are read into separate dicts. If they aren't
    merged onto the Covariates records, the summary reports 100% missing for
    bmi/smoking/parity even though cohort.tsv has them populated."""
    _stage(tmp_path)
    out = tmp_path / "cohort_out3"
    build(_cfg(tmp_path), out)

    rows = {r["sample_name"]: r for r in _read(out / "cohort.tsv")}
    # These ARE populated in the table...
    assert rows["TEST_001"]["bmi"] == "31.2"
    assert rows["TEST_001"]["smoking_status"] == "current"
    assert rows["TEST_002"]["parity"] == "3"

    # ...so the summary must not claim they are fully missing.
    summary = (out / "cohort_summary.md").read_text()
    for field in ["bmi", "smoking_status", "parity"]:
        assert f"| {field} | 100.0% |" not in summary, \
            f"{field} reported 100% missing despite being populated"


def test_icd_screen_can_be_disabled(tmp_path: Path):
    """With the screen off, cases are still cases (they're in the case file),
    but nobody gets the excluded_icd_cancer label."""
    _stage(tmp_path)
    cfg = _cfg(tmp_path)
    cfg["cohort"]["screen_icd_for_cancer"] = False
    out = tmp_path / "cohort_out2"
    build(cfg, out)
    groups = {r["sample_name"]: r["group"] for r in _read(out / "cohort.tsv")}
    assert "excluded_icd_cancer" not in groups.values()
    assert groups["TEST_001"] == "case"
