"""Tests pinned to the ACTUAL MSM value vocabularies.

Each block below encodes a trap that would otherwise produce a silently wrong
number rather than an error:

  "Never Assessed"  contains "never" but means NOT ASKED. Scored as
                    never-smoker it would bias the smoking exposure.
  "Not Currently"   is a FORMER drinker, not a never-drinker.
  "Para"/"Gravida"  are summary counts, not pregnancy outcomes; counting them
                    inflates parity.
  no BMI rows       Vitals has Height and Weight but no BMI measure_type, so
                    every BMI is derived and the weight unit decides whether
                    the column is usable at all.
  `Group`           is a curated vocabulary; pattern-matching it loses detail
                    (Endometrium vs Uterus, Multiple Myeloma vs Leukemia).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.phenotype import covariates as cov
from src.phenotype.cancer_cases import slugify_cancer_group

# --- Real vocabularies, transcribed from value counts on the cluster ---------

TOBACCO_USER_VALUES = [
    "Never", "Former", "Never Smoker", "Former Smoker", "Every Day",
    "Current Every Day Smoker", "Never Assessed", "Some Days",
    "Current Some Day Smoker",
]
IS_ALCOHOL_USER_VALUES = ["No", "Yes", "Not Currently", "NA", "Not Asked", "Never", ""]

SMOKING_MAP = {
    "unknown": ["never assessed", "not asked", "not assessed", "na", "unknown"],
    "never": ["never", "never smoker", "non-smoker", "nonsmoker"],
    "former": ["former", "former smoker", "quit", "ex-smoker", "previous"],
    "current": ["every day", "some days", "current every day smoker",
                "current some day smoker", "yes"],
}
ALCOHOL_MAP = {
    "unknown": ["na", "not asked", "not assessed", "unknown"],
    "never": ["no", "never", "none", "non-drinker"],
    "former": ["not currently", "former", "quit", "previous"],
    "current": ["yes", "current", "occasional", "social", "moderate", "heavy", "daily"],
}


# ---------- smoking: the "Never Assessed" trap ----------

def test_never_assessed_is_unknown_not_never():
    """96,717 rows say "Never Assessed", meaning not asked. Substring matching
    on "never" would score them as non-smokers."""
    assert cov.normalize_category("Never Assessed", SMOKING_MAP) == "unknown"
    # And a genuine never-smoker is still never.
    assert cov.normalize_category("Never", SMOKING_MAP) == "never"
    assert cov.normalize_category("Never Smoker", SMOKING_MAP) == "never"


def test_every_real_tobacco_value_maps_sensibly():
    expected = {
        "Never": "never",
        "Never Smoker": "never",
        "Former": "former",
        "Former Smoker": "former",
        "Every Day": "current",
        "Some Days": "current",
        "Current Every Day Smoker": "current",
        "Current Some Day Smoker": "current",
        "Never Assessed": "unknown",
    }
    for raw in TOBACCO_USER_VALUES:
        got = cov.normalize_category(raw, SMOKING_MAP)
        assert got == expected[raw], f"{raw!r} -> {got!r}, expected {expected[raw]!r}"


# ---------- alcohol: "Not Currently" is former, not never ----------

def test_not_currently_is_former_not_never():
    assert cov.normalize_category("Not Currently", ALCOHOL_MAP) == "former"
    assert cov.normalize_category("No", ALCOHOL_MAP) == "never"
    assert cov.normalize_category("Never", ALCOHOL_MAP) == "never"


def test_every_real_alcohol_value_maps_sensibly():
    expected = {
        "No": "never", "Never": "never",
        "Yes": "current",
        "Not Currently": "former",
        "NA": "unknown", "Not Asked": "unknown",
        "": "unknown",                      # blank -> nullish -> unknown
    }
    for raw in IS_ALCOHOL_USER_VALUES:
        got = cov.normalize_category(raw, ALCOHOL_MAP)
        assert got == expected[raw], f"{raw!r} -> {got!r}, expected {expected[raw]!r}"


# ---------- BMI: no BMI rows, so units decide everything ----------

def test_infer_weight_unit_separates_lb_kg_oz():
    # Adult medians: ~77 kg, ~170 lb, ~2,700 oz.
    assert cov.infer_weight_unit([70.0, 80.0, 77.0]) == "kg"
    assert cov.infer_weight_unit([160.0, 180.0, 170.0]) == "lb"
    assert cov.infer_weight_unit([2600.0, 2800.0, 2700.0]) == "oz"
    assert cov.infer_weight_unit([]) == "lb"          # documented default


def test_infer_height_unit_separates_in_cm():
    assert cov.infer_height_unit([64.0, 66.0, 70.0]) == "in"
    assert cov.infer_height_unit([160.0, 168.0, 175.0]) == "cm"
    assert cov.infer_height_unit([]) == "in"


def test_bmi_from_is_consistent_across_unit_systems():
    """The same person expressed in three unit systems must give one BMI."""
    metric = cov.bmi_from(170.0, 70.0, "cm", "kg")
    imperial = cov.bmi_from(66.93, 154.32, "in", "lb")
    ounces = cov.bmi_from(66.93, 2469.1, "in", "oz")
    assert round(metric, 1) == 24.2
    assert abs(metric - imperial) < 0.2
    assert abs(metric - ounces) < 0.2


def test_ounce_weights_would_be_unusable_without_unit_inference():
    """Regression: treating ounces as pounds makes every BMI implausible, so the
    column silently comes back empty rather than obviously wrong."""
    # 2700 oz (~168 lb) at 66 in. Misread as pounds -> absurd -> rejected.
    assert cov.bmi_from(66.0, 2700.0, "in", "lb") is None
    # With the unit correctly inferred, it is a normal BMI.
    ok = cov.bmi_from(66.0, 2700.0, "in", "oz")
    assert ok is not None and 20.0 <= ok <= 35.0


def test_bmi_report_surfaces_inferred_units(tmp_path: Path):
    """Height/Weight-only long vitals must derive BMI and report the inference."""
    p = tmp_path / "vitals.tsv"
    # Weight in pounds, height in inches, no BMI rows at all.
    lines = []
    for i in range(1, 6):
        pid = f"TEST_{i:03d}"
        lines.append(f"{pid}|Height|66|in|01/05/2019|05-JAN-19")
        lines.append(f"{pid}|Weight|{150 + i}|lb|01/05/2019|05-JAN-19")
        lines.append(f"{pid}|Pulse|72|bpm|01/05/2019|05-JAN-19")
    p.write_text("\n".join(lines) + "\n")

    schema = {
        "has_header": False,
        "columns": ["sample_name", "measure_type", "value", "unit",
                    "encounter_date", "entry_date"],
        "person_id_col": "sample_name", "measure_type_col": "measure_type",
        "value_col": "value", "date_col": "encounter_date",
        "bmi_labels": [], "height_labels": ["Height"], "weight_labels": ["Weight"],
        "date_format": "%m/%d/%Y", "sep": "|",
    }
    report: dict = {}
    bmis = cov.read_bmi(p, schema, report=report)
    assert len(bmis) == 5
    assert report["height_unit_inferred"] == "in"
    assert report["weight_unit_inferred"] == "lb"
    assert report["n_bmi_derived"] == 5
    assert 22.0 <= report["bmi_median"] <= 28.0
    # Non-anthropometric measure types are ignored.
    assert all(10.0 <= v <= 80.0 for v, _ in bmis.values())


# ---------- Group: curated passthrough ----------

def test_slugify_cancer_group_real_values():
    expected = {
        "Breast": "breast",
        "Colorectum": "colorectum",
        "Multiple Myeloma": "multiple_myeloma",
        "Prostate": "prostate",
        "Head & Neck": "head_and_neck",
        "Neuroendocrine Tumors": "neuroendocrine_tumors",
        "Endometrium": "endometrium",
        "Uterus": "uterus",
        "Skin/Mucosa": "skin_mucosa",
        "Urinary Tract": "urinary_tract",
        "Soft Tissue": "soft_tissue",
    }
    for raw, want in expected.items():
        assert slugify_cancer_group(raw) == want, f"{raw!r}"


def test_site_less_groups_return_none_so_they_are_not_a_spurious_site():
    """These people ARE cases, but the group carries no primary site."""
    assert slugify_cancer_group("Cancer of Unknown Primary (CUP)") is None
    assert slugify_cancer_group("Other/Ill-Defined Primary Sites") is None
    assert slugify_cancer_group(
        "Hematologic/Lymphatic - Other/Unspecified B-cell") is None
    assert slugify_cancer_group("") is None
    assert slugify_cancer_group("NA") is None


def test_passthrough_keeps_distinctions_that_pattern_matching_would_merge():
    """Endometrium and Uterus are distinct groups; a uterine pattern would
    collapse them. Multiple Myeloma must not become 'myeloma' alongside
    Leukemia/Lymphoma."""
    assert slugify_cancer_group("Endometrium") != slugify_cancer_group("Uterus")
    assert slugify_cancer_group("Multiple Myeloma") == "multiple_myeloma"
    assert slugify_cancer_group("Leukemia") == "leukemia"
    assert slugify_cancer_group("Lymphoma") == "lymphoma"


# ---------- OB: Para/Gravida are not outcomes ----------

def test_live_birth_filter_excludes_summary_and_nonbirth_codes(tmp_path: Path):
    p = tmp_path / "ob.tsv"
    p.write_text(
        "sample_name|ob_hx_outcome_dt|ob_hx_outcome_c|encounter_date\n"
        "TEST_002|06/15/2000|Term|01/01/2015\n"
        "TEST_002|09/02/2003|Preterm|01/01/2015\n"
        "TEST_002|04/10/2006|Abortion|01/01/2015\n"
        "TEST_002|04/10/2006|Spontaneous Abortion|01/01/2015\n"
        "TEST_002|04/10/2006|Therapeutic Abortion|01/01/2015\n"
        "TEST_002|01/01/2008|Ectopic|01/01/2015\n"
        "TEST_002||Para|01/01/2015\n"
        "TEST_002||Gravida|01/01/2015\n"
        "TEST_002||Current|01/01/2015\n"
    )
    schema = {
        "person_id_col": "sample_name",
        "outcome_date_col": "ob_hx_outcome_dt",
        "outcome_code_col": "ob_hx_outcome_c",
        "live_birth_codes": ["Term", "Preterm"],
        "date_format": "%m/%d/%Y", "sep": "|",
    }
    obs = cov.read_ob_history(p, schema, {"TEST_002": 1972})
    parity, afb = obs["TEST_002"]
    assert parity == 2                       # Term + Preterm only
    assert afb == 28.0                       # 2000 - 1972

    # Without the filter, all 9 rows count -> parity 9, which is the bug.
    unfiltered = cov.read_ob_history(p, {**schema, "live_birth_codes": []},
                                     {"TEST_002": 1972})
    assert unfiltered["TEST_002"][0] == 9
