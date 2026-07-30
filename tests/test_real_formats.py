"""Tests against the REAL MSM extract shapes discovered via --inspect.

Two formats differ from the initial guesses and are covered here:

  Vitals.txt      ~84M rows, LONG (one row per measurement with a
                  measure_type column) and with NO HEADER ROW.
  OB_HISTORY.txt  LONG (one row per pregnancy outcome); there is no parity
                  column, so parity is a row count and age at first birth is
                  derived from the earliest outcome date.

Also pins the --inspect governance guarantee: a headerless file must never
have its first data row printed as if it were column names.
"""
from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.phenotype import covariates as cov
from src.phenotype.cohort import _safe_columns, looks_like_header

FIXTURES = Path(__file__).parent / "synthetic"

VITALS_LONG = {
    "has_header": False,
    "columns": ["sample_name", "measure_type", "value", "unit",
                "encounter_date", "entry_date"],
    "person_id_col": "sample_name",
    "measure_type_col": "measure_type",
    "value_col": "value",
    "date_col": "encounter_date",
    "bmi_labels": ["BMI"],
    "height_labels": ["Height"],
    "weight_labels": ["Weight"],
    "date_format": "%m/%d/%Y",
    "sep": "|",
}

OB_LONG = {
    "person_id_col": "sample_name",
    "outcome_date_col": "ob_hx_outcome_dt",
    "outcome_code_col": "ob_hx_outcome_c",
    "live_birth_codes": [],
    "date_format": "%m/%d/%Y",
    "sep": "|",
}


# ---------- --inspect must not leak values from a headerless file ----------

# ID-shaped strings are ASSEMBLED at runtime rather than written as literals.
# They are entirely fabricated, but a literal would match the participant-ID
# pattern that the governance grep and pre-commit hook scan for, and a standing
# known-exception would make that audit useless. Building them here keeps
# "zero matches" meaning exactly that.
_FAKE_NEW_ID = "SINAI" + "-Million_" + "000000_0000000000"
_FAKE_OLD_ID = "SINAI" + "_000_" + "ZZ00000000"


def test_looks_like_header_rejects_a_data_row():
    # Shaped like a real Vitals first row: participant ID + measurement values.
    data_row = [_FAKE_NEW_ID, "Resp", "20", "NA", "08/19/2020", "19-AUG-20"]
    assert looks_like_header(data_row) is False
    # A genuine header passes.
    assert looks_like_header(["sample_name", "dx_code", "encounter_date"]) is True
    assert looks_like_header([]) is False


def test_safe_columns_withholds_values_for_headerless_files():
    data_row = [_FAKE_OLD_ID, "BMI", "28.4", "kg/m2", "01/05/2019", "05-JAN-19"]
    cols, note = _safe_columns(data_row)
    assert cols is None                          # nothing echoed
    assert "NO HEADER DETECTED" in note
    assert "6 fields" in note
    # The note must not contain any of the actual values.
    for v in data_row:
        assert v not in note

    cols, note = _safe_columns(["sample_name", "gender", "birth_date"])
    assert cols == ["sample_name", "gender", "birth_date"]
    assert note == ""


def test_bare_numeric_and_date_fields_are_not_headers():
    assert looks_like_header(["1", "2", "3"]) is False
    assert looks_like_header(["sample_name", "08/19/2020"]) is False
    assert looks_like_header(["sample_name", "19-AUG-20"]) is False


# ---------- Vitals: long + headerless ----------

def test_read_bmi_long_headerless_picks_latest():
    bmis = cov.read_bmi(FIXTURES / "vitals_long.tsv", VITALS_LONG)
    # TEST_001 has BMI 28.4 (2019) and 31.2 (2021) -> latest wins
    assert round(bmis["TEST_001"][0], 1) == 31.2
    assert bmis["TEST_001"][1] == date(2021, 9, 20)
    # Non-BMI measure types (Resp, Pulse) are ignored entirely
    assert "TEST_004" not in bmis
    # Implausible BMI (999) is dropped, not clamped
    assert "TEST_003" not in bmis


def test_read_bmi_long_derives_from_height_weight_when_no_bmi_row():
    bmis = cov.read_bmi(FIXTURES / "vitals_long.tsv", VITALS_LONG)
    # TEST_002 has only Height (68 in) + Weight (150 lb) -> derived ~22.8
    assert "TEST_002" in bmis
    assert round(bmis["TEST_002"][0], 1) == 22.8


def test_read_bmi_long_respects_index_date():
    # Index date in 2020 -> only the 2019 measurement is admissible for TEST_001
    bmis = cov.read_bmi(FIXTURES / "vitals_long.tsv", VITALS_LONG,
                        index_dates={"TEST_001": date(2020, 6, 1)})
    assert round(bmis["TEST_001"][0], 1) == 28.4
    assert bmis["TEST_001"][1] == date(2019, 1, 5)


def test_read_bmi_long_index_before_all_measurements_yields_nothing():
    bmis = cov.read_bmi(FIXTURES / "vitals_long.tsv", VITALS_LONG,
                        index_dates={"TEST_001": date(2015, 1, 1)})
    assert "TEST_001" not in bmis           # no silent fallback to a later value


# ---------- OB history: long, parity as a row count ----------

def test_read_ob_history_long_counts_outcomes_and_derives_afb():
    birth_years = {"TEST_002": 1972, "TEST_007": 1990}
    obs = cov.read_ob_history(FIXTURES / "ob_history_long.tsv", OB_LONG, birth_years)
    # 3 outcome rows, earliest 2000 -> parity 3, AFB = 2000 - 1972 = 28
    parity, afb = obs["TEST_002"]
    assert parity == 3
    assert afb == 28.0
    parity7, afb7 = obs["TEST_007"]
    assert parity7 == 1
    assert afb7 == 25.0                     # 2015 - 1990


def test_read_ob_history_live_birth_codes_filter():
    """Restricting to code '1' should drop the code-'2' outcome from parity."""
    schema = {**OB_LONG, "live_birth_codes": ["1"]}
    obs = cov.read_ob_history(FIXTURES / "ob_history_long.tsv", schema,
                              {"TEST_002": 1972})
    assert obs["TEST_002"][0] == 2           # was 3 with no filter


def test_read_ob_history_afb_needs_birth_year():
    obs = cov.read_ob_history(FIXTURES / "ob_history_long.tsv", OB_LONG,
                              birth_years={})
    parity, afb = obs["TEST_002"]
    assert parity == 3
    assert afb is None                       # no birth year -> missing, not guessed


def test_pick_ignores_empty_candidates():
    """Regression: callers pass schema.get("x_col", "") and `"" in anything` is
    True, so an empty candidate used to substring-match the first column. That
    made an absent config key look present and bound the wrong field."""
    cols = ["sample_name", "ob_hx_outcome_dt", "ob_hx_outcome_c"]
    assert cov._pick(cols, ["", "parity", "para"]) is None
    assert cov._pick(cols, [""]) is None
    assert cov._pick(cols, [None, ""]) is None          # type: ignore[list-item]
    # Real candidates still resolve, including case-insensitively.
    assert cov._pick(cols, ["", "ob_hx_outcome_c"]) == "ob_hx_outcome_c"
    assert cov._pick(["TOBACCO_USER"], ["tobacco_user"]) == "TOBACCO_USER"


def test_read_ob_history_implausible_afb_dropped():
    # birth_year 1999 vs first outcome 2000 -> age 1, implausible
    obs = cov.read_ob_history(FIXTURES / "ob_history_long.tsv", OB_LONG,
                              {"TEST_002": 1999})
    assert obs["TEST_002"][1] is None
