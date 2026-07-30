"""Tests for config validation and the case-file-driven cancer code screen.

Both cover defects the first real run exposed:

  * Inserting a top-level `cohort:` section above `phenotype.roster/famhx`
    reparented those blocks, surfacing much later as a bare
    KeyError: 'famhx'. validate_config now names the actual location.

  * The built-in ICD->phecodeX map covers ~6 cancers while the case file spans
    54 Group categories, so screening controls with it missed ~45 categories
    and admitted people with undocumented cancers as "non-cancer" controls.
    The code vocabulary is now taken from the case file itself.
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.phenotype.cancer_cases import load_cancer_code_prefixes
from src.phenotype.cohort import validate_config
from src.phenotype.controls import people_with_cancer_codes

REPO = Path(__file__).resolve().parents[1]


# ---------- the shipped config must be correctly nested ----------

def test_shipped_config_passes_validation():
    cfg = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    assert validate_config(cfg) == []


def test_shipped_config_keeps_famhx_and_roster_under_phenotype():
    """Regression: a top-level `cohort:` inserted above these blocks silently
    reparented them into `cohort:`, and the build died on KeyError: 'famhx'."""
    cfg = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    assert "famhx" in cfg["phenotype"]
    assert "roster" in cfg["phenotype"]
    assert "famhx" not in cfg.get("cohort", {})
    assert "roster" not in cfg.get("cohort", {})
    # famhx must still carry the keys the readers dereference.
    fam = cfg["phenotype"]["famhx"]
    assert "cancer_free_text_patterns" in fam
    assert "cancer_type_patterns" in fam
    assert "relative_degrees" in fam


def test_validate_config_points_at_the_wrong_parent():
    cfg = {
        "inputs": {"cancer_cases": {}, "phenotypes": {}},
        "phenotype": {"ehr_schema": {}, "cancer_phecodes": {}},
        "cohort": {"famhx": {}},          # mis-nested, as happened for real
    }
    problems = validate_config(cfg)
    assert len(problems) == 1
    assert "phenotype.famhx" in problems[0]
    assert "under `cohort:`" in problems[0]      # names where it actually is


def test_validate_config_reports_missing_section():
    problems = validate_config({})
    assert any("top-level `inputs:`" in p for p in problems)
    assert any("top-level `phenotype:`" in p for p in problems)


def test_build_exits_nonzero_on_bad_config(tmp_path: Path, capsys):
    from src.phenotype.cohort import build
    rc = build({"phenotype": {}}, tmp_path / "out")
    assert rc == 2                                # no traceback, actionable rc
    err = capsys.readouterr().err
    assert "CONFIG ERROR" in err


# ---------- cancer code vocabulary from the case file ----------

def _case_file(tmp_path: Path) -> Path:
    """Case file spanning categories the built-in map does NOT cover:
    C22 liver, C90 myeloma, C91 leukemia, C73 thyroid, alongside C18/C50."""
    p = tmp_path / "cancers.tsv"
    p.write_text(
        "sample_name\tdx_code\tdx_code_type\tdx_name\tGroup\n"
        "TEST_001\tC18.7\tICD10\tSigmoid colon ca\tColorectum\n"
        "TEST_002\tC50.911\tICD10\tBreast ca\tBreast\n"
        "TEST_003\tC22.0\tICD10\tLiver cell ca\tLiver\n"
        "TEST_004\tC90.00\tICD10\tMultiple myeloma\tMultiple Myeloma\n"
        "TEST_005\tC91.00\tICD10\tALL\tLeukemia\n"
        "TEST_006\tC73\tICD10\tThyroid ca\tThyroid\n"
    )
    return p


def test_code_prefixes_cover_every_category_in_the_case_file(tmp_path: Path):
    prefixes = load_cancer_code_prefixes(_case_file(tmp_path))
    assert prefixes == {"C18", "C50", "C22", "C90", "C91", "C73"}


def test_screen_catches_cancers_the_builtin_map_would_miss(tmp_path: Path):
    """Liver / myeloma / leukemia / thyroid are absent from the built-in
    ICD->phecodeX map. Those people must still be kept out of the controls."""
    prefixes = load_cancer_code_prefixes(_case_file(tmp_path))

    enc = tmp_path / "enc.tsv"
    enc.write_text(
        "sample_name|dx_code|dx_code_type|encounter_date\n"
        "TEST_100|C22.1|ICD10|01/01/2020\n"        # liver — builtin misses
        "TEST_101|C90.01|ICD10|01/01/2020\n"       # myeloma — builtin misses
        "TEST_102|C91.10|ICD10|01/01/2020\n"       # leukemia — builtin misses
        "TEST_103|C73|ICD10|01/01/2020\n"          # thyroid — builtin misses
        "TEST_104|C18.9|ICD10|01/01/2020\n"        # colorectal — builtin covers
        "TEST_105|I10|ICD10|01/01/2020\n"          # hypertension — not cancer
        "TEST_106|Z00.00|ICD10|01/01/2020\n"       # well visit — not cancer
    )
    schema = {"person_id_col": "sample_name", "icd_col": "dx_code", "sep": "|"}
    flagged = people_with_cancer_codes([(enc, schema)], prefixes)
    assert flagged == {"TEST_100", "TEST_101", "TEST_102", "TEST_103", "TEST_104"}
    assert "TEST_105" not in flagged
    assert "TEST_106" not in flagged


def test_prefix_matching_absorbs_subcode_variation(tmp_path: Path):
    """Case file has C18.7; the roster may record C18.70 or bare C18."""
    prefixes = load_cancer_code_prefixes(_case_file(tmp_path))
    enc = tmp_path / "enc2.tsv"
    enc.write_text(
        "sample_name|dx_code|sep_pad\n"
        "TEST_200|C18|x\n"
        "TEST_201|C18.70|x\n"
        "TEST_202|C1|x\n"                          # too short to match
    )
    schema = {"person_id_col": "sample_name", "icd_col": "dx_code", "sep": "|"}
    flagged = people_with_cancer_codes([(enc, schema)], prefixes)
    assert flagged == {"TEST_200", "TEST_201"}


def test_empty_prefix_set_flags_nobody(tmp_path: Path):
    """A missing case file must not silently flag everyone (or nobody
    unnoticed) — it yields an empty set, and the caller's counts show it."""
    assert load_cancer_code_prefixes(tmp_path / "nope.tsv") == set()
    enc = tmp_path / "enc3.tsv"
    enc.write_text("sample_name|dx_code\nTEST_1|C18.7\n")
    schema = {"person_id_col": "sample_name", "icd_col": "dx_code", "sep": "|"}
    assert people_with_cancer_codes([(enc, schema)], set()) == set()
