"""Unit tests for the pure-Python classification rules (no data touched)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.classification.clinvar import (
    clinvar_stars,
    is_clinvar_plp,
    is_clinvar_blb,
)
from src.classification.alphamissense import call_carrier, promote


PLP_TERMS = ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"]
BLB_TERMS = ["Benign", "Likely_benign", "Benign/Likely_benign"]


def test_clinvar_stars_map():
    assert clinvar_stars("reviewed_by_expert_panel") == 3
    assert clinvar_stars("criteria_provided,_multiple_submitters,_no_conflicts") == 2
    assert clinvar_stars("no_assertion_criteria_provided") == 0
    assert clinvar_stars(None) == 0


def test_is_clinvar_plp_requires_stars_and_no_conflict():
    ok = is_clinvar_plp("Pathogenic",
                       "criteria_provided,_multiple_submitters,_no_conflicts",
                       PLP_TERMS, min_stars=2)
    assert ok is True
    # Conflicting classifications disqualify.
    conflict = is_clinvar_plp("Pathogenic",
                              "criteria_provided,_conflicting_classifications",
                              PLP_TERMS, min_stars=2)
    assert conflict is False
    # Below star threshold.
    low = is_clinvar_plp("Pathogenic", "criteria_provided,_single_submitter",
                        PLP_TERMS, min_stars=2)
    assert low is False


def test_is_clinvar_blb_used_by_acmg_removal():
    assert is_clinvar_blb("Benign", BLB_TERMS) is True
    assert is_clinvar_blb("Likely_benign", BLB_TERMS) is True
    assert is_clinvar_blb("Pathogenic", BLB_TERMS) is False
    assert is_clinvar_blb(None, BLB_TERMS) is False


def test_alphamissense_gene_specific_promotion():
    # Fake per-gene calibration: BRCA1 with two cutoffs.
    gene_thresholds = {
        "BRCA1": {"cutoffs": [
            {"threshold": 0.564, "evidence_label": "PP3_Supporting"},
            {"threshold": 0.792, "evidence_label": "PP3_Moderate"},
            {"threshold": 0.972, "evidence_label": "PP3_Strong"},
        ]},
    }
    assert call_carrier(0.99, "BRCA1", gene_thresholds) == "PP3_Strong"
    assert call_carrier(0.80, "BRCA1", gene_thresholds) == "PP3_Moderate"
    assert call_carrier(0.60, "BRCA1", gene_thresholds) == "PP3_Supporting"
    assert call_carrier(0.10, "BRCA1", gene_thresholds) == "Indeterminate"
    # Unknown gene -> Indeterminate (no domain-aggregate promotion allowed).
    assert call_carrier(0.99, "MADE_UP", gene_thresholds) == "Indeterminate"
    # Promotion threshold: Moderate default -> Supporting does not promote.
    assert promote("PP3_Moderate") is True
    assert promote("PP3_Supporting") is False
