"""Unit tests for QC helpers (pure Python; no data touched)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.qc.sample_qc import (
    SampleQcRow,
    sex_mismatch,
    contamination_fail,
    is_duplicate,
    pick_from_duplicate_cluster,
)
from src.qc.site_gt_qc import (
    site_filter_expr,
    genotype_mask_include_expr,
    site_drop_expr,
)


def _row(sid: str, cr: float) -> SampleQcRow:
    return SampleQcRow(sid, cr, None, None, False, None, "", True, "")


def test_sex_mismatch_only_when_both_known_and_differ():
    assert sex_mismatch("M", "F") is True
    assert sex_mismatch("Female", "Male") is True
    assert sex_mismatch("M", "M") is False
    assert sex_mismatch(None, "M") is False
    assert sex_mismatch("U", "F") is False


def test_contamination_threshold():
    assert contamination_fail(0.05, 0.03) is True
    assert contamination_fail(0.02, 0.03) is False
    assert contamination_fail(None, 0.03) is False


def test_kinship_duplicate_call():
    assert is_duplicate(0.45, 0.354) is True
    assert is_duplicate(0.35, 0.354) is False


def test_duplicate_cluster_keeps_highest_call_rate():
    cluster = [_row("A", 0.90), _row("B", 0.99), _row("C", 0.95)]
    assert pick_from_duplicate_cluster(cluster, "highest_call_rate") == "B"


def test_site_gt_expr_pieces():
    assert site_filter_expr("PASS") == "PASS"
    expr = genotype_mask_include_expr(10, 20, 0.20, 0.80)
    assert "FMT/DP < 10" in expr
    assert "FMT/GQ < 20" in expr
    assert "0.2" in expr and "0.8" in expr
    drop = site_drop_expr(0.10, drop_monoallelic=True)
    assert "F_MISSING > 0.1" in drop
    assert "AC == 0" in drop and "AC == AN" in drop
    # mono-allelic drop can be disabled
    drop_off = site_drop_expr(0.10, drop_monoallelic=False)
    assert "AC == 0" not in drop_off
