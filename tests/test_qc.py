"""Unit tests for QC helpers (pure Python; no data touched)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.qc.sample_qc import (
    SampleQcRow,
    build_duplicate_clusters,
    het_hom_outlier,
    is_duplicate,
    missingness_fail,
    pick_from_duplicate_cluster,
    sex_from_f_stat,
    sex_mismatch,
)
from src.qc.site_gt_qc import (
    genotype_mask_include_expr,
    site_drop_expr,
    site_filter_expr,
)
from src.qc.vcf_stats import PSC, merge_psc, parse_psc_lines


def _row(sid: str, cr: float) -> SampleQcRow:
    return SampleQcRow(sid, cr, None, None, None, False, "", True, "")


# --- site / GT QC expression builders ---------------------------------------

def test_site_gt_expr_pieces():
    assert site_filter_expr("PASS") == "PASS"
    expr = genotype_mask_include_expr(10, 20, 0.20, 0.80)
    assert "FMT/DP < 10" in expr
    assert "FMT/GQ < 20" in expr
    assert "0.2" in expr and "0.8" in expr
    drop = site_drop_expr(0.10, drop_monoallelic=True)
    assert "F_MISSING > 0.1" in drop
    assert "AC == 0" in drop and "AC == AN" in drop
    drop_off = site_drop_expr(0.10, drop_monoallelic=False)
    assert "AC == 0" not in drop_off


# --- sample QC building blocks ---------------------------------------------

def test_missingness_threshold():
    assert missingness_fail(0.06, 0.05) is True
    assert missingness_fail(0.04, 0.05) is False


def test_sex_mismatch_only_when_both_known_and_differ():
    assert sex_mismatch("M", "F") is True
    assert sex_mismatch("Female", "Male") is True
    assert sex_mismatch("M", "M") is False
    assert sex_mismatch(None, "M") is False


def test_sex_from_f_stat_zones():
    assert sex_from_f_stat(0.10, 0.20, 0.80) == "F"
    assert sex_from_f_stat(0.90, 0.20, 0.80) == "M"
    assert sex_from_f_stat(0.50, 0.20, 0.80) is None


def test_het_hom_outlier_z_threshold():
    # Sample well within cohort -> not an outlier.
    assert het_hom_outlier(1.6, 1.5, 0.1, 4.0) is False
    # 5 SDs above the mean -> outlier.
    assert het_hom_outlier(2.0, 1.5, 0.1, 4.0) is True
    # sd=0 (degenerate cohort) -> never outlier.
    assert het_hom_outlier(2.0, 1.5, 0.0, 4.0) is False


def test_kinship_duplicate_call():
    assert is_duplicate(0.45, 0.354) is True
    assert is_duplicate(0.35, 0.354) is False


def test_duplicate_cluster_keeps_highest_call_rate():
    cluster = [_row("A", 0.90), _row("B", 0.99), _row("C", 0.95)]
    assert pick_from_duplicate_cluster(cluster, "highest_call_rate") == "B"


def test_build_duplicate_clusters_union_find():
    # A-B are dups, B-C are dups, A-C would be dups by transitivity -> one cluster.
    pairs = [("A", "B", 0.45), ("B", "C", 0.40), ("D", "E", 0.10)]
    clusters = build_duplicate_clusters(pairs, duplicate_min=0.354)
    assert len(clusters) == 1
    assert clusters[0] == {"A", "B", "C"}


# --- PSC parsing / merging --------------------------------------------------

_STATS = """\
# PSC, Per-sample counts
# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTs\t[8]nTv\t[9]nIndels\t[10]avgDP\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing
PSC\t0\tS1\t1000\t100\t150\t0\t0\t0\t35.0\t20\t0\t0\t50
PSC\t0\tS2\t900\t80\t120\t0\t0\t0\t32.0\t15\t0\t0\t500
"""


def test_parse_psc_and_derived():
    rows = parse_psc_lines(_STATS)
    assert len(rows) == 2
    s1 = next(r for r in rows if r.sample == "S1")
    assert s1.n_called == 1000 + 100 + 150
    assert 0.0 < s1.missingness < 1.0
    # S2 has 10x more missing -> higher missingness
    s2 = next(r for r in rows if r.sample == "S2")
    assert s2.missingness > s1.missingness


def test_merge_psc_sums_across_chunks():
    a = PSC("S1", 1000, 100, 150, 20, 50)
    b = PSC("S1", 500, 60, 80, 10, 25)
    agg = merge_psc([a, b])
    assert agg["S1"].n_ref_hom == 1500
    assert agg["S1"].n_missing == 75
    assert agg["S1"].n_singletons == 30
