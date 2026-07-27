"""Tests for results-* discovery and ranking.

Builds fake NF output trees (synthetic IDs only) and asserts:
  - batch vs pilot classification from chromosome coverage
  - latest complete batch wins the auto-selection
  - tolerance of missing / malformed manifests and carrier matrices
  - n_samples is NEVER backfilled with the carrier count
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.utils.discover_results import (
    classify,
    discover,
    inspect,
    normalize_chrom,
    parse_batch,
    rank_batches,
    read_manifest,
    scan_carrier_matrix,
    select_best,
    summarize_chroms,
    to_yaml_dict,
)

MANIFEST = {
    "reference_build": "GRCh38",
    "vep_version": "112",
    "clinvar_release": "2026-07-16",
    "gnomad_version": "v4.1",
    "pipeline_sha": "abc1234",
    "tool_versions": {"bcftools": "1.19"},
}


def _write_carriers(d: Path, chroms: list[str], people: list[str]) -> None:
    """One carrier row per (chrom, person)."""
    (d / "carriers").mkdir(parents=True, exist_ok=True)
    lines = ["chr\tpos\tref\talt\tgene\tperson_id\tis_clinvar_PLP\tis_acmg_PLP\tis_AM_PLP"]
    pos = 1000
    for c in chroms:
        for p in people:
            pos += 10
            lines.append(f"{c}\t{pos}\tA\tG\tGENE{c}\t{p}\t1\t0\t0")
    (d / "carriers" / "carrier_matrix.tsv").write_text("\n".join(lines) + "\n")


def _write_variants(d: Path) -> None:
    (d / "variants").mkdir(parents=True, exist_ok=True)
    for f in ["clinvar_plp.tsv", "acmg_plp.tsv", "am_plp.tsv", "qc_per_gene.tsv"]:
        (d / "variants" / f).write_text("chr\tpos\tref\talt\tgene\n")


def _make_run(root: Path, name: str, chroms: list[str],
              people: list[str] | None = None,
              manifest: dict | None = MANIFEST,
              mtime_offset: int = 0) -> Path:
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    _write_carriers(d, chroms, people or ["TEST_001", "TEST_002", "TEST_003"])
    _write_variants(d)
    if manifest is not None:
        (d / "manifest.json").write_text(json.dumps(manifest))
        if mtime_offset:
            ts = time.time() + mtime_offset
            os.utime(d / "manifest.json", (ts, ts))
    return d


def _full_autosomes() -> list[str]:
    return [str(i) for i in range(1, 23)]


# ---------- small helpers ----------

def test_normalize_chrom():
    assert normalize_chrom("chr1") == "1"
    assert normalize_chrom("1") == "1"
    assert normalize_chrom("chrX") == "X"
    assert normalize_chrom("CHRY") == "Y"
    assert normalize_chrom("") == ""


def test_summarize_chroms_collapses_ranges():
    assert summarize_chroms(["1", "2", "3", "X"]) == "1-3,X"
    assert summarize_chroms(["1", "3", "5"]) == "1,3,5"
    assert summarize_chroms(["17"]) == "17"
    assert summarize_chroms([]) == ""


def test_parse_batch():
    assert parse_batch("results-batch2") == ("batch2", 2)
    assert parse_batch("results-chr17") == ("chr17", None)
    assert parse_batch("results-pilot") == ("pilot", None)


def test_classify_prefers_chrom_coverage_over_name():
    # Coverage wins: a multi-chrom run is a batch even with an odd folder name.
    assert classify("results-pilot", ["1", "2", "3"]) == "batch"
    # Single chromosome -> pilot.
    assert classify("results-batch9", ["17"]) == "pilot"
    # No coverage info -> fall back to the folder name.
    assert classify("results-chr17", []) == "pilot"
    assert classify("results-chr11-single", []) == "pilot"
    assert classify("results-batch1", []) == "batch"


# ---------- manifest tolerance ----------

def test_read_manifest_flattens_nested_and_tolerates_junk(tmp_path: Path):
    p = tmp_path / "manifest.json"
    p.write_text(json.dumps(MANIFEST))
    flat, err = read_manifest(p)
    assert err is None
    assert flat["reference_build"] == "GRCh38"
    assert flat["tool_versions.bcftools"] == "1.19"      # flattened one level

    missing, err = read_manifest(tmp_path / "nope.json")
    assert missing == {} and "missing" in err

    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    flat2, err2 = read_manifest(bad)
    assert flat2 == {} and "unreadable" in err2

    notobj = tmp_path / "arr.json"
    notobj.write_text("[1,2,3]")
    flat3, err3 = read_manifest(notobj)
    assert flat3 == {} and "not a JSON object" in err3


# ---------- carrier matrix scan ----------

def test_scan_carrier_matrix_counts(tmp_path: Path):
    d = _make_run(tmp_path, "results-batch1", ["1", "2"],
                  people=["TEST_001", "TEST_002"])
    cm = scan_carrier_matrix(d / "carriers" / "carrier_matrix.tsv")
    assert cm["error"] is None
    assert cm["chrs"] == ["1", "2"]
    assert cm["n_carriers"] == 2          # distinct people
    assert cm["n_rows"] == 4              # 2 chroms x 2 people
    assert cm["n_variants"] == 4          # distinct (chr,pos,ref,alt)


def test_scan_carrier_matrix_unexpected_columns_reports(tmp_path: Path):
    d = tmp_path / "results-weird" / "carriers"
    d.mkdir(parents=True)
    (d / "carrier_matrix.tsv").write_text("foo\tbar\n1\t2\n")
    cm = scan_carrier_matrix(d / "carrier_matrix.tsv")
    assert cm["error"] is not None and "unexpected columns" in cm["error"]


# ---------- the n_samples guard ----------

def test_n_samples_is_not_the_carrier_count(tmp_path: Path):
    """Regression guard: distinct person_id in the carrier matrix is the number
    of CARRIERS. Using it as cohort size would understate every rate."""
    d = _make_run(tmp_path, "results-batch1", _full_autosomes(),
                  people=["TEST_001", "TEST_002", "TEST_003"])
    info = inspect(d)
    assert info.n_carriers == 3
    assert info.n_samples is None                        # left unknown, not 3
    assert any("NOT cohort size" in n for n in info.notes)


def test_n_samples_read_from_explicit_sample_list(tmp_path: Path):
    d = _make_run(tmp_path, "results-batch1", _full_autosomes(),
                  people=["TEST_001", "TEST_002"])
    (d / "qc").mkdir()
    (d / "qc" / "sample_keep_list.tsv").write_text(
        "sample_id\nTEST_001\nTEST_002\nTEST_003\nTEST_004\n")
    info = inspect(d)
    assert info.n_carriers == 2
    assert info.n_samples == 4                           # header dropped
    assert any("n_samples from qc/sample_keep_list.tsv" in n for n in info.notes)


# ---------- completeness ----------

def test_complete_requires_all_files_and_full_autosomes(tmp_path: Path):
    full = _make_run(tmp_path, "results-batch1", _full_autosomes())
    assert inspect(full).complete is True

    partial = _make_run(tmp_path, "results-batch2", ["1", "2", "3"])
    pi = inspect(partial)
    assert pi.complete is False
    assert any("autosomes not fully covered" in n for n in pi.notes)

    nomanifest = _make_run(tmp_path, "results-batch3", _full_autosomes(), manifest=None)
    ni = inspect(nomanifest)
    assert ni.complete is False
    assert "manifest.json" in ni.missing_files


# ---------- ranking / selection ----------

def test_latest_complete_batch_is_selected(tmp_path: Path):
    # batch1 complete but older; batch3 complete and newer; chr17 is a pilot.
    _make_run(tmp_path, "results-batch1", _full_autosomes(), mtime_offset=-86400 * 10)
    _make_run(tmp_path, "results-batch3", _full_autosomes(), mtime_offset=-60)
    _make_run(tmp_path, "results-chr17", ["17"])

    infos = discover(tmp_path)
    assert len(infos) == 3
    best = select_best(infos)
    assert best is not None
    assert best.folder == "results-batch3"
    # Pilots are never selected.
    assert all(i.kind == "pilot" for i in infos if i.folder == "results-chr17")


def test_complete_batch_beats_incomplete_newer_one(tmp_path: Path):
    """A newer but incomplete run must not outrank a complete one."""
    _make_run(tmp_path, "results-batch1", _full_autosomes(), mtime_offset=-86400)
    _make_run(tmp_path, "results-batch2", ["1", "2"], mtime_offset=0)   # newer, partial
    infos = discover(tmp_path)
    best = select_best(infos)
    assert best.folder == "results-batch1"


def test_ranked_batches_exclude_pilots(tmp_path: Path):
    _make_run(tmp_path, "results-batch1", _full_autosomes())
    _make_run(tmp_path, "results-chr13", ["13"])
    _make_run(tmp_path, "results-chr11-single", ["11"])
    ranked = rank_batches(discover(tmp_path))
    assert [r.folder for r in ranked] == ["results-batch1"]


def test_no_batch_runs_yields_no_selection(tmp_path: Path):
    _make_run(tmp_path, "results-chr17", ["17"])
    infos = discover(tmp_path)
    assert select_best(infos) is None


# ---------- yaml payload ----------

def test_yaml_payload_records_manifest_keys_and_warning(tmp_path: Path):
    _make_run(tmp_path, "results-batch1", _full_autosomes())
    infos = discover(tmp_path)
    payload = to_yaml_dict(infos, select_best(infos), tmp_path)
    assert payload["selected"]["folder"] == "results-batch1"
    assert payload["selected"]["clinvar_release"] == "2026-07-16"
    assert payload["selected"]["vep_version"] == "112"
    # Real manifest schema is captured so we learn it from the first real run.
    assert "reference_build" in payload["manifest_keys_seen"]
    assert "tool_versions.bcftools" in payload["manifest_keys_seen"]
    assert "NOT the cohort size" in payload["_warning"]
