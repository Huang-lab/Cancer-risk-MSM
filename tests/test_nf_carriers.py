"""Tests for reading germline-plp-carrier-nf outputs (the read-boundary)."""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.carriers.nf_carriers import (
    assert_consistent_builds,
    filter_to_keep_list,
    panel_carrier_flags,
    per_gene_flags,
    read_carriers,
    read_manifests,
    read_sample_keep_list,
)

FIXTURES = Path(__file__).parent / "synthetic"
PANEL = ["BRCA1", "BRCA2", "MLH1", "MSH2", "MUTYH", "PMS2", "TP53"]


def _cfg(root: Path, results_dirs: list[str]) -> dict:
    return {
        "inputs": {
            "carrier_source": {
                "nf_results_root": str(root),
                "results_dirs": results_dirs,
                "carrier_matrix": "carriers/carrier_matrix.tsv",
                "variants_dir": "variants",
                "manifest": "manifest.json",
                "columns": {
                    "person_id": "person_id", "chrom": "chr", "pos": "pos",
                    "ref": "ref", "alt": "alt", "gene": "gene",
                    "is_clinvar_plp": "is_clinvar_PLP",
                    "is_acmg_plp": "is_acmg_PLP",
                    "is_am_plp": "is_AM_PLP",
                },
                "sep": "\t",
            }
        }
    }


def _stage_nf_run(root: Path, name: str, build: str = "GRCh38") -> Path:
    d = root / name
    (d / "carriers").mkdir(parents=True)
    shutil.copy(FIXTURES / "nf_carrier_matrix.tsv", d / "carriers" / "carrier_matrix.tsv")
    (d / "manifest.json").write_text(json.dumps({
        "reference_build": build, "pipeline_sha": "abc123",
        "tool_versions": {"vep": "112", "bcftools": "1.19"},
    }))
    return d


def test_read_carriers_parses_long_format(tmp_path: Path):
    _stage_nf_run(tmp_path, "results-batch1")
    rows = read_carriers(_cfg(tmp_path, ["results-batch1"]))
    assert len(rows) == 6
    by_person = {}
    for r in rows:
        by_person.setdefault(r.person_id, []).append(r)
    # TEST_001 carries two variants
    assert len(by_person["TEST_001"]) == 2
    mutyh = next(r for r in by_person["TEST_001"] if r.gene == "MUTYH")
    assert mutyh.is_clinvar_plp is True
    assert mutyh.is_acmg_plp is True
    assert mutyh.is_am_plp is False
    assert mutyh.pos == 45331558


def test_pooling_two_runs_ors_framework_flags(tmp_path: Path):
    # Same variants in both runs, but the second marks MUTYH as AM-positive.
    _stage_nf_run(tmp_path, "results-batch1")
    d2 = _stage_nf_run(tmp_path, "results-batch2")
    p = d2 / "carriers" / "carrier_matrix.tsv"
    p.write_text(p.read_text().replace(
        "MUTYH\tTEST_001\t1\t1\t0", "MUTYH\tTEST_001\t0\t0\t1"))

    rows = read_carriers(_cfg(tmp_path, ["results-batch1", "results-batch2"]))
    # Deduplicated on (person, variant) — still 6 unique rows, not 12.
    assert len(rows) == 6
    mutyh = next(r for r in rows if r.gene == "MUTYH" and r.person_id == "TEST_001")
    # Flags OR'd across runs: clinvar/acmg from run1, AM from run2.
    assert mutyh.is_clinvar_plp is True
    assert mutyh.is_acmg_plp is True
    assert mutyh.is_am_plp is True


def test_missing_results_dir_is_skipped_not_fatal(tmp_path: Path):
    _stage_nf_run(tmp_path, "results-batch1")
    rows = read_carriers(_cfg(tmp_path, ["results-batch1", "results-nonexistent"]))
    assert len(rows) == 6


def test_assert_consistent_builds_rejects_mismatch(tmp_path: Path):
    _stage_nf_run(tmp_path, "results-batch1", build="GRCh38")
    _stage_nf_run(tmp_path, "results-batch2", build="GRCh37")
    manifests = read_manifests(_cfg(tmp_path, ["results-batch1", "results-batch2"]))
    assert len(manifests) == 2
    with pytest.raises(ValueError, match="disagree on reference build"):
        assert_consistent_builds(manifests)


def test_assert_consistent_builds_accepts_match(tmp_path: Path):
    _stage_nf_run(tmp_path, "results-batch1")
    _stage_nf_run(tmp_path, "results-batch2")
    manifests = read_manifests(_cfg(tmp_path, ["results-batch1", "results-batch2"]))
    assert assert_consistent_builds(manifests) == "GRCh38"


def test_keep_list_filter(tmp_path: Path):
    _stage_nf_run(tmp_path, "results-batch1")
    rows = read_carriers(_cfg(tmp_path, ["results-batch1"]))
    kept = filter_to_keep_list(rows, {"TEST_001", "TEST_002"})
    assert {r.person_id for r in kept} == {"TEST_001", "TEST_002"}
    # Empty keep-list is a no-op (QC not run yet), not a wipe.
    assert len(filter_to_keep_list(rows, set())) == len(rows)


def test_read_sample_keep_list_with_header(tmp_path: Path):
    d = _stage_nf_run(tmp_path, "results-batch1")
    (d / "qc").mkdir()
    (d / "qc" / "sample_keep_list.tsv").write_text(
        "sample_id\tcall_rate\nTEST_001\t0.99\nTEST_002\t0.98\n")
    cfg = _cfg(tmp_path, ["results-batch1"])
    cfg["inputs"]["carrier_source"]["sample_keep_list"] = "qc/sample_keep_list.tsv"
    assert read_sample_keep_list(cfg) == {"TEST_001", "TEST_002"}


def test_read_sample_keep_list_headerless(tmp_path: Path):
    d = _stage_nf_run(tmp_path, "results-batch1")
    (d / "qc").mkdir()
    (d / "qc" / "keep.txt").write_text("TEST_001\nTEST_003\n")
    cfg = _cfg(tmp_path, ["results-batch1"])
    cfg["inputs"]["carrier_source"]["sample_keep_list"] = "qc/keep.txt"
    assert read_sample_keep_list(cfg) == {"TEST_001", "TEST_003"}


def test_read_sample_keep_list_unconfigured_is_empty(tmp_path: Path):
    _stage_nf_run(tmp_path, "results-batch1")
    cfg = _cfg(tmp_path, ["results-batch1"])   # sample_keep_list not set
    keep = read_sample_keep_list(cfg)
    assert keep == set()
    # ...and an empty keep-set must not wipe the carrier rows.
    rows = read_carriers(cfg)
    assert len(filter_to_keep_list(rows, keep)) == len(rows)


def test_per_gene_and_panel_flags(tmp_path: Path):
    _stage_nf_run(tmp_path, "results-batch1")
    rows = read_carriers(_cfg(tmp_path, ["results-batch1"]))

    clinvar = per_gene_flags(rows, "clinvar")
    assert clinvar["TEST_001"] == {"MUTYH"}          # MSH2 row is AM-only
    assert clinvar["TEST_002"] == {"BRCA1"}
    assert "TEST_004" not in clinvar                 # PMS2 row is AM-only

    am = per_gene_flags(rows, "am")
    assert am["TEST_001"] == {"MSH2"}
    assert am["TEST_004"] == {"PMS2"}

    panel = panel_carrier_flags(rows, "acmg", PANEL)
    assert panel["TEST_001"] is True                 # MUTYH acmg+
    assert panel["TEST_003"] is True                 # BRCA2 acmg+
    assert panel.get("TEST_006", False) is False     # TP53 is clinvar-only
