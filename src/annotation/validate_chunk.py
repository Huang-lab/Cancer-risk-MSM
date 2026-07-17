"""Validate an annotated per-chunk VCF.

Checks (called by 04_validate_chunk.sh and by the worker's skip logic):
  1. File exists, non-empty, tabix index present.
  2. bcftools view exits 0 on the file.
  3. Expected INFO fields present at least once (ClinVar, AlphaMissense, LOFTEE, gnomAD).
  4. Variant count within a reasonable ratio of the input chunk.
  5. Spot-check: at least one variant on a known cancer-predisposition gene appears.

Exits 0 on pass; non-zero with a short reason on failure.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REQUIRED_INFO = ["CLNSIG", "am_pathogenicity", "LoF", "gnomAD_AF"]
SPOTCHECK_GENES = {"BRCA1", "BRCA2", "MLH1", "APC", "TP53"}


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def _need(tool: str) -> None:
    if shutil.which(tool) is None:
        print(f"FAIL: {tool} not on PATH", file=sys.stderr)
        sys.exit(3)


def validate(input_vcf: Path, annotated_vcf: Path, min_ratio: float) -> tuple[bool, str]:
    if not annotated_vcf.exists() or annotated_vcf.stat().st_size == 0:
        return False, "annotated file missing or empty"
    tbi = annotated_vcf.with_suffix(annotated_vcf.suffix + ".tbi")
    if not tbi.exists():
        return False, "tabix index missing"

    _need("bcftools")
    r = _run(["bcftools", "view", "-h", str(annotated_vcf)])
    if r.returncode != 0:
        return False, f"bcftools view failed: {r.stderr.strip()[:200]}"

    header = r.stdout
    missing = [f for f in REQUIRED_INFO if f"ID={f}" not in header]
    if missing:
        return False, f"missing INFO fields: {missing}"

    # Variant counts (approximate)
    in_n = int(_run(["bcftools", "view", "-H", str(input_vcf)]).stdout.count("\n") or 0)
    out_n = int(_run(["bcftools", "view", "-H", str(annotated_vcf)]).stdout.count("\n") or 0)
    if in_n and out_n / max(in_n, 1) < min_ratio:
        return False, f"variant count too low ({out_n}/{in_n})"

    # Spot-check: any of our known genes appears in CSQ
    csq = _run(["bcftools", "view", "-H", str(annotated_vcf)])
    if not any(g in csq.stdout for g in SPOTCHECK_GENES):
        return False, "no spot-check genes appear in output"

    return True, f"OK ({out_n} variants, {in_n} input)"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--min-ratio", type=float, default=0.5)
    args = ap.parse_args()
    ok, msg = validate(Path(args.input), Path(args.annotated), args.min_ratio)
    print(msg)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
