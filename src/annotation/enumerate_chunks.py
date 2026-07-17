"""Enumerate WES pVCF Target chunks for the LSF array.

Usage (on Minerva):
    python -m src.annotation.enumerate_chunks --config config/config.yaml \
        --out $OUTPUT_ROOT/logs/chunks.list

Emits one absolute path per line (LSF array index i -> line i, 1-based).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from src.common.config import load, resolve
from src.common.logging_utils import get_logger


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load(args.config)
    log = get_logger("enumerate_chunks")

    r = resolve(cfg)
    if not r.wes_chunks:
        log.error("No chunks found under %s", r.wes_target_dir)
        return 2

    expected = cfg["inputs"]["wes"].get("expected_chunks")
    if expected and abs(len(r.wes_chunks) - expected) > 5:
        log.warning("Chunk count %d differs from expected %d", len(r.wes_chunks), expected)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(str(p) for p in r.wes_chunks) + "\n")
    log.info("Wrote %d chunk paths to %s", len(r.wes_chunks), out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
