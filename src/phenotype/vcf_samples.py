"""Load the WES sample list (produced by `bcftools query -l`).

Used to intersect the EHR-derived phenotype roster with samples that actually
have genotype data. Missing file -> empty set + warning; the pipeline still
runs but the "with WES" count columns will all be zero.
"""
from __future__ import annotations

from pathlib import Path


def load_wes_samples(cfg: dict) -> set[str]:
    """Return a set of sample IDs present in the WES data. Empty if unavailable."""
    path = (cfg.get("inputs", {}).get("wes", {}) or {}).get("sample_list_file")
    if not path:
        return set()
    p = Path(path)
    if not p.exists():
        return set()
    out: set[str] = set()
    with p.open() as fh:
        for line in fh:
            sid = line.strip()
            if sid:
                out.add(sid)
    return out
