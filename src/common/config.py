"""Config loader with path auto-resolution (latest batch / date / release).

No participant data is read here — only directory listings on Minerva.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass
class Resolved:
    """Concrete Minerva paths after auto-resolution.

    Note: reference resources for annotation (VEP cache, ClinVar, gnomAD,
    AlphaMissense) are NOT resolved here — germline-plp-carrier-nf owns
    annotation. We only resolve inputs this repo actually reads.
    """

    wes_target_dir: Path
    wes_chunks: list[Path]
    common_variants_dir: Path
    phenotypes_dir: Path


def load(config_path: str | os.PathLike) -> dict[str, Any]:
    with open(config_path) as fh:
        return yaml.safe_load(fh)


# --- Auto-resolvers ---------------------------------------------------------


_BATCH_RE = re.compile(r"batch_(\d+)$")
_DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def _newest_batch(parent_glob: Path) -> Path:
    parent = parent_glob.parent
    pattern = parent_glob.name
    candidates = [p for p in parent.glob(pattern) if _BATCH_RE.search(p.name)]
    if not candidates:
        raise FileNotFoundError(f"No batch_### directory matching {parent_glob}")
    return max(candidates, key=lambda p: int(_BATCH_RE.search(p.name).group(1)))


def _newest_dated_dir(parent: Path) -> Path:
    candidates = [p for p in parent.iterdir() if p.is_dir() and _DATE_RE.match(p.name)]
    if not candidates:
        raise FileNotFoundError(f"No yyyy-mm-dd directory in {parent}")
    return max(candidates, key=lambda p: p.name)


def resolve(cfg: dict) -> Resolved:
    """Turn config globs into concrete Minerva paths.

    Call this on-cluster; from a non-Minerva host, paths will not exist.
    """
    root = Path(cfg["inputs"]["msm_data_root"])

    # WES pVCF Target chunks
    wes_glob = root / cfg["inputs"]["wes"]["pvcf_target_glob"]
    if cfg["inputs"]["wes"].get("batch"):
        wes_target = root / cfg["inputs"]["wes"]["pvcf_target_glob"].replace(
            "batch_*", cfg["inputs"]["wes"]["batch"]
        )
    else:
        wes_target = _newest_batch(wes_glob)
    chunks = sorted(wes_target.glob(cfg["inputs"]["wes"]["chunk_glob"]))

    # Common variants
    cv_glob = root / cfg["inputs"]["common_variants"]["root_glob"]
    cv_dir = (
        root / cv_glob.name.replace("batch_*", cfg["inputs"]["common_variants"]["batch"])
        if cfg["inputs"]["common_variants"].get("batch")
        else _newest_batch(cv_glob)
    )

    # Phenotypes (latest dated folder)
    ph_root = root / "phenotypes"
    ph_dir = (
        ph_root / cfg["inputs"]["phenotypes"]["date"]
        if cfg["inputs"]["phenotypes"].get("date")
        else _newest_dated_dir(ph_root)
    )

    return Resolved(
        wes_target_dir=wes_target,
        wes_chunks=chunks,
        common_variants_dir=cv_dir,
        phenotypes_dir=ph_dir,
    )
