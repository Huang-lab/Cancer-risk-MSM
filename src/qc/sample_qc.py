"""Sample-level QC — produces $OUTPUT_ROOT/results/qc/sample_keep_list.tsv.

Runs once, before annotation. Uses the biobank's existing common-variants
kinship + PCA + PLINK sex-check outputs (paths in config.inputs.common_variants
and config.qc.sample.*) so we do not have to re-compute anything expensive.

Checks
------
1. Missingness      -> per-sample variant missingness <= qc.sample.missingness_max.
2. Sex check        -> reported sex matches PLINK --check-sex inferred sex.
3. Contamination    -> VerifyBamID2 FREEMIX <= qc.sample.contamination.freemix_max.
4. Kinship          -> for pairs with kinship >= duplicate_kin_min, keep one
                       (highest call rate); flag but keep <2nd-degree relatives.

Outputs
-------
- results/qc/sample_keep_list.tsv        : one row per KEPT sample_id.
- results/qc/sample_qc_report.tsv        : per-sample flags (kept + dropped + reasons).
- results/qc/sample_qc_summary.md        : counts (never per-sample values).

PHI: sample IDs live in the tsvs, which are gitignored (results/**). Never
logged at INFO — only counts and reasons.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class SampleQcRow:
    sample_id: str
    call_rate: float | None
    sex_reported: str | None
    sex_inferred: str | None
    sex_mismatch: bool
    freemix: float | None
    kinship_flag: str          # "" | "duplicate_dropped" | "duplicate_kept" | "related_kept"
    kept: bool
    reasons: str               # comma-separated failed checks


# --- Building blocks --------------------------------------------------------


def sex_mismatch(reported: str | None, inferred: str | None) -> bool:
    if not reported or not inferred:
        return False
    r = reported.strip().upper()[:1]      # "M"/"F"
    i = inferred.strip().upper()[:1]
    return r in ("M", "F") and i in ("M", "F") and r != i


def contamination_fail(freemix: float | None, freemix_max: float) -> bool:
    return freemix is not None and freemix > freemix_max


def is_duplicate(kinship: float, duplicate_min: float) -> bool:
    return kinship >= duplicate_min


def pick_from_duplicate_cluster(cluster: list[SampleQcRow], strategy: str) -> str:
    """Return the sample_id to KEEP from a set of duplicate/MZ-twin samples."""
    if strategy == "highest_call_rate":
        return max(cluster, key=lambda r: (r.call_rate or 0.0)).sample_id
    return cluster[0].sample_id


# TODO(minerva): implement:
#   load_missingness(cfg) -> dict[sample_id -> call_rate]
#   load_sex_check(cfg)   -> dict[sample_id -> (reported, inferred)]
#   load_kinship(cfg)     -> list[(id_a, id_b, kinship)]
#   run(cfg) -> writes results/qc/sample_keep_list.tsv and reports.
