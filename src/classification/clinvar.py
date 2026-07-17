"""ClinVar P/LP classification.

Rule: variant is ClinVar P/LP if CLNSIG matches config plp_terms AND CLNREVSTAT
implies >= min_review_stars gold stars (default 2 -> "criteria_provided_multiple_submitters"
or better, no conflicts).

Reads per-chunk annotated VCFs; writes a variant-level table under
$OUTPUT_ROOT/results/variants/clinvar_plp.tsv.

Scaffolded today; runs on Minerva after annotation completes.
"""
from __future__ import annotations

from typing import Iterable

# Gold-star mapping: CLNREVSTAT -> stars. Reference:
# https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
_STAR_MAP = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,  # conflicts -> not clean P/LP
    "criteria_provided,_single_submitter": 1,
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_classification_provided": 0,
    "no_classifications_from_unflagged_records": 0,
}


def clinvar_stars(clnrevstat: str | None) -> int:
    if not clnrevstat:
        return 0
    return _STAR_MAP.get(clnrevstat.strip().lower(), 0)


def is_clinvar_plp(clnsig: str | None, clnrevstat: str | None,
                   plp_terms: Iterable[str], min_stars: int = 2) -> bool:
    if not clnsig:
        return False
    tokens = {t.strip() for t in clnsig.replace("|", ",").split(",")}
    if not any(t in tokens for t in plp_terms):
        return False
    # A "conflicting classifications" review-status disqualifies even if a P/LP term is present.
    if clnrevstat and "conflicting" in clnrevstat.lower():
        return False
    return clinvar_stars(clnrevstat) >= min_stars


def is_clinvar_blb(clnsig: str | None, blb_terms: Iterable[str]) -> bool:
    """Used by ACMG QC rule to remove ACMG-P/LP variants ClinVar calls B/LB."""
    if not clnsig:
        return False
    tokens = {t.strip() for t in clnsig.replace("|", ",").split(",")}
    return any(t in tokens for t in blb_terms)


# TODO(minerva): implement classify_chunk(annot_vcf, cfg) -> DataFrame with
#   chr, pos, ref, alt, gene, csq, clnsig, clnrevstat, stars, is_clinvar_PLP
# then classify_all(cfg) which concatenates per-chunk tables and writes
# results/variants/clinvar_plp.tsv (long, one row per variant).
