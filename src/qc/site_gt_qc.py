"""Site + genotype QC command builders.

Called by workflow/01_norm_vep_chunk.sh between `bcftools norm` and `vep`.
The heavy lifting is bcftools; this module just centralizes the expressions
and thresholds so tests can pin them.
"""
from __future__ import annotations


def site_filter_expr(keep_filter: str) -> str:
    """`bcftools view -f <expr>` — keep only VQSR/GATK-PASS sites."""
    return keep_filter or "PASS"


def genotype_mask_include_expr(dp_min: int, gq_min: int,
                               ab_het_min: float, ab_het_max: float) -> str:
    """`bcftools +setGT` -include expression: mark genotypes to be set to `./.`.

    Nulls genotypes when ANY of:
      - FMT/DP < dp_min
      - FMT/GQ < gq_min
      - het call with allele balance outside [ab_het_min, ab_het_max]
    """
    return (
        f"FMT/DP < {dp_min}"
        f" | FMT/GQ < {gq_min}"
        f" | (GT=\"het\" & (FMT/AD[:1]/(FMT/AD[:0]+FMT/AD[:1])) < {ab_het_min})"
        f" | (GT=\"het\" & (FMT/AD[:1]/(FMT/AD[:0]+FMT/AD[:1])) > {ab_het_max})"
    )


def site_drop_expr(site_missing_max: float, drop_monoallelic: bool) -> str:
    """`bcftools view -e <expr>` — drop sites after genotype masking."""
    parts = [f"F_MISSING > {site_missing_max}"]
    if drop_monoallelic:
        parts += ["AC == 0", "AC == AN"]     # mono-ref or mono-alt after masking
    return " | ".join(parts)
