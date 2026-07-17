"""Build the per-person carrier table in LONG format.

Per user requirement:
    One row per (variant x person). Columns:
        chr  pos  ref  alt  gene  person_id  is_clinvar_PLP  is_acmg_PLP  is_AM_PLP
    Carrier = at least one qualifying alternate allele (0/1 or 1/1).

Pipeline
--------
1. Build a `qualifying_variants` set from the three classification outputs:
     - results/variants/clinvar_plp.tsv   -> is_clinvar_PLP flags
     - results/variants/acmg_plp.tsv      -> is_acmg_PLP flags (already ClinVar-B/LB-filtered)
     - results/variants/am_plp.tsv        -> is_AM_PLP flags
   Union of variants across the three frameworks; per-framework booleans on each row.
2. For each qualifying variant, run:
       bcftools query -r "$CHR:$POS-$POS" -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' <annotated.vcf.gz>
   parse genotypes, emit one row per (variant, carrier-person).
3. Write long TSV to results/carriers/carriers.long.tsv.
4. Optional wide pivot (config.carriers.emit_wide_pivot) -> sample x gene matrix
   with per-framework aggregation (e.g. is_clinvar_PLP_gene = any is_clinvar_PLP on that gene).

Scaffolded; runs on Minerva after classification.
"""
from __future__ import annotations

from dataclasses import dataclass

LONG_COLUMNS = [
    "chr", "pos", "ref", "alt", "gene", "person_id",
    "is_clinvar_PLP", "is_acmg_PLP", "is_AM_PLP",
]


@dataclass
class QualifyingVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    is_clinvar_plp: bool
    is_acmg_plp: bool
    is_am_plp: bool


def bcftools_query_cmd(annotated_vcf: str, chrom: str, pos: int) -> list[str]:
    """Emit sample-level genotypes at a single position."""
    return [
        "bcftools", "query",
        "-r", f"{chrom}:{pos}-{pos}",
        "-f", r"%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n",
        annotated_vcf,
    ]


def is_carrier(gt: str) -> bool:
    """0/0 / ./. -> False; anything with an alt allele -> True."""
    if not gt or gt in (".", "./.", ".|."):
        return False
    alleles = gt.replace("|", "/").split("/")
    return any(a not in ("0", ".") for a in alleles)


# TODO(minerva): implement:
#   union_qualifying_variants(cfg) -> Iterator[QualifyingVariant]
#   emit_long_table(cfg) -> Path(results/carriers/carriers.long.tsv)
#   emit_wide_pivot(long_df, cfg) -> Path(results/carriers/carriers.wide.tsv)
