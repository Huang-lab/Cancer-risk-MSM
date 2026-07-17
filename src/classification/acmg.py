"""ACMG/AMP P/LP classification via ANNOVAR + InterVar.

Approach (per user request; mirrors the BioMe ACMG track):
  1. Convert per-chunk annotated VCFs -> ANNOVAR input (avinput) via
     convert2annovar.pl, then annotate with table_annovar.pl using humandb.
  2. Run InterVar (autoInterVar) over the ANNOVAR output to emit
     rule-based ACMG P/LP calls (PVS1 hc-pLoF in LoF-intolerant genes,
     PM2 gnomAD rarity via popmax AF <= pm2_af_max, PP3 from calibrated
     in-silico, plus the standard InterVar rule set).
  3. Post-process:
       - Drop variants that ClinVar classifies B/LB (config: remove_if_clinvar_blb).
       - Emit a variant-level TSV: chr, pos, ref, alt, gene, is_acmg_PLP, rules_fired.

Scaffolded today; runs on Minerva after annotation.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class AcmgCall:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    is_plp: bool
    rules: str          # e.g. "PVS1,PM2,PP3"
    clinvar_clnsig: str | None = None
    clinvar_removed: bool = False


# --- Command builders (executed on Minerva) ---------------------------------

def annovar_convert_cmd(annotated_vcf: str, avinput: str, annovar_root: str) -> list[str]:
    return [
        f"{annovar_root}/convert2annovar.pl",
        "-format", "vcf4",
        "-allsample", "-withfreq",
        annotated_vcf,
        "-outfile", avinput,
    ]


def annovar_table_cmd(avinput: str, out_prefix: str, humandb: str,
                      buildver: str = "hg38") -> list[str]:
    protocols = "refGene,gnomad41_exome,clinvar_latest,dbnsfp47a"
    ops = "g,f,f,f"
    return [
        "table_annovar.pl", avinput, humandb,
        "-buildver", buildver,
        "-out", out_prefix,
        "-remove",
        "-protocol", protocols,
        "-operation", ops,
        "-nastring", ".",
        "-vcfinput",
    ]


def intervar_cmd(annovar_txt: str, out_prefix: str, intervar_root: str,
                 intervar_config: str, buildver: str = "hg38") -> list[str]:
    return [
        f"{intervar_root}/Intervar.py",
        "-b", buildver,
        "-i", annovar_txt,
        "--input_type=AVinput",
        "-o", out_prefix,
        "-c", intervar_config,
    ]


# --- Post-processing --------------------------------------------------------

def apply_clinvar_blb_removal(calls: list[AcmgCall],
                              clinvar_lookup,
                              blb_terms: list[str]) -> list[AcmgCall]:
    """Per user rule: remove ACMG-P/LP variants that ClinVar classifies B/LB.

    `clinvar_lookup(chrom, pos, ref, alt) -> Optional[str]` returns CLNSIG or None.
    """
    from .clinvar import is_clinvar_blb
    out: list[AcmgCall] = []
    for c in calls:
        clnsig = clinvar_lookup(c.chrom, c.pos, c.ref, c.alt)
        c.clinvar_clnsig = clnsig
        if c.is_plp and is_clinvar_blb(clnsig, blb_terms):
            c.is_plp = False
            c.clinvar_removed = True
        out.append(c)
    return out


# TODO(minerva): parse InterVar output into AcmgCall list; wire clinvar_lookup to
# the annotated-chunk ClinVar CSQ; emit results/variants/acmg_plp.tsv.
