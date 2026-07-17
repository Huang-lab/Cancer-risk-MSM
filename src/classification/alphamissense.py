"""AlphaMissense-calibrated P/LP classification (gene-specific).

Mirrors Huang-lab/Biome-alphamissense-calibration:
- Use Chen/Pejaver per-variant `evidence` labels rather than raw scores.
- Primary call uses GENE-SPECIFIC calibrated thresholds only.
- Domain-aggregate thresholds are RECORDED but NEVER promote a carrier.
- Carrier requires evidence >= configured min strength (default Moderate).

Scaffolded today; runs on Minerva.
"""
from __future__ import annotations

from dataclasses import dataclass

# Evidence-strength ranks (Chen/Pejaver ACMG calibration).
_RANK = {
    "BP4_Very_Strong": -4, "BP4_Strong": -3, "BP4_Moderate": -2, "BP4_Supporting": -1,
    "Indeterminate": 0,
    "PP3_Supporting": 1, "PP3_Moderate": 2, "PP3_Strong": 3, "PP3_Very_Strong": 4,
}


@dataclass
class CalibratedCall:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    am_score: float | None
    evidence_gene: str          # gene-specific evidence label (or "Indeterminate")
    evidence_domain: str        # domain-aggregate (recorded only)
    is_am_plp: bool             # gene-specific evidence >= min strength


def _rank(label: str | None) -> int:
    return _RANK.get(label or "Indeterminate", 0)


def call_carrier(am_score: float | None, gene: str, gene_thresholds: dict[str, dict],
                 min_evidence: str = "PP3_Moderate") -> str:
    """Look up the gene-specific threshold and return the evidence label."""
    row = gene_thresholds.get(gene)
    if row is None or am_score is None:
        return "Indeterminate"
    # gene_thresholds[gene] holds a list of {threshold: float, evidence_label: str}
    # sorted by rank descending; pick the highest whose threshold is met.
    best = "Indeterminate"
    for cutoff in row.get("cutoffs", []):
        if am_score >= cutoff["threshold"]:
            if _rank(cutoff["evidence_label"]) > _rank(best):
                best = cutoff["evidence_label"]
    return best


def promote(label_gene: str, min_evidence: str = "PP3_Moderate") -> bool:
    return _rank(label_gene) >= _rank(min_evidence)


# TODO(minerva): load calibration table -> gene_thresholds dict; iterate over
# per-chunk annotated VCFs (VEP CSQ carries am_pathogenicity); emit
# results/variants/am_plp.tsv with columns chr,pos,ref,alt,gene,am_score,
# evidence_gene,evidence_domain,is_AM_PLP.
