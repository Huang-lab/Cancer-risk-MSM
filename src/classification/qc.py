"""Per-gene QC for P/LP call sets.

Checks:
  1. Paralog/pseudogene artifact scrutiny (PMS2/PMS2CL, CHEK2, NBN, SDHA, ...):
     - Flag variants in high-identity paralog regions unless supported by
       paralog-aware read alignment (marker in header).
  2. Cohort vs gnomAD carrier-frequency comparison:
     - For each gene, compare cohort P/LP carrier frequency vs gnomAD popmax P/LP
       carrier frequency. If ratio > cfg.classification.qc.gnomad_inflation_ratio_max,
       raise a warning row in the QC report.

Scaffolded today; runs on Minerva after classification.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class GeneQcRow:
    gene: str
    n_variants: int
    n_carriers: int
    cohort_freq: float
    gnomad_freq: float
    ratio: float
    flag: str          # "OK" / "INFLATED" / "PARALOG_SUSPECT"


# TODO(minerva): implement compare_freqs(cohort_df, gnomad_df, cfg) -> [GeneQcRow];
# write results/variants/qc_per_gene.tsv and a short SUMMARY.md.
