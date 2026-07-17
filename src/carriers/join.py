"""Join the long carrier table with harmonized IDs, ancestry PCs, and phenotypes.

Emits an analysis-ready table under results/analysis/ for the ML pipeline.

Scaffolded; runs on Minerva.
"""
from __future__ import annotations

# TODO(minerva): implement build_analysis_table(cfg) that:
#   1. Loads carriers.long.tsv.
#   2. Maps sample IDs to canonical person_id via src.phenotype.id_harmonize.
#   3. Left-joins PCs (from common_variants batch) and phenotype rows.
#   4. Writes results/analysis/analysis_ready.tsv (long) and (optional) wide.
