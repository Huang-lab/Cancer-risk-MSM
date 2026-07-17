"""Family history of cancer from Family_History.txt.

Filter to cancer conditions; keep first-degree relatives by default (config).
Expect messy relationship coding and missing onset age. Emit:
  person_id, cancer, has_famhx_1deg, n_affected_1deg, has_famhx_any

Scaffolded; runs on Minerva.
"""
from __future__ import annotations

# TODO(minerva): implement build_famhx_table(cfg) -> DataFrame with the columns above.
