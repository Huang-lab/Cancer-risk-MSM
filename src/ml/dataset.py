"""Assemble the modeling dataset for a given cancer.

Loads config/<cancer>.yaml on top of config/config.yaml, joins the analysis-ready
carrier + covariate table, applies temporal-leakage windows, and emits X, y, groups.

Scaffolded; runs on Minerva.
"""
from __future__ import annotations

# TODO(minerva): build_dataset(cancer: str, cfg) -> (X, y, groups, meta)
