"""Propensity-score k:1 caliper matching for risk-set sampling.

Controls are drawn from the clean roster (no cancer family history unless
overridden per cancer YAML). Mirrors the approach in
Huang-lab/BioMe-cancer-risk-prediction.
"""
from __future__ import annotations

# TODO(minerva): match_cases_to_controls(df, ratio, caliper, covariates) -> matched_df
