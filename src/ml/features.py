"""Feature engineering with strict temporal-leakage control.

Feature ascertainment window (per config): [index - 730d, index - 182d].
Only records in this window contribute to features; dx date must postdate the window.
"""
from __future__ import annotations

# TODO(minerva): build_features(cases_df, covariates_df, cfg) -> DataFrame
