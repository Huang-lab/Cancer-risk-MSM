"""Covariate extraction for ML.

Sources: Demographics, Social_History (smoking/alcohol), Vitals (BMI),
OB_HISTORY (parity/age at first birth), HEALTH_MAINTENANCE_* (screening).

Longitudinal fields are collapsed to a single per-person value using the
temporal-leakage window from config.ml.temporal_leakage.feature_window_days
(pre-index [-730, -182] by default), joined against the case index date.

Scaffolded; runs on Minerva.
"""
from __future__ import annotations

# TODO(minerva): implement build_covariates(cfg, cases_df) -> DataFrame keyed by person_id.
