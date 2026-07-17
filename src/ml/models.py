"""Model zoo. Interpretable baselines first; hooks for richer models.

- Penalized logistic regression (sklearn.linear_model.LogisticRegression, L1/L2)
- Cox proportional hazards (lifelines.CoxPHFitter) with time-to-diagnosis
- Random Forest (sklearn.ensemble.RandomForestClassifier)
- XGBoost (xgboost.XGBClassifier)

All models exposed through a single `fit_and_score(model_name, X, y, groups, cfg)` API.
"""
from __future__ import annotations

# TODO(minerva): fit_and_score(...)
