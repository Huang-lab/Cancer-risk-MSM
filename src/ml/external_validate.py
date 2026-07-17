"""BioMe cross-cohort external validation hook.

Score an MSM-trained model on BioMe held-out data (and vice versa) using the
same temporal-leakage window and covariate schema. Requires the BioMe carrier
matrix + phenotype tables to be present on Minerva.
"""
from __future__ import annotations

# TODO(minerva): external_validate(model, external_cfg) -> metrics dict
