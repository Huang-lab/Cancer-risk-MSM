"""Harmonize MSM participant IDs across pre-2022 and post-2022 conventions.

Pre-2022 form:  SINAI_###_AB########
Post-2022 form: SINAI-Million_#####_##########

Downstream joins (WES <-> PCA <-> phenotypes) use a canonical `person_id`
column. This module owns the mapping.

No participant data lives here in the repo — only regex + join utilities that
run on Minerva.
"""
from __future__ import annotations

import re

_SINAI_OLD = re.compile(r"^SINAI_\d+_[A-Z]{2}\d+$")
_SINAI_NEW = re.compile(r"^SINAI-Million_\d+_\d+$")


def id_form(raw: str) -> str:
    if _SINAI_OLD.match(raw):
        return "sinai_old"
    if _SINAI_NEW.match(raw):
        return "sinai_new"
    return "unknown"


# TODO(minerva): implement build_id_map(cfg) -> DataFrame[canonical_id, sinai_old, sinai_new]
# using the crosswalk file the biobank provides (path in config). All joins
# downstream import canonical_id.
