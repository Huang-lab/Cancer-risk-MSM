"""Cancer case/control identification.

Sources (per user):
  - `PheWas_MSM_phecodeX.tsv`   — standardized phecodeX roster (case/control).
  - `Encounter_Diagnosis.txt`   — ICD-coded; provides first-diagnosis dates.
  - `Problem_List.txt`          — ICD-coded; cross-check for prevalent status.
  - `Medical_History.txt`, `Surgical_History.txt` — cross-check.

Strategy
--------
1. Load the phecodeX roster (already standardized) -> primary case/control label.
2. Map ICD codes in Encounter_Diagnosis + Problem_List to phecodeX via
   `phenotype.icd_to_phecodex_map`.
3. For each cancer defined in `phenotype.cancer_phecodes`:
     - Case set  = roster cases  OR (>=1 ICD->phecodeX match for that cancer).
     - Case date = min(encounter_date) with a matching ICD -> phecodeX code.
     - Cross-check Problem_List / Medical_History / Surgical_History to
       distinguish incident vs prevalent cases (e.g. prior colectomy/mastectomy).
4. Emit an ICD-mapping audit table so downstream ML can debug leakage.

Scaffolded; runs on Minerva.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class CaseRow:
    person_id: str
    cancer: str
    is_case: bool
    is_incident: bool
    first_dx_date: str | None   # yyyy-mm-dd
    source: str                  # "phecodeX_roster" | "encounter_dx" | "problem_list"


# TODO(minerva): implement:
#   load_phecodex_roster(cfg) -> DataFrame
#   map_icd_to_phecodex(encounter_df, problem_df, cfg) -> DataFrame  # incl. first-dx date
#   build_case_table(cfg) -> list[CaseRow]  # merges roster + ICD-derived
