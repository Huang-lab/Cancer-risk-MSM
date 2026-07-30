"""Define non-cancer controls.

Controls = (cohort roster) − (people in the cancer case file) − (people with a
cancer ICD code anywhere in the EHR).

That last term is deliberate. The case file is authoritative for *cases*, but it
was built for a specific purpose and may not enumerate every malignancy in the
EHR. Anyone carrying a cancer ICD code that the case file omits would otherwise
be silently admitted as a "non-cancer" control and bias every comparison, so
they are excluded and counted separately as `excluded_icd_cancer` rather than
being quietly folded into either group.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from src.phenotype.icd_mapping import IcdToPhecodeX, infer_icd_version

# phecodeX / phecode prefixes that denote a malignancy. Kept broad on purpose:
# for control screening, a false positive costs one control while a false
# negative contaminates the control group.
_CANCER_PHECODE_PREFIXES = ("CA_", "153", "174", "185", "162", "157", "155", "183")


@dataclass
class RosterReport:
    n_roster: int = 0
    n_cases: int = 0
    n_controls: int = 0
    n_excluded_icd_cancer: int = 0
    n_cases_not_in_roster: int = 0
    n_too_few_encounters: int = 0
    warnings: list[str] = field(default_factory=list)


def read_person_ids(path: str | Path, person_id_col: str,
                    sep: str = "|") -> set[str]:
    """Collect distinct person IDs from any EHR table."""
    p = Path(path)
    if not p.exists():
        return set()
    out: set[str] = set()
    with p.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=sep)
        cols = rdr.fieldnames or []
        col = person_id_col if person_id_col in cols else (cols[0] if cols else None)
        if col is None:
            return set()
        for row in rdr:
            pid = (row.get(col) or "").strip()
            if pid:
                out.add(pid)
    return out


def encounter_counts(path: str | Path, person_id_col: str,
                     sep: str = "|") -> dict[str, int]:
    p = Path(path)
    counts: dict[str, int] = {}
    if not p.exists():
        return counts
    with p.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=sep)
        cols = rdr.fieldnames or []
        col = person_id_col if person_id_col in cols else (cols[0] if cols else None)
        if col is None:
            return counts
        for row in rdr:
            pid = (row.get(col) or "").strip()
            if pid:
                counts[pid] = counts.get(pid, 0) + 1
    return counts


def is_cancer_phecode(code: str) -> bool:
    c = (code or "").strip().upper()
    return bool(c) and c.startswith(_CANCER_PHECODE_PREFIXES)


def people_with_cancer_codes(paths_and_schemas: list[tuple[Path, dict]],
                             code_prefixes: set[str],
                             prefix_len: int = 3) -> set[str]:
    """Flag anyone whose ICD code matches a cancer code prefix.

    Preferred over `people_with_cancer_icd`, which is limited to whatever the
    built-in ICD->phecodeX map covers. `code_prefixes` comes from the case file
    (see cancer_cases.load_cancer_code_prefixes) and therefore spans every
    cancer category present in this cohort.
    """
    flagged: set[str] = set()
    if not code_prefixes:
        return flagged
    for path, schema in paths_and_schemas:
        if not Path(path).exists():
            continue
        sep = schema.get("sep", "|")
        pid_col = schema["person_id_col"]
        icd_col = schema["icd_col"]
        with Path(path).open(newline="") as fh:
            rdr = csv.DictReader(fh, delimiter=sep)
            for row in rdr:
                pid = (row.get(pid_col) or "").strip()
                if not pid or pid in flagged:
                    continue
                code = (row.get(icd_col) or "").strip().upper().replace(" ", "")
                if len(code) >= prefix_len and code[:prefix_len] in code_prefixes:
                    flagged.add(pid)
    return flagged


def people_with_cancer_icd(paths_and_schemas: list[tuple[Path, dict]],
                           mapper: IcdToPhecodeX) -> set[str]:
    """Scan ICD-coded EHR tables for anyone with a malignancy code.

    `paths_and_schemas` is [(path, ehr_schema_dict), ...] for
    Encounter_Diagnosis and Problem_List.
    """
    flagged: set[str] = set()
    for path, schema in paths_and_schemas:
        if not Path(path).exists():
            continue
        sep = schema.get("sep", "|")
        pid_col = schema["person_id_col"]
        icd_col = schema["icd_col"]
        ver_col = schema.get("icd_version_col")
        with Path(path).open(newline="") as fh:
            rdr = csv.DictReader(fh, delimiter=sep)
            for row in rdr:
                pid = (row.get(pid_col) or "").strip()
                code = (row.get(icd_col) or "").strip()
                if not pid or not code or pid in flagged:
                    continue
                ver = None
                if ver_col and row.get(ver_col):
                    v = str(row[ver_col]).strip().upper().replace("ICD-", "").replace("ICD", "")
                    ver = "9" if v.startswith("9") else "10" if v.startswith("10") else None
                phe = mapper.map_code(code, ver or infer_icd_version(code))
                if phe and is_cancer_phecode(phe):
                    flagged.add(pid)
    return flagged


def assign_groups(roster: set[str], case_ids: set[str],
                  icd_cancer_ids: set[str] | None = None,
                  enc_counts: dict[str, int] | None = None,
                  min_encounters: int = 0) -> tuple[dict[str, str], RosterReport]:
    """Return ({person_id -> group}, report).

    Groups: `case`, `control`, `excluded_icd_cancer`, `excluded_few_encounters`.
    Cases are kept even if absent from the roster (counted in the report), since
    dropping a known case would understate the numerator.
    """
    icd_cancer_ids = icd_cancer_ids or set()
    rep = RosterReport()
    rep.n_roster = len(roster)

    groups: dict[str, str] = {}

    for pid in case_ids:
        groups[pid] = "case"
        if pid not in roster:
            rep.n_cases_not_in_roster += 1
    rep.n_cases = len(case_ids)

    for pid in roster:
        if pid in groups:            # already a case
            continue
        if pid in icd_cancer_ids:
            groups[pid] = "excluded_icd_cancer"
            rep.n_excluded_icd_cancer += 1
            continue
        if min_encounters and (enc_counts or {}).get(pid, 0) < min_encounters:
            groups[pid] = "excluded_few_encounters"
            rep.n_too_few_encounters += 1
            continue
        groups[pid] = "control"
        rep.n_controls += 1

    if rep.n_cases_not_in_roster:
        rep.warnings.append(
            f"{rep.n_cases_not_in_roster} case(s) absent from the roster source; "
            "kept as cases. Check the roster covers the same batch as the case file.")
    if rep.n_controls == 0 and rep.n_roster:
        rep.warnings.append("zero controls — check the roster source and ID forms match")
    return groups, rep
