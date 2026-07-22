"""Cancer case/control identification from ICD-coded EHR.

Scans Encounter_Diagnosis + Problem_List, maps each ICD code to a phecodeX,
and for each cancer defined in `config.phenotype.cancer_phecodes` emits:

    person_id, cancer, is_case, first_dx_date, incident_flag,
    n_dx_encounters, source

Cross-checks with the pre-computed `PheWas_MSM_phecodeX.tsv` roster when the
file is present.

Sources are read directly from files pointed to by
`config.inputs.phenotypes.files.*`; column names come from
`config.phenotype.ehr_schema.*`.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from src.phenotype.icd_mapping import IcdToPhecodeX, infer_icd_version


@dataclass
class DxRow:
    person_id: str
    phecodex: str
    date: str | None    # ISO yyyy-mm-dd if parseable, else None
    source: str          # "encounter_diagnosis" | "problem_list" | "phecodex_roster"


@dataclass
class CaseRow:
    person_id: str
    cancer: str
    is_case: bool
    first_dx_date: str | None
    incident_flag: bool | None      # None when enrollment date is unknown
    n_dx_encounters: int
    sources: str                    # comma-separated set of sources


def _parse_date(raw: str, fmt: str) -> str | None:
    from datetime import datetime
    if not raw or raw in ("NA", "NULL", "."):
        return None
    try:
        return datetime.strptime(raw.strip(), fmt).date().isoformat()
    except ValueError:
        return None


def read_dx_ehr(path: str | Path, schema: dict, mapper: IcdToPhecodeX,
                source: str) -> list[DxRow]:
    """Read Encounter_Diagnosis or Problem_List into DxRow list (ICD -> phecodeX)."""
    out: list[DxRow] = []
    if not Path(path).exists():
        return out
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=schema.get("sep", "\t"))
        pid_col = schema["person_id_col"]
        icd_col = schema["icd_col"]
        ver_col = schema.get("icd_version_col")
        date_col = schema.get("date_col") or schema.get("noted_date_col")
        date_fmt = schema.get("date_format", "%Y-%m-%d")
        for row in rdr:
            pid = (row.get(pid_col) or "").strip()
            code = (row.get(icd_col) or "").strip()
            if not pid or not code:
                continue
            ver = None
            if ver_col and row.get(ver_col):
                ver = str(row[ver_col]).strip().replace("ICD-", "").replace("ICD", "").strip()
                ver = "9" if ver.startswith("9") else "10" if ver.startswith("10") else None
            phe = mapper.map_code(code, ver or infer_icd_version(code))
            if phe is None:
                continue
            dt = _parse_date(row.get(date_col, ""), date_fmt) if date_col else None
            out.append(DxRow(pid, phe, dt, source))
    return out


def read_phecodex_roster(path: str | Path, schema: dict,
                         cancer_phecodes: dict[str, list[str]]) -> list[DxRow]:
    """Read the pre-computed phecodeX roster. Supports wide OR long formats."""
    out: list[DxRow] = []
    if not Path(path).exists():
        return out
    interesting = {p for pcs in cancer_phecodes.values() for p in pcs}
    sep = schema.get("sep", "\t")
    fmt = schema.get("format", "auto")
    pid_col = schema.get("person_id_col", "MRN")

    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=sep)
        cols = rdr.fieldnames or []
        # Auto-detect: wide if any interesting phecodeX appears as a column name.
        is_wide = fmt == "wide" or (fmt == "auto" and any(c in interesting for c in cols))
        if is_wide:
            for row in rdr:
                pid = (row.get(pid_col) or "").strip()
                if not pid:
                    continue
                for phe in interesting:
                    if phe in row and str(row[phe]).strip() in ("1", "True", "true", "case"):
                        out.append(DxRow(pid, phe, None, "phecodex_roster"))
        else:
            phe_col = schema.get("phecodex_col", "phecodeX")
            val_col = schema.get("value_col", "case_control")
            for row in rdr:
                pid = (row.get(pid_col) or "").strip()
                phe = (row.get(phe_col) or "").strip()
                if not pid or phe not in interesting:
                    continue
                if str(row.get(val_col, "")).strip() in ("1", "True", "true", "case"):
                    out.append(DxRow(pid, phe, None, "phecodex_roster"))
    return out


def build_cases(dx_rows: list[DxRow], cancer_phecodes: dict[str, list[str]],
                enrollment_by_pid: dict[str, str] | None = None) -> list[CaseRow]:
    """Aggregate per (person, cancer) -> CaseRow. is_case iff >=1 matching dx."""
    # Index rows: (pid, cancer) -> [DxRow, ...]
    idx: dict[tuple[str, str], list[DxRow]] = {}
    for r in dx_rows:
        for cancer, codes in cancer_phecodes.items():
            if r.phecodex in codes:
                idx.setdefault((r.person_id, cancer), []).append(r)

    out: list[CaseRow] = []
    for (pid, cancer), rows in idx.items():
        dates = sorted([r.date for r in rows if r.date])
        first_dx = dates[0] if dates else None
        incident = None
        if enrollment_by_pid and first_dx:
            enr = enrollment_by_pid.get(pid)
            incident = (first_dx > enr) if enr else None
        out.append(CaseRow(
            person_id=pid,
            cancer=cancer,
            is_case=True,
            first_dx_date=first_dx,
            incident_flag=incident,
            n_dx_encounters=len(rows),
            sources=",".join(sorted({r.source for r in rows})),
        ))
    return out
