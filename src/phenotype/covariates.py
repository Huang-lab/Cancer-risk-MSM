"""Lifestyle + demographic covariate extraction.

Sources (column names all config-driven under phenotype.ehr_schema.*):
  Demographics    -> age / birth_year, sex, race-ethnicity
  Social_History  -> smoking status, alcohol use
  Vitals          -> BMI (most recent, or closest to an index date)
  OB_HISTORY      -> parity, age at first birth

Longitudinal sources carry many rows per person. Two selection rules:
  * no index date  -> most recent non-null measurement
  * index date     -> closest measurement STRICTLY BEFORE the index date, so a
                      value recorded at or after diagnosis can never leak in

Free-text lifestyle fields are normalized through the config maps
(`cohort.smoking_map` / `cohort.alcohol_map`) to never/former/current; anything
unmatched becomes `unknown` and is counted, never silently coerced.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from datetime import date, datetime
from pathlib import Path

_NULLISH = {"", "na", "n/a", "null", ".", "-", "unknown", "not recorded", "none recorded"}

_SMOKING_CANDIDATES = ["smoking_status", "smoking", "tobacco_use", "smoking_use", "tobacco"]
_ALCOHOL_CANDIDATES = ["alcohol_use", "alcohol", "alcohol_status", "etoh"]
_BMI_CANDIDATES = ["bmi", "body_mass_index"]
_HEIGHT_CANDIDATES = ["height", "height_cm", "height_in"]
_WEIGHT_CANDIDATES = ["weight", "weight_kg", "weight_lb"]
_SEX_CANDIDATES = ["sex", "gender", "sex_at_birth", "legal_sex"]
_RACE_CANDIDATES = ["race_ethnicity", "race", "ethnicity", "self_reported_race"]
_BIRTH_YEAR_CANDIDATES = ["birth_year", "year_of_birth", "yob", "dob", "date_of_birth", "birth_date"]
_PARITY_CANDIDATES = ["parity", "num_births", "gravida_para", "para"]
_AFB_CANDIDATES = ["age_at_first_birth", "age_first_birth", "afb"]


@dataclass
class Covariates:
    person_id: str
    birth_year: int | None = None
    age: float | None = None
    sex: str | None = None
    race_ethnicity: str | None = None
    smoking_status: str = "unknown"
    alcohol_use: str = "unknown"
    bmi: float | None = None
    bmi_date: str | None = None
    parity: int | None = None
    age_at_first_birth: float | None = None
    notes: list[str] = field(default_factory=list)


def _nullish(v: str | None) -> bool:
    return (v or "").strip().lower() in _NULLISH


def _pick(cols: list[str], candidates: list[str]) -> str | None:
    """First matching column: exact, then case-insensitive, then substring.

    Falsy candidates are dropped first. Callers routinely pass
    `schema.get("x_col", "")`, and an empty string would substring-match the
    very first column (`"" in anything` is True), silently binding the wrong
    field and, worse, making an absent config key look present.
    """
    cands = [c for c in candidates if c]
    lower = {c.lower(): c for c in cols}
    for cand in cands:
        if cand in cols:
            return cand
        if cand.lower() in lower:
            return lower[cand.lower()]
    for cand in cands:
        for c in cols:
            if cand.lower() in c.lower():
                return c
    return None


def parse_date(raw: str | None, fmt: str) -> date | None:
    if _nullish(raw):
        return None
    raw = raw.strip()
    for f in (fmt, "%m/%d/%Y", "%Y-%m-%d", "%m/%d/%y", "%d/%m/%Y"):
        try:
            return datetime.strptime(raw, f).date()
        except ValueError:
            continue
    return None


def _float(raw: str | None) -> float | None:
    if _nullish(raw):
        return None
    try:
        return float(str(raw).strip())
    except ValueError:
        return None


def normalize_category(raw: str | None, mapping: dict[str, list[str]]) -> str:
    """Map free text through {canonical: [patterns]}. Unmatched -> 'unknown'.

    Exact matches win over substring matches so 'never smoker' can't be
    captured by a 'smoker' -> current pattern.
    """
    if _nullish(raw):
        return "unknown"
    t = raw.strip().lower()
    for canon, pats in mapping.items():
        if any(t == p.lower() for p in pats):
            return canon
    for canon, pats in mapping.items():
        if any(p.lower() in t for p in pats):
            return canon
    return "unknown"


def _birth_year(raw: str | None) -> int | None:
    if _nullish(raw):
        return None
    s = str(raw).strip()
    if len(s) == 4 and s.isdigit():
        return int(s)
    d = parse_date(s, "%Y-%m-%d")
    return d.year if d else None


# --- Demographics -----------------------------------------------------------

def read_demographics(path: str | Path, schema: dict,
                      as_of_year: int | None = None) -> dict[str, Covariates]:
    out: dict[str, Covariates] = {}
    p = Path(path)
    if not p.exists():
        return out
    with p.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=schema.get("sep", "|"))
        cols = list(rdr.fieldnames or [])
        pid_col = schema.get("person_id_col", "sample_name")
        if pid_col not in cols:
            pid_col = _pick(cols, ["sample_name", "person_id"]) or (cols[0] if cols else "")
        by_col = _pick(cols, [schema.get("birth_year_col", "")] + _BIRTH_YEAR_CANDIDATES)
        sex_col = _pick(cols, [schema.get("sex_col", "")] + _SEX_CANDIDATES)
        race_col = _pick(cols, _RACE_CANDIDATES)
        pc_cols = [c for c in cols if c.upper().startswith("PC") and c[2:].isdigit()]

        for row in rdr:
            pid = (row.get(pid_col) or "").strip()
            if not pid:
                continue
            rec = out.setdefault(pid, Covariates(person_id=pid))
            if by_col:
                rec.birth_year = _birth_year(row.get(by_col))
            if sex_col and not _nullish(row.get(sex_col)):
                rec.sex = row[sex_col].strip().upper()[:1]
            if race_col and not _nullish(row.get(race_col)):
                rec.race_ethnicity = row[race_col].strip()
            if pc_cols:
                rec.notes.append(f"ancestry PCs present ({len(pc_cols)})")
            if rec.birth_year and as_of_year:
                rec.age = float(as_of_year - rec.birth_year)
    return out


# --- Longitudinal selection -------------------------------------------------

def pick_measurement(rows: list[tuple[date | None, float]],
                     index_date: date | None) -> tuple[float | None, date | None]:
    """Closest value strictly BEFORE index_date; else most recent overall.

    Rows with no date are usable only when there is no index date — otherwise
    they can't be proven pre-index and are skipped.
    """
    if not rows:
        return None, None
    if index_date is not None:
        eligible = [(d, v) for d, v in rows if d is not None and d < index_date]
        if not eligible:
            return None, None
        d, v = max(eligible, key=lambda t: t[0])
        return v, d
    dated = [(d, v) for d, v in rows if d is not None]
    if dated:
        d, v = max(dated, key=lambda t: t[0])
        return v, d
    return rows[-1][1], None


# --- Social history (smoking / alcohol) -------------------------------------

def read_social_history(path: str | Path, schema: dict,
                        smoking_map: dict, alcohol_map: dict,
                        index_dates: dict[str, date] | None = None,
                        ) -> dict[str, tuple[str, str]]:
    """{person_id -> (smoking_status, alcohol_use)}, latest pre-index value."""
    out: dict[str, tuple[str, str]] = {}
    p = Path(path)
    if not p.exists():
        return out
    index_dates = index_dates or {}
    smoke_rows: dict[str, list[tuple[date | None, str]]] = {}
    alc_rows: dict[str, list[tuple[date | None, str]]] = {}

    with p.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=schema.get("sep", "|"))
        cols = list(rdr.fieldnames or [])
        pid_col = schema.get("person_id_col", "sample_name")
        if pid_col not in cols:
            pid_col = _pick(cols, ["sample_name", "person_id"]) or (cols[0] if cols else "")
        smoke_col = _pick(cols, [schema.get("smoking_status_col", "")] + _SMOKING_CANDIDATES)
        alc_col = _pick(cols, [schema.get("alcohol_col", "")] + _ALCOHOL_CANDIDATES)
        date_col = schema.get("date_col")
        fmt = schema.get("date_format", "%m/%d/%Y")

        for row in rdr:
            pid = (row.get(pid_col) or "").strip()
            if not pid:
                continue
            d = parse_date(row.get(date_col), fmt) if date_col else None
            if smoke_col and not _nullish(row.get(smoke_col)):
                smoke_rows.setdefault(pid, []).append((d, row[smoke_col]))
            if alc_col and not _nullish(row.get(alc_col)):
                alc_rows.setdefault(pid, []).append((d, row[alc_col]))

    for pid in set(smoke_rows) | set(alc_rows):
        idx = index_dates.get(pid)
        smoke = _latest_category(smoke_rows.get(pid, []), idx, smoking_map)
        alc = _latest_category(alc_rows.get(pid, []), idx, alcohol_map)
        out[pid] = (smoke, alc)
    return out


def _latest_category(rows: list[tuple[date | None, str]], idx: date | None,
                     mapping: dict) -> str:
    if not rows:
        return "unknown"
    if idx is not None:
        eligible = [(d, v) for d, v in rows if d is not None and d < idx]
        if not eligible:
            return "unknown"
        raw = max(eligible, key=lambda t: t[0])[1]
    else:
        dated = [(d, v) for d, v in rows if d is not None]
        raw = max(dated, key=lambda t: t[0])[1] if dated else rows[-1][1]
    return normalize_category(raw, mapping)


# --- Vitals (BMI) -----------------------------------------------------------

def _matches(value: str, labels: list[str]) -> bool:
    v = (value or "").strip().lower()
    return bool(v) and any(lab.strip().lower() in v for lab in labels if lab)


def open_rows(path: Path, schema: dict):
    """Yield dict rows, supporting headerless files via positional `columns`.

    Some MSM extracts (Vitals) ship with no header. Setting
    `has_header: false` plus a `columns:` list addresses fields by position
    instead of letting csv.DictReader consume row 1 as the header.
    """
    sep = schema.get("sep", "|")
    if schema.get("has_header", True):
        with path.open(newline="") as fh:
            yield from csv.DictReader(fh, delimiter=sep)
    else:
        names = schema.get("columns") or []
        with path.open(newline="") as fh:
            for rec in csv.reader(fh, delimiter=sep):
                yield {names[i]: rec[i] for i in range(min(len(names), len(rec)))}


def read_bmi(path: str | Path, schema: dict,
             index_dates: dict[str, date] | None = None,
             ) -> dict[str, tuple[float | None, date | None]]:
    """{person_id -> (bmi, measurement_date)} from a LONG vitals table.

    The real Vitals extract is ~84M rows in long form (one row per
    measurement, with a `measure_type` column) and has no header. Rather than
    accumulate every row, this keeps only the current best candidate per
    person while streaming, so memory scales with the number of people rather
    than the number of rows.

    Height/weight rows are retained only when no direct BMI row exists for
    that person, since a recorded BMI beats a derived one.
    """
    out: dict[str, tuple[float | None, date | None]] = {}
    p = Path(path)
    if not p.exists():
        return out
    index_dates = index_dates or {}

    type_col = schema.get("measure_type_col")
    if not type_col:                      # wide layout — fall back to columns
        return _read_bmi_wide(p, schema, index_dates)

    pid_col = schema.get("person_id_col", "sample_name")
    val_col = schema.get("value_col", "value")
    date_col = schema.get("date_col")
    fmt = schema.get("date_format", "%m/%d/%Y")
    bmi_labels = schema.get("bmi_labels", ["BMI"])
    h_labels = schema.get("height_labels", ["Height"])
    w_labels = schema.get("weight_labels", ["Weight"])

    # best[pid] = (date, bmi); hw[pid] = {"h": (date, val), "w": (date, val)}
    best: dict[str, tuple[date | None, float]] = {}
    hw: dict[str, dict[str, tuple[date | None, float]]] = {}

    def _better(pid: str, d: date | None, v: float) -> bool:
        """Keep the latest admissible measurement (strictly pre-index if known)."""
        idx = index_dates.get(pid)
        if idx is not None and (d is None or d >= idx):
            return False
        cur = best.get(pid)
        if cur is None:
            return True
        cur_d = cur[0]
        if d is None:
            return False
        return cur_d is None or d > cur_d

    for row in open_rows(p, schema):
        pid = (row.get(pid_col) or "").strip()
        if not pid:
            continue
        mtype = row.get(type_col) or ""
        is_bmi = _matches(mtype, bmi_labels)
        is_h = (not is_bmi) and _matches(mtype, h_labels)
        is_w = (not is_bmi) and (not is_h) and _matches(mtype, w_labels)
        if not (is_bmi or is_h or is_w):
            continue
        val = _float(row.get(val_col))
        if val is None:
            continue
        d = parse_date(row.get(date_col), fmt) if date_col else None

        if is_bmi:
            if 10.0 <= val <= 80.0 and _better(pid, d, val):
                best[pid] = (d, val)
        else:
            slot = hw.setdefault(pid, {})
            key = "h" if is_h else "w"
            prev = slot.get(key)
            if prev is None or (d is not None and (prev[0] is None or d > prev[0])):
                slot[key] = (d, val)

    # Derive BMI only for people with no directly recorded value.
    for pid, slot in hw.items():
        if pid in best or "h" not in slot or "w" not in slot:
            continue
        bmi = derive_bmi(slot["h"][1], slot["w"][1])
        if bmi is None:
            continue
        d = slot["w"][0] or slot["h"][0]
        if _better(pid, d, bmi):
            best[pid] = (d, bmi)

    for pid, (d, v) in best.items():
        out[pid] = (v, d)
    return out


def _read_bmi_wide(p: Path, schema: dict,
                   index_dates: dict[str, date]
                   ) -> dict[str, tuple[float | None, date | None]]:
    """Wide-layout fallback: bmi (or height+weight) as their own columns."""
    out: dict[str, tuple[float | None, date | None]] = {}
    rows_by_pid: dict[str, list[tuple[date | None, float]]] = {}
    with p.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=schema.get("sep", "|"))
        cols = list(rdr.fieldnames or [])
        pid_col = schema.get("person_id_col", "sample_name")
        if pid_col not in cols:
            pid_col = _pick(cols, ["sample_name", "person_id"]) or (cols[0] if cols else "")
        bmi_col = _pick(cols, [schema.get("bmi_col", "")] + _BMI_CANDIDATES)
        h_col = _pick(cols, [schema.get("height_col", "")] + _HEIGHT_CANDIDATES)
        w_col = _pick(cols, [schema.get("weight_col", "")] + _WEIGHT_CANDIDATES)
        date_col = schema.get("date_col")
        fmt = schema.get("date_format", "%m/%d/%Y")
        for row in rdr:
            pid = (row.get(pid_col) or "").strip()
            if not pid:
                continue
            bmi = _float(row.get(bmi_col)) if bmi_col else None
            if bmi is None and h_col and w_col:
                bmi = derive_bmi(_float(row.get(h_col)), _float(row.get(w_col)))
            if bmi is None or not (10.0 <= bmi <= 80.0):
                continue
            d = parse_date(row.get(date_col), fmt) if date_col else None
            rows_by_pid.setdefault(pid, []).append((d, bmi))
    for pid, rows in rows_by_pid.items():
        out[pid] = pick_measurement(rows, index_dates.get(pid))
    return out


def derive_bmi(height: float | None, weight: float | None) -> float | None:
    """BMI from height/weight, guessing units. cm+kg or in+lb.

    Returns None unless the result is physiologically plausible, so a unit
    misread produces a missing value rather than a wrong number.
    """
    if not height or not weight or height <= 0:
        return None
    if height > 100:                      # cm + kg
        bmi = weight / ((height / 100.0) ** 2)
    else:                                 # inches + pounds
        bmi = (weight * 703.0) / (height ** 2)
    return bmi if 10.0 <= bmi <= 80.0 else None


# --- OB history -------------------------------------------------------------

def read_ob_history(path: str | Path, schema: dict,
                    birth_years: dict[str, int] | None = None,
                    ) -> dict[str, tuple[int | None, float | None]]:
    """{person_id -> (parity, age_at_first_birth)}.

    The real OB_HISTORY extract is LONG — one row per pregnancy outcome
    (`ob_hx_outcome_dt`, `ob_hx_outcome_c`) with no parity column. So:

      parity              = number of qualifying outcome rows
      age_at_first_birth  = year(earliest outcome date) - birth_year

    `live_birth_codes` restricts which outcome codes count; when empty, every
    outcome row counts, which over-counts parity if the file also records
    miscarriages or terminations. That is why the config carries an explicit
    TODO to confirm the code vocabulary — the fallback is deliberately visible
    rather than silently wrong.

    A wide layout with real parity / age columns is still honored if present.
    """
    out: dict[str, tuple[int | None, float | None]] = {}
    p = Path(path)
    if not p.exists():
        return out
    birth_years = birth_years or {}

    with p.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=schema.get("sep", "|"))
        cols = list(rdr.fieldnames or [])
        pid_col = schema.get("person_id_col", "sample_name")
        if pid_col not in cols:
            pid_col = _pick(cols, ["sample_name", "person_id"]) or (cols[0] if cols else "")

        par_col = _pick(cols, [schema.get("parity_col", "")] + _PARITY_CANDIDATES)
        afb_col = _pick(cols, [schema.get("age_at_first_birth_col", "")] + _AFB_CANDIDATES)
        out_date_col = schema.get("outcome_date_col") or _pick(cols, ["ob_hx_outcome_dt"])
        out_code_col = schema.get("outcome_code_col") or _pick(cols, ["ob_hx_outcome_c"])
        live_codes = [str(c).strip().lower() for c in (schema.get("live_birth_codes") or [])]
        fmt = schema.get("date_format", "%m/%d/%Y")

        # Wide layout: explicit parity / AFB columns.
        if par_col or afb_col:
            for row in rdr:
                pid = (row.get(pid_col) or "").strip()
                if not pid:
                    continue
                par = _float(row.get(par_col)) if par_col else None
                afb = _float(row.get(afb_col)) if afb_col else None
                prev_par, prev_afb = out.get(pid, (None, None))
                new_par = int(max(par, prev_par)) if par is not None and prev_par is not None \
                    else (int(par) if par is not None else prev_par)
                new_afb = min(afb, prev_afb) if afb is not None and prev_afb is not None \
                    else (afb if afb is not None else prev_afb)
                out[pid] = (new_par, new_afb)
            return out

        # Long layout: count outcome rows, track the earliest date.
        counts: dict[str, int] = {}
        first: dict[str, date] = {}
        for row in rdr:
            pid = (row.get(pid_col) or "").strip()
            if not pid:
                continue
            if live_codes:
                code = str(row.get(out_code_col) or "").strip().lower()
                if code not in live_codes:
                    continue
            counts[pid] = counts.get(pid, 0) + 1
            d = parse_date(row.get(out_date_col), fmt) if out_date_col else None
            if d is not None and (pid not in first or d < first[pid]):
                first[pid] = d

        for pid, n in counts.items():
            afb = None
            by = birth_years.get(pid)
            d = first.get(pid)
            if by and d:
                age = d.year - by
                if 10 <= age <= 60:          # implausible -> leave missing
                    afb = float(age)
            out[pid] = (n, afb)
    return out


def missingness(records: dict[str, Covariates]) -> dict[str, float]:
    """Fraction missing per covariate — surfaces unusable columns early."""
    n = len(records) or 1
    fields = ["birth_year", "age", "sex", "race_ethnicity", "bmi",
              "parity", "age_at_first_birth"]
    out = {f: sum(1 for r in records.values() if getattr(r, f) is None) / n
           for f in fields}
    for f, unk in (("smoking_status", "unknown"), ("alcohol_use", "unknown")):
        out[f] = sum(1 for r in records.values() if getattr(r, f) == unk) / n
    return out
