"""Load the pre-built cancer case list (EHR + genetics combined).

This file is the AUTHORITATIVE case definition — we do not re-derive cases from
ICD codes. Controls are a set-difference against it (see `controls.py`).

Only the ID column is known up front (`sample_name`); everything else is
auto-detected, and `describe_schema()` reports what was found so the real
layout can be confirmed from a live run without anyone pasting PHI.

A person may legitimately appear on multiple rows (several cancers, or one row
per cancer x data-source), so cases are aggregated to one record per person
with the set of cancer types attached.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

# Candidate columns naming the cancer type, best first.
_CANCER_TYPE_CANDIDATES = [
    "cancer_type", "cancer", "cancer_site", "primary_site", "site",
    "phenotype", "diagnosis", "dx_name", "tumor_type", "phecode_description",
]

# Candidate columns flagging genetic-data availability.
_GENETICS_CANDIDATES = [
    "has_genetics", "has_wes", "genotyped", "in_wes", "has_genotype",
]

# Cancer-type normalization: free text -> canonical key. Order matters;
# longer/more specific patterns are checked before generic ones.
_TYPE_PATTERNS: list[tuple[str, list[str]]] = [
    ("colorectal",  ["colorect", "colon", "rectal", "rectum"]),
    ("breast",      ["breast"]),
    ("prostate",    ["prostate"]),
    ("lung",        ["lung", "bronchus"]),
    ("ovarian",     ["ovary", "ovarian"]),
    ("pancreatic",  ["pancrea"]),
    ("endometrial", ["endometri", "uterine", "uterus"]),
    ("gastric",     ["gastric", "stomach"]),
    ("kidney",      ["kidney", "renal"]),
    ("bladder",     ["bladder"]),
    ("melanoma",    ["melanoma"]),
    ("thyroid",     ["thyroid"]),
    ("liver",       ["liver", "hepatocellular", "hepatic"]),
    ("esophageal",  ["esophag", "oesophag"]),
    ("lymphoma",    ["lymphoma"]),
    ("leukemia",    ["leukemia", "leukaemia"]),
    ("myeloma",     ["myeloma"]),
    ("brain",       ["glioma", "glioblastoma", "brain"]),
]

_TRUTHY = {"1", "true", "t", "yes", "y"}


@dataclass
class CancerCase:
    person_id: str
    cancer_types: set[str] = field(default_factory=set)
    raw_types: set[str] = field(default_factory=set)
    n_rows: int = 0
    has_genetics: bool | None = None

    @property
    def n_cancers(self) -> int:
        return len(self.cancer_types) or (1 if self.n_rows else 0)


def normalize_cancer_type(raw: str) -> str | None:
    """'Colon adenocarcinoma' -> 'colorectal'. None when unclassifiable.

    For free-text diagnosis names (`dx_name`). When the source has a curated
    category column (`Group`), prefer `slugify_cancer_group` instead -- pattern
    matching a curated vocabulary only loses information.
    """
    if not raw:
        return None
    t = raw.strip().lower()
    if not t or t in ("na", "null", ".", "-", "unknown"):
        return None
    for key, pats in _TYPE_PATTERNS:
        if any(p in t for p in pats):
            return key
    return None


# Curated categories that carry no usable primary site. Kept as cases (the
# person does have cancer) but excluded from per-site counts.
_UNSPECIFIED_GROUPS = {
    "cancer of unknown primary (cup)",
    "other/ill-defined primary sites",
    "hematologic/lymphatic - other/unspecified b-cell",
}


def slugify_cancer_group(raw: str) -> str | None:
    """Curated `Group` value -> a stable key. 'Head & Neck' -> 'head_neck'.

    Passthrough by design: the Group column is already a clean vocabulary
    (Breast, Colorectum, Multiple Myeloma, ...), so it is slugified rather than
    re-classified. Site-less buckets such as 'Cancer of Unknown Primary (CUP)'
    return None so they do not become a spurious site, while the person still
    counts as a case.
    """
    if not raw:
        return None
    t = raw.strip().lower()
    if not t or t in ("na", "null", ".", "-", "unknown"):
        return None
    if t in _UNSPECIFIED_GROUPS:
        return None
    t = t.replace("&", " and ")
    t = re.sub(r"[^a-z0-9]+", "_", t).strip("_")
    return t or None


def _pick(cols: list[str], candidates: list[str]) -> str | None:
    """Exact match first, then case-insensitive, then substring.

    Falsy candidates are dropped: an empty string would substring-match the
    first column, since `"" in anything` is True.
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


def describe_schema(path: str | Path, sep: str = "\t") -> dict:
    """Report the file's structure WITHOUT emitting any cell values.

    Safe to print or paste: returns column names, row count, and which columns
    were auto-detected — never participant data.
    """
    p = Path(path)
    if not p.exists():
        return {"error": f"file not found: {p}"}
    with p.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=sep)
        cols = list(rdr.fieldnames or [])
        n_rows = sum(1 for _ in rdr)
    return {
        "path": str(p),
        "n_columns": len(cols),
        "n_rows": n_rows,
        "columns": cols,
        "detected_cancer_type_col": _pick(cols, _CANCER_TYPE_CANDIDATES),
        "detected_genetics_col": _pick(cols, _GENETICS_CANDIDATES),
    }


def load_cancer_cases(path: str | Path, person_id_col: str = "sample_name",
                      cancer_type_col: str | None = None,
                      sep: str = "\t",
                      cancer_type_mode: str = "normalize",
                      ) -> tuple[dict[str, CancerCase], dict]:
    """Return ({person_id -> CancerCase}, report).

    The report records the detected schema, unclassifiable type strings (as
    counts, not values), and rows dropped for a blank ID — so nothing is
    silently discarded.
    """
    p = Path(path)
    report: dict = {
        "path": str(p), "n_rows": 0, "n_people": 0,
        "cancer_type_col": None, "genetics_col": None,
        "n_rows_blank_id": 0, "n_rows_unclassified_type": 0,
        "unclassified_type_examples": [], "types_seen": {},
        "warnings": [],
    }
    if not p.exists():
        report["warnings"].append(f"file not found: {p}")
        return {}, report

    cases: dict[str, CancerCase] = {}
    type_counts: dict[str, int] = defaultdict(int)
    unclassified: set[str] = set()

    with p.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=sep)
        cols = list(rdr.fieldnames or [])
        if person_id_col not in cols:
            alt = _pick(cols, [person_id_col, "sample_name", "person_id", "sample_id"])
            if alt is None:
                report["warnings"].append(
                    f"ID column {person_id_col!r} not found; columns are {cols[:10]}")
                return {}, report
            report["warnings"].append(
                f"ID column {person_id_col!r} absent; using {alt!r}")
            person_id_col = alt

        type_col = cancer_type_col or _pick(cols, _CANCER_TYPE_CANDIDATES)
        gen_col = _pick(cols, _GENETICS_CANDIDATES)
        report["cancer_type_col"] = type_col
        report["genetics_col"] = gen_col
        if type_col is None:
            report["warnings"].append(
                "no cancer-type column detected; all cases treated as "
                "cancer-type-agnostic (group still 'case')")

        for row in rdr:
            report["n_rows"] += 1
            pid = (row.get(person_id_col) or "").strip()
            if not pid:
                report["n_rows_blank_id"] += 1
                continue
            rec = cases.setdefault(pid, CancerCase(person_id=pid))
            rec.n_rows += 1

            if type_col:
                raw = (row.get(type_col) or "").strip()
                if raw:
                    rec.raw_types.add(raw)
                    norm = (slugify_cancer_group(raw)
                            if cancer_type_mode == "passthrough"
                            else normalize_cancer_type(raw))
                    if norm:
                        rec.cancer_types.add(norm)
                        type_counts[norm] += 1
                    else:
                        report["n_rows_unclassified_type"] += 1
                        if len(unclassified) < 20:
                            unclassified.add(raw)

            if gen_col:
                v = str(row.get(gen_col) or "").strip().lower()
                if v:
                    rec.has_genetics = v in _TRUTHY

    report["n_people"] = len(cases)
    report["types_seen"] = dict(sorted(type_counts.items(), key=lambda kv: -kv[1]))
    # Type strings are clinical vocabulary, not participant data — safe to show,
    # and necessary for extending _TYPE_PATTERNS.
    report["unclassified_type_examples"] = sorted(unclassified)
    return cases, report


def load_cancer_code_prefixes(path: str | Path, code_col: str = "dx_code",
                              sep: str = "\t", prefix_len: int = 3) -> set[str]:
    """Cancer ICD code prefixes taken from the case file itself.

    The case file is already cancer-filtered and carries the diagnosis codes, so
    its distinct `dx_code` values ARE the cancer vocabulary for this cohort --
    all 54 Group categories, not just the handful the built-in ICD->phecodeX map
    happens to cover. Screening controls with the built-in map misses ~45
    categories (liver, leukemia, lymphoma, thyroid, ...), which lets people with
    an undocumented cancer be admitted as "non-cancer" controls.

    Truncating to a 3-character category ("C18.70" -> "C18") absorbs subcode
    variation between the case file and the roster tables.
    """
    p = Path(path)
    out: set[str] = set()
    if not p.exists():
        return out
    with p.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=sep)
        cols = list(rdr.fieldnames or [])
        col = _pick(cols, [code_col, "dx_code", "Code", "icd_code", "code"])
        if col is None:
            return out
        for row in rdr:
            code = (row.get(col) or "").strip().upper().replace(" ", "")
            if len(code) >= prefix_len:
                out.add(code[:prefix_len])
    return out


def cancer_types_present(cases: dict[str, CancerCase]) -> list[str]:
    seen: set[str] = set()
    for c in cases.values():
        seen |= c.cancer_types
    return sorted(seen)
