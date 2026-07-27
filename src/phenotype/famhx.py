"""Family history of cancer from Family_History.txt.

For each person we emit:
    person_id
    has_famhx_1deg_any_cancer   bool
    has_famhx_2deg_any_cancer   bool
    has_famhx_1deg_<cancer>     bool     # per cancer in config.phenotype.cancer_phecodes
    has_famhx_2deg_<cancer>     bool
    n_affected_1deg             int
    n_affected_2deg             int

Relationship classification
---------------------------
1st degree: parent (mother, father), sibling (brother, sister, full-sib),
            child (son, daughter).
2nd degree: grandparent (maternal/paternal), aunt/uncle, niece/nephew,
            half-sibling, grandchild.

We normalize free-text relationship strings to these buckets. If the source
provides a `degree_col` (numeric 1 or 2), we use that directly.

Condition classification
------------------------
Preferred: coded `condition_code_col` (ICD or SNOMED) mapped to phecodeX.
Fallback: free-text `condition_col` matched against config
`phenotype.famhx.cancer_free_text_patterns` (case-insensitive).
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

from src.phenotype.icd_mapping import IcdToPhecodeX


# --- Relationship -> degree normalization -----------------------------------

_FIRST_DEGREE = {
    "MOTHER", "FATHER", "PARENT",
    "SISTER", "BROTHER", "SIBLING", "FULL_SIBLING", "FULL SIB", "FULLSIB",
    "SON", "DAUGHTER", "CHILD",
}
_SECOND_DEGREE = {
    "GRANDMOTHER", "GRANDFATHER", "GRANDPARENT",
    "MATERNAL_GRANDMOTHER", "MATERNAL GRANDMOTHER", "MGM",
    "MATERNAL_GRANDFATHER", "MATERNAL GRANDFATHER", "MGF",
    "PATERNAL_GRANDMOTHER", "PATERNAL GRANDMOTHER", "PGM",
    "PATERNAL_GRANDFATHER", "PATERNAL GRANDFATHER", "PGF",
    "AUNT", "MATERNAL_AUNT", "PATERNAL_AUNT",
    "UNCLE", "MATERNAL_UNCLE", "PATERNAL_UNCLE",
    "NIECE", "NEPHEW",
    "HALF_SISTER", "HALF SISTER", "HALFSISTER",
    "HALF_BROTHER", "HALF BROTHER", "HALFBROTHER",
    "HALF_SIBLING", "HALF SIBLING", "HALFSIB",
    "GRANDSON", "GRANDDAUGHTER", "GRANDCHILD",
}


def _normalize_rel(raw: str) -> str:
    return re.sub(r"[\s\-]+", "_", (raw or "").strip().upper())


def relationship_degree(raw: str) -> int | None:
    """Return 1 for 1st-degree, 2 for 2nd-degree, None for other."""
    if not raw:
        return None
    norm = _normalize_rel(raw)
    # Accept variants with spaces too
    canon = norm.replace("_", " ")
    if norm in _FIRST_DEGREE or canon in _FIRST_DEGREE:
        return 1
    if norm in _SECOND_DEGREE or canon in _SECOND_DEGREE:
        return 2
    return None


# --- Cancer condition detection --------------------------------------------

def is_cancer_free_text(text: str, patterns: list[str]) -> bool:
    t = (text or "").lower()
    return any(p.lower() in t for p in patterns)


def is_cancer_coded(code: str, mapper: IcdToPhecodeX,
                    cancer_phecodes_all: set[str]) -> str | None:
    """Return the phecodeX if `code` (ICD) maps to any cancer phecodeX; else None."""
    phe = mapper.map_code(code)
    return phe if phe in cancer_phecodes_all else None


def classify_cancer_type(text: str, type_patterns: dict[str, list[str]]) -> str | None:
    """Map free-text condition ('Prostate Cancer') to a cancer-type key ('prostate').

    Case-insensitive substring match; first cancer whose patterns hit wins.
    Returns None if the text is not classifiable (still counts as any-cancer
    via `is_cancer_free_text` upstream).
    """
    if not text:
        return None
    t = text.lower()
    for cancer, patterns in type_patterns.items():
        if any(p.lower() in t for p in patterns):
            return cancer
    return None


# --- Aggregation ------------------------------------------------------------

@dataclass
class PersonFamHx:
    """`cancer_by_degree[deg]` holds cancer-type keys, e.g. {'colorectal','breast'},
    plus the sentinel '_ANY_CANCER' when a row was cancer-positive but not
    classified to a specific type."""
    person_id: str
    n_affected_1deg: int = 0
    n_affected_2deg: int = 0
    cancer_by_degree: dict[int, set[str]] = field(default_factory=lambda: {1: set(), 2: set()})

    @property
    def has_1deg_any(self) -> bool:
        return self.n_affected_1deg > 0

    @property
    def has_2deg_any(self) -> bool:
        return self.n_affected_2deg > 0


def read_family_history(
    path: str | Path,
    schema: dict,
    cancer_free_text_patterns: list[str],
    cancer_type_patterns: dict[str, list[str]],
    cancer_phecodes: dict[str, list[str]],
    mapper: IcdToPhecodeX,
    relative_degrees: list[int],
) -> dict[str, PersonFamHx]:
    """Parse Family_History.txt -> {person_id -> PersonFamHx}.

    Detection order per row:
      1. If a coded condition column is present and maps to a cancer phecodeX,
         classify by phecodeX -> cancer key.
      2. Else if free text matches a cancer_type_patterns entry, use that key.
      3. Else if free text matches the generic cancer_free_text_patterns,
         mark the row as `_ANY_CANCER` only.
    """
    out: dict[str, PersonFamHx] = {}
    if not Path(path).exists():
        return out
    cancer_phecodes_all = {p for pcs in cancer_phecodes.values() for p in pcs}
    phecodex_to_cancer = {phe: cancer for cancer, phes in cancer_phecodes.items() for phe in phes}
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=schema.get("sep", "\t"))
        pid_col = schema["person_id_col"]
        rel_col = schema["relationship_col"]
        cond_col = schema["condition_col"]
        code_col = schema.get("condition_code_col")
        degree_col = schema.get("degree_col")
        for row in rdr:
            pid = (row.get(pid_col) or "").strip()
            if not pid:
                continue

            # Degree
            deg: int | None
            if degree_col and row.get(degree_col):
                try:
                    deg = int(str(row[degree_col]).strip())
                except ValueError:
                    deg = relationship_degree(row.get(rel_col, ""))
            else:
                deg = relationship_degree(row.get(rel_col, ""))
            if deg not in relative_degrees:
                continue

            # Cancer detection & classification
            matched: set[str] = set()
            if code_col and row.get(code_col):
                phe = is_cancer_coded(row[code_col], mapper, cancer_phecodes_all)
                if phe and phe in phecodex_to_cancer:
                    matched.add(phecodex_to_cancer[phe])
            if not matched:
                cond_text = row.get(cond_col, "")
                typed = classify_cancer_type(cond_text, cancer_type_patterns)
                if typed:
                    matched.add(typed)
                elif is_cancer_free_text(cond_text, cancer_free_text_patterns):
                    matched.add("_ANY_CANCER")
            if not matched:
                continue

            rec = out.setdefault(pid, PersonFamHx(pid))
            if deg == 1:
                rec.n_affected_1deg += 1
            else:
                rec.n_affected_2deg += 1
            rec.cancer_by_degree[deg].update(matched)
    return out


def famhx_flags(rec: PersonFamHx | None, cancer_phecodes: dict[str, list[str]]
                ) -> dict[str, bool | int]:
    """Convert a PersonFamHx into the flat flag dict for the roster."""
    cancers = list(cancer_phecodes.keys())
    if rec is None:
        d = {"has_famhx_1deg_any_cancer": False, "has_famhx_2deg_any_cancer": False,
             "n_affected_1deg": 0, "n_affected_2deg": 0}
        for cancer in cancers:
            d[f"has_famhx_1deg_{cancer}"] = False
            d[f"has_famhx_2deg_{cancer}"] = False
        return d
    d = {
        "has_famhx_1deg_any_cancer": rec.has_1deg_any,
        "has_famhx_2deg_any_cancer": rec.has_2deg_any,
        "n_affected_1deg": rec.n_affected_1deg,
        "n_affected_2deg": rec.n_affected_2deg,
    }
    for cancer in cancers:
        d[f"has_famhx_1deg_{cancer}"] = cancer in rec.cancer_by_degree[1]
        d[f"has_famhx_2deg_{cancer}"] = cancer in rec.cancer_by_degree[2]
    return d
