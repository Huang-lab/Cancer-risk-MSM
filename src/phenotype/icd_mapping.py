"""ICD-9-CM / ICD-10-CM  ->  phecodeX mapping.

If `phenotype.icd_to_phecodex_map` is set, loads that TSV (columns:
`icd, phecodeX, icd_version`). Otherwise falls back to a built-in map that
covers the cancer phecodeX codes referenced in `config.phenotype.cancer_phecodes`.
"""
from __future__ import annotations

import re
from pathlib import Path

# --- Built-in cancer map ----------------------------------------------------
# Only what we need for CRC + breast. Extend on cluster with the full map file.
# Keys are ICD prefixes (checked longest-first); values are phecodeX codes.
# References: WHO ICD-10, CMS ICD-10-CM, PheWAS-X (Wei et al.) v1.2.
_BUILTIN_MAP: dict[str, list[tuple[str, str]]] = {
    "10": [
        # Colorectal
        ("C18",   "CA_101.1"),
        ("C19",   "CA_101.11"),
        ("C20",   "CA_101.11"),
        ("C21",   "CA_101.11"),
        ("D01.0", "CA_101.1"),      # carcinoma in situ colon
        ("D01.1", "CA_101.11"),     # rectosigmoid
        ("D01.2", "CA_101.11"),     # rectum
        # Breast
        ("C50",   "CA_105.1"),
        ("D05",   "CA_105.1"),      # carcinoma in situ breast
        # Prostate
        ("C61",   "CA_105.2"),
        ("D07.5", "CA_105.2"),
        # Lung
        ("C33",   "CA_101.6"),
        ("C34",   "CA_101.6"),
        # Ovarian
        ("C56",   "CA_105.3"),
        ("D07.39", "CA_105.3"),
        # Pancreatic
        ("C25",   "CA_101.2"),
    ],
    "9": [
        # Colorectal
        ("153",   "CA_101.1"),
        ("154",   "CA_101.11"),
        ("230.3", "CA_101.1"),
        ("230.4", "CA_101.11"),
        # Breast
        ("174",   "CA_105.1"),
        ("175",   "CA_105.1"),
        ("233.0", "CA_105.1"),
        # Prostate
        ("185",   "CA_105.2"),
        ("233.4", "CA_105.2"),
        # Lung
        ("162",   "CA_101.6"),
        # Ovarian
        ("183.0", "CA_105.3"),
        # Pancreatic
        ("157",   "CA_101.2"),
    ],
}

_ICD10_RE = re.compile(r"^[A-TV-Z][0-9][0-9AB]")
_ICD9_RE = re.compile(r"^\d{3}(\.\d+)?$|^V\d{2}")


def infer_icd_version(code: str) -> str:
    """Return '9' or '10'. Best-effort; use icd_version_col from EHR if present."""
    code = (code or "").strip().upper().replace(" ", "")
    if _ICD10_RE.match(code):
        return "10"
    if _ICD9_RE.match(code):
        return "9"
    # Default to 10 (post-2015 US clinical data is overwhelmingly ICD-10-CM).
    return "10"


class IcdToPhecodeX:
    """Prefix-matching lookup. `map_code(code, version) -> phecodeX | None`.

    Longest-prefix match wins so "C18.4" -> CA_101.1 (matches "C18") and
    "C19.0" -> CA_101.11 (matches "C19").
    """

    def __init__(self, entries: dict[str, list[tuple[str, str]]]):
        # entries[version] = list of (prefix, phecodeX), sorted long-first.
        self._by_version: dict[str, list[tuple[str, str]]] = {
            v: sorted(pairs, key=lambda p: -len(p[0])) for v, pairs in entries.items()
        }

    @classmethod
    def load(cls, path: str | None) -> "IcdToPhecodeX":
        if not path:
            return cls(_BUILTIN_MAP)
        entries: dict[str, list[tuple[str, str]]] = {"9": [], "10": []}
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i_icd = header.index("icd")
            i_phe = header.index("phecodeX")
            i_ver = header.index("icd_version") if "icd_version" in header else None
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(i_icd, i_phe):
                    continue
                icd = f[i_icd].strip().upper()
                phe = f[i_phe].strip()
                ver = f[i_ver].strip() if i_ver is not None else infer_icd_version(icd)
                entries.setdefault(ver, []).append((icd, phe))
        return cls(entries)

    def map_code(self, code: str, version: str | None = None) -> str | None:
        code = (code or "").strip().upper().replace(" ", "")
        if not code:
            return None
        v = version or infer_icd_version(code)
        for prefix, phe in self._by_version.get(v, ()):
            if code.startswith(prefix):
                return phe
        return None
