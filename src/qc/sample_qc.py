"""Sample-level QC — build $OUTPUT_ROOT/results/qc/sample_keep_list.tsv from
stats computed directly on the annotated WES VCFs.

Inputs (all produced by workflow/00_wes_qc_*.sh, from VCFs):
    - $OUTPUT_ROOT/results/qc/psc_agg.tsv        (aggregated bcftools PSC stats)
    - $OUTPUT_ROOT/results/qc/king.kin0          (plink2 --make-king output)
    - $OUTPUT_ROOT/results/qc/sexcheck.sexcheck  (plink2 --check-sex output)
    - $OUTPUT_ROOT/results/qc/reported_sex.tsv   (person_id, reported_sex)

Checks (thresholds from config.qc.sample):
  1. Missingness   -> from PSC (n_missing / n_total).
  2. Het/hom ratio -> from PSC; flag samples with |z-score| > het_hom_z_max
                       (weak contamination proxy usable without BAM/CRAM).
  3. Sex-check     -> from plink2 --check-sex F-stat; mismatch_action = flag|exclude.
  4. Kinship       -> from plink2 --make-king. For pairs with kinship >=
                       duplicate_kin_min, keep the sample with the lower
                       missingness; related pairs (>= pi_hat_related_max) kept
                       but flagged; ML uses group-aware CV downstream.

Outputs (all gitignored; on Minerva only):
    - sample_keep_list.tsv        one row per KEPT sample_id
    - sample_qc_report.tsv        per-sample flags + reasons
    - sample_qc_summary.md        aggregate counts (no per-sample values)
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class SampleQcRow:
    sample_id: str
    call_rate: float | None
    het_hom: float | None
    sex_reported: str | None
    sex_inferred: str | None
    sex_mismatch: bool
    kinship_flag: str          # "" | "duplicate_dropped" | "duplicate_kept" | "related_kept"
    kept: bool
    reasons: str               # comma-separated failed checks


# --- Building blocks (unit-testable) ---------------------------------------

def missingness_fail(missingness: float, missingness_max: float) -> bool:
    return missingness > missingness_max


def sex_mismatch(reported: str | None, inferred: str | None) -> bool:
    if not reported or not inferred:
        return False
    r = reported.strip().upper()[:1]      # "M"/"F"
    i = inferred.strip().upper()[:1]
    return r in ("M", "F") and i in ("M", "F") and r != i


def sex_from_f_stat(f_stat: float, female_max: float, male_min: float) -> str | None:
    """plink2 F-stat convention: F<female_max -> Female; F>male_min -> Male."""
    if f_stat < female_max:
        return "F"
    if f_stat > male_min:
        return "M"
    return None                            # ambiguous


def het_hom_outlier(sample_ratio: float, cohort_mean: float,
                    cohort_sd: float, z_max: float) -> bool:
    """Contamination proxy: excess het vs cohort. Requires cohort_sd > 0."""
    if cohort_sd <= 0:
        return False
    z = (sample_ratio - cohort_mean) / cohort_sd
    return abs(z) > z_max


def is_duplicate(kinship: float, duplicate_min: float) -> bool:
    return kinship >= duplicate_min


def pick_from_duplicate_cluster(cluster: list[SampleQcRow], strategy: str) -> str:
    """Return sample_id to KEEP within a duplicate/MZ-twin cluster."""
    if strategy == "highest_call_rate":
        return max(cluster, key=lambda r: (r.call_rate or 0.0)).sample_id
    return cluster[0].sample_id


# --- Cluster building from pairwise kinship ---------------------------------

def build_duplicate_clusters(pairs: list[tuple[str, str, float]],
                             duplicate_min: float) -> list[set[str]]:
    """Union-find over duplicate pairs (kinship >= duplicate_min).

    Returns connected components (clusters of samples that are all mutually
    or transitively duplicates).
    """
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent.get(x, x), parent.get(x, x))
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    nodes: set[str] = set()
    for a, b, k in pairs:
        if k < duplicate_min:
            continue
        parent.setdefault(a, a)
        parent.setdefault(b, b)
        nodes.add(a); nodes.add(b)
        union(a, b)

    comps: dict[str, set[str]] = {}
    for n in nodes:
        r = find(n)
        comps.setdefault(r, set()).add(n)
    return [c for c in comps.values() if len(c) > 1]


# TODO(minerva): run(cfg) — reads psc_agg.tsv, king.kin0, sexcheck.sexcheck,
# reported_sex.tsv; applies the checks above; writes the three output files.
