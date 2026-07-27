"""Build the phenotype roster: case / famhx / control per (person, cancer).

Usage:
    python -m src.phenotype.roster --config config/config.yaml \
        --out-dir $OUTPUT_ROOT/results/phenotype

Outputs (long format, gitignored):
    cases.tsv        one row per (person, cancer) where is_case=True
    famhx.tsv        one row per person with famhx flags + counts
    roster.tsv       one row per (person, cancer): group in {case,famhx,control}
    roster_summary.md  aggregate counts (no per-sample values)
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import yaml

from src.phenotype.cases import (
    build_cases,
    read_dx_ehr,
    read_phecodex_roster,
)
from src.phenotype.famhx import famhx_flags, read_family_history
from src.phenotype.icd_mapping import IcdToPhecodeX
from src.phenotype.vcf_samples import load_wes_samples


def _phenotypes_dir(cfg: dict) -> Path:
    """Resolve `inputs.phenotypes` dir; on Minerva uses `resolve()`, otherwise best-effort."""
    from datetime import datetime
    root = Path(cfg["inputs"]["msm_data_root"]) / "phenotypes"
    if not root.exists():
        return root
    date = cfg["inputs"]["phenotypes"].get("date")
    if date:
        return root / date
    dated = sorted(
        [p for p in root.iterdir() if p.is_dir() and len(p.name) == 10 and p.name[4] == "-"],
        key=lambda p: p.name,
    )
    return dated[-1] if dated else root


def _read_enrollment(demog_path: Path, schema: dict) -> dict[str, str]:
    if not demog_path.exists() or not schema.get("enrollment_date_col"):
        return {}
    out: dict[str, str] = {}
    with open(demog_path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=schema.get("sep", "\t"))
        pid_col = schema["person_id_col"]
        enr_col = schema["enrollment_date_col"]
        for row in rdr:
            pid = (row.get(pid_col) or "").strip()
            enr = (row.get(enr_col) or "").strip()
            if pid and enr:
                out[pid] = enr[:10]      # ISO date; comparable lexicographically
    return out


def _read_encounter_counts(enc_path: Path, schema: dict) -> Counter:
    """Count all encounter rows per person (any diagnosis)."""
    n = Counter()
    if not enc_path.exists():
        return n
    with open(enc_path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=schema.get("sep", "\t"))
        pid_col = schema["person_id_col"]
        for row in rdr:
            pid = (row.get(pid_col) or "").strip()
            if pid:
                n[pid] += 1
    return n


def assign_group(is_case: bool, has_famhx_any: bool) -> str:
    if is_case:
        return "case"
    return "famhx" if has_famhx_any else "control"


def build(cfg: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    ph_cfg = cfg["phenotype"]
    schema = ph_cfg["ehr_schema"]
    ph_dir = _phenotypes_dir(cfg)
    files = cfg["inputs"]["phenotypes"]["files"]

    mapper = IcdToPhecodeX.load(ph_cfg.get("icd_to_phecodex_map"))
    cancer_phecodes = ph_cfg["cancer_phecodes"]

    # --- Cases ---
    dx_rows = []
    dx_rows += read_dx_ehr(ph_dir / files["encounter_diagnosis"],
                           schema["encounter_diagnosis"], mapper,
                           source="encounter_diagnosis")
    dx_rows += read_dx_ehr(ph_dir / files["problem_list"],
                           schema["problem_list"], mapper,
                           source="problem_list")
    dx_rows += read_phecodex_roster(ph_dir / files["phecodex"],
                                    schema["phecodex_roster"], cancer_phecodes)

    enrollment = _read_enrollment(ph_dir / files["demographics"], schema["demographics"])
    cases = build_cases(dx_rows, cancer_phecodes, enrollment)
    case_pids_any_cancer = {c.person_id for c in cases if c.is_case}
    case_pids_per_cancer: dict[str, set[str]] = {
        cancer: {c.person_id for c in cases if c.is_case and c.cancer == cancer}
        for cancer in cancer_phecodes
    }

    with (out_dir / "cases.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["person_id", "cancer", "is_case", "first_dx_date",
                    "incident_flag", "n_dx_encounters", "sources"])
        for c in cases:
            w.writerow([c.person_id, c.cancer, int(c.is_case),
                        c.first_dx_date or "", "" if c.incident_flag is None else int(c.incident_flag),
                        c.n_dx_encounters, c.sources])

    # --- Family history ---
    famhx_by_pid = read_family_history(
        ph_dir / files["family_history"],
        schema["family_history"],
        ph_cfg["famhx"]["cancer_free_text_patterns"],
        ph_cfg["famhx"].get("cancer_type_patterns", {}),
        cancer_phecodes,
        mapper,
        ph_cfg["famhx"]["relative_degrees"],
    )

    # --- WES sample list (for the "with genotype" overlap counts) ---
    wes_samples = load_wes_samples(cfg)

    with (out_dir / "famhx.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        cancer_cols = []
        for cancer in cancer_phecodes:
            cancer_cols += [f"has_famhx_1deg_{cancer}", f"has_famhx_2deg_{cancer}"]
        header = ["person_id", "has_famhx_1deg_any_cancer", "has_famhx_2deg_any_cancer",
                  "n_affected_1deg", "n_affected_2deg"] + cancer_cols
        w.writerow(header)
        for pid, rec in famhx_by_pid.items():
            flags = famhx_flags(rec, cancer_phecodes)
            w.writerow([pid,
                        int(flags["has_famhx_1deg_any_cancer"]),
                        int(flags["has_famhx_2deg_any_cancer"]),
                        flags["n_affected_1deg"], flags["n_affected_2deg"],
                        *(int(flags[c]) for c in cancer_cols)])

    # --- Roster (per person x cancer) ---
    all_pids = set(enrollment.keys()) | case_pids_any_cancer | set(famhx_by_pid.keys())
    if not all_pids:
        # Fall back to whoever appears in encounter data
        all_pids = {r.person_id for r in dx_rows}

    enc_counts = _read_encounter_counts(
        ph_dir / files["encounter_diagnosis"], schema["encounter_diagnosis"]
    )
    min_enc = ph_cfg["roster"].get("min_control_encounters", 0)
    control_exclusion = ph_cfg["roster"].get("control_exclusion", "any_cancer")
    famhx_scope = ph_cfg["roster"].get("famhx_scope", "any_cancer")

    case_lookup = {(c.person_id, c.cancer): c for c in cases}
    summary = Counter()

    # Track counts per (cancer, group, subset) for cancer_counts.tsv.
    # subset "" = all; "wes" = intersected with wes_samples.
    counts: Counter = Counter()

    with (out_dir / "roster.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "person_id", "cancer", "group", "first_dx_date", "incident_flag",
            "has_wes",
            "has_famhx_1deg_any_cancer", "has_famhx_2deg_any_cancer",
            "has_famhx_1deg_this_cancer", "has_famhx_2deg_this_cancer",
            "n_encounters", "excluded_from_control_reason",
        ])
        for pid in sorted(all_pids):
            fam = famhx_by_pid.get(pid)
            fam_flags = famhx_flags(fam, cancer_phecodes)
            for cancer in cancer_phecodes:
                c = case_lookup.get((pid, cancer))
                is_case = c is not None and c.is_case
                # Determine "has famhx" per scope
                if famhx_scope == "cancer_specific":
                    has_fam = bool(fam_flags.get(f"has_famhx_1deg_{cancer}") or
                                   fam_flags.get(f"has_famhx_2deg_{cancer}"))
                else:
                    has_fam = fam_flags["has_famhx_1deg_any_cancer"] or fam_flags["has_famhx_2deg_any_cancer"]

                group = assign_group(is_case, has_fam)

                # Cross-cancer exclusion: a person with a DIFFERENT cancer isn't
                # a case for this cancer, and shouldn't be counted as famhx or
                # control here either — they're excluded from this cancer's
                # roster entirely (control_exclusion=any_cancer).
                reason = ""
                if (not is_case
                        and control_exclusion == "any_cancer"
                        and pid in case_pids_any_cancer):
                    group = "excluded"
                    reason = "has_other_cancer"
                elif group == "control" and min_enc and enc_counts.get(pid, 0) < min_enc:
                    group = "excluded"
                    reason = "too_few_encounters"

                has_wes = pid in wes_samples
                w.writerow([
                    pid, cancer, group,
                    (c.first_dx_date or "") if c else "",
                    "" if not c or c.incident_flag is None else int(c.incident_flag),
                    int(has_wes),
                    int(fam_flags["has_famhx_1deg_any_cancer"]),
                    int(fam_flags["has_famhx_2deg_any_cancer"]),
                    int(fam_flags.get(f"has_famhx_1deg_{cancer}", False)),
                    int(fam_flags.get(f"has_famhx_2deg_{cancer}", False)),
                    enc_counts.get(pid, 0),
                    reason,
                ])
                summary[(cancer, group)] += 1

                # Cancer-counts aggregation — dedupe per (person, cancer) automatically
                # (one row per (pid, cancer)).
                if is_case:
                    counts[(cancer, "ehr_cases", "all")] += 1
                    if has_wes:
                        counts[(cancer, "ehr_cases", "wes")] += 1
                if fam_flags["has_famhx_1deg_any_cancer"]:
                    counts[(cancer, "famhx_1deg_any", "all")] += 1
                    if has_wes:
                        counts[(cancer, "famhx_1deg_any", "wes")] += 1
                if fam_flags.get(f"has_famhx_1deg_{cancer}", False):
                    counts[(cancer, "famhx_1deg_this", "all")] += 1
                    if has_wes:
                        counts[(cancer, "famhx_1deg_this", "wes")] += 1
                if group == "control":
                    counts[(cancer, "controls", "all")] += 1
                    if has_wes:
                        counts[(cancer, "controls", "wes")] += 1

    # --- Summary (aggregate counts only) ---
    with (out_dir / "roster_summary.md").open("w") as fh:
        fh.write("# Phenotype roster — summary\n\n")
        fh.write("| cancer | case | famhx | control | excluded |\n")
        fh.write("|---|---:|---:|---:|---:|\n")
        for cancer in cancer_phecodes:
            fh.write(
                f"| {cancer} "
                f"| {summary[(cancer, 'case')]} "
                f"| {summary[(cancer, 'famhx')]} "
                f"| {summary[(cancer, 'control')]} "
                f"| {summary[(cancer, 'excluded')]} |\n"
            )

    # --- Cancer counts (the "clean table with the numbers") ---
    # One row per cancer; columns dedup per participant.
    with (out_dir / "cancer_counts.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "cancer",
            "n_ehr_cases",
            "n_ehr_cases_with_wes",
            "n_famhx_1deg_any_cancer",
            "n_famhx_1deg_any_cancer_with_wes",
            "n_famhx_1deg_this_cancer",
            "n_famhx_1deg_this_cancer_with_wes",
            "n_controls",
            "n_controls_with_wes",
        ])
        for cancer in cancer_phecodes:
            w.writerow([
                cancer,
                counts[(cancer, "ehr_cases", "all")],
                counts[(cancer, "ehr_cases", "wes")],
                counts[(cancer, "famhx_1deg_any", "all")],
                counts[(cancer, "famhx_1deg_any", "wes")],
                counts[(cancer, "famhx_1deg_this", "all")],
                counts[(cancer, "famhx_1deg_this", "wes")],
                counts[(cancer, "controls", "all")],
                counts[(cancer, "controls", "wes")],
            ])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    build(cfg, Path(args.out_dir))
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
