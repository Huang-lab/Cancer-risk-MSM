"""Build the EHR cohort table: case/control + family history + lifestyle.

    # 1. Confirm the file schemas before trusting any numbers
    python -m src.phenotype.cohort --config config/config.yaml --inspect

    # 2. Build the cohort table
    python -m src.phenotype.cohort --config config/config.yaml \
        --out-dir $OUTPUT_ROOT/results/phenotype

Outputs (gitignored, Minerva only):
    cohort.tsv          one row per person: group, cancers, famhx, lifestyle
    cohort_summary.md   aggregate counts + covariate missingness (no per-person values)

`--inspect` prints column names and row counts only — never cell values — so its
output is safe to paste back for schema confirmation.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import date
from pathlib import Path

import yaml

from src.phenotype import covariates as cov
from src.phenotype.cancer_cases import (
    cancer_types_present,
    describe_schema,
    load_cancer_cases,
)
from src.phenotype.controls import (
    assign_groups,
    encounter_counts,
    people_with_cancer_icd,
    read_person_ids,
)
from src.phenotype.famhx import famhx_flags, read_family_history
from src.phenotype.icd_mapping import IcdToPhecodeX


def phenotypes_dir(cfg: dict) -> Path:
    """Latest date-stamped phenotypes subfolder (auto-advances on new extracts)."""
    root = Path(cfg["inputs"]["msm_data_root"]) / "phenotypes"
    date_pin = cfg["inputs"]["phenotypes"].get("date")
    if date_pin:
        return root / date_pin
    if not root.exists():
        return root
    dated = sorted(p for p in root.iterdir()
                   if p.is_dir() and len(p.name) == 10 and p.name[4] == "-")
    return dated[-1] if dated else root


def _paths(cfg: dict) -> dict[str, Path]:
    ph = phenotypes_dir(cfg)
    f = cfg["inputs"]["phenotypes"]["files"]
    return {
        "cancer_cases": Path(cfg["inputs"]["cancer_cases"]["path"]),
        "demographics": ph / f["demographics"],
        "family_history": ph / f["family_history"],
        "social_history": ph / f["social_history"],
        "vitals": ph / f["vitals"],
        "ob_history": ph / f["ob_history"],
        "encounter_diagnosis": ph / f["encounter_diagnosis"],
        "problem_list": ph / f["problem_list"],
    }


def inspect(cfg: dict) -> int:
    """Report file presence + schemas. No cell values are emitted."""
    paths = _paths(cfg)
    print(f"phenotypes dir: {phenotypes_dir(cfg)}\n")

    cc = cfg["inputs"]["cancer_cases"]
    print("=== cancer case file ===")
    info = describe_schema(cc["path"], cc.get("sep", "\t"))
    if "error" in info:
        print(f"  ERROR: {info['error']}")
    else:
        print(f"  rows: {info['n_rows']:,}   columns: {info['n_columns']}")
        print(f"  columns: {info['columns']}")
        print(f"  detected cancer-type column: {info['detected_cancer_type_col']}")
        print(f"  detected genetics column:    {info['detected_genetics_col']}")
        if cc["person_id_col"] not in info["columns"]:
            print(f"  WARNING: configured ID column {cc['person_id_col']!r} not present")

    print("\n=== EHR files ===")
    schemas = cfg["phenotype"]["ehr_schema"]
    for key in ["demographics", "family_history", "social_history", "vitals",
                "ob_history", "encounter_diagnosis", "problem_list"]:
        p = paths[key]
        sep = schemas.get(key, {}).get("sep", "|")
        if not p.exists():
            print(f"  {key:22s} MISSING  {p}")
            continue
        with p.open(newline="") as fh:
            rdr = csv.DictReader(fh, delimiter=sep)
            cols = list(rdr.fieldnames or [])
            n = sum(1 for _ in rdr)
        print(f"  {key:22s} rows={n:<10,} cols={len(cols)}")
        print(f"    {cols}")
        expected = schemas.get(key, {})
        missing = [v for k, v in expected.items()
                   if k.endswith("_col") and v and v not in cols]
        if missing:
            print(f"    WARNING: configured columns absent: {missing}")
    return 0


def build(cfg: dict, out_dir: Path) -> int:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = _paths(cfg)
    schemas = cfg["phenotype"]["ehr_schema"]
    coh = cfg.get("cohort", {})
    cc_cfg = cfg["inputs"]["cancer_cases"]

    # --- 1. Cancer cases (authoritative) ---
    cases, cc_report = load_cancer_cases(
        cc_cfg["path"], cc_cfg.get("person_id_col", "sample_name"),
        cc_cfg.get("cancer_type_col"), cc_cfg.get("sep", "\t"))
    print(f"cancer cases: {cc_report['n_people']:,} people "
          f"from {cc_report['n_rows']:,} rows")
    if cc_report["cancer_type_col"]:
        print(f"  cancer-type column: {cc_report['cancer_type_col']}")
        print(f"  types: {cc_report['types_seen']}")
    for w in cc_report["warnings"]:
        print(f"  WARNING: {w}")
    if cc_report["unclassified_type_examples"]:
        print(f"  unclassified type strings ({cc_report['n_rows_unclassified_type']} rows): "
              f"{cc_report['unclassified_type_examples']}")

    # --- 2. Roster + control assignment ---
    demog_schema = schemas["demographics"]
    pid_col = demog_schema.get("person_id_col", "sample_name")
    if coh.get("roster_source", "demographics") == "union":
        roster: set[str] = set()
        for key in ["demographics", "encounter_diagnosis", "problem_list"]:
            roster |= read_person_ids(paths[key], schemas[key].get("person_id_col", pid_col),
                                      schemas[key].get("sep", "|"))
    else:
        roster = read_person_ids(paths["demographics"], pid_col,
                                 demog_schema.get("sep", "|"))
    print(f"roster ({coh.get('roster_source','demographics')}): {len(roster):,} people")

    icd_cancer: set[str] = set()
    if coh.get("screen_icd_for_cancer", True):
        mapper = IcdToPhecodeX.load(cfg["phenotype"].get("icd_to_phecodex_map"))
        icd_cancer = people_with_cancer_icd(
            [(paths["encounter_diagnosis"], schemas["encounter_diagnosis"]),
             (paths["problem_list"], schemas["problem_list"])], mapper)
        print(f"people with a cancer ICD code: {len(icd_cancer):,}")

    enc_counts = {}
    min_enc = coh.get("min_control_encounters", 0)
    if min_enc:
        enc_counts = encounter_counts(
            paths["encounter_diagnosis"],
            schemas["encounter_diagnosis"].get("person_id_col", pid_col),
            schemas["encounter_diagnosis"].get("sep", "|"))

    groups, roster_rep = assign_groups(roster, set(cases), icd_cancer,
                                       enc_counts, min_enc)
    print(f"  cases={roster_rep.n_cases:,}  controls={roster_rep.n_controls:,}  "
          f"excluded_icd_cancer={roster_rep.n_excluded_icd_cancer:,}  "
          f"excluded_few_encounters={roster_rep.n_too_few_encounters:,}")
    for w in roster_rep.warnings:
        print(f"  WARNING: {w}")

    # --- 3. Family history (first-degree, everyone) ---
    ph_cfg = cfg["phenotype"]
    cancer_keys = sorted(set(cancer_types_present(cases))
                         | set(ph_cfg.get("cancer_phecodes", {})))
    fam_key_map = {k: [k] for k in cancer_keys}   # famhx_flags enumerates by key
    fam = read_family_history(
        paths["family_history"], schemas["family_history"],
        ph_cfg["famhx"]["cancer_free_text_patterns"],
        ph_cfg["famhx"].get("cancer_type_patterns", {}),
        fam_key_map, IcdToPhecodeX.load(None),
        relative_degrees=[1],                     # first-degree per request
    )
    print(f"family history: {len(fam):,} people with >=1 first-degree cancer relative")

    # --- 4. Lifestyle + demographics ---
    as_of = date.today().year
    covs = cov.read_demographics(paths["demographics"], demog_schema, as_of_year=as_of)
    social = cov.read_social_history(paths["social_history"], schemas["social_history"],
                                    coh.get("smoking_map", {}), coh.get("alcohol_map", {}))
    bmis = cov.read_bmi(paths["vitals"], schemas["vitals"])
    obs = cov.read_ob_history(paths["ob_history"], schemas["ob_history"])
    print(f"covariates: demographics={len(covs):,} social={len(social):,} "
          f"bmi={len(bmis):,} ob={len(obs):,}")

    # Merge the per-source lifestyle values back onto the Covariates records so
    # missingness() reports the assembled table rather than demographics alone.
    for pid in set(covs) | set(social) | set(bmis) | set(obs):
        rec = covs.setdefault(pid, cov.Covariates(person_id=pid))
        if pid in social:
            rec.smoking_status, rec.alcohol_use = social[pid]
        if pid in bmis:
            rec.bmi, _bmi_d = bmis[pid]
            rec.bmi_date = _bmi_d.isoformat() if _bmi_d else None
        if pid in obs:
            rec.parity, rec.age_at_first_birth = obs[pid]


    # --- 5. Emit cohort.tsv ---
    all_pids = sorted(set(groups) | set(covs) | set(fam))
    fam_cols = [f"has_famhx_1deg_{c}" for c in cancer_keys]
    header = (["sample_name", "group", "cancer_types", "n_cancers",
               "has_famhx_1deg_any_cancer", "n_affected_1deg"] + fam_cols +
              ["age", "sex", "race_ethnicity", "smoking_status", "alcohol_use",
               "bmi", "bmi_date", "parity", "age_at_first_birth"])

    group_counts: dict[str, int] = {}
    famhx_any_by_group: dict[str, int] = {}

    with (out_dir / "cohort.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for pid in all_pids:
            grp = groups.get(pid, "not_in_roster")
            group_counts[grp] = group_counts.get(grp, 0) + 1
            case = cases.get(pid)
            flags = famhx_flags(fam.get(pid), fam_key_map)
            if flags["has_famhx_1deg_any_cancer"]:
                famhx_any_by_group[grp] = famhx_any_by_group.get(grp, 0) + 1
            c = covs.get(pid, cov.Covariates(person_id=pid))
            smoke, alc = social.get(pid, ("unknown", "unknown"))
            bmi, bmi_d = bmis.get(pid, (None, None))
            parity, afb = obs.get(pid, (None, None))
            w.writerow([
                pid, grp,
                ";".join(sorted(case.cancer_types)) if case else "",
                case.n_cancers if case else 0,
                int(flags["has_famhx_1deg_any_cancer"]), flags["n_affected_1deg"],
                *(int(flags.get(fc, False)) for fc in fam_cols),
                "" if c.age is None else f"{c.age:.0f}",
                c.sex or "", c.race_ethnicity or "",
                smoke, alc,
                "" if bmi is None else f"{bmi:.1f}",
                bmi_d.isoformat() if bmi_d else "",
                "" if parity is None else parity,
                "" if afb is None else f"{afb:.0f}",
            ])

    # --- 6. Summary (aggregate only) ---
    miss = cov.missingness(covs) if covs else {}
    with (out_dir / "cohort_summary.md").open("w") as fh:
        fh.write("# EHR cohort summary\n\n")
        fh.write(f"- cancer case file: `{cc_cfg['path']}`\n")
        fh.write(f"- phenotypes dir: `{phenotypes_dir(cfg)}`\n")
        fh.write(f"- roster source: {coh.get('roster_source','demographics')} "
                 f"({roster_rep.n_roster:,} people)\n")
        fh.write(f"- ICD cancer screen: {'on' if coh.get('screen_icd_for_cancer', True) else 'off'}\n\n")
        fh.write("## Groups\n\n| group | n | with 1st-deg cancer famhx |\n|---|---:|---:|\n")
        for g, n in sorted(group_counts.items(), key=lambda kv: -kv[1]):
            fh.write(f"| {g} | {n:,} | {famhx_any_by_group.get(g, 0):,} |\n")
        fh.write("\n## Cancer types among cases\n\n| type | n rows |\n|---|---:|\n")
        for t, n in cc_report["types_seen"].items():
            fh.write(f"| {t} | {n:,} |\n")
        if miss:
            fh.write("\n## Covariate missingness\n\n| covariate | fraction missing |\n|---|---:|\n")
            for k, v in sorted(miss.items(), key=lambda kv: -kv[1]):
                fh.write(f"| {k} | {v:.1%} |\n")
        if roster_rep.warnings or cc_report["warnings"]:
            fh.write("\n## Warnings\n\n")
            for wmsg in roster_rep.warnings + cc_report["warnings"]:
                fh.write(f"- {wmsg}\n")

    print(f"\nWrote {out_dir/'cohort.tsv'} ({len(all_pids):,} rows)")
    print(f"Wrote {out_dir/'cohort_summary.md'}")
    if miss:
        worst = sorted(miss.items(), key=lambda kv: -kv[1])[:4]
        print("highest covariate missingness: " +
              ", ".join(f"{k}={v:.0%}" for k, v in worst))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--inspect", action="store_true",
                    help="report file schemas (columns + row counts only) and exit")
    args = ap.parse_args(argv)

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    if args.inspect:
        return inspect(cfg)

    out = Path(args.out_dir) if args.out_dir else \
        Path(cfg["project"]["output_root"]) / "results" / "phenotype"
    return build(cfg, out)


if __name__ == "__main__":
    sys.exit(main())
