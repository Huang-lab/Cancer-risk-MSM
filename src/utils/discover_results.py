"""Discover and rank germline-plp-carrier-nf output directories.

Scans a project root for `results-*/` dirs, reads each `manifest.json`, derives
coverage/count metrics from the carrier matrix, classifies each run as a
genome-wide `batch` or a single-chromosome `pilot`, and auto-selects the latest
complete batch run as the default carrier source.

Usage
-----
    python -m src.utils.discover_results \
        --project-root /sc/arion/projects/rg_huangk06/variants_PLP_MSM \
        --out config/discovered_results.yaml

What comes from where
---------------------
`manifest.json` carries PROVENANCE ONLY (reference build, VEP + cache version,
ClinVar release, AlphaMissense version, gnomAD version, pipeline git SHA,
container digests). It does NOT record run date, batch name, chromosomes
covered, or any counts. Those are derived here:

    run_date      mtime of manifest.json (fallback: carrier_matrix.tsv)
    batch         parsed from the folder name
    chrs_covered  distinct chr column in carrier_matrix.tsv
    n_variants    distinct (chr,pos,ref,alt)
    n_carriers    distinct person_id
    n_samples     NOT derivable from the carrier matrix -- see below

IMPORTANT — n_carriers is not n_samples. The carrier matrix contains only
people with at least one P/LP call, so its distinct person_id count is the
number of *carriers*. Reporting that as cohort size would understate the
denominator of every downstream rate. `n_samples` is therefore left blank
unless an explicit sample list is found, and the table says so.

The manifest reader is deliberately tolerant: unknown keys are recorded rather
than rejected, so a schema change surfaces as data in the report instead of a
crash.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

# Canonical chromosome sets. NF emits Ensembl-style names ("1"), but we accept
# "chr1" too so a naming change doesn't silently classify everything as pilot.
AUTOSOMES = {str(i) for i in range(1, 23)}
ALL_CHROMS = AUTOSOMES | {"X", "Y", "M", "MT"}

# Files we expect a complete run to have produced.
EXPECTED_FILES = [
    "carriers/carrier_matrix.tsv",
    "variants/clinvar_plp.tsv",
    "variants/acmg_plp.tsv",
    "variants/am_plp.tsv",
    "manifest.json",
]

# Candidate locations for an explicit sample list (cohort denominator).
SAMPLE_LIST_CANDIDATES = [
    "qc/sample_keep_list.tsv",
    "qc/samples.tsv",
    "qc/sample_list.txt",
    "samples.txt",
]

_PILOT_NAME_RE = re.compile(r"(results-chr|[-_]single|pilot)", re.I)
_BATCH_NUM_RE = re.compile(r"batch[-_]?(\d+)", re.I)


def normalize_chrom(raw: str) -> str:
    """'chr1' -> '1', 'CHRX' -> 'X'. Leaves unrecognized values as-is."""
    c = (raw or "").strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c.upper() if c.upper() in ("X", "Y", "M", "MT") else c


@dataclass
class ResultInfo:
    folder: str
    path: Path
    kind: str = "unknown"                  # "batch" | "pilot" | "unknown"
    batch: str = ""
    batch_num: int | None = None
    chrs_covered: list[str] = field(default_factory=list)
    n_variants: int | None = None
    n_carriers: int | None = None
    n_samples: int | None = None           # None => not determinable
    n_carrier_rows: int | None = None
    run_date: str = ""                     # ISO date, derived from mtime
    manifest: dict = field(default_factory=dict)
    missing_files: list[str] = field(default_factory=list)
    complete: bool = False
    notes: list[str] = field(default_factory=list)

    # --- convenience accessors over the tolerant manifest -------------------

    def _m(self, *keys: str) -> str:
        """First present manifest key among aliases, stringified."""
        for k in keys:
            v = self.manifest.get(k)
            if v not in (None, ""):
                return str(v)
        return ""

    @property
    def clinvar_release(self) -> str:
        return self._m("clinvar_release", "clinvar_version", "clinvar", "clinvar_date")

    @property
    def vep_version(self) -> str:
        return self._m("vep_version", "vep", "vep_cache_version")

    @property
    def reference_build(self) -> str:
        return self._m("reference_build", "genome", "assembly", "build")

    @property
    def pipeline_sha(self) -> str:
        return self._m("pipeline_sha", "pipeline_git_sha", "git_sha", "commit")

    @property
    def chrs_display(self) -> str:
        if not self.chrs_covered:
            return "-"
        return summarize_chroms(self.chrs_covered)


def summarize_chroms(chroms: list[str]) -> str:
    """['1','2','3','X'] -> '1-3,X'. Keeps the table narrow."""
    nums = sorted((int(c) for c in chroms if c.isdigit()))
    others = sorted(c for c in chroms if not c.isdigit())
    parts: list[str] = []
    start = prev = None
    for n in nums:
        if start is None:
            start = prev = n
            continue
        if n == prev + 1:
            prev = n
            continue
        parts.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = n
    if start is not None:
        parts.append(str(start) if start == prev else f"{start}-{prev}")
    return ",".join(parts + others)


def read_manifest(path: Path) -> tuple[dict, str | None]:
    """Tolerant manifest read. Returns (flattened_dict, error_or_None).

    Nested objects are flattened one level with dotted keys so provenance like
    {"tool_versions": {"vep": "112"}} becomes {"tool_versions.vep": "112"} and
    still shows up in the report.
    """
    if not path.exists():
        return {}, "manifest.json missing"
    try:
        raw = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as e:
        return {}, f"manifest.json unreadable ({type(e).__name__})"
    if not isinstance(raw, dict):
        return {}, "manifest.json is not a JSON object"

    flat: dict = {}
    for k, v in raw.items():
        if isinstance(v, dict):
            for k2, v2 in v.items():
                flat[f"{k}.{k2}"] = v2
            flat[k] = v
        else:
            flat[k] = v
    return flat, None


def scan_carrier_matrix(path: Path) -> dict:
    """Single streaming pass over the carrier matrix.

    Returns chrs / n_variants / n_carriers / n_rows. Does not load the file into
    memory (these can be large), and tolerates a missing chr/person_id column by
    reporting what it could not find rather than raising.
    """
    out: dict = {"chrs": [], "n_variants": None, "n_carriers": None,
                 "n_rows": None, "error": None}
    if not path.exists():
        out["error"] = "carrier_matrix.tsv missing"
        return out

    chroms: set[str] = set()
    variants: set[tuple[str, str, str, str]] = set()
    people: set[str] = set()
    rows = 0
    try:
        with path.open(newline="") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            cols = rdr.fieldnames or []
            c_chr = _pick(cols, ["chr", "chrom", "CHROM", "#CHROM"])
            c_pos = _pick(cols, ["pos", "POS", "position"])
            c_ref = _pick(cols, ["ref", "REF"])
            c_alt = _pick(cols, ["alt", "ALT"])
            c_pid = _pick(cols, ["person_id", "sample_id", "sample", "IID"])
            if c_chr is None or c_pid is None:
                out["error"] = f"unexpected columns: {cols[:8]}"
                return out
            for row in rdr:
                rows += 1
                ch = normalize_chrom(row.get(c_chr, ""))
                if ch:
                    chroms.add(ch)
                pid = (row.get(c_pid) or "").strip()
                if pid:
                    people.add(pid)
                if c_pos and c_ref and c_alt:
                    variants.add((ch, (row.get(c_pos) or "").strip(),
                                  (row.get(c_ref) or "").strip(),
                                  (row.get(c_alt) or "").strip()))
    except OSError as e:
        out["error"] = f"carrier_matrix.tsv unreadable ({type(e).__name__})"
        return out

    out["chrs"] = sorted(chroms)
    out["n_rows"] = rows
    out["n_carriers"] = len(people)
    out["n_variants"] = len(variants) if variants else None
    return out


def _pick(cols: list[str], names: list[str]) -> str | None:
    for n in names:
        if n in cols:
            return n
    return None


def find_sample_count(result_dir: Path) -> tuple[int | None, str]:
    """Look for an explicit sample list to use as the cohort denominator.

    Returns (count, source). (None, "") when no list exists — we do NOT fall
    back to the carrier count, which would be a different quantity.
    """
    for rel in SAMPLE_LIST_CANDIDATES:
        p = result_dir / rel
        if not p.exists():
            continue
        try:
            lines = [ln.strip() for ln in p.read_text().splitlines() if ln.strip()]
        except OSError:
            continue
        if not lines:
            continue
        # Drop a header row if the first line looks like one.
        if re.search(r"sample|id", lines[0], re.I) and len(lines) > 1:
            lines = lines[1:]
        return len(lines), rel
    return None, ""


def classify(folder: str, chrs: list[str]) -> str:
    """`pilot` for single-chromosome / pilot runs, else `batch`.

    Chromosome coverage is the primary signal; the folder name is the fallback
    for runs whose carrier matrix could not be read.
    """
    if chrs:
        return "pilot" if len(chrs) <= 1 else "batch"
    return "pilot" if _PILOT_NAME_RE.search(folder) else "batch"


def parse_batch(folder: str) -> tuple[str, int | None]:
    m = _BATCH_NUM_RE.search(folder)
    if m:
        return f"batch{m.group(1)}", int(m.group(1))
    label = folder[len("results-"):] if folder.startswith("results-") else folder
    return label, None


def derive_run_date(result_dir: Path) -> str:
    for rel in ("manifest.json", "carriers/carrier_matrix.tsv"):
        p = result_dir / rel
        if p.exists():
            try:
                ts = p.stat().st_mtime
            except OSError:
                continue
            return datetime.fromtimestamp(ts, tz=timezone.utc).date().isoformat()
    return ""


def inspect(result_dir: Path) -> ResultInfo:
    folder = result_dir.name
    info = ResultInfo(folder=folder, path=result_dir)

    info.manifest, man_err = read_manifest(result_dir / "manifest.json")
    if man_err:
        info.notes.append(man_err)

    cm = scan_carrier_matrix(result_dir / "carriers" / "carrier_matrix.tsv")
    if cm["error"]:
        info.notes.append(cm["error"])
    info.chrs_covered = cm["chrs"]
    info.n_variants = cm["n_variants"]
    info.n_carriers = cm["n_carriers"]
    info.n_carrier_rows = cm["n_rows"]

    n_samp, src = find_sample_count(result_dir)
    info.n_samples = n_samp
    if n_samp is None:
        info.notes.append("no sample list found; n_samples unknown (n_carriers is NOT cohort size)")
    else:
        info.notes.append(f"n_samples from {src}")

    info.kind = classify(folder, info.chrs_covered)
    info.batch, info.batch_num = parse_batch(folder)
    info.run_date = derive_run_date(result_dir)

    info.missing_files = [f for f in EXPECTED_FILES if not (result_dir / f).exists()]
    covers_autosomes = AUTOSOMES.issubset(set(info.chrs_covered))
    info.complete = (not info.missing_files) and covers_autosomes
    if info.missing_files:
        info.notes.append(f"missing: {','.join(info.missing_files)}")
    elif not covers_autosomes and info.kind == "batch":
        missing_chr = sorted(AUTOSOMES - set(info.chrs_covered), key=int)
        info.notes.append(f"autosomes not fully covered (missing {summarize_chroms(missing_chr)})")
    return info


def discover(project_root: Path) -> list[ResultInfo]:
    dirs = sorted(p for p in project_root.glob("results-*") if p.is_dir())
    return [inspect(d) for d in dirs]


def rank_batches(infos: list[ResultInfo]) -> list[ResultInfo]:
    """Complete batch runs first, then by run date desc, then batch number desc."""
    batches = [i for i in infos if i.kind == "batch"]
    return sorted(
        batches,
        key=lambda i: (i.complete, i.run_date, i.batch_num or -1),
        reverse=True,
    )


def select_best(infos: list[ResultInfo]) -> ResultInfo | None:
    ranked = rank_batches(infos)
    return ranked[0] if ranked else None


# --- Reporting --------------------------------------------------------------

_COLS = [
    ("folder", 24), ("type", 6), ("batch", 8), ("chrs", 12),
    ("n_variants", 11), ("n_carriers", 11), ("n_samples", 10),
    ("clinvar", 12), ("vep", 6), ("run_date", 10), ("complete", 8),
]


def _fmt(v) -> str:
    if v is None:
        return "?"
    if isinstance(v, int):
        return f"{v:,}"
    return str(v) if str(v) else "-"


def render_table(infos: list[ResultInfo], best: ResultInfo | None) -> str:
    lines = []
    header = "  ".join(name.ljust(w) for name, w in _COLS)
    lines.append(header)
    lines.append("  ".join("-" * w for _, w in _COLS))
    for i in infos:
        marker = " *" if best is not None and i.folder == best.folder else ""
        cells = [
            i.folder + marker, i.kind, i.batch, i.chrs_display,
            _fmt(i.n_variants), _fmt(i.n_carriers), _fmt(i.n_samples),
            i.clinvar_release or "-", i.vep_version or "-",
            i.run_date or "-", "Y" if i.complete else "N",
        ]
        lines.append("  ".join(str(c).ljust(w)[:w] for c, (_, w) in zip(cells, _COLS)))
    return "\n".join(lines)


def render_notes(infos: list[ResultInfo]) -> str:
    out = []
    for i in infos:
        if i.notes:
            out.append(f"  {i.folder}:")
            out += [f"    - {n}" for n in i.notes]
    return "\n".join(out)


def manifest_keys_seen(infos: list[ResultInfo]) -> list[str]:
    keys: set[str] = set()
    for i in infos:
        keys.update(str(k) for k in i.manifest)
    return sorted(keys)


def to_yaml_dict(infos: list[ResultInfo], best: ResultInfo | None,
                 project_root: Path) -> dict:
    def entry(i: ResultInfo) -> dict:
        return {
            "folder": i.folder,
            "path": str(i.path),
            "type": i.kind,
            "batch": i.batch,
            "chrs_covered": i.chrs_covered,
            "n_variants": i.n_variants,
            "n_carriers": i.n_carriers,
            "n_samples": i.n_samples,
            "n_carrier_rows": i.n_carrier_rows,
            "run_date": i.run_date,
            "reference_build": i.reference_build,
            "clinvar_release": i.clinvar_release,
            "vep_version": i.vep_version,
            "pipeline_sha": i.pipeline_sha,
            "complete": i.complete,
            "missing_files": i.missing_files,
            "notes": i.notes,
        }

    return {
        "_generated_by": "src/utils/discover_results.py",
        "_warning": (
            "Auto-generated; gitignored. n_carriers is the number of people with "
            ">=1 P/LP call, NOT the cohort size. n_samples is null unless an "
            "explicit sample list was found."
        ),
        "project_root": str(project_root),
        "selected": entry(best) if best else None,
        "candidates": [entry(i) for i in infos],
        "manifest_keys_seen": manifest_keys_seen(infos),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--project-root",
                    default="/sc/arion/projects/rg_huangk06/variants_PLP_MSM")
    ap.add_argument("--out", default="config/discovered_results.yaml")
    ap.add_argument("--json", action="store_true", help="emit JSON to stdout instead of a table")
    args = ap.parse_args(argv)

    root = Path(args.project_root)
    if not root.exists():
        print(f"ERROR: project root does not exist: {root}", file=sys.stderr)
        print("Are you running this on Minerva?", file=sys.stderr)
        return 2

    infos = discover(root)
    if not infos:
        print(f"ERROR: no results-*/ directories under {root}", file=sys.stderr)
        return 3

    best = select_best(infos)
    payload = to_yaml_dict(infos, best, root)

    if args.json:
        print(json.dumps(payload, indent=2, default=str))
    else:
        print(f"Scanned: {root}")
        print(f"Found {len(infos)} results-*/ dirs "
              f"({sum(1 for i in infos if i.kind == 'batch')} batch, "
              f"{sum(1 for i in infos if i.kind == 'pilot')} pilot)\n")
        print(render_table(infos, best))
        notes = render_notes(infos)
        if notes:
            print("\nNotes:")
            print(notes)
        print("\nManifest keys seen across runs:")
        print("  " + (", ".join(manifest_keys_seen(infos)) or "(none)"))
        print()
        if best:
            print(f"AUTO-SELECTED (*): {best.folder}   "
                  f"complete={'Y' if best.complete else 'N'}  run_date={best.run_date or '?'}")
            print(f"  carrier matrix: {best.path / 'carriers' / 'carrier_matrix.tsv'}")
        else:
            print("AUTO-SELECTED: none — no batch-type run found.")
        print("\nReminder: n_carriers = people with >=1 P/LP call, not cohort size.")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        import yaml
        out.write_text(yaml.safe_dump(payload, sort_keys=False, default_flow_style=False))
    except ImportError:
        out = out.with_suffix(".json")
        out.write_text(json.dumps(payload, indent=2, default=str))
    if not args.json:
        print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
