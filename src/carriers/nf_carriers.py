"""Read carrier calls produced by `Huang-lab/germline-plp-carrier-nf`.

We do NOT extract carriers ourselves — the NF pipeline owns annotation,
P/LP classification, and carrier extraction. This module is the read-boundary:
it pools one or more NF run-output dirs into a single long-format table and
aggregates to per-person / per-gene features for ML.

NF long-format schema (config: inputs.carrier_source.columns):
    chr  pos  ref  alt  gene  person_id  is_clinvar_PLP  is_acmg_PLP  is_AM_PLP
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass
class CarrierRow:
    person_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    is_clinvar_plp: bool
    is_acmg_plp: bool
    is_am_plp: bool

    def framework_flag(self, framework: str) -> bool:
        return {
            "clinvar": self.is_clinvar_plp,
            "acmg": self.is_acmg_plp,
            "am": self.is_am_plp,
        }[framework]


def _truthy(v: str | None) -> bool:
    return str(v or "").strip().lower() in ("1", "true", "t", "yes", "y")


def nf_result_dirs(cfg: dict) -> list[Path]:
    """Absolute paths of the NF run-output dirs we pool."""
    src = cfg["inputs"]["carrier_source"]
    root = Path(src["nf_results_root"])
    return [root / d for d in src["results_dirs"]]


def read_manifests(cfg: dict) -> list[dict]:
    """Read each NF run's manifest.json (ref build, tool versions, pipeline SHA).

    Used to assert all pooled runs share a reference build before we merge them.
    """
    src = cfg["inputs"]["carrier_source"]
    out: list[dict] = []
    for d in nf_result_dirs(cfg):
        p = d / src.get("manifest", "manifest.json")
        if p.exists():
            try:
                m = json.loads(p.read_text())
                m["_source_dir"] = d.name
                out.append(m)
            except json.JSONDecodeError:
                continue
    return out


def assert_consistent_builds(manifests: list[dict]) -> str | None:
    """Return the shared reference build, or raise if runs disagree.

    Pooling carrier calls across different reference builds would silently
    mis-join coordinates, so this is a hard check rather than a warning.
    """
    builds = {m.get("reference_build") or m.get("genome") for m in manifests}
    builds.discard(None)
    if len(builds) > 1:
        raise ValueError(
            f"NF runs disagree on reference build: {sorted(builds)}. "
            "Refusing to pool carrier calls across builds."
        )
    return next(iter(builds)) if builds else None


def read_carriers(cfg: dict) -> list[CarrierRow]:
    """Pool carrier_matrix.tsv across all configured NF result dirs.

    Duplicate (person, variant) rows across runs are collapsed, OR-ing the
    per-framework flags — a variant called P/LP in any pooled run counts.
    """
    src = cfg["inputs"]["carrier_source"]
    cols = src["columns"]
    sep = src.get("sep", "\t")

    merged: dict[tuple[str, str, int, str, str], CarrierRow] = {}
    for d in nf_result_dirs(cfg):
        path = d / src["carrier_matrix"]
        if not path.exists():
            continue
        with path.open(newline="") as fh:
            rdr = csv.DictReader(fh, delimiter=sep)
            for row in rdr:
                pid = (row.get(cols["person_id"]) or "").strip()
                gene = (row.get(cols["gene"]) or "").strip()
                chrom = (row.get(cols["chrom"]) or "").strip()
                if not pid or not chrom:
                    continue
                try:
                    pos = int(str(row.get(cols["pos"], "")).strip())
                except ValueError:
                    continue
                ref = (row.get(cols["ref"]) or "").strip()
                alt = (row.get(cols["alt"]) or "").strip()
                key = (pid, chrom, pos, ref, alt)
                new = CarrierRow(
                    person_id=pid, chrom=chrom, pos=pos, ref=ref, alt=alt, gene=gene,
                    is_clinvar_plp=_truthy(row.get(cols["is_clinvar_plp"])),
                    is_acmg_plp=_truthy(row.get(cols["is_acmg_plp"])),
                    is_am_plp=_truthy(row.get(cols["is_am_plp"])),
                )
                prev = merged.get(key)
                if prev is None:
                    merged[key] = new
                else:
                    # OR the framework flags across pooled runs.
                    prev.is_clinvar_plp |= new.is_clinvar_plp
                    prev.is_acmg_plp |= new.is_acmg_plp
                    prev.is_am_plp |= new.is_am_plp
    return list(merged.values())


def read_sample_keep_list(cfg: dict) -> set[str]:
    """Read the NF pipeline's sample-QC keep-list (sample-level QC is NF's job).

    Returns an empty set when the path isn't configured, which
    `filter_to_keep_list` treats as a no-op rather than excluding everyone.
    """
    src = cfg["inputs"]["carrier_source"]
    rel = src.get("sample_keep_list")
    if not rel:
        return set()
    col = src.get("sample_keep_list_col", "sample_id")
    keep: set[str] = set()
    for d in nf_result_dirs(cfg):
        p = d / rel
        if not p.exists():
            continue
        with p.open(newline="") as fh:
            first = fh.readline()
            has_header = col in first
            if not has_header and first.strip():
                keep.add(first.strip().split("\t")[0])
            rdr = csv.DictReader(fh, delimiter="\t", fieldnames=None) if has_header else None
            if has_header:
                header = [h.strip() for h in first.rstrip("\n").split("\t")]
                idx = header.index(col) if col in header else 0
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) > idx and f[idx].strip():
                        keep.add(f[idx].strip())
            else:
                for line in fh:
                    sid = line.strip().split("\t")[0]
                    if sid:
                        keep.add(sid)
    return keep


def filter_to_keep_list(rows: list[CarrierRow], keep: set[str]) -> list[CarrierRow]:
    """Restrict to QC-passing samples. Empty keep-set is a no-op, not a wipe."""
    if not keep:
        return rows
    return [r for r in rows if r.person_id in keep]


def per_gene_flags(rows: list[CarrierRow], framework: str
                   ) -> dict[str, set[str]]:
    """{person_id -> set(genes)} where the person carries a P/LP call for `framework`."""
    out: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        if r.framework_flag(framework) and r.gene:
            out[r.person_id].add(r.gene)
    return dict(out)


def panel_carrier_flags(rows: list[CarrierRow], framework: str,
                        panel_genes: list[str]) -> dict[str, bool]:
    """{person_id -> bool} aggregate: carries any P/LP in the focus panel."""
    panel = set(panel_genes)
    by_person = per_gene_flags(rows, framework)
    return {pid: bool(genes & panel) for pid, genes in by_person.items()}
