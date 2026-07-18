"""Command builders + parsers for VCF-derived QC stats.

We compute site- and sample-level QC directly from the annotated WES VCFs
(after site/GT masking), not from the biobank's common-variants outputs.

Stats sources
-------------
- `bcftools stats -s -` per chunk
    PSC lines carry per-sample counts (nRefHom, nNonRefHom, nHets, nSingletons,
    nMissing, etc.). Aggregated across chunks -> per-sample call rate,
    het/hom ratio (weak contamination proxy), singleton count.
- `plink2 --make-king` on autosomal MAF>0.01, LD-pruned variants extracted
  from the WES data -> pairwise kinship.
- `plink2 --check-sex` on chrX non-PAR -> F-stat for reported-vs-inferred sex.

WES-derived kinship is inherently rougher than array-common-variants kinship
(smaller marker count, more depth-dependent noise), but it is self-contained
and does not depend on external artifacts. Duplicate/MZ-twin detection at
kinship >= 0.354 remains robust; 2nd-degree calls are less certain.
"""
from __future__ import annotations

from dataclasses import dataclass


# --- bcftools stats ---------------------------------------------------------

def bcftools_stats_cmd(vcf: str, out: str, fasta: str | None = None) -> list[str]:
    """`bcftools stats -s - -F fasta <vcf> > out.stats.txt`"""
    cmd = ["bcftools", "stats", "-s", "-"]
    if fasta:
        cmd += ["-F", fasta]
    cmd += [vcf]
    return cmd  # caller redirects stdout to `out`


@dataclass
class PSC:
    """One row from `bcftools stats -s -` PSC section."""
    sample: str
    n_ref_hom: int
    n_non_ref_hom: int
    n_hets: int
    n_singletons: int
    n_missing: int

    @property
    def n_called(self) -> int:
        return self.n_ref_hom + self.n_non_ref_hom + self.n_hets

    @property
    def n_total(self) -> int:
        return self.n_called + self.n_missing

    @property
    def missingness(self) -> float:
        t = self.n_total
        return (self.n_missing / t) if t else 1.0

    @property
    def het_hom_ratio(self) -> float:
        return (self.n_hets / self.n_non_ref_hom) if self.n_non_ref_hom else float("nan")


def parse_psc_lines(stats_text: str) -> list[PSC]:
    """Parse PSC section of bcftools stats output.

    Format:
      # PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs.
      # PSC   [2]id  [3]sample  [4]nRefHom  [5]nNonRefHom  [6]nHets  [7]nTransitions
      #       [8]nTransversions  [9]nIndels  [10]average depth  [11]nSingletons
      #       [12]nHapRef  [13]nHapAlt  [14]nMissing
      PSC     0      sample_x   1000       200            300      ...
    """
    out: list[PSC] = []
    for line in stats_text.splitlines():
        if not line.startswith("PSC"):
            continue
        f = line.split("\t")
        if len(f) < 14:
            continue
        try:
            out.append(PSC(
                sample=f[2],
                n_ref_hom=int(f[3]),
                n_non_ref_hom=int(f[4]),
                n_hets=int(f[5]),
                n_singletons=int(f[10]) if f[10] not in (".", "") else 0,
                n_missing=int(f[13]) if f[13] not in (".", "") else 0,
            ))
        except (ValueError, IndexError):
            continue
    return out


def merge_psc(rows: list[PSC]) -> dict[str, PSC]:
    """Aggregate across chunks: sum counts per sample."""
    agg: dict[str, PSC] = {}
    for r in rows:
        if r.sample not in agg:
            agg[r.sample] = PSC(r.sample, 0, 0, 0, 0, 0)
        a = agg[r.sample]
        agg[r.sample] = PSC(
            r.sample,
            a.n_ref_hom + r.n_ref_hom,
            a.n_non_ref_hom + r.n_non_ref_hom,
            a.n_hets + r.n_hets,
            a.n_singletons + r.n_singletons,
            a.n_missing + r.n_missing,
        )
    return agg


# --- plink2 commands --------------------------------------------------------

def plink2_extract_common_ld_pruned_cmd(
    vcf_glob: str, out_prefix: str, maf_min: float, ld_window_kb: int,
    ld_step: int, ld_r2: float, threads: int,
) -> list[str]:
    """Build a common-variant, LD-pruned bfile from many annotated chunks."""
    return [
        "plink2",
        "--vcf", vcf_glob, "dosage=DS",
        "--maf", str(maf_min),
        "--geno", "0.05",
        "--hwe", "1e-10",
        "--autosome",
        "--indep-pairwise", f"{ld_window_kb}kb", str(ld_step), str(ld_r2),
        "--threads", str(threads),
        "--out", out_prefix,
    ]


def plink2_make_king_cmd(bfile_prefix: str, out_prefix: str, threads: int) -> list[str]:
    """King-robust kinship on the LD-pruned bfile."""
    return [
        "plink2",
        "--bfile", bfile_prefix,
        "--make-king", "triangle", "bin",
        "--threads", str(threads),
        "--out", out_prefix,
    ]


def plink2_check_sex_cmd(bfile_prefix: str, out_prefix: str,
                        female_max: float, male_min: float) -> list[str]:
    """chrX heterozygosity F-stat for sex-check."""
    return [
        "plink2",
        "--bfile", bfile_prefix,
        "--check-sex", str(female_max), str(male_min),
        "--out", out_prefix,
    ]


# --- kinship parsing --------------------------------------------------------

@dataclass
class KinshipPair:
    id_a: str
    id_b: str
    kinship: float


def parse_king_kin0(king_kin0_path: str) -> list[KinshipPair]:
    """Parse `plink2 --make-king` `.king.kin0` output (tab-separated)."""
    out: list[KinshipPair] = []
    with open(king_kin0_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_a, i_b, i_k = header.index("ID1"), header.index("ID2"), header.index("KINSHIP")
        except ValueError:
            return out
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                out.append(KinshipPair(f[i_a], f[i_b], float(f[i_k])))
            except (ValueError, IndexError):
                continue
    return out
