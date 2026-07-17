# Annotation runbook

Everything below runs **on Minerva**. This session (in a cloud container) built
the scripts; execution happens on-cluster by the user.

## 0. One-time setup

```bash
cd <repo>
mamba env create -f envs/annotation.yml
mamba activate msm-annotation
export CONFIG=config/config.yaml
export OUTPUT_ROOT=$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")
bash workflow/00_setup_dirs.sh
bash workflow/00_confirm_build.sh          # asserts GRCh38 + resource presence
```

## 1. Step-1 gate: single chunk

Runs ONE chunk end-to-end, then STOPS for review.

```bash
bash workflow/02_submit_single.lsf         # picks first Target chunk
bjobs -J vep_gate
```

**Review** (all must pass; ping the ML lead before scaling):

- Exit code 0 (see `$OUTPUT_ROOT/logs/annotation/lsf/gate.*.out`).
- `python -m src.annotation.validate_chunk --input <chunk> --annotated <annot>` prints `OK ...`.
- Header carries INFO IDs for `CLNSIG`, `am_pathogenicity`, `LoF`, `gnomAD_AF`.
- Variant count sane vs input (default: ≥ 50%).
- Spot-check: a few records on `BRCA1`, `MLH1`, `TP53` look plausible.

## 2. Step-2: full array

Only after Step-1 passes review:

```bash
bash workflow/03_submit_array.lsf
```

- Concurrency capped at `config.lsf.array_concurrency_cap` (default `%50`).
- Each worker: `bcftools norm -m -any -f $FASTA` → `vep --fork N` with LOFTEE,
  AlphaMissense, ClinVar (`--custom`), gnomAD AF (`--custom`), dbNSFP if present
  → bgzip + tabix into `$OUTPUT_ROOT/data/annotated/`.
- **Parallelism at two levels:**
  1. *Within a chunk:* `vep --fork N` runs N worker threads inside the single
     worker. LSF `-n N` matches so the slot reserves N cores.
  2. *Across chunks:* the LSF array runs many single-chunk workers concurrently,
     capped at `%50`.
- **Resume-safe:** re-submitting skips any chunk whose annotated output already
  exists and passes `validate_chunk`.

Progress:

```bash
bjobs -A                       # array status
grep -c '^DONE' $OUTPUT_ROOT/logs/annotation/*.log
ls $OUTPUT_ROOT/data/annotated/*.annot.vcf.gz | wc -l
```

## Resources referenced (from config.resources)

- FASTA (GRCh38 primary assembly)
- VEP cache (assembly=GRCh38, species=homo_sapiens, cache_version matching install)
- LOFTEE plugin dir + ancestor FA + conservation SQLite
- AlphaMissense plugin dir + scores TSV.gz + gene-specific calibration table
  (mirrors Huang-lab/Biome-alphamissense-calibration)
- ClinVar VCF (latest dated release; `resources.clinvar.release: latest`)
- gnomAD v4 exomes VCF (AF, AF_popmax)
- dbNSFP (optional)
- ANNOVAR humandb + InterVar (for the ACMG classification track, downstream)
