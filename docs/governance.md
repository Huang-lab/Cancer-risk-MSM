# Data governance

**IRB-restricted human genetic + clinical data.** These rules are absolute.

## Rules

1. **Everything protected stays on Minerva.** Never upload, push, email, or move
   off-cluster: VCFs, genotypes, phenotype/EHR files, carrier matrices,
   participant IDs, or any per-participant record.
2. **The git repository is code only.** Tracked: `src/`, `workflow/`, `config/`,
   `envs/`, `docs/`, `README`, `.gitignore`, `tests/synthetic/`.
3. **No participant data in code/configs/comments/README/tests.** If you need a
   fixture, generate synthetic data.
4. **Logs may contain PHI.** They live under `$OUTPUT_ROOT/logs/` on Minerva
   (gitignored). Never paste raw record contents anywhere outside the cluster.
5. **Remote pushes are code-only.** Run the audit below before any commit.

## `.gitignore` (see repo root)

Excludes `data/`, `results/`, `resources/`, `logs/`, `intermediate/`, all common
bio file extensions (`*.vcf*`, `*.bcf`, `*.bam`, `*.bed`, `*.fam`, ...), tabular
data extensions, and MSM ID patterns (`SINAI_*`, `SINAI-Million_*`).

## Pre-commit audit (run before every commit)

```bash
git ls-files | grep -Ei '\.(vcf|bcf|bgen|bam|cram|bim|fam|tsv|csv|txt)(\.gz)?$' \
  | grep -v '^tests/synthetic/'                # must be empty

git ls-files -z | xargs -0 grep -lE 'SINAI_[0-9]+_[A-Z]{2}[0-9]+|SINAI-Million_[0-9]+_[0-9]+' \
  || true                                      # must be empty

git status --short | grep -E '^\?\? .*/(data|results|resources|logs)/' \
  || true                                      # must be empty
```

If any check produces output, do not commit until it is resolved.

## What the code-review skill checks

- No `/sc/arion/.../data/**` literal record content in source or configs.
- No participant-ID regex matches in tracked files.
- Test fixtures are clearly synthetic and marked so.
