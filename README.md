# Cancer-risk-MSM

Pipeline for the Mount Sinai Million (MSM) biobank:
1. Annotate WES pVCFs (bcftools norm + VEP with LOFTEE, AlphaMissense, ClinVar, gnomAD).
2. Classify P/LP variants under three frameworks: ClinVar, ACMG (ANNOVAR+InterVar,
   with ClinVar-B/LB removal), AlphaMissense gene-specific calibrated.
3. Phenotype cancer cases + family history + covariates.
4. Build per-person carrier tables (long format: one row per variant × person).
5. Predict cancer risk (colorectal + breast first; add cancers by adding a YAML).

Mirrors [`Huang-lab/Biome-alphamissense-calibration`](https://github.com/Huang-lab/Biome-alphamissense-calibration)
and [`Huang-lab/BioMe-cancer-risk-prediction`](https://github.com/Huang-lab/BioMe-cancer-risk-prediction)
for methodology and config style.

## Data governance (read first)

This repository contains **code only**. IRB-restricted human genetic and clinical
data lives on Minerva and never enters git. See `docs/governance.md`.

Runtime tree (gitignored, Minerva-only), created by `workflow/00_setup_dirs.sh`:

```
$OUTPUT_ROOT = /sc/arion/projects/rg_huangk06/variants_PLP_MSM
  data/annotated/          per-chunk annotated VCFs + tabix indexes
  resources/               local mirrors/symlinks to references
  results/variants/        classified variant tables per framework
  results/carriers/        per-person carrier tables (long + optional wide)
  results/analysis/        analysis-ready joined tables
  results/models/          ML outputs
  logs/                    per-stage logs (may contain PHI; never leaves cluster)
```

## Quick start (on Minerva)

```bash
git clone <this repo> && cd Cancer-risk-MSM
mamba env create -f envs/annotation.yml && mamba activate msm-annotation

# 1. Set up the runtime tree
bash workflow/00_setup_dirs.sh

# 2. Confirm GRCh38 + resource paths + tool versions
bash workflow/00_confirm_build.sh

# 3. STEP-1 GATE: run one chunk, review before scaling
bash workflow/02_submit_single.lsf

# 4. STEP-2: capped LSF array over all chunks (resume-safe)
bash workflow/03_submit_array.lsf
```

See `docs/annotation_runbook.md` for details, and `docs/governance.md` for
governance/audit steps.

## Layout

```
config/            Single-source-of-truth YAML (paths, thresholds, gene panels)
src/annotation/    Chunk enumeration, per-chunk validation
src/classification/  clinvar.py, acmg.py, alphamissense.py, qc.py
src/phenotype/     id_harmonize, cases (ICD -> phecodeX), famhx, covariates
src/carriers/      Long-format carrier table + optional wide pivot
src/ml/            Dataset/match/features/models/evaluate/external_validate
src/common/        Config loader, PHI-safe logger
workflow/          LSF submit scripts (00-04)
envs/              Conda specs for annotation + analysis
docs/              Runbook, data dictionary, governance
tests/synthetic/   Synthetic-only unit fixtures (no real data)
```
