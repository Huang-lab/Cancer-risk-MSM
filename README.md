# Cancer-risk-MSM

End-to-end pipeline for the Mount Sinai Million (MSM) biobank: annotate WES
pVCFs, classify P/LP variants under three frameworks, build per-person carrier
tables, phenotype cancer cases + family history + covariates, and model cancer
risk (colorectal + breast first).

Mirrors [`Huang-lab/Biome-alphamissense-calibration`](https://github.com/Huang-lab/Biome-alphamissense-calibration)
and [`Huang-lab/BioMe-cancer-risk-prediction`](https://github.com/Huang-lab/BioMe-cancer-risk-prediction)
for methodology and config style.

> **Governance first.** The repo contains **code only**. IRB-restricted data
> lives on Minerva under `$OUTPUT_ROOT` and is gitignored. See
> [`docs/governance.md`](docs/governance.md).

---

## Workflow overview

```
                       ┌─────────────────────────────────┐
Pre-flight  ─────────► │ 00_setup_dirs  00_confirm_build │
                       └───────────────┬─────────────────┘
                                       │
Step 0: Sample QC  ──────►  00_run_sample_qc.sh  ──►  sample_keep_list.tsv
                                       │
Step 1: Annotation ──────► [ 02_submit_single.lsf → REVIEW → 03_submit_array.lsf ]
        per chunk:        bcftools norm → SITE+GT QC → VEP (all annotations)
                                       │
Step 2: Classification (3 in parallel)  ──►  clinvar_plp / acmg_plp / am_plp
                                       │
Step 3: Carrier matrix (long format)   ──►  carriers.long.tsv
                                       │
Step 4: Phenotype + covariates         ──►  cases / famhx / covariates
                                       │
Step 5: Analysis-ready join            ──►  analysis_ready.tsv
                                       │
Step 6: ML per cancer                  ──►  model / metrics / SHAP / calibration
```

The **only phase runnable today** is Steps 0 + 1 (sample QC → annotation) via
the LSF scripts under `workflow/`. Steps 2–6 are scaffolded (structure and
config-plumbing in place, per-cluster implementation pending). This is
intentional: run and review annotation first, then unlock downstream stages.

---

## Step-by-step (on Minerva)

### Pre-flight — one-time setup

```bash
git clone <this repo> && cd Cancer-risk-MSM
mamba env create -f envs/annotation.yml && mamba activate msm-annotation
export CONFIG=config/config.yaml

bash workflow/00_setup_dirs.sh        # creates $OUTPUT_ROOT/{data,resources,results,logs}
bash workflow/00_confirm_build.sh     # asserts GRCh38 + resource paths + tool versions
```

`00_confirm_build.sh` writes a version log to `$OUTPUT_ROOT/logs/annotation/`
and fails fast if any resource path in `config/config.yaml` is missing.

### Step 0 — Sample QC (produces the keep-list)

Runs **once, before annotation**. Uses the biobank's already-computed common-
variants outputs (kinship, PLINK sex-check, PCA) so we don't re-genotype.

```bash
bash workflow/00_run_sample_qc.sh
```

**What it checks** (thresholds in `config.qc.sample`):
- **Missingness** — drop samples with variant missingness above `missingness_max` (default 5%).
- **Sex check** — reported sex vs PLINK `--check-sex` inferred sex; `mismatch_action` = `flag` or `exclude`.
- **Contamination** — VerifyBamID2 `FREEMIX > freemix_max` (default 3%) → exclude.
  Set `qc.sample.contamination.enabled: true` once the summary is on cluster.
- **Kinship** — from `common_variants/.../kinship/kinship.tsv`. For pairs with
  KING kinship ≥ `duplicate_kin_min` (0.354, MZ-twin / duplicate), keep the sample
  with the highest call rate. Related pairs (≥ 2nd degree) are **kept but flagged**;
  ML uses group-aware cross-validation.

**Outputs** (all gitignored on Minerva):
- `results/qc/sample_keep_list.tsv`   — one row per KEPT `sample_id`.
- `results/qc/sample_qc_report.tsv`   — per-sample kept/dropped + reasons.
- `results/qc/sample_qc_summary.md`   — aggregate counts (no per-sample values).

### Step 1 — Annotation (gated: single chunk, then array)

Per-chunk worker (`workflow/01_norm_vep_chunk.sh`) does, in order:

1. **`bcftools norm`** — split multi-allelics, left-align to `$FASTA`.
2. **Site + genotype QC** (`config.qc.site_gt`, expressions in `src/qc/site_gt_qc.py`):
   - `bcftools view -f PASS` (keep VQSR/GATK PASS sites).
   - `bcftools +setGT` null low-quality genotypes:
     `FMT/DP<10 | FMT/GQ<20 | het-AB outside [0.20, 0.80]`.
   - `bcftools +fill-tags` recomputes AC/AN/F_MISSING.
   - `bcftools view -e ...` drops sites with `F_MISSING > 10%` or that
     became mono-allelic after masking.
3. **VEP** — one pass, attaches **all** raw fields the downstream classifiers need:
   MANE Select / canonical / HGVS, **LOFTEE** (`hc` pLoF), **AlphaMissense**
   scores, **ClinVar** via `--custom` (CLNSIG/CLNREVSTAT/CLNDN/CLNVI),
   **gnomAD** AF/popmax via `--custom`, **dbNSFP** if present. Output bgzip + tabix
   into `$OUTPUT_ROOT/data/annotated/`.

Parallelism runs at two levels:
- *Within a chunk:* `vep --fork N` runs N threads (LSF `-n N` matches).
- *Across chunks:* the LSF array runs many workers concurrently, capped at
  `config.lsf.array_concurrency_cap` (default `%50`).

**Submit — gated:**

```bash
# STEP 1a: single-chunk GATE — run one chunk, review, do NOT continue on failure.
bash workflow/02_submit_single.lsf
# → review $OUTPUT_ROOT/data/annotated/ + logs; validate_chunk must print OK.

# STEP 1b: only after 1a passes review, capped LSF array over all chunks.
bash workflow/03_submit_array.lsf     # bsub -J "vep[1-N]%50"
```

The workers are **resume-safe**: rerunning skips any chunk whose output already
exists and passes `src.annotation.validate_chunk`.

### Step 2 — Classification (three frameworks, scaffolded)

Per-chunk classifiers read the annotated VCFs and emit long variant tables
under `results/variants/`:

| Module                              | Framework        | Output                 |
|-------------------------------------|------------------|------------------------|
| `src/classification/clinvar.py`     | ClinVar P/LP     | `clinvar_plp.tsv`      |
| `src/classification/acmg.py`        | ACMG via ANNOVAR + InterVar (auto), then **drop variants ClinVar calls B/LB** | `acmg_plp.tsv` |
| `src/classification/alphamissense.py` | AlphaMissense gene-specific evidence-label calibration (Chen/Pejaver; min-strength `Moderate`; domain-aggregate recorded, never promotes) | `am_plp.tsv` |
| `src/classification/qc.py`          | Per-gene QC: paralog artifacts (PMS2/PMS2CL, CHEK2, NBN, ...) + cohort-vs-gnomAD carrier-freq inflation | `qc_per_gene.tsv` |

### Step 3 — Carrier matrix (long format; scaffolded)

`src/carriers/carrier_matrix.py` unions the three variant sets, runs
`bcftools query` per qualifying variant, and emits **one row per (variant × person)**:

```
chr  pos  ref  alt  gene  person_id  is_clinvar_PLP  is_acmg_PLP  is_AM_PLP
```

The genotype-parse filters to `results/qc/sample_keep_list.tsv` from Step 0.
Optional wide pivot (`config.carriers.emit_wide_pivot: true`) produces the
sample × gene matrix for ML.

### Step 4 — Phenotype + covariates (scaffolded)

- `src/phenotype/id_harmonize.py` — canonical `person_id` across
  `SINAI_*` (pre-2022) and `SINAI-Million_*` (post-2022) forms.
- `src/phenotype/cases.py` — `PheWas_MSM_phecodeX.tsv` roster, plus
  **ICD → phecodeX** mapping from the ICD-coded `Encounter_Diagnosis.txt`
  and `Problem_List.txt` for first-diagnosis dates; incident vs prevalent
  distinction using Medical/Surgical History.
- `src/phenotype/famhx.py` — first-degree cancer family history from `Family_History.txt`.
- `src/phenotype/covariates.py` — demographics, smoking/alcohol, BMI,
  parity/age-at-first-birth, screening flags — collapsed to the
  temporal-leakage window in `config.ml.temporal_leakage.feature_window_days`.

### Step 5 — Analysis-ready join (scaffolded)

`src/carriers/join.py` joins the long carrier table with harmonized IDs +
ancestry PCs (from `common_variants/.../pca/pcs.tsv`) + phenotype rows,
writing `results/analysis/analysis_ready.tsv`.

### Step 6 — ML per cancer (scaffolded)

Adding a cancer = adding a `config/<cancer>.yaml`. Provided today: `crc.yaml`,
`breast.yaml`. Pipeline modules under `src/ml/` mirror the BioMe repo:

```
dataset → match (propensity k:1 caliper) → features (window [-730d, -182d]) →
models (penalized logistic / Cox / RF / XGBoost) → evaluate (AUC, PR-AUC,
ancestry-stratified calibration, SHAP) → external_validate (cross-cohort with
BioMe)
```

---

## Repo layout

```
config/            config.yaml (single source of truth) + crc.yaml + breast.yaml
src/qc/            sample_qc, site_gt_qc (helpers imported by workflow scripts)
src/annotation/    chunk enumeration + per-chunk output validation
src/classification/  clinvar.py, acmg.py, alphamissense.py, qc.py
src/phenotype/     id_harmonize, cases (ICD→phecodeX), famhx, covariates
src/carriers/      long-format carrier table + optional wide pivot
src/ml/            dataset, match, features, models, evaluate, external_validate
src/common/        config loader, PHI-safe logger
workflow/          numbered scripts (00_setup, 00_confirm_build, 00_run_sample_qc,
                   01_norm_vep_chunk, 02_submit_single, 03_submit_array, 04_validate)
envs/              conda specs: annotation.yml, analysis.yml
docs/              annotation_runbook.md, data_dictionary.md, governance.md
tests/             synthetic fixtures + unit tests for rules and QC helpers
```

## Runtime tree (Minerva-only, gitignored)

```
$OUTPUT_ROOT = /sc/arion/projects/rg_huangk06/variants_PLP_MSM
  data/annotated/           per-chunk annotated VCFs + tabix indexes
  resources/                symlinks/mirrors to references
  results/qc/               sample_keep_list.tsv + reports
  results/variants/         classified variant tables per framework
  results/carriers/         long-format carriers + optional wide pivot
  results/analysis/         analysis-ready joined tables
  results/models/           ML outputs
  logs/                     per-stage logs (may contain PHI; never leaves cluster)
```

## What the config controls (single source of truth)

- Input auto-resolution: latest WES batch, latest phenotypes dated folder, latest ClinVar release.
- Reference paths: FASTA, VEP cache/plugins, LOFTEE, AlphaMissense scores + calibration table, ClinVar, gnomAD, dbNSFP, ANNOVAR humandb, InterVar.
- **QC thresholds** (`qc.site_gt` and `qc.sample`) — everything QC touches lives here.
- Classification thresholds: ClinVar stars, ACMG "remove if ClinVar B/LB" flag, AlphaMissense min-evidence.
- LSF: queue, project, walltime, memory, `--fork`, array cap, modules to load.
- Phenotype: ICD→phecodeX table, cancer phecodes, case-definition rules.
- ML: temporal-leakage window, matching ratio/caliper, covariates, models, CV.

See [`docs/annotation_runbook.md`](docs/annotation_runbook.md) for the runbook,
[`docs/data_dictionary.md`](docs/data_dictionary.md) for output schemas, and
[`docs/governance.md`](docs/governance.md) for the pre-commit audit steps.
