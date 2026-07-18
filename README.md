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
                       ┌──────────────────────────────────────────┐
Pre-flight  ─────────► │ 00_setup_dirs  00_confirm_build          │
                       │ 00_download_refs  (fresh ClinVar VCF)    │
                       └──────────────────┬───────────────────────┘
                                          │
Step 1: Annotation (chunk parallel)       │
   per chunk: bcftools norm → SITE+GT QC → VEP
              (VEP --custom attaches CLNSIG/CLNREVSTAT/CLNDN/CLNVI
               from the downloaded ClinVar VCF, plus AlphaMissense,
               LOFTEE, gnomAD)
                                          │
Step 2: Sample QC FROM VCFs               │
   bcftools stats per chunk → aggregate PSC
   plink2 --make-king (kinship) + --check-sex on WES
   → sample_keep_list.tsv
                                          │
Step 3: Classification (3 frameworks)     │
   ClinVar P/LP · ACMG (ANNOVAR+InterVar, drop ClinVar B/LB) · AlphaMissense
                                          │
Step 4: Carrier matrix (long format)      │
   one row per variant × person; is_clinvar_PLP · is_acmg_PLP · is_AM_PLP
                                          │
Step 5: Phenotype + covariates            │
   phecodeX roster + ICD→phecodeX first-dx dates · famhx · covariates
                                          │
Step 6: Analysis-ready join               │
                                          │
Step 7: ML per cancer                     │
   propensity match · temporal-leakage window · penalized logistic / Cox /
   RF / XGBoost · ancestry-stratified evaluation · cross-cohort validation
```

Steps **runnable today**: 1 (annotation) and 2 (VCF-based sample QC). Steps
3–7 are scaffolded and unlock after annotation completes.

---

## QC design

QC runs at **two levels**, both computed from our data (no reliance on the
biobank's common-variants outputs):

### Site + genotype QC — per chunk, between `bcftools norm` and VEP

Applied inside `workflow/01_norm_vep_chunk.sh` so annotation only sees clean
data. Config keys under `qc.site_gt`:

- `bcftools view -f PASS` — keep VQSR/GATK PASS sites only.
- `bcftools +setGT` nulls genotypes when any of:
  - `FMT/DP < genotype_dp_min` (default 10)
  - `FMT/GQ < genotype_gq_min` (default 20)
  - het call with allele balance outside `[genotype_ab_het_min, genotype_ab_het_max]` (default `[0.20, 0.80]`)
- `bcftools +fill-tags` recomputes AC, AN, F_MISSING.
- `bcftools view -e ...` drops sites with `F_MISSING > site_missing_max` (default 10%)
  or sites that became mono-allelic after masking.

Helper expressions live in `src/qc/site_gt_qc.py` and are unit-tested.

### Sample QC — computed *from the annotated WES VCFs* (post site/GT QC)

Runs after Step 1. Two-phase, both under `workflow/`:

**Phase A — per-chunk stats (LSF array):** `workflow/00_wes_qc_chunks.lsf`
runs `bcftools stats -s -` on each annotated chunk, writing
`results/qc/stats/<chunk>.stats.txt`. Resume-safe.

**Phase B — aggregate + kinship + sex-check + keep-list:**
`workflow/00_wes_qc_aggregate.sh` chains three LSF jobs:
1. Aggregate `PSC` lines across chunks → `results/qc/psc_agg.tsv`
   (per-sample n_called, n_missing, missingness, het/hom, singletons).
2. **Kinship**: `bcftools concat` autosomal chunks → `plink2 --maf 0.01
   --geno 0.05 --hwe 1e-10 --indep-pairwise 500kb 50 0.2` → `plink2 --make-king
   triangle bin` on the LD-pruned bfile → `results/qc/king.king` (+ `.kin0`).
   Bonus: `plink2 --pca 20` on the same set gives WES-derived ancestry PCs.
3. **Sex-check**: `bcftools view -r chrX` → concat → `plink2 --split-par b38
   --check-sex 0.20 0.80` → `results/qc/sexcheck.sexcheck` (X het F-stat).
4. Once (2) and (3) finish (`bsub -w done(...)`), `src.qc.sample_qc` reads
   all four tables and emits:
   - `results/qc/sample_keep_list.tsv` — rows to KEEP (feeds Step 4 carrier).
   - `results/qc/sample_qc_report.tsv` — per-sample flags + reasons.
   - `results/qc/sample_qc_summary.md` — aggregate counts (no per-sample values).

**Config keys** (`qc.sample`):

| Check                 | Threshold                                | Source (VCF-derived)                    |
|-----------------------|------------------------------------------|-----------------------------------------|
| Missingness           | `missingness_max` (default 0.05)         | `bcftools stats -s -` PSC               |
| Contamination (proxy) | het/hom `|z| > het_hom_z_max` (default 4)| `bcftools stats -s -` PSC               |
| Sex-check             | F<`f_stat_female_max`, F>`f_stat_male_min`; `mismatch_action` = flag \| exclude | `plink2 --check-sex` on chrX |
| Kinship duplicate     | `duplicate_kin_min` (default 0.354)      | `plink2 --make-king`                    |
| Kinship related       | `pi_hat_related_max` (default 0.1875) — flagged, kept | `plink2 --make-king`         |

**Trade-off worth naming:** WES-based kinship is inherently rougher than
array-common-variants kinship (fewer LD-pruned markers, more depth-dependent
noise). Duplicate/MZ-twin calls at kinship ≥ 0.354 remain robust; 2nd-degree
calls are less certain but only used to flag+keep with group-aware CV in the
ML stage, not to exclude. True contamination scoring (VerifyBamID2) needs
BAM/CRAM and is optional (`qc.sample.contamination_verifybamid.enabled`).

---

## ClinVar handling

VEP annotation uses the **freshly-downloaded weekly ClinVar VCF**. There is no
"look up in a static bundled table" — every run against a new download re-
annotates against the current release.

**Fetch:**
```bash
bash workflow/00_download_refs.sh
```
This wgets `clinvar.vcf.gz` from `ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/`,
verifies md5, parses `##fileDate=` from the header, and stages the file into
`$OUTPUT_ROOT/resources/clinvar/YYYY-MM-DD/`. It also brings down the
NCBI-published tabix index (or builds one if that is missing).

**Annotate:** VEP's `--custom` mechanism, invoked in `01_norm_vep_chunk.sh`,
reads the downloaded VCF at each variant site and attaches
`CLNSIG, CLNREVSTAT, CLNDN, CLNVI` to the CSQ record. That is the annotation —
`--custom` is exactly "for each site, look it up in this VCF and copy these
fields", so the output annotated VCFs carry the ClinVar release we downloaded.

**Config:** `config.resources.clinvar.release: latest` — the loader picks the
newest dated folder under `resources/clinvar/`, so a fresh download becomes
the default with no config edit.

**ACMG track (parallel):** for `src/classification/acmg.py`, ANNOVAR uses its
own `clinvar_latest` protocol against `resources.annovar.humandb`. Refresh that
periodically with ANNOVAR's `annotate_variation.pl -downdb clinvar_latest`.

---

## Step-by-step (on Minerva)

### Pre-flight

```bash
git clone <this repo> && cd Cancer-risk-MSM
mamba env create -f envs/annotation.yml && mamba activate msm-annotation
export CONFIG=config/config.yaml

bash workflow/00_setup_dirs.sh        # creates $OUTPUT_ROOT/{data,resources,results,logs}
bash workflow/00_confirm_build.sh     # asserts GRCh38 + resource paths + tool versions
bash workflow/00_download_refs.sh     # fetches the current ClinVar release
```

### Step 1 — Annotation (gated: single chunk, then array)

Per-chunk worker: `bcftools norm` → **site + GT QC** → `vep` with LOFTEE,
AlphaMissense, ClinVar via `--custom`, gnomAD AF via `--custom`, dbNSFP if present.

```bash
bash workflow/02_submit_single.lsf     # STEP-1 GATE: one chunk, review before scaling
bash workflow/03_submit_array.lsf      # capped LSF array (resume-safe)
```

Parallelism at two levels: `vep --fork N` threads within a chunk (LSF `-n N`
matches), and many chunks concurrently at `%config.lsf.array_concurrency_cap`
(default 50). Workers skip chunks whose validated output already exists.

### Step 2 — Sample QC (from VCFs)

Runs after Step 1 finishes.

```bash
bash workflow/00_run_sample_qc.sh      # submits per-chunk stats array
# after the array is DONE:
bash workflow/00_wes_qc_aggregate.sh   # chains kinship + sex-check + keep-list
```

Outputs: `sample_keep_list.tsv`, `sample_qc_report.tsv`, `sample_qc_summary.md`
under `$OUTPUT_ROOT/results/qc/`. Downstream carrier extraction filters to
this keep-list.

### Step 3 — Classification (three frameworks, scaffolded)

| Module                                | Framework                                                                                                     | Output              |
|---------------------------------------|---------------------------------------------------------------------------------------------------------------|---------------------|
| `src/classification/clinvar.py`       | ClinVar P/LP; `min_review_stars: 2`, no conflicts                                                              | `clinvar_plp.tsv`   |
| `src/classification/acmg.py`          | ANNOVAR + InterVar rule-based; **drop variants ClinVar calls B/LB**                                            | `acmg_plp.tsv`      |
| `src/classification/alphamissense.py` | Chen/Pejaver gene-specific evidence-label calibration; min-strength Moderate; domain-aggregate recorded, never promotes | `am_plp.tsv`        |
| `src/classification/qc.py`            | Per-gene QC: paralog artifacts (PMS2/PMS2CL...) + cohort-vs-gnomAD carrier-freq inflation                      | `qc_per_gene.tsv`   |

### Step 4 — Carrier matrix (long format)

`src/carriers/carrier_matrix.py` unions the three variant sets, runs
`bcftools query` per qualifying variant, and emits **one row per (variant × person)**:

```
chr  pos  ref  alt  gene  person_id  is_clinvar_PLP  is_acmg_PLP  is_AM_PLP
```

Filters to `results/qc/sample_keep_list.tsv`. Optional wide pivot for ML.

### Step 5 — Phenotype + covariates

- `id_harmonize.py` — canonical `person_id` across `SINAI_*` and `SINAI-Million_*` forms.
- `cases.py` — phecodeX roster + **ICD → phecodeX** from ICD-coded
  `Encounter_Diagnosis.txt` and `Problem_List.txt` for first-dx dates.
- `famhx.py` — first-degree cancer family history.
- `covariates.py` — demographics, smoking, BMI, parity, screening; collapsed to
  the temporal-leakage window in `config.ml.temporal_leakage.feature_window_days`.

### Step 6 — Analysis-ready join

`src/carriers/join.py` joins the long carrier table with harmonized IDs, PCs
(from Step 2's WES-derived PCA or the biobank's common-variants PCs), and
phenotype rows → `results/analysis/analysis_ready.tsv`.

### Step 7 — ML per cancer

Add a cancer by adding `config/<cancer>.yaml` (provided: `crc.yaml`, `breast.yaml`).
Pipeline: dataset → match (propensity k:1 caliper) → features (window `[-730d, -182d]`)
→ models (penalized logistic / Cox / RF / XGBoost) → evaluate (AUC, PR-AUC,
ancestry-stratified calibration, SHAP) → external_validate (cross-cohort with BioMe).

---

## Repo layout

```
config/            config.yaml (single source of truth) + crc.yaml + breast.yaml
src/qc/            sample_qc, site_gt_qc, vcf_stats (bcftools + plink2 helpers)
src/annotation/    chunk enumeration + per-chunk output validation
src/classification/  clinvar.py, acmg.py, alphamissense.py, qc.py
src/phenotype/     id_harmonize, cases (ICD→phecodeX), famhx, covariates
src/carriers/      long-format carrier table + optional wide pivot
src/ml/            dataset, match, features, models, evaluate, external_validate
src/common/        config loader, PHI-safe logger
workflow/          numbered scripts:
                     00_setup_dirs, 00_confirm_build, 00_download_refs,
                     01_norm_vep_chunk, 02_submit_single, 03_submit_array,
                     04_validate_chunk,
                     00_run_sample_qc → 00_wes_qc_chunks → 00_wes_qc_aggregate
envs/              conda specs: annotation.yml, analysis.yml
docs/              annotation_runbook.md, data_dictionary.md, governance.md
tests/             synthetic fixtures + unit tests for classification + QC helpers
```

## Runtime tree (Minerva-only, gitignored)

```
$OUTPUT_ROOT = /sc/arion/projects/rg_huangk06/variants_PLP_MSM
  data/annotated/           per-chunk annotated VCFs + tabix indexes
  resources/clinvar/YYYY-MM-DD/  weekly ClinVar releases (downloaded on demand)
  resources/                other references (VEP cache, LOFTEE, AlphaMissense, ...)
  results/qc/               sample_keep_list + reports + stats/ + king + sexcheck
  results/variants/         classified variant tables per framework
  results/carriers/         long-format carriers + optional wide pivot
  results/analysis/         analysis-ready joined tables
  results/models/           ML outputs
  logs/                     per-stage logs (may contain PHI; never leaves cluster)
```

See [`docs/annotation_runbook.md`](docs/annotation_runbook.md) for the runbook,
[`docs/data_dictionary.md`](docs/data_dictionary.md) for output schemas, and
[`docs/governance.md`](docs/governance.md) for the pre-commit audit steps.
