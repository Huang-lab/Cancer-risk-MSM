# Cancer-risk-MSM

Cancer-risk modeling for the Mount Sinai Million (MSM) biobank: sample QC,
EHR phenotyping, and germline-carrier → cancer-risk prediction.

> **Scope.** This repo is the **downstream** layer. Variant annotation, P/LP
> classification, and carrier extraction are owned by
> [`Huang-lab/germline-plp-carrier-nf`](https://github.com/Huang-lab/germline-plp-carrier-nf).
> We consume its outputs read-only and never re-annotate.

Methodology and config style mirror
[`Huang-lab/BioMe-cancer-risk-prediction`](https://github.com/Huang-lab/BioMe-cancer-risk-prediction)
and [`Huang-lab/Biome-alphamissense-calibration`](https://github.com/Huang-lab/Biome-alphamissense-calibration).

> **Governance.** Code only. IRB-restricted data stays on Minerva under
> `$OUTPUT_ROOT`, gitignored. See [`docs/governance.md`](docs/governance.md).

---

## Division of labor

```
┌─────────────────────────────────────────────────────────────┐
│ germline-plp-carrier-nf   (Nextflow — NOT this repo)        │
│   pVCF → NORM_QC → VEP (ClinVar/AlphaMissense/LOFTEE/gnomAD)│
│        → ANNOVAR+InterVar → 3 P/LP classifiers → carriers    │
│   writes: results-*/carriers/carrier_matrix.tsv              │
│           results-*/variants/{clinvar,acmg,am}_plp.tsv       │
│           results-*/manifest.json                            │
└──────────────────────────┬──────────────────────────────────┘
                           │  read-only
                           ▼
┌─────────────────────────────────────────────────────────────┐
│ Cancer-risk-MSM  (THIS repo)                                │
│   Step 1  Sample QC      from raw pVCFs → sample_keep_list   │
│   Step 2  Phenotyping    ICD → cases / famhx / controls      │
│   Step 3  Analysis join  carriers ⋈ QC ⋈ PCs ⋈ phenotypes    │
│   Step 4  Cancer-risk ML matching, leakage control, models   │
└─────────────────────────────────────────────────────────────┘
```

**Runnable today:** Steps 1 and 2. Steps 3–4 are scaffolded.

## Paths on Minerva

```
/sc/arion/projects/rg_huangk06/
  variants_PLP_MSM/                  # NF pipeline workspace — we only READ
    refs/                            #   staged references (VEP cache, ClinVar, ...)
    results-batch1/                  #   carriers/carrier_matrix.tsv
    results-batch2/                  #   variants/*_plp.tsv, manifest.json
    results-batch3/
    cancer_risk/                     # ── $OUTPUT_ROOT for THIS repo ──
      resources/                     #   wes_samples.txt
      results/qc/                    #   sample_keep_list.tsv, psc_agg, king, sexcheck, wes_pca
      results/phenotype/             #   cancer_counts.tsv, roster.tsv, cases.tsv, famhx.tsv
      results/analysis/              #   analysis_ready.tsv
      results/models/                #   ML outputs
      logs/
  code/Cancer-risk-MSM/              # this git clone
```

Which `results-batch*` dirs are authoritative is set in
`config.inputs.carrier_source.results_dirs`. Pooling multiple runs OR-s the
per-framework P/LP flags per (person, variant), and refuses to pool runs whose
`manifest.json` disagrees on reference build.

## Quick start (on Minerva)

```bash
cd /sc/arion/projects/rg_huangk06/code
git clone https://github.com/Huang-lab/Cancer-risk-MSM.git
cd Cancer-risk-MSM

module load anaconda3
mamba env create -f envs/analysis.yml && mamba activate msm-analysis
export CONFIG=config/config.yaml

bash workflow/00_setup_dirs.sh          # create the runtime tree
```

### Step 1 — Sample QC (from raw pVCFs)

Sample QC needs genotypes only, so it reads the immutable source pVCFs directly
and does not wait on the NF pipeline.

```bash
bash workflow/00_run_sample_qc.sh       # submits per-chunk `bcftools stats` array
# when the array is DONE:
bash workflow/00_wes_qc_aggregate.sh    # PSC aggregate + plink2 kinship/sex-check → keep-list
```

| Check | Threshold (`qc.sample`) | Source |
|---|---|---|
| Missingness | `missingness_max` (0.05) | `bcftools stats -s -` PSC |
| Contamination proxy | het/hom `\|z\| > het_hom_z_max` (4) | `bcftools stats -s -` PSC |
| Sex-check | F < 0.20 female, F > 0.80 male; `mismatch_action` | `plink2 --check-sex` on chrX |
| Duplicate / MZ twin | kinship ≥ `duplicate_kin_min` (0.354) → keep highest call rate | `plink2 --make-king` |
| Related (≥2nd degree) | ≥ `pi_hat_related_max` (0.1875) → flagged, kept | `plink2 --make-king` |

Bonus: `plink2 --pca 20` on the same LD-pruned substrate gives WES-derived
ancestry PCs (`results/qc/wes_pca.eigenvec`) for use as ML covariates.

*Caveat:* WES-based kinship is rougher than array-common-variant kinship (fewer
LD-pruned markers, more depth noise). Duplicate calls at ≥0.354 are robust;
2nd-degree calls are less certain and only flag+keep, never exclude.

### Step 2 — Phenotype roster (case / famhx / control)

```bash
bash workflow/09_export_wes_samples.sh   # bcftools query -l → resources/wes_samples.txt
# set config.inputs.wes.sample_list_file to that path, then:
bash workflow/10_build_phenotype_roster.sh
column -t -s $'\t' $OUTPUT_ROOT/results/phenotype/cancer_counts.tsv
```

Groups, per (person, cancer):

| Group | Definition |
|---|---|
| `case` | Cancer diagnosis (ICD → phecodeX) for this cancer |
| `famhx` | No cancer diagnosis, but a 1st- or 2nd-degree relative had cancer |
| `control` | No cancer diagnosis, no family history of cancer |
| `excluded` | Has a *different* cancer (`control_exclusion: any_cancer`) |

**Deliverable — `cancer_counts.tsv`:** one row per cancer with
`n_ehr_cases`, `n_ehr_cases_with_wes`, `n_famhx_1deg_any_cancer`,
`n_famhx_1deg_this_cancer`, `n_controls`, each also intersected with WES.
Counts are deduplicated per participant (many ICD rows → one case).

Modules:
- `src/phenotype/icd_mapping.py` — ICD-9/ICD-10 → phecodeX, longest-prefix match.
  Built-in map covers colorectal, breast, prostate, lung, ovarian, pancreatic;
  `phenotype.icd_to_phecodex_map` overrides with a full crosswalk.
- `src/phenotype/cases.py` — `Encounter_Diagnosis` + `Problem_List` (both
  ICD-coded), first-dx dates, incident vs prevalent.
- `src/phenotype/famhx.py` — relationship → 1st/2nd degree; free-text condition
  → cancer type (BioMe famhx has no ICD codes).
- `src/phenotype/roster.py` — group assignment + `cancer_counts.tsv`.

All EHR column names live in `config.phenotype.ehr_schema.*`, matching the real
BioMe phenofile layouts (pipe-delimited, `MM/DD/YYYY`).

### Step 3 — Analysis-ready join (scaffolded)

`src/carriers/nf_carriers.py` is the read-boundary to the NF pipeline:
`read_carriers` pools `carrier_matrix.tsv` across runs, `read_manifests` +
`assert_consistent_builds` fail fast on reference-build mismatch,
`per_gene_flags` / `panel_carrier_flags` build per-framework features.
`src/carriers/join.py` then joins carriers ⋈ keep-list ⋈ PCs ⋈ roster into
`results/analysis/analysis_ready.tsv`.

### Step 4 — Cancer-risk ML (scaffolded)

Add a cancer by adding `config/<cancer>.yaml` (provided: `crc.yaml`,
`breast.yaml`). Pipeline under `src/ml/`: dataset → match (propensity k:1
caliper) → features (window `[-730d, -182d]`) → models (penalized logistic /
Cox / RF / XGBoost) → evaluate (AUC, PR-AUC, ancestry-stratified calibration,
SHAP) → external_validate (cross-cohort with BioMe).

## Repo layout

```
config/          config.yaml (single source of truth) + crc.yaml + breast.yaml
src/qc/          sample_qc, vcf_stats (bcftools/plink2 helpers), site_gt_qc
src/phenotype/   icd_mapping, cases, famhx, roster, covariates, id_harmonize, vcf_samples
src/carriers/    nf_carriers (read germline-plp-carrier-nf), join
src/ml/          dataset, match, features, models, evaluate, external_validate
src/common/      config loader (input auto-resolution), PHI-safe logger
workflow/        00_setup_dirs, 00_run_sample_qc → 00_wes_qc_chunks →
                 00_wes_qc_aggregate, 09_export_wes_samples,
                 10_build_phenotype_roster
envs/            analysis.yml (python + bcftools/plink2 + ML stack)
docs/            data_dictionary.md, governance.md
tests/           synthetic fixtures + 28 unit tests
```

## Tests

```bash
python -m pytest tests -q      # 28 tests, synthetic fixtures only, no PHI
```

See [`docs/data_dictionary.md`](docs/data_dictionary.md) for output schemas and
[`docs/governance.md`](docs/governance.md) for the pre-commit audit.
