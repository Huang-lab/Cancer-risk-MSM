# Cancer-risk-MSM

EHR phenotyping and cancer-risk modeling for the Mount Sinai Million (MSM)
biobank.

> **Scope.** This repo is the **downstream** layer. Variant annotation, P/LP
> classification, carrier extraction, and sample-level QC are all owned by
> [`Huang-lab/germline-plp-carrier-nf`](https://github.com/Huang-lab/germline-plp-carrier-nf).
> We consume its outputs read-only.

Methodology and config style mirror
[`Huang-lab/BioMe-cancer-risk-prediction`](https://github.com/Huang-lab/BioMe-cancer-risk-prediction).

> **Governance.** Code only. IRB-restricted data stays on Minerva, gitignored.
> See [`docs/governance.md`](docs/governance.md).

---

## Division of labor

| | germline-plp-carrier-nf | Cancer-risk-MSM (this repo) |
|---|---|---|
| Annotation (VEP, ClinVar, AlphaMissense, LOFTEE, gnomAD) | ✅ | — |
| P/LP classification (ClinVar / ACMG / AlphaMissense) | ✅ | — |
| Carrier extraction | ✅ | — |
| Sample-level QC (missingness, kinship, sex-check) | ✅ | — |
| **EHR phenotyping** (cases / famhx / controls) | — | ✅ |
| **Analysis-ready joins** | — | ✅ |
| **Cancer-risk ML** | — | ✅ |

## Directory layout on Minerva

Both repos live inside the shared project dir, as siblings. Code and outputs
are separated by naming, not by a separate tree:

```
/sc/arion/projects/rg_huangk06/variants_PLP_MSM/
│
├── germline-plp-carrier-nf/      # REPO 1 (git) — upstream pipeline
├── Cancer-risk-MSM/              # REPO 2 (git) — this repo
│
├── refs/                         # shared reference data (VEP cache, ClinVar, ...)
│
├── nf-home/                      # NF runtime state
├── batch1.config  batch2.config  batch3.config  chr13.config  ...
├── input-batch1/   input-batch2/   input-chr13/   ...
├── launch-batch1/  launch-batch2/  launch-chr13/  ...
├── results-batch1/ results-batch2/ results-batch3/         ← NF outputs we READ
├── results-pilot/  results-chr1/   results-chr13/  ...        (pilot / per-chr runs)
├── logs/
│
└── cancer_risk/                  # ← OUR outputs ($OUTPUT_ROOT). We're the only writer.
    ├── resources/                #   wes_samples.txt
    ├── results/phenotype/        #   cancer_counts.tsv, roster.tsv, cases.tsv, famhx.tsv
    ├── results/analysis/         #   analysis_ready.tsv
    ├── results/models/           #   ML outputs
    └── logs/
```

**Two rules that keep this sane as it grows:**

1. **One writer per directory.** NF writes `results-*`, `input-*`, `launch-*`,
   `nf-home`. We write only `cancer_risk/`. Neither pipeline ever writes into
   the other's tree — so a rerun of either can never clobber the other's state.
2. **Code is a git clone, never a data directory.** Both repo folders are
   `git status`-clean checkouts. Nothing is generated inside them; every output
   goes to `cancer_risk/` (us) or `results-*` (NF).

### Optional future tidy-up

The top level currently interleaves four things per NF run (`batchN.config`,
`input-batchN/`, `launch-batchN/`, `results-batchN/`). If that gets unwieldy,
a run-scoped convention groups them:

```
runs/batch1/{config, input/, launch/, results/}
runs/batch2/...
```

This requires updating path references inside NF's `*.config` and `launch-*`
scripts, so it's worth doing only while no NF run is in flight — and it's
entirely optional. Nothing in this repo depends on it: we resolve NF outputs
through `config.inputs.carrier_source.results_dirs`, so a reorganization is a
one-line config edit on our side.

## Quick start

```bash
cd /sc/arion/projects/rg_huangk06/variants_PLP_MSM
git clone https://github.com/Huang-lab/Cancer-risk-MSM.git
cd Cancer-risk-MSM

module load anaconda3
mamba env create -f envs/analysis.yml && mamba activate msm-analysis
export CONFIG=config/config.yaml

bash workflow/00_setup_dirs.sh          # creates ../cancer_risk/{resources,results,logs}
```

### Step 1 — Phenotype roster (runnable now)

```bash
bash workflow/09_export_wes_samples.sh   # bcftools query -l → resources/wes_samples.txt
# set config.inputs.wes.sample_list_file to that path, then:
bash workflow/10_build_phenotype_roster.sh
column -t -s $'\t' ../cancer_risk/results/phenotype/cancer_counts.tsv
```

Groups assigned per (person, cancer):

| Group | Definition |
|---|---|
| `case` | Cancer diagnosis (ICD → phecodeX) for this cancer |
| `famhx` | No cancer diagnosis, but a 1st- or 2nd-degree relative had cancer |
| `control` | No cancer diagnosis, no family history of cancer |
| `excluded` | Has a *different* cancer (`control_exclusion: any_cancer`) |

**Deliverable — `cancer_counts.tsv`:** one row per cancer with `n_ehr_cases`,
`n_famhx_1deg_any_cancer`, `n_famhx_1deg_this_cancer`, `n_controls`, each also
intersected with genotyped samples (`*_with_wes`). Counts are deduplicated per
participant, so a person with 12 encounter rows for one cancer counts once.

Modules:
- `src/phenotype/icd_mapping.py` — ICD-9/ICD-10 → phecodeX, longest-prefix match.
  Built-in map covers colorectal, breast, prostate, lung, ovarian, pancreatic;
  `phenotype.icd_to_phecodex_map` overrides with a full crosswalk.
- `src/phenotype/cases.py` — `Encounter_Diagnosis` + `Problem_List` (both
  ICD-coded), first-dx dates, incident vs prevalent.
- `src/phenotype/famhx.py` — relationship → 1st/2nd degree; free-text condition
  → cancer type (BioMe famhx has no ICD codes, only `problem_description`).
- `src/phenotype/roster.py` — group assignment + `cancer_counts.tsv`.

EHR column names live in `config.phenotype.ehr_schema.*`, matching the real
BioMe phenofile layouts (pipe-delimited, `MM/DD/YYYY` dates).

### Step 2 — Analysis-ready join (scaffolded)

`src/carriers/nf_carriers.py` is the read-boundary to NF:

- `read_carriers()` pools `carrier_matrix.tsv` across the configured
  `results-*` dirs, deduplicating on (person, variant) and OR-ing the
  per-framework P/LP flags.
- `read_manifests()` + `assert_consistent_builds()` **hard-fail** if pooled runs
  disagree on reference build, which would otherwise silently mis-join coordinates.
- `read_sample_keep_list()` reads NF's sample-QC keep-list so we restrict to
  QC-passing samples without recomputing QC.
- `per_gene_flags()` / `panel_carrier_flags()` build per-framework ML features.

`src/carriers/join.py` then joins carriers ⋈ keep-list ⋈ PCs ⋈ roster into
`cancer_risk/results/analysis/analysis_ready.tsv`.

### Step 3 — Cancer-risk ML (scaffolded)

Add a cancer by adding `config/<cancer>.yaml` (provided: `crc.yaml`,
`breast.yaml`). Pipeline under `src/ml/`: dataset → match (propensity k:1
caliper) → features (window `[-730d, -182d]`) → models (penalized logistic /
Cox / RF / XGBoost) → evaluate (AUC, PR-AUC, ancestry-stratified calibration,
SHAP) → external_validate (cross-cohort with BioMe).

## Repo layout

```
config/          config.yaml (single source of truth) + crc.yaml + breast.yaml
src/phenotype/   icd_mapping, cases, famhx, roster, covariates, id_harmonize, vcf_samples
src/carriers/    nf_carriers (read germline-plp-carrier-nf), join
src/ml/          dataset, match, features, models, evaluate, external_validate
src/common/      config loader (input auto-resolution), PHI-safe logger
workflow/        00_setup_dirs, 09_export_wes_samples, 10_build_phenotype_roster
envs/            analysis.yml
docs/            data_dictionary.md, governance.md
tests/           synthetic fixtures + 21 unit tests
```

## Tests

```bash
python -m pytest tests -q      # 21 tests, synthetic fixtures only, no PHI
```

See [`docs/data_dictionary.md`](docs/data_dictionary.md) for schemas.
