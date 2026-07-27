# Data dictionary (schemas only — no records)

## Inputs we consume (produced by `germline-plp-carrier-nf`)

### `<nf_results_root>/<results-dir>/carriers/carrier_matrix.tsv`
Long format — one row per (variant × person). Column names are configurable
under `inputs.carrier_source.columns`.

| column         | type | note                                          |
|----------------|------|-----------------------------------------------|
| chr            | str  | Ensembl-style contig (NF normalizes `chr1`→`1`)|
| pos            | int  | 1-based                                       |
| ref            | str  |                                               |
| alt            | str  | single allele (post-norm)                     |
| gene           | str  | HGNC symbol                                   |
| person_id      | str  | sample ID as it appears in the pVCF           |
| is_clinvar_PLP | bool | ClinVar framework P/LP                        |
| is_acmg_PLP    | bool | ACMG framework (ANNOVAR+InterVar)             |
| is_AM_PLP      | bool | AlphaMissense gene-specific calibrated P/LP   |

### `<results-dir>/manifest.json`
Reference build, tool versions, pipeline SHA. We read `reference_build` and
**refuse to pool runs that disagree** (`assert_consistent_builds`).

### `<results-dir>/variants/{clinvar_plp,acmg_plp,am_plp,qc_per_gene}.tsv`
Variant-level classification tables. Available for per-gene QC review; the
carrier matrix is the primary join key for ML features.

### NF sample QC outputs
Sample-level QC (missingness, kinship, sex-check) is performed by the NF
pipeline. Configure where to find its outputs under
`inputs.carrier_source`:

| config key | expects | used for |
|---|---|---|
| `sample_keep_list` | one QC-passing sample per row (with or without header) | restrict carriers + roster to QC-passing samples |
| `sample_keep_list_col` | column name if the file has a header (default `sample_id`) | — |
| `relatedness_table` | pairwise kinship (e.g. KING `.kin0`) | group-aware cross-validation so relatives never split across CV folds |

Both are `null` until confirmed. With `sample_keep_list` unset,
`filter_to_keep_list` is a no-op rather than excluding every sample.

---

## Outputs we produce

### `results/phenotype/cases.tsv`
| column | type | note |
|---|---|---|
| person_id | str | |
| cancer | str | key from `phenotype.cancer_phecodes` |
| is_case | bool | |
| first_dx_date | date | earliest matching ICD encounter, ISO |
| incident_flag | bool | first_dx_date > enrollment_date; blank if unknown |
| n_dx_encounters | int | count of matching rows (pre-dedup) |
| sources | str | `encounter_diagnosis` / `problem_list` / `phecodex_roster` |

### `results/phenotype/famhx.tsv`
One row per person: `has_famhx_1deg_any_cancer`, `has_famhx_2deg_any_cancer`,
`n_affected_1deg`, `n_affected_2deg`, plus `has_famhx_{1,2}deg_<cancer>` per
cancer key.

### `results/phenotype/roster.tsv`
One row per (person × cancer).

| column | type | note |
|---|---|---|
| person_id | str | |
| cancer | str | |
| group | str | `case` / `famhx` / `control` / `excluded` |
| first_dx_date | date | cases only |
| incident_flag | bool | |
| has_wes | bool | person appears in `inputs.wes.sample_list_file` |
| has_famhx_{1,2}deg_any_cancer | bool | |
| has_famhx_{1,2}deg_this_cancer | bool | cancer-type-specific famhx |
| n_encounters | int | total EHR encounter rows |
| excluded_from_control_reason | str | `has_other_cancer` / `too_few_encounters` |

### `results/phenotype/cancer_counts.tsv` — primary deliverable
One row per cancer; every count deduplicated per participant.

| column | note |
|---|---|
| cancer | |
| n_ehr_cases | ICD-defined cases |
| n_ehr_cases_with_wes | ∩ genotyped samples |
| n_famhx_1deg_any_cancer | 1st-degree relative with any cancer |
| n_famhx_1deg_any_cancer_with_wes | ∩ genotyped |
| n_famhx_1deg_this_cancer | 1st-degree relative with *this* cancer |
| n_famhx_1deg_this_cancer_with_wes | ∩ genotyped |
| n_controls | no cancer dx, no cancer famhx |
| n_controls_with_wes | ∩ genotyped |

### `results/analysis/analysis_ready.tsv`
Roster ⋈ carrier features ⋈ PCs ⋈ covariates — one row per (person, cancer)
with per-framework per-gene carrier flags and the focus-panel aggregate.
