# Data dictionary (schemas only — no records)

## `results/variants/clinvar_plp.tsv` (long)
| column          | type   | note                                          |
|-----------------|--------|-----------------------------------------------|
| chr             | str    | GRCh38, "chr1".."chrY"                        |
| pos             | int    | 1-based                                       |
| ref             | str    |                                               |
| alt             | str    | single allele (post-`bcftools norm -m -any`)  |
| gene            | str    | HGNC symbol (MANE Select / canonical)         |
| clnsig          | str    | ClinVar CLNSIG                                |
| clnrevstat      | str    | ClinVar review status                         |
| stars           | int    | 0..4 gold stars                               |
| is_clinvar_PLP  | bool   | P/LP and >= min_review_stars, no conflicts    |

## `results/variants/acmg_plp.tsv` (long)
| column           | type | note                                            |
|------------------|------|-------------------------------------------------|
| chr,pos,ref,alt,gene | ...  | as above                                        |
| rules            | str  | ACMG rules fired, comma-separated (e.g. PVS1,PM2,PP3) |
| is_acmg_PLP      | bool | InterVar P/LP call, post ClinVar-B/LB filter    |
| clinvar_removed  | bool | true if the row was ACMG-P/LP but ClinVar B/LB  |

## `results/variants/am_plp.tsv` (long)
| column               | type  | note                                     |
|----------------------|-------|------------------------------------------|
| chr,pos,ref,alt,gene | ...   |                                          |
| am_score             | float | AlphaMissense pathogenicity              |
| evidence_gene        | str   | gene-specific Chen/Pejaver label         |
| evidence_domain      | str   | recorded only; never promotes            |
| is_AM_PLP            | bool  | evidence_gene >= min_evidence_strength   |

## `results/carriers/carriers.long.tsv`  (**one row per variant × person**)
| column         | type | note                                                     |
|----------------|------|----------------------------------------------------------|
| chr            | str  |                                                          |
| pos            | int  |                                                          |
| ref            | str  |                                                          |
| alt            | str  |                                                          |
| gene           | str  | HGNC symbol                                              |
| person_id      | str  | canonical (harmonized old/new SINAI IDs)                 |
| is_clinvar_PLP | bool | this variant is ClinVar-framework P/LP                   |
| is_acmg_PLP    | bool | ACMG framework (ClinVar-B/LB removed)                    |
| is_AM_PLP      | bool | AlphaMissense gene-specific calibrated P/LP              |

## `results/carriers/carriers.wide.tsv` (optional pivot)
Sample × gene matrix with per-framework aggregation
(e.g. `<gene>__is_clinvar_PLP = any(is_clinvar_PLP over that gene for the person)`).

## `results/analysis/analysis_ready.tsv`
Long carrier table joined with:
- canonical `person_id`
- ancestry PCs (from `common_variants/.../pca/pcs.tsv`)
- phenotype: case/control, first_dx_date, incident flag, family-hx flags,
  demographics + covariates (age, sex, BMI, smoking, alcohol, parity, screening).
