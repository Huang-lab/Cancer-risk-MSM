"""Join NF carrier calls with sample QC, harmonized IDs, PCs, and phenotypes.

Produces the analysis-ready table the ML stage consumes.

Inputs
------
- germline-plp-carrier-nf carrier calls  (via src.carriers.nf_carriers)
- results/qc/sample_keep_list.tsv        (our sample QC)
- results/qc/wes_pca.eigenvec            (our WES-derived PCs) OR the biobank's
                                          common-variants PCs
- results/phenotype/roster.tsv           (our case/famhx/control roster)

Output
------
results/analysis/analysis_ready.tsv — one row per (person, cancer) with the
roster group, covariates, and per-framework carrier features.
"""
from __future__ import annotations

# TODO(minerva): implement build_analysis_table(cfg) that:
#   1. rows = nf_carriers.read_carriers(cfg)
#      manifests = nf_carriers.read_manifests(cfg)
#      nf_carriers.assert_consistent_builds(manifests)   # fail fast on build mismatch
#   2. keep = load sample_keep_list.tsv; rows = filter_to_keep_list(rows, keep)
#   3. per-framework features:
#        nf_carriers.per_gene_flags(rows, framework)       -> per-gene carrier sets
#        nf_carriers.panel_carrier_flags(rows, framework, panel)  -> aggregate flag
#   4. left-join roster.tsv (person_id, cancer, group, first_dx_date, famhx flags)
#   5. left-join PCs + covariates
#   6. write results/analysis/analysis_ready.tsv
