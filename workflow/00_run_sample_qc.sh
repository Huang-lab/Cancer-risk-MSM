#!/usr/bin/env bash
# =============================================================================
# Sample QC from VCFs (runs AFTER Step 1 annotation is complete).
#
# Pipeline:
#   00_wes_qc_chunks.lsf        per-chunk `bcftools stats -s -` (LSF array)
#   00_wes_qc_aggregate.sh      aggregate PSC + plink2 kinship/sex-check + build keep-list
#
# This dispatcher submits the array and prints the aggregate step to run next.
# We do not chain the two via bsub -w because the aggregate step is its own
# script (it submits child jobs and chains them internally).
# =============================================================================
set -euo pipefail

bash workflow/00_wes_qc_chunks.lsf
echo ""
echo "When the wes_qc_stats array is DONE, run:"
echo "  bash workflow/00_wes_qc_aggregate.sh"
