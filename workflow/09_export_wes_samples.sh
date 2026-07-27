#!/usr/bin/env bash
# =============================================================================
# Export the WES sample list from a raw pVCF chunk via `bcftools query -l`.
# All chunks share the same sample set (joint-called), so one chunk suffices.
#
# Output: <output_root>/resources/wes_samples.txt (one sample per line).
# Then set config.inputs.wes.sample_list_file to that path so the phenotype
# roster can report EHR-case counts intersected with genotyped samples.
# =============================================================================
set -euo pipefail

CONFIG="${CONFIG:-config/config.yaml}"
OUTPUT_ROOT="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")"

if command -v module >/dev/null 2>&1; then
  module load bcftools/1.19 || true
fi

mkdir -p "$OUTPUT_ROOT/resources"
OUT="$OUTPUT_ROOT/resources/wes_samples.txt"

FIRST=$(python -c "from src.common.config import load, resolve; print(resolve(load('$CONFIG')).wes_chunks[0])")
echo "Reading sample list from: $FIRST"
bcftools query -l "$FIRST" > "$OUT"
echo "Wrote $(wc -l < "$OUT") sample IDs to $OUT"
echo "Now set config.inputs.wes.sample_list_file: $OUT"
