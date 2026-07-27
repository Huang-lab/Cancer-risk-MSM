#!/usr/bin/env bash
# =============================================================================
# Export the WES sample list from any pVCF chunk via `bcftools query -l`.
# All chunks share the same sample set (joint-called), so one chunk is enough.
#
# Output: $OUTPUT_ROOT/resources/wes_samples.txt (one sample per line).
# Set config.inputs.wes.sample_list_file to this path for the roster step
# to intersect EHR cases with WES samples.
# =============================================================================
set -euo pipefail

CONFIG="${CONFIG:-config/config.yaml}"
OUTPUT_ROOT="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")"

if command -v module >/dev/null 2>&1; then
  module load bcftools/1.19 || true
fi

mkdir -p "$OUTPUT_ROOT/resources"
OUT="$OUTPUT_ROOT/resources/wes_samples.txt"

# Prefer any annotated chunk; fall back to raw pVCF chunk if annotation not done yet.
FIRST=$(ls "$OUTPUT_ROOT/data/annotated"/*.annot.vcf.gz 2>/dev/null | head -1 || true)
if [[ -z "$FIRST" ]]; then
  FIRST=$(python -c "from src.common.config import load, resolve; print(resolve(load('$CONFIG')).wes_chunks[0])")
fi
echo "Reading sample list from: $FIRST"
bcftools query -l "$FIRST" > "$OUT"
echo "Wrote $(wc -l < "$OUT") sample IDs to $OUT"
