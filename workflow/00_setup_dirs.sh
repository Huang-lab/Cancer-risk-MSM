#!/usr/bin/env bash
# Create the gitignored runtime directory tree on Minerva.
#
# NOTE: we do NOT create data/annotated/ — variant annotation is owned by
# germline-plp-carrier-nf, whose outputs we read from
# inputs.carrier_source.nf_results_root.
set -euo pipefail

CONFIG="${CONFIG:-config/config.yaml}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")}"

echo "Setting up runtime tree under: $OUTPUT_ROOT"
mkdir -p \
  "$OUTPUT_ROOT/resources" \
  "$OUTPUT_ROOT/results/qc" \
  "$OUTPUT_ROOT/results/qc/stats" \
  "$OUTPUT_ROOT/results/phenotype" \
  "$OUTPUT_ROOT/results/analysis" \
  "$OUTPUT_ROOT/results/models" \
  "$OUTPUT_ROOT/logs/qc" \
  "$OUTPUT_ROOT/logs/phenotype" \
  "$OUTPUT_ROOT/logs/analysis" \
  "$OUTPUT_ROOT/logs/ml"

echo "OK"
