#!/usr/bin/env bash
# Create the gitignored runtime directory tree on Minerva.
# Reads $OUTPUT_ROOT from config.yaml if not overridden.
set -euo pipefail

CONFIG="${CONFIG:-config/config.yaml}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$(python -c "import yaml,sys; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")}"

echo "Setting up runtime tree under: $OUTPUT_ROOT"
mkdir -p \
  "$OUTPUT_ROOT/data/annotated" \
  "$OUTPUT_ROOT/resources" \
  "$OUTPUT_ROOT/results/variants" \
  "$OUTPUT_ROOT/results/carriers" \
  "$OUTPUT_ROOT/results/analysis" \
  "$OUTPUT_ROOT/results/models" \
  "$OUTPUT_ROOT/logs/annotation" \
  "$OUTPUT_ROOT/logs/classification" \
  "$OUTPUT_ROOT/logs/phenotype" \
  "$OUTPUT_ROOT/logs/carriers" \
  "$OUTPUT_ROOT/logs/ml"

echo "OK"
