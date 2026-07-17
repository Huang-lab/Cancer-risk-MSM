#!/usr/bin/env bash
# =============================================================================
# Sample QC (one-shot pre-annotation).
# Emits $OUTPUT_ROOT/results/qc/{sample_keep_list.tsv, sample_qc_report.tsv,
# sample_qc_summary.md}. Downstream carrier extraction filters to this keep-list.
# =============================================================================
set -euo pipefail

CONFIG="${CONFIG:-config/config.yaml}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")}"
QUEUE="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['lsf']['queue'])")"
PROJECT="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['lsf']['project'])")"

LSF_LOG_DIR="$OUTPUT_ROOT/logs/qc/lsf"
mkdir -p "$LSF_LOG_DIR" "$OUTPUT_ROOT/results/qc"

bsub \
  -J "sample_qc" \
  -q "$QUEUE" -P "$PROJECT" \
  -n 4 -R "rusage[mem=8000]" \
  -W "02:00" \
  -o "$LSF_LOG_DIR/sample_qc.%J.out" \
  -e "$LSF_LOG_DIR/sample_qc.%J.err" \
  bash -lc "cd $(pwd) && CONFIG='$CONFIG' OUTPUT_ROOT='$OUTPUT_ROOT' \
    python -m src.qc.sample_qc --config '$CONFIG' \
      --out-dir '$OUTPUT_ROOT/results/qc'"

echo "Submitted. Monitor: bjobs -J sample_qc"
echo "When it finishes, review $OUTPUT_ROOT/results/qc/sample_qc_summary.md"
