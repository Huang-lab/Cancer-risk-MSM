#!/usr/bin/env bash
# =============================================================================
# Build phenotype roster (case / famhx / control per person x cancer).
# Independent of the annotation phase — can run any time on Minerva.
# =============================================================================
set -euo pipefail

CONFIG="${CONFIG:-config/config.yaml}"
OUTPUT_ROOT="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")"
QUEUE="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['lsf']['queue'])")"
PROJECT="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['lsf']['project'])")"

OUT_DIR="$OUTPUT_ROOT/results/phenotype"
LSF_LOG_DIR="$OUTPUT_ROOT/logs/phenotype/lsf"
mkdir -p "$OUT_DIR" "$LSF_LOG_DIR"

bsub \
  -J "phenotype_roster" \
  -q "$QUEUE" -P "$PROJECT" \
  -n 2 -R "rusage[mem=8000]" -W "01:00" \
  -o "$LSF_LOG_DIR/roster.%J.out" \
  -e "$LSF_LOG_DIR/roster.%J.err" \
  bash -lc "cd $(pwd) && CONFIG='$CONFIG' \
    python -m src.phenotype.roster --config '$CONFIG' --out-dir '$OUT_DIR'"

echo "Submitted. Monitor: bjobs -J phenotype_roster"
echo "When it finishes, review $OUT_DIR/roster_summary.md"
