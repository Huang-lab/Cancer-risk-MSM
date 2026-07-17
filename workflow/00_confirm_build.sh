#!/usr/bin/env bash
# Confirm GRCh38 from a WES VCF header before annotation. Log tool versions.
# Run once, on Minerva, before submitting anything.
set -euo pipefail

CONFIG="${CONFIG:-config/config.yaml}"
OUTPUT_ROOT="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")"
LOG="$OUTPUT_ROOT/logs/annotation/confirm_build.$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"

# Load Minerva modules if available
if command -v module >/dev/null 2>&1; then
  while read -r mod; do module load "$mod" || true; done < <(
    python -c "import yaml; [print(m) for m in yaml.safe_load(open('$CONFIG'))['lsf']['modules']]"
  )
fi

exec > >(tee -a "$LOG") 2>&1

echo "=== Tool versions ==="
bcftools --version | head -2
vep --help 2>&1 | grep -E 'ensembl-vep' | head -1 || vep 2>&1 | head -5
tabix --version 2>&1 | head -1
samtools --version | head -1

echo "=== Resolve one WES chunk ==="
FIRST_CHUNK="$(python -c "
from src.common.config import load, resolve
r = resolve(load('$CONFIG'))
print(r.wes_chunks[0])
")"
echo "First chunk: $FIRST_CHUNK"

echo "=== Reference build assertion ==="
if bcftools view -h "$FIRST_CHUNK" | grep -Ei 'GRCh38|assembly=hg38' >/dev/null; then
  echo "OK: header advertises GRCh38"
else
  if bcftools view -h "$FIRST_CHUNK" | grep -E '^##contig=<ID=(chr)?1,length=248956422' >/dev/null; then
    echo "OK: chr1 length matches GRCh38 (248,956,422)"
  else
    echo "FAIL: could not confirm GRCh38 from header"
    exit 2
  fi
fi

echo "=== Resources present ==="
python -c "
from pathlib import Path
import yaml
cfg = yaml.safe_load(open('$CONFIG'))
res = cfg['resources']
missing = []
for label, p in [
    ('fasta',                res['fasta']),
    ('vep.cache_dir',        res['vep']['cache_dir']),
    ('loftee.plugin_dir',    res['loftee']['plugin_dir']),
    ('alphamissense.scores', res['alphamissense']['scores_tsv_gz']),
    ('gnomad.vcf',           res['gnomad']['vcf']),
    ('clinvar.root',         res['clinvar']['root']),
]:
    if not Path(p).exists():
        missing.append((label, p))
for label, p in missing:
    print(f'MISSING: {label} -> {p}')
raise SystemExit(1 if missing else 0)
"

echo "=== Build confirmation complete ==="
