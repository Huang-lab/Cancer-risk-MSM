#!/usr/bin/env bash
# =============================================================================
# Download / refresh reference resources into $OUTPUT_ROOT/resources/.
# Currently: ClinVar weekly release (GRCh38). Extensible to gnomAD, AlphaMissense,
# ANNOVAR humandb clinvar_latest, etc.
#
# ClinVar is downloaded from NCBI FTP, md5-verified, placed in a dated folder
# `resources/clinvar/YYYY-MM-DD/clinvar.vcf.gz`, and tabix-indexed. VEP's
# `--custom` step in 01_norm_vep_chunk.sh then annotates against this file
# (that's the annotation — --custom is exactly "read this VCF and attach these
# fields per matching site"). config.resources.clinvar.release: "latest" makes
# the resolver pick the newest dated folder automatically.
# =============================================================================
set -euo pipefail

CONFIG="${CONFIG:-config/config.yaml}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")}"

CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_MD5="${CLINVAR_URL}.md5"
CLINVAR_TBI="${CLINVAR_URL}.tbi"

# Load modules for htslib (tabix) + curl if needed
if command -v module >/dev/null 2>&1; then
  module load htslib/1.19 2>/dev/null || true
fi

STAGING="$(mktemp -d -t clinvar.XXXXXX)"
trap "rm -rf '$STAGING'" EXIT

echo "[$(date -Iseconds)] Downloading ClinVar from NCBI"
curl -sSfL "$CLINVAR_URL"     -o "$STAGING/clinvar.vcf.gz"
curl -sSfL "$CLINVAR_MD5"     -o "$STAGING/clinvar.vcf.gz.md5" || true
curl -sSfL "$CLINVAR_TBI"     -o "$STAGING/clinvar.vcf.gz.tbi" || true

# Verify md5 if NCBI provided one
if [[ -s "$STAGING/clinvar.vcf.gz.md5" ]]; then
  echo "Verifying md5"
  ( cd "$STAGING" && md5sum -c clinvar.vcf.gz.md5 )
else
  echo "WARN: no md5 sidecar from NCBI; skipping verification"
fi

# Release date from ##fileDate= header, e.g. ##fileDate=2026-07-16
REL_DATE="$(zcat "$STAGING/clinvar.vcf.gz" | head -400 | awk -F= '/^##fileDate=/ {print $2; exit}')"
if [[ -z "$REL_DATE" ]]; then
  echo "FAIL: could not parse ##fileDate= from ClinVar header" >&2
  exit 2
fi

DEST_DIR="$OUTPUT_ROOT/resources/clinvar/$REL_DATE"
mkdir -p "$DEST_DIR"
mv "$STAGING/clinvar.vcf.gz" "$DEST_DIR/clinvar.vcf.gz"

# Trust the FTP .tbi if present; otherwise build one
if [[ -s "$STAGING/clinvar.vcf.gz.tbi" ]]; then
  mv "$STAGING/clinvar.vcf.gz.tbi" "$DEST_DIR/clinvar.vcf.gz.tbi"
else
  tabix -p vcf -f "$DEST_DIR/clinvar.vcf.gz"
fi

echo "$REL_DATE" > "$DEST_DIR/RELEASE_DATE"
echo "[$(date -Iseconds)] ClinVar release $REL_DATE staged at $DEST_DIR"
echo "Set config.resources.clinvar.release: latest to auto-select it."
