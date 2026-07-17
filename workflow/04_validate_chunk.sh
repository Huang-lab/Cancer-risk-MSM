#!/usr/bin/env bash
# Thin shell wrapper over src.annotation.validate_chunk for ad-hoc checks.
# Usage: bash workflow/04_validate_chunk.sh <input.vcf.gz> <annotated.vcf.gz>
set -euo pipefail
python -m src.annotation.validate_chunk --input "$1" --annotated "$2" --min-ratio 0.5
