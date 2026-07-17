#!/usr/bin/env bash
# =============================================================================
# Single-chunk annotation worker: bcftools norm -> VEP -> bgzip + tabix.
# Called by 02_submit_single.lsf and by each array element of 03_submit_array.lsf.
# Resume-safe: exits 0 immediately if a validated output already exists.
#
# Args:
#   $1 (INPUT_VCF)  absolute path to a bgzipped WES Target chunk .vcf.gz
# Env (from LSF submit):
#   CONFIG          path to config.yaml
#   OUTPUT_ROOT     runtime tree root
# =============================================================================
set -euo pipefail

INPUT_VCF="${1:?usage: 01_norm_vep_chunk.sh <input.vcf.gz>}"
CONFIG="${CONFIG:-config/config.yaml}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")}"

chunk="$(basename "$INPUT_VCF" .vcf.gz)"
OUT_DIR="$OUTPUT_ROOT/data/annotated"
LOG_DIR="$OUTPUT_ROOT/logs/annotation"
NORM_VCF="$OUT_DIR/${chunk}.norm.vcf.gz"
QC_VCF="$OUT_DIR/${chunk}.norm.qc.vcf.gz"
OUT_VCF="$OUT_DIR/${chunk}.annot.vcf.gz"
LOG="$LOG_DIR/${chunk}.log"
mkdir -p "$OUT_DIR" "$LOG_DIR"

exec > >(tee -a "$LOG") 2>&1
echo "[$(date -Iseconds)] START $chunk on $(hostname)"

# --- Resume-safe skip ---
if [[ -s "$OUT_VCF" && -s "${OUT_VCF}.tbi" ]]; then
  if python -m src.annotation.validate_chunk --input "$INPUT_VCF" --annotated "$OUT_VCF" --min-ratio 0.5; then
    echo "SKIP: validated output already exists"
    exit 0
  else
    echo "Existing output failed validation; re-running"
    rm -f "$OUT_VCF" "${OUT_VCF}.tbi" "$NORM_VCF" "${NORM_VCF}.tbi"
  fi
fi

# --- Load modules (LSF worker) ---
if command -v module >/dev/null 2>&1; then
  while read -r mod; do module load "$mod" || true; done < <(
    python -c "import yaml; [print(m) for m in yaml.safe_load(open('$CONFIG'))['lsf']['modules']]"
  )
fi

# --- Read config values ---
read FASTA VEP_CACHE VEP_ASSEMBLY VEP_SPECIES VEP_FORK \
     LOFTEE_DIR LOFTEE_ANC LOFTEE_CONS \
     AM_SCORES AM_PLUGIN_DIR \
     GNOMAD_VCF DBNSFP CLINVAR_VCF \
     <<< "$(python -c "
import yaml
from src.common.config import load, resolve
cfg = load('$CONFIG'); r = resolve(cfg)
V = cfg['resources']
print(V['fasta'], V['vep']['cache_dir'], V['vep']['assembly'], V['vep']['species'], V['vep']['fork'],
      V['loftee']['plugin_dir'], V['loftee']['human_ancestor_fa'], V['loftee']['conservation_file'],
      V['alphamissense']['scores_tsv_gz'], V['alphamissense']['plugin_dir'],
      V['gnomad']['vcf'], V.get('dbnsfp',{}).get('tsv_gz','NONE'), r.clinvar_vcf)
")"

# --- 1. bcftools norm: split multiallelics, left-align to reference ---
echo "[$(date -Iseconds)] bcftools norm"
bcftools norm \
  -m -any \
  -f "$FASTA" \
  --check-ref w \
  -Oz -o "$NORM_VCF" \
  "$INPUT_VCF"
tabix -p vcf -f "$NORM_VCF"

# --- 1b. Site + genotype QC (config-driven; between norm and VEP) ---
echo "[$(date -Iseconds)] site + genotype QC"
read KEEP_FILTER DP_MIN GQ_MIN AB_MIN AB_MAX SITE_MISS DROP_MONO \
     <<< "$(python -c "
import yaml
q = yaml.safe_load(open('$CONFIG'))['qc']['site_gt']
print(q['keep_filter'], q['genotype_dp_min'], q['genotype_gq_min'],
      q['genotype_ab_het_min'], q['genotype_ab_het_max'],
      q['site_missing_max'], str(q['drop_monoallelic_after_mask']).lower())
")"
GT_INCLUDE="$(python -c "
from src.qc.site_gt_qc import genotype_mask_include_expr
print(genotype_mask_include_expr($DP_MIN, $GQ_MIN, $AB_MIN, $AB_MAX))
")"
SITE_DROP="$(python -c "
from src.qc.site_gt_qc import site_drop_expr
print(site_drop_expr($SITE_MISS, ${DROP_MONO@Q} == 'true'))
")"

# view -f PASS  ->  +setGT to null low-quality genotypes  ->  drop sites failing site filters
bcftools view -f "$KEEP_FILTER" "$NORM_VCF" \
  | bcftools +setGT -- -t q -n . -i "$GT_INCLUDE" \
  | bcftools +fill-tags -- -t AC,AN,F_MISSING \
  | bcftools view -e "$SITE_DROP" -Oz -o "$QC_VCF"
tabix -p vcf -f "$QC_VCF"

# --- 2. VEP annotation ---
echo "[$(date -Iseconds)] vep --fork $VEP_FORK"
VEP_CMD=(vep
  --input_file "$QC_VCF"
  --output_file "$OUT_VCF"
  --vcf --compress_output bgzip
  --cache --offline
  --dir_cache "$VEP_CACHE"
  --dir_plugins "$LOFTEE_DIR:$AM_PLUGIN_DIR"
  --fasta "$FASTA"
  --assembly "$VEP_ASSEMBLY"
  --species "$VEP_SPECIES"
  --fork "$VEP_FORK"
  --mane_select --pick_allele_gene --canonical --biotype --symbol --hgvs --numbers
  --plugin "LoF,loftee_path:$LOFTEE_DIR,human_ancestor_fa:$LOFTEE_ANC,conservation_file:$LOFTEE_CONS"
  --plugin "AlphaMissense,file=$AM_SCORES"
  --custom "$CLINVAR_VCF,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNVI"
  --custom "$GNOMAD_VCF,gnomAD,vcf,exact,0,AF,AF_popmax"
  --stats_file "$LOG_DIR/${chunk}.vep_stats.html"
  --warning_file "$LOG_DIR/${chunk}.vep_warnings.txt"
  --force_overwrite
)
if [[ "$DBNSFP" != "NONE" && -s "$DBNSFP" ]]; then
  VEP_CMD+=(--plugin "dbNSFP,$DBNSFP,SIFT_score,Polyphen2_HDIV_score,REVEL_score,CADD_phred")
fi
"${VEP_CMD[@]}"

tabix -p vcf -f "$OUT_VCF"

# --- 3. Validate ---
python -m src.annotation.validate_chunk --input "$INPUT_VCF" --annotated "$OUT_VCF" --min-ratio 0.5

# --- 4. Clean intermediate norm + QC outputs (keep only annotated) ---
rm -f "$NORM_VCF" "${NORM_VCF}.tbi" "$QC_VCF" "${QC_VCF}.tbi"

echo "[$(date -Iseconds)] DONE $chunk"
