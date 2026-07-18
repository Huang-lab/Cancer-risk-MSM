#!/usr/bin/env bash
# =============================================================================
# Aggregate per-chunk bcftools stats into a per-sample PSC table, then run
# plink2 --make-king + --check-sex over autosomal / chrX WES data, then build
# results/qc/sample_keep_list.tsv via src.qc.sample_qc.
#
# Runs on Minerva after 00_wes_qc_chunks.lsf array finishes.
# =============================================================================
set -euo pipefail

CONFIG="${CONFIG:-config/config.yaml}"
OUTPUT_ROOT="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['project']['output_root'])")"
QUEUE="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['lsf']['queue'])")"
PROJECT="$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['lsf']['project'])")"

QC_DIR="$OUTPUT_ROOT/results/qc"
LSF_LOG_DIR="$OUTPUT_ROOT/logs/qc/lsf"
mkdir -p "$QC_DIR" "$LSF_LOG_DIR"

# --- 1. Aggregate PSC across chunks ---
echo "[$(date -Iseconds)] Aggregating per-chunk PSC stats"
python - <<'PY'
import glob, os, yaml, sys
from pathlib import Path
sys.path.insert(0, ".")
from src.qc.vcf_stats import parse_psc_lines, merge_psc

cfg = yaml.safe_load(open(os.environ["CONFIG"]))
root = cfg["project"]["output_root"]
files = sorted(glob.glob(f"{root}/results/qc/stats/*.stats.txt"))
rows = []
for p in files:
    rows.extend(parse_psc_lines(Path(p).read_text()))
agg = merge_psc(rows)

out = Path(f"{root}/results/qc/psc_agg.tsv")
with out.open("w") as fh:
    fh.write("sample_id\tn_called\tn_missing\tmissingness\thet_hom_ratio\tn_singletons\n")
    for s, r in sorted(agg.items()):
        fh.write(f"{s}\t{r.n_called}\t{r.n_missing}\t{r.missingness:.6f}\t{r.het_hom_ratio:.6f}\t{r.n_singletons}\n")
print(f"Wrote {out} ({len(agg)} samples)")
PY

# --- 2. Kinship: build LD-pruned common-variant bfile from WES, run --make-king ---
KIN_JOB=$(bsub \
  -J "wes_kinship" \
  -q "$QUEUE" -P "$PROJECT" \
  -n 16 -R "rusage[mem=64000]" -W "12:00" \
  -o "$LSF_LOG_DIR/kinship.%J.out" \
  -e "$LSF_LOG_DIR/kinship.%J.err" \
  bash -lc '
    set -euo pipefail
    module load plink2 bcftools/1.19 htslib/1.19 || true

    QCDIR='"$QC_DIR"'
    ROOT='"$OUTPUT_ROOT"'
    TMP=$QCDIR/kinship_wrk; mkdir -p "$TMP"

    # Merge autosomal chunks into a single site list for pruning.
    ls $ROOT/data/annotated/*.annot.vcf.gz > $TMP/vcfs.list
    bcftools concat -f $TMP/vcfs.list -a --regions $(seq -s, 1 22 | sed "s/,/,chr/g; s/^/chr/") \
      -Oz -o $TMP/autosomal.vcf.gz
    tabix -p vcf -f $TMP/autosomal.vcf.gz

    plink2 --vcf $TMP/autosomal.vcf.gz \
           --maf 0.01 --geno 0.05 --hwe 1e-10 \
           --indep-pairwise 500kb 50 0.2 \
           --threads 16 --out $TMP/prune
    plink2 --vcf $TMP/autosomal.vcf.gz \
           --extract $TMP/prune.prune.in \
           --make-bed --threads 16 --out $TMP/pruned
    plink2 --bfile $TMP/pruned \
           --make-king triangle bin \
           --king-table-filter 0.0442 \
           --threads 16 --out $QCDIR/king

    # Optional: PCA on the same pruned set (self-contained ancestry PCs)
    plink2 --bfile $TMP/pruned --pca 20 --threads 16 --out $QCDIR/wes_pca
  ' | awk "/is submitted/ {gsub(/[<>]/, \"\", \$2); print \$2}")

echo "kinship job: $KIN_JOB"

# --- 3. Sex-check: chrX non-PAR heterozygosity via plink2 --check-sex ---
SEX_JOB=$(bsub \
  -J "wes_sexcheck" \
  -q "$QUEUE" -P "$PROJECT" \
  -n 8 -R "rusage[mem=32000]" -W "06:00" \
  -o "$LSF_LOG_DIR/sexcheck.%J.out" \
  -e "$LSF_LOG_DIR/sexcheck.%J.err" \
  bash -lc '
    set -euo pipefail
    module load plink2 bcftools/1.19 htslib/1.19 || true
    QCDIR='"$QC_DIR"'; ROOT='"$OUTPUT_ROOT"'
    TMP=$QCDIR/sex_wrk; mkdir -p "$TMP"
    ls $ROOT/data/annotated/*.annot.vcf.gz \
      | xargs -I{} bcftools view -r chrX {} -Oz -o $TMP/$(basename {} .annot.vcf.gz).chrX.vcf.gz \
        --threads 2
    ls $TMP/*.chrX.vcf.gz > $TMP/chrX.list
    bcftools concat -f $TMP/chrX.list -a -Oz -o $TMP/chrX.vcf.gz
    tabix -p vcf -f $TMP/chrX.vcf.gz
    plink2 --vcf $TMP/chrX.vcf.gz \
           --split-par b38 \
           --check-sex 0.20 0.80 \
           --threads 8 --out $QCDIR/sexcheck
  ' | awk "/is submitted/ {gsub(/[<>]/, \"\", \$2); print \$2}")

echo "sexcheck job: $SEX_JOB"

# --- 4. Once both finish, build the keep-list ---
bsub \
  -J "sample_keep_list" \
  -w "done($KIN_JOB) && done($SEX_JOB)" \
  -q "$QUEUE" -P "$PROJECT" \
  -n 2 -R "rusage[mem=8000]" -W "01:00" \
  -o "$LSF_LOG_DIR/keep_list.%J.out" \
  -e "$LSF_LOG_DIR/keep_list.%J.err" \
  bash -lc "cd $(pwd) && CONFIG='$CONFIG' OUTPUT_ROOT='$OUTPUT_ROOT' \
    python -m src.qc.sample_qc --config '$CONFIG' --out-dir '$QC_DIR'"

echo "Chained keep-list job. When done: $QC_DIR/sample_keep_list.tsv + sample_qc_summary.md"
