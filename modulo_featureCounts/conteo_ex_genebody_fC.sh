#!/usr/bin/env bash
set -euo pipefail

# ================== Configuración ==================

# ===== Entradas =====
# BAMs de entrada
BAMS=(
    /mnt/d/work_dir/memoria_trials/data/KO/control/CTRL{1..4}.sorted.bam
    /mnt/d/work_dir/memoria_trials/data/KO/RNU5A/A{1..4}.sorted.bam
    /mnt/d/work_dir/memoria_trials/data/KO/RNU5B/B{1..4}.sorted.bam
    /mnt/d/work_dir/memoria_trials/data/KO/RNU5D/D{1..4}.sorted.bam
    /mnt/d/work_dir/memoria_trials/data/KO/RNU5E/E{1..4}.sorted.bam
    /mnt/d/work_dir/memoria_trials/data/KO/RNU5F/F{1..4}.sorted.bam
)

# Archivos SAF
EXON_SAF="/mnt/d/work_dir/memoria/obj2/anotacion/table_genebodies_1_fixed_overlap/exons_extended_final.saf"
GENEBODY_SAF="/mnt/d/work_dir/memoria/obj2/anotacion/table_genebodies_1_fixed_overlap/genebodies_final.saf"

# Salidas
OUTDIR="/mnt/d/work_dir/memoria/obj2/featureCounts"
LOGDIR="/mnt/d/work_dir/memoria/obj2/featureCounts/logs"
mkdir -p "$OUTDIR" "$LOGDIR"

# ===== Parámetros =====
THREADS=8
STRAND=0              # 0 = unstranded, 1 = forward, 2 = reverse
MIN_MAPQ=10           # filtra mapeos débiles

# Algunas opciones
COMMON_OPTS=(
  -T "$THREADS"
  -s "$STRAND"
  -p
  --countReadPairs
  -B
  -C
  -Q "$MIN_MAPQ"
  -F SAF
)

# ================== Conteo EXÓNICO ==================
echo "[INFO] Contando EXONES..."
featureCounts \
  "${COMMON_OPTS[@]}" \
  -a "$EXON_SAF" \
  -o "$OUTDIR/exon_counts.txt" \
  "${BAMS[@]}" \
  > "$LOGDIR/log_exonic.txt" 2>&1

# ================== Conteo GENE BODY ==================
echo "[INFO] Contando GENE BODY..."
featureCounts \
  "${COMMON_OPTS[@]}" \
  -a "$GENEBODY_SAF" \
  -o "$OUTDIR/genebody_counts.txt" \
  "${BAMS[@]}" \
  > "$LOGDIR/log_genebody.txt" 2>&1

echo "[OK] Listo:"
echo "     $OUTDIR/exon_counts.txt (+ .summary)"
echo "     $OUTDIR/genebody_counts.txt (+ .summary)"
echo "[OK] Logs:"
echo "     $LOGDIR/log_exonic.txt"
echo "     $LOGDIR/log_genebody.txt"