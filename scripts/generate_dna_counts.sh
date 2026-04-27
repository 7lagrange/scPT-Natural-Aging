#!/bin/bash
set -euo pipefail

conda activate scDNA

tissues=(
  "skin" "brainFC" "brainHip" "brainCB"
  "kidney" "lung" "heart"
  "BAT" "colon" "muscle"
  "liver" "testis" "spleen"
)


histones=(
  "H3K27me3" "H3K36me3" "H3K9me3"
  "H3K4me1" "H3K27ac" "H3K4me3"
)

BASE_DIR="/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag"
DATA_DIR="$BASE_DIR/data/merge"
PIPELINE_DIR="/storage/zhangyanxiaoLab/qihongjian/github/Paired-Tag-Pipeline-modified/pipelines"
BED_DIR="/storage/zhangyanxiaoLab/qihongjian/datasets/genomic_regions/mm10_bed"
OUT_DIR="/storage/zhangyanxiaoLab/qihongjian/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/DNA/counts/data"
THREADS=64

get_bin_size() {
  local histone="$1"

  case "$histone" in
    H3K27me3|H3K36me3|H3K9me3)
      echo "10kb"
      ;;
    H3K4me1|H3K27ac|H3K4me3)
      echo "1kb"
      ;;
    *)
      echo "[ERROR] Unsupported histone mark: $histone" >&2
      return 1
      ;;
  esac
}

for histone in "${histones[@]}"; do
  bin_size="$(get_bin_size "$histone")"
  bed_file="$BED_DIR/mm10_${bin_size}_bins.bed"
  output_prefix="${histone}_${bin_size}_bins"

  if [[ ! -f "$bed_file" ]]; then
    echo "[WARN] BED file not found for $histone: $bed_file, skipping this histone."
    continue
  fi

  for tissue in "${tissues[@]}"; do
    echo "======================================="
    echo "[INFO] Processing $tissue - $histone"
    echo "======================================="

    mkdir -p "$OUT_DIR"
    cd "$OUT_DIR"

    merged_bam="${tissue}-${histone}-merged.bam"
    matrix_name="${tissue}-${histone}-${output_prefix}"

    echo "[INFO] Merging BAM..."
    files=("$DATA_DIR"/*"${tissue}"*"${histone}"*/dna_bam/histone_dna_bam/*"${histone}"*.bam)

    n=${#files[@]}
    if [ $n -eq 0 ]; then
        echo "[WARN] No BAM files found for $tissue $histone, skipping..."
        continue
    elif [ $n -eq 1 ]; then
        echo "[INFO] Only one BAM file found, copying..."
        cp "${files[@]}" "$merged_bam"
        samtools index -@ $THREADS "$merged_bam"
    else
        echo "[INFO] Merging $n BAM files..."
        samtools merge -f -@ $THREADS "$merged_bam" "${files[@]}"
        samtools index -@ $THREADS "$merged_bam"
    fi


    bash "$PIPELINE_DIR/bam2mtx.sh" \
      -n "$matrix_name" \
      -g "$bed_file" \
      -i "$merged_bam" \
      -l DNA \
      -o "$matrix_name" \
      --thread "$THREADS"

    # counts in K9me3 common peaks
    bash "$PIPELINE_DIR/bam2mtx.sh" \
      -n "${tissue}-${histone}-H3K9me3_zhuojie" \
      -g /storage/zhangyanxiaoLab/qihongjian/github/zhanglab-code/projects/paired_seq_tag/code/zhuojie/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.bed  \
      -i "$merged_bam" \
      -l DNA \
      -o "${tissue}-${histone}-H3K9me3_zhuojie" \
      --thread $THREADS
    
    rm -f "$merged_bam" "${merged_bam}.bai" "${matrix_name}/CountMatrix.txt"
    echo "[INFO] Done: $tissue - $histone"
  done
done

echo "[INFO] ALL DONE"
