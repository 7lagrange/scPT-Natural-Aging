#!/bin/bash
set -euo pipefail

conda activate sicer_env

# ---------- 1. 参数 ----------
tissues=(
  "brainFC" "brainHip" "brainCB"
  "kidney" "lung" "heart"
  "BAT" "colon" "muscle"
  "liver" "testis" "spleen"
)

histones=(
  "H3K27me3" "H3K9me3"
  "H3K4me1" "H3K27ac" 
  "H3K4me3"  "H3K36me3"
  
)

tissues=(
  "BAT"  "muscle" "brainFC" "brainHip" "brainCB"
  "kidney" "lung" "heart" "colon" 
  "liver" "testis" "spleen"
)

histones=("H3K9me3"
  
)



BASE_DIR="/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag"
DATA_DIR="$BASE_DIR/data/merge"

THREADS=64

# ---------- 2. 主循环 ----------
for tissue in "${tissues[@]}"; do
    for histone in "${histones[@]}"; do
    
    echo "======================================="
    echo "[INFO] Processing $tissue - $histone"
    echo "======================================="

    OUT_DIR="/storage/zhangyanxiaoLab/qihongjian/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/DNA/peaks/data"
    mkdir -p "$OUT_DIR"
    cd "$OUT_DIR"

    MERGED_BAM="${tissue}-${histone}-merged.bam"

    # ---------- merge ----------
    echo "[INFO] Merging BAM..."

    files=($DATA_DIR/*${tissue}*${histone}*/dna_bam/histone_dna_bam/*${histone}*.bam)
    n=${#files[@]}
    
    if [ $n -eq 0 ]; then
        echo "[WARN] No BAM files found for $tissue $histone, skipping..."
        continue
    elif [ $n -eq 1 ]; then
        echo "[INFO] Only one BAM file found, copying..."
        cp "${files[@]}" "$MERGED_BAM"
        samtools index -@ $THREADS "$MERGED_BAM"
    else
        echo "[INFO] Merging $n BAM files..."
        samtools merge -@ $THREADS "$MERGED_BAM" "${files[@]}"
        samtools index -@ $THREADS "$MERGED_BAM"
    fi
    

    # ---------- barcode ----------
    echo "[INFO] Extracting barcodes..."
    
    awk -v h="$histone" 'NR==1 || $0 ~ h' "$BASE_DIR/results/${tissue}-merged/celltype_histone.txt" > celltype_histone.txt
    
    python ~/github/Paired-Tag-Pipeline-modified/scripts/extract_sample_barcode_from_bam.py \
        --bam "$MERGED_BAM" \
        --statH celltype_histone.txt \
        --outPrefix DNA

    rm celltype_histone.txt
    
    for BAM_FILE in DNA.*.bam; do
        echo "[INFO] Processing $BAM_FILE"

        samtools index -@ $THREADS "$BAM_FILE"

        mark=$(basename "$BAM_FILE" .bam | sed 's/^DNA\.//')
        BED_FILE="${tissue}.${mark}.bed"

        echo "[INFO] Mark = $mark"

        rm -rf peaks
        mkdir peaks

        # ---------- MACS3 ----------
        if [[ "$mark" == *"H3K4me1" || "$mark" == *"H3K4me3" || "$mark" == *"H3K27ac" ]]; then

            macs3 callpeak \
                --nomodel \
                --format BAM \
                --treatment "$BAM_FILE" \
                --outdir peaks \
                -g mm \
                -q 0.05 \
                --keep-dup all \
                -n "$mark"

            bedtools sort -i "peaks/${mark}_peaks.narrowPeak" > "$BED_FILE"
        fi

        # ---------- H3K27me3 ----------
        if [[ "$mark" == *"H3K27me3" ]]; then
            sicer \
                --treatment_file "$BAM_FILE" \
                -s mm10 \
                --gap_size 2000 \
                --window_size 200 \
                --output_directory peaks \
                -fdr 0.01

            FILE=$(ls peaks/*W200-G2000.scoreisland | head -n 1)
            bedtools sort -i "$FILE" > "$BED_FILE"
        fi

        # ---------- H3K36me3 ----------
        # https://doi.org/10.1016/j.celrep.2019.12.060
        # Peak calls were created with SICER (Zang et al., 2009) using redundancy threshold = 1, window size = 500, fragment size = 250, effective genome fraction = 0.79, gap size = 1500 and an FDR threshold of 0.01
        if [[ "$mark" == *"H3K36me3" ]]; then

            sicer \
                --treatment_file "$BAM_FILE" \
                -s mm10 \
                --window_size 500 \
                --gap_size 1500 \
                -fdr 0.01 \
                --output_directory peaks

            FILE=$(ls peaks/*W500-G1500.scoreisland | head -n 1)
            bedtools sort -i "$FILE" > "$BED_FILE"
        fi

        # ---------- H3K9me3 ----------
        # https://www.nature.com/articles/s41422-022-00719-6
        # H3K9me3 associated domains (peaks) were called using SICER93 on aggregated H3K9me3 signals from Paired-tag (without input). All default parameters were used, except that window size parameter was set to 5000 and gap size was set to 10000 to detect large peaks.
        if [[ "$mark" == *"H3K9me3" ]]; then

            sicer \
                --treatment_file "$BAM_FILE" \
                -s mm10 \
                --window_size 5000 \
                --gap_size 10000 \
                -fdr 0.01 \
                --output_directory peaks

            FILE=$(ls peaks/*W5000-G10000.scoreisland | head -n 1)
            bedtools sort -i "$FILE" > "$BED_FILE"
        fi

        echo "[INFO] Output BED: $BED_FILE"
    done
    rm *bam
    rm *bam.bai
    rm DNA*bed
    echo "[INFO] Done: $tissue - $histone"

  done
done

echo "[INFO] ALL DONE"