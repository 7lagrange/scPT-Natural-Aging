#!/bin/bash

# sh generate_bam.sh config/20240808-scPT-Natural-Aging-liver-H3K4me1-H3K4me3.sh /path/to/xx_pipeline /path/to/merge
 
CONFIG_FILE="$1"
PIPELINE_BASE="$2"
MERGE_BASE="$3"

source "$CONFIG_FILE"

cd "$MERGE_BASE/$proj"

################################################ RNA part ###########################
rm -r rna_bam
mkdir -p rna_bam
cd rna_bam

# filter cells dna: 300 rna: 300 
echo "RNA_final.bam" > merge_list_rna_bam.txt
for pair in "${dna_rna[@]}"
do
    echo $pair
    dna=$(echo "$pair" | cut -d "_" -f 1)
    rna=$(echo "$pair" | cut -d "_" -f 2)
    rna_bam="$PIPELINE_BASE/$rna/$depth/$rna/bam/${rna}_UsefulReads.bam"

    bc_counts="../qc/${dna}_${rna}_bc_counts.txt"
    awk -v dna_min_reads=$dna_min_reads -v rna_min_reads=$rna_min_reads 'NR == 1 || ($2 >= dna_min_reads && $3 >= rna_min_reads)' $bc_counts > filtered_bc.txt

    python ~/github/Paired-Tag-Pipeline-modified/scripts/filter_bam_based_on_bc.py --bam $rna_bam --statH filtered_bc.txt
    echo "${rna}_UsefulReads_filtered_bc.bam ${rna}" >> merge_list_rna_bam.txt
done


echo "start to merge and add lib for RNA"
perl ~/github/Paired-Tag/perlscripts/mergeBam.pl merge_list_rna_bam.txt

echo "start to generate RNA matrix"
bash ~/github/Paired-Tag-Pipeline-modified/pipelines/bam2mtx.sh -n merged_RNA -g mm10 -i RNA_final_sorted.bam -l RNA -o ./RNA_matrix 

bash ~/github/Paired-Tag-Pipeline-modified/pipelines/bam2mtx-TEsubfamily.sh -n merged_RNA -g mm10 -i RNA_final_sorted.bam -l RNA -o ./RNA_TE_matrix 


rm *_filtered_bc.bam

cd ..


###################### DNA filter and merge ##################################################
rm -r dna_bam
mkdir -p dna_bam
cd dna_bam

echo "DNA_merge.bam" > merge_list_dna_bam.txt
for pair in "${dna_rna[@]}"
do
    echo $pair
    dna=$(echo "$pair" | cut -d "_" -f 1)
    rna=$(echo "$pair" | cut -d "_" -f 2)

    dna_bam="$PIPELINE_BASE/$dna/$depth/$dna/bam/${dna}_UsefulReads.bam"
  
    # filter cells dna: 300 rna: 300 
    bc_counts="../qc/${dna}_${rna}_bc_counts.txt"
    awk -v dna_min_reads=$dna_min_reads -v rna_min_reads=$rna_min_reads 'NR == 1 || ($2 >= dna_min_reads && $3 >= rna_min_reads)' $bc_counts > filtered_bc.txt

    python ~/github/Paired-Tag-Pipeline-modified/scripts/filter_bam_based_on_bc.py --bam $dna_bam --statH filtered_bc.txt
    echo "${dna}_UsefulReads_filtered_bc.bam ${rna}" >> merge_list_dna_bam.txt

done

echo "start to merge and add lib for DNA"
perl ~/github/Paired-Tag/perlscripts/mergeBam.pl merge_list_dna_bam.txt
samtools index DNA_merge_sorted.bam

echo "filter high pileups for DNA"  
python ~/github/Paired-Tag/remove_pileup/count_pileups.py DNA_merge_sorted.bam DNA_merge_sorted.tally.txt
python ~/github/Paired-Tag/remove_pileup/remove_pileups.py DNA_merge_sorted.bam DNA_merge_sorted.tally.txt DNA_final_sorted.bam 100
samtools index DNA_final_sorted.bam

#rm DNA*merge*
rm *filter*


###################### split DNA bam into histone mark ############

mkdir -p histone_dna_bam
cd histone_dna_bam

echo "barcode sample" > barcode.txt
for line in "${barcode_samples[@]}"; do
    echo "$line" >> barcode.txt
done

python ~/github/Paired-Tag-Pipeline/scripts/extract_sample_barcode_from_bam.py --bam ../DNA_final_sorted.bam --statH barcode.txt  --outPrefix 'dna'


for bam_name in *bam; do
    echo "Processing $bam_name"
    
    samtools view -@ 16 "$bam_name" | awk '{lib=substr($1,length($1)-25,6); cb=substr($1,length($1)-18,8); print lib, cb, $1}' OFS='\t' > tmp.library_CB.txt
    sort --parallel 16 tmp.library_CB.txt | bedtools groupby -g 1,2 -c 3 -o count -i - > library_CB_counts.txt
    
    # 输出格式：文件名 平均值
    awk -v bam="$bam_name" '{total+=$3; count++} END {if(count>0) print bam, total/count}' library_CB_counts.txt >> bam.fragments.txt
    
    awk '{total[$1]+=$3; count[$1]++} END {for (lib in total) print lib, total[lib]/count[lib]}' library_CB_counts.txt  >> bam.library.fragments.txt
    
    rm -f tmp.library_CB.txt library_CB_counts.txt
done
