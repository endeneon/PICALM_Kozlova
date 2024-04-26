#! /bin/bash

# Siwei 19 May 2021

mkdir merged

## list all .bam files
shopt -s nullglob
bam_array=(*.bam)
echo "${bam_array[@]}"

date
echo "merging..."
samtools merge \
	-l 9 -@ 40 -f \
	-h ../DN-19_S11_WASPed.bam \
	merged/merged_unsorted.bam \
	${bam_array[@]}

# samtools index -@ 40 merged/merged_unsorted.bam

date
echo "sorting..."
samtools sort -l 9 -m 5G -@ 40 \
	-o merged/DN_20_lines_30M_merged_sorted_4_macs2.bam \
	merged/merged_unsorted.bam

samtools index -@ 40 merged/DN_20_lines_30M_merged_sorted_4_macs2.bam


