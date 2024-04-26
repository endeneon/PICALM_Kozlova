#!/bin/bash

# Siwei 09 Dec 2022
# Count BAMs at exon level

featureCounts_path="/home/zhangs3/Data/Tools/subread-2.0.3-Linux-x86_64/bin/featureCounts"
star_ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
# gencode_v35_gtf="/home/zhangs3/Data/Databases/Genomes/hg38/gencode.v35.annotation.sorted.gtf"
saf_table="GSE143666_sorted_MG_Ast_unified_peaks_to_count.saf"

mkdir -p output

#$featureCounts_path \
#	-p \
#	-O \
#	-a $gencode_v35_gtf \
#	-t exon \
#	-g exon_id \
#	-o VPS45_3xKD_gencode_v35_by_exon_09Dec2022.txt \
#	-C \
#	-T 32 \
#	*.bam

$featureCounts_path \
	-F SAF \
	-p \
	-B \
	-O \
	--countReadPairs \
	-P \
	-D 2000 \
	-G $star_ref_genome \
	-J \
	-a $saf_table \
	-o output/all_Ast_MG_use_GSE143666_interval_09Jun2023.txt \
	-T 16 \
	*.bam

