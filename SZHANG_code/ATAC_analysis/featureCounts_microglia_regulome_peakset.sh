#!/bin/bash

# Siwei 01 Jun 2023
# Count BAMs use the unified peakset used in Kosoy et al., Nat Genet 2022
# use an SAF file

featureCounts_path="/home/zhangs3/Data/Tools/subread-2.0.3-Linux-x86_64/bin/featureCounts"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
peak_saf="peak_file/OCRs_microglia_regulome.saf"

#for eachfile in *.bam
#do
#	echo $eachfile
#	samtools index -@ 20 $eachfile
#done



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
	-a $peak_saf \
	-F SAF \
	-p \
	-B \
	-O \
	--countReadPairs \
	--primary \
	-P \
	-D 2000 \
	-G $ref_genome \
	-o "all_MGs_Asts_OCR_microglia_regulome_count_matrix_featureCounts_01Jun2023.txt" \
	-T 32 \
	*.bam

