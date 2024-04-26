#!/bin/bash

# Siwei 11 Aug 2023
# revised, first down-sample every sample to 20M to ensure all samples have the same weight

# Siwei 04 Jul 2023
# merge and downsample rs1532278 het bams to 100M reads for peak plotting

cell_type="NGN2"
temp_folder="/home/zhangs3/NVME/package_temp/downsample_100M_"$cell_type

ref_fasta="/home/zhangs3/Data/Databases/Genomes/hg38/INDEX/Homo_sapiens_assembly38.fasta"

if [[ -d $temp_folder ]]
then
        rm -r $temp_folder
fi

mkdir -p $temp_folder
mkdir -p output

samtools merge \
	-o $temp_folder"/merged_unsorted.bam" \
	-s 42 \
	-l 5 \
	-@ 20 \
	down_20M/*.bam

samtools sort \
	-m 5G \
	-l 5 \
	-@ 20 \
	-o $temp_folder/merged_sorted.bam \
	-T $temp_folder/temp \
	$temp_folder/merged_unsorted.bam

samtools index \
	-@ 30 \
	$temp_folder/merged_sorted.bam

reformat.sh \
	in=$temp_folder/merged_sorted.bam \
	out=$temp_folder/down_sampled.bam \
	sampleseed=42 \
	ref=$ref_fasta \
	overwrite=true \
	samplereadstarget=50000000

samtools index \
	-@ 30 \
	$temp_folder/down_sampled.bam

# cut +/- 1MB of rs1532278
samtools view \
	-b -h \
	-@ 20 \
	$temp_folder/down_sampled.bam \
	"chr8:26608798-28608798" \
	| samtools sort \
	-l 9 -m 5G \
	-@ 20 \
	-o "output/"$cell_type"_downsampled_50M_het_rs1532278_2MB.bam"

rm -r $temp_folder
