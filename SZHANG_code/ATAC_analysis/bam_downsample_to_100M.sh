#!/bin/bash

# Siwei 04 Jul 2023
# merge and downsample rs10792832 het bams to 100M reads for peak plotting

cell_type="DN"
temp_folder="/home/zhangs3/NVME/package_temp/downsample_100M_"$cell_type

if [[ -d $temp_folder ]]
then
	rm -r $temp_folder
fi

mkdir -p $temp_folder
mkdir -p output

echo $temp_folder


samtools merge \
	-o $temp_folder"/merged_unsorted.bam" \
	-s 42 \
	-l 5 \
	-@ 20 \
	*.bam

samtools sort \
	-m 5G \
	-l 5 \
	-@ 20 \
	-o $temp_folder"/merged_sorted.bam" \
	-T $temp_folder"/temp" \
	$temp_folder"/merged_unsorted.bam"

samtools index \
	-@ 30 \
	$temp_folder/merged_sorted.bam

reformat.sh \
	in=$temp_folder"/merged_sorted.bam" \
	out="output/"$cell_type"_downsampled_100M_het_rs10792832.bam" \
	sampleseed=42 \
	samplereadstarget=100000000

rm -r $temp_folder
