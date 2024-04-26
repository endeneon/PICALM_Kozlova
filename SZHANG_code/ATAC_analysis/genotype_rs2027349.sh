#!/bin/bash


# Siwei 14 Feb 2018
# Siwei 22 Dec 2019
# Siwei 04 Jul 2023

# Genotype PICALM rs10792832

output_file_name=$1"_rs2027349_VPS45_genotyping.txt"

date > $output_file_name

echo "All genomic coordinates are based on GRCh38p7." >> $output_file_name

for EACHFILE in *.bam
do
	echo $EACHFILE >> $output_file_name
	echo $EACHFILE
	printf "rs2027349\t" >> $output_file_name
	samtools mpileup -r "chr1:150067621-150067621" $EACHFILE \
		| grep -i "ag" \
		>> $output_file_name

done


