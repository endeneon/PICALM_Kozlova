#!/bin/bash


# Siwei 14 Feb 2018
# Siwei 22 Dec 2019
# Siwei 04 Jul 2023

# Genotype CLU

output_file_name=$1"_rs1532278_CLU_genotyping.txt"

date > $output_file_name

echo "All genomic coordinates are based on GRCh38p7." >> $output_file_name

for EACHFILE in *.bam
do
	echo $EACHFILE >> $output_file_name
	echo $EACHFILE
	printf "rs1532278\t" >> $output_file_name
	samtools mpileup -r chr8:27608798-27608798 $EACHFILE >> $output_file_name

done


