#!/bin/bash


#Siwei 14 Feb 2018
# Siwei 22 Dec 2019


# Genotype RERE rs301791

date > RERE_genotyping.txt

echo "All genomic coordinates are based on GRCh38p7." >> RERE_genotyping.txt

for EACHFILE in *.bam
do
	echo $EACHFILE >> RERE_genotyping.txt
	echo $EACHFILE
	printf "rs301791\t" >> RERE_genotyping.txt
	samtools mpileup -r chr1:8408312-8408312 $EACHFILE >> RERE_genotyping.txt

done


