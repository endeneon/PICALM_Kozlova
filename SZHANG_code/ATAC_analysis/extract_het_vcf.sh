#!/bin/bash

# Siwei 13 May 2021
# extract vcf files generated from WASPed bams and HaplotypeCaller,
# retain GT of 0/1 (isHet == 1) only

# init
gatk4="/home/zhangs3/Data/Tools/gatk-4.1.8.1/gatk"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/hg38/INDEX/Homo_sapiens_assembly38.fasta"

rm *.txt

# check if het_vcf directory exists
if [ ! -d "het_vcf" ]; then
	mkdir het_vcf
fi


for eachfile in *.vcf
do
	echo $eachfile

	cat $eachfile | grep "^#" > head.txt
	cat $eachfile | grep -v "^#" | grep '0/1:' > body.txt
	cat head.txt body.txt > het_vcf/${eachfile/%.vcf/_het.vcf}


#        $gatk4 --java-options "-Xmx50g" \
#                VariantFiltration \
#                -R $ref_genome \
#                -V $eachfile \
#		-G-filter-name "het_filter" -G-filter "isHet == 1" \
#                -O het_vcf/${eachfile/%.vcf/_het.vcf}
done



