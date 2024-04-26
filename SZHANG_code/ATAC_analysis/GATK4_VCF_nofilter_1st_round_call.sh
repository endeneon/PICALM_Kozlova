#!/bin/bash

# 26 Sept 2017 Siwei
# modified 26 Oct 2017 Siwei
# modified 1 Nov 2017 Siwei
# modified 9 Nov 2017 Siwei
# modified 24 Jan 2020 Siwei
# modified 05 May 2021 Siwei
# updated to GATK version 4

# modified 26 May 2021 Siwei
# updated to dbsnp version 154

gatk="/home/zhangs3/Data/Tools/gatk-4.1.8.1/gatk"
ref_genome="/home/zhangs3/Data/Databases/Genomes/hg38/INDEX/Homo_sapiens_assembly38.fasta"
dbsnp="/home/zhangs3/Data/Databases/Genomes/hg38/dbsnp_154.hg38.vcf"
vcf_suffix='_unfiltered.vcf'

# check if output directory exists
if [ ! -d "unfiltered_vcf" ]; then
        mkdir unfiltered_vcf
fi

date > log.txt

for eachfile in *.bam
do
	echo $eachfile
	echo $eachfile >> log.txt
	echo $date >> log.txt
#############call variants##############
	
#	samtools index -@ 23 $EACHFILE

	$gatk --java-options "-Xmx200g" \
		HaplotypeCaller \
		--native-pair-hmm-threads 20 \
		-R $ref_genome \
		--dbsnp $dbsnp \
		-I $eachfile \
		-O unfiltered_vcf/${eachfile/%bam/unfiltered.vcf}

#	~/2TB/Tools/oracle_JDK/bin/java -Xms10G -Xmx60G -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
#		-T HaplotypeCaller \
#		-R ~/2TB/Databases/hg38/Homo_sapiens_assembly38.fasta \
#		--dbsnp ~/2TB/Databases/hg38/dbsnp_146.hg38.vcf \
#		-gt_mode DISCOVERY \
#		-nct 23 \
#		-stand_call_conf 30 \
#		-I $EACHFILE \
#		-o $EACHFILE$vcf_suffix
done

