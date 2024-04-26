#! /bin/bash

# 07 May 2021 Siwei

## use conda aligners environment

## Set these environment vars to point to
## your local installation

gatk37_path="/home/zhangs3/Data/Tools/GATK37"
gatk36_path="/home/zhangs3/Data/Tools/gatk_36/opt/gatk-3.6"
gatk4_path="/home/zhangs3/Data/Tools/gatk-4.1.8.1"

jre_8="/home/zhangs3/Data/Tools/jre1.8.0_291/bin"
openjdk_16="/home/zhangs3/Data/Tools/jdk-16.0.1/bin/java"

hg38_index="/home/zhangs3/Data/Databases/Genomes/hg38/INDEX/Homo_sapiens_assembly38.fasta"

## list all .vcf files
#shopt -s nullglob
#vcf_array=(*.vcf)
#echo "${vcf_array[@]}"

## send all .vcf files to a .list file
ls *.vcf > vcf_list.list

## use picard to merge all vcfs
## note that picard is located under the gatk4 folder
#$openjdk_16 -Xmx200G -jar $gatk4_path/picard2_25_4.jar MergeVcfs \
#	-I vcf_list.list \
#	-O raw_merged_vcfs.vcf

## !! Note!! Use GATK3.6 and jre 1.8 CombineVariants to merge all VCFs
## The new Picard from GATK4 only accepts GVCF files therefore 
## cannot merge vcfs with different sample identities

$jre_8/java -Xmx200G -jar $gatk36_path/GenomeAnalysisTK.jar \
	-T CombineVariants \
	-R $hg38_index \
	--variant vcf_list.list \
	-genotypeMergeOptions UNIQUIFY \
	-o output/raw_merged_vcfs.vcf



## select variants that only includes SNP
$gatk4_path/gatk SelectVariants \
	-R $hg38_index \
	-V output/raw_merged_vcfs.vcf \
	--select-type-to-include SNP \
	-O output/DN_20_lines_4_WASP.vcf


