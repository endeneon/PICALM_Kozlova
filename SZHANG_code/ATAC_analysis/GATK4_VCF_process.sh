#!/bin/bash

# 29 May 2021
# Siwei update dnsnp version to v154

# 11 May 2021
# Siwei rewrite in GATK4

# 11 Jun 2020

vcf_suffix="_995.vcf"
gatk4="/home/zhangs3/Data/Tools/gatk-4.1.8.1/gatk"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/hg38/INDEX/Homo_sapiens_assembly38.fasta"
dbsnp="/home/zhangs3/Data/Databases/Genomes/hg38/dbsnp_154.hg38.vcf"

# check if output directory exists
if [ ! -d "vcf_output" ]; then
        mkdir vcf_output
fi


date > log.txt

rm *.vcf
rm *.vcf.idx
rm *.recal
rm *.tranches
rm *.R


for EACHFILE in *.bam
do
	rm *.vcf
	rm *.vcf.idx
	rm *.recal
	rm *.tranches
	rm *.R


	echo $EACHFILE
	echo $EACHFILE >> log.txt
	date >> log.txt
############## call variants ##############
	
	samtools index -@ 20 $EACHFILE
	
	$gatk4 --java-options "-Xmx150g" \
		HaplotypeCaller \
		--native-pair-hmm-threads 20 \
		-R $ref_genome \
		--dbsnp $dbsnp \
		-I $EACHFILE \
		-O raw_variants.vcf

############variants filtering##########
############SNP#########################
############Remove previous calibration files######

############ SNP Recalibration ###########

        $gatk4 --java-options "-Xmx150g" \
	 	VariantRecalibrator \
		-R $ref_genome \
		-V raw_variants.vcf \
		--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $ref_path/hapmap_3.3.hg38.vcf \
		--resource:omni,known=false,training=true,truth=true,prior=12.0 $ref_path/1000G_omni2.5.hg38.vcf \
		--resource:1000G,known=false,training=true,truth=false,prior=10.0 $ref_path/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
		--max-gaussians 4 \
		-an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum \
		-mode SNP \
		-AS \
		-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 90.0 \
		-O raw_recalibrate_SNP.recal \
		--tranches-file raw_recalibrate_SNP.tranches \
		--rscript-file raw_recalibrate_SNP_plots.R  

########### apply SNP recalib###########

        $gatk4 --java-options "-Xmx150g" \
		ApplyVQSR \
		-R $ref_genome \
		-V raw_variants.vcf \
		-mode SNP \
		-AS \
		--truth-sensitivity-filter-level 99.5 \
		--recal-file raw_recalibrate_SNP.recal \
		--tranches-file raw_recalibrate_SNP.tranches \
		-O variants_calibrated.vcf

########## filter SNPs located in the hg38 blacklisted region #######

	$gatk4 --java-options "-Xmx150g" \
		VariantFiltration \
		-R $ref_genome \
		-V variants_calibrated.vcf \
		--mask $ref_path/hg38_blacklisted_regions.bed --mask-name "blacklisted_regions" \
		-O variants_calibrated_blacklisted.vcf

########## select SNP for output from autosomes only #######

	$gatk4 --java-options "-Xmx150g" \
		SelectVariants \
		-R $ref_genome \
		-V variants_calibrated_blacklisted.vcf \
		--select-type-to-include SNP \
		--restrict-alleles-to BIALLELIC \
		-L $ref_path/autosomes.list \
		-O vcf_output/${EACHFILE/%.bam/}$vcf_suffix

done

