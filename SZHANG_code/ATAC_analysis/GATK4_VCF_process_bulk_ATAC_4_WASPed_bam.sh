#!/bin/bash

# 03 May 2023
# Siwei

# 11 May 2021
# Siwei rewrite in GATK4

# 11 Jun 2020

vcf_suffix="_995.vcf"
gatk4="/home/zhangs3/Data/Tools/gatk-4.2.6.1/gatk"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
dbsnp="/home/zhangs3/Data/Databases/Genomes/hg38/dbsnp_154.hg38.vcf"

mkdir -p vcf_output 

date > log.log

rm *.vcf
rm *.vcf.idx
rm *.recal
rm *.tranches
rm *.R


for EACHFILE in *.bam
do

#	rm *.vcf
#	rm *.vcf.idx
#	rm *.recal
#	rm *.tranches
#	rm *.R


	echo $EACHFILE
	echo $EACHFILE >> log.log
	date >> log.log
############## call variants ##############
	
	samtools index -@ 20 $EACHFILE
	
	$gatk4 --java-options "-Xmx150g" \
		HaplotypeCaller \
		--native-pair-hmm-threads 15 \
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
		--filter-expression "AC > 0" --filter-name "PASS" \
		-O variants_calibrated_blacklisted.vcf

########## select SNP for output from autosomes only #######

	$gatk4 --java-options "-Xmx150g" \
		SelectVariants \
		-R $ref_genome \
		-V variants_calibrated_blacklisted.vcf \
		--select-type-to-include SNP \
		--restrict-alleles-to BIALLELIC \
		-L $ref_path/autosomes.list \
		--exclude-filtered true \
		-O pre_output.vcf

######## extract het only

        cat pre_output.vcf \
                | grep "^#" \
                > header.txt

        cat pre_output.vcf \
                | grep -v "^#" \
                | grep "0/1" \
                > body.txt

        cat header.txt body.txt \
                > vcf_output/${EACHFILE/%.bam/}$vcf_suffix

###### cleanup

	rm *.vcf
        rm *.idx
        rm *.recal
        rm *.tranches
        rm *.R
	rm *.txt
	rm *.pdf
	rm *.idx
	rm *.log
done

