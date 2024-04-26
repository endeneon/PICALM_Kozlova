#!/bin/bash

# Siwei 21 Jan 2019

for EACHFILE in *.fastq.gz
do
	echo $EACHFILE
	####align for unpaired reads
	bowtie2 -p 23 -X 2000 --mm --qc-filter --met 1 -t --sensitive -x ~/1TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta -U $EACHFILE --met-file align_metrics/${EACHFILE/%.fastq.gz/}_U.txt | samtools view -b -u -h -@ 23 | samtools sort -l 9 -@ 23 -m 2G -o merged_sorted.bam
	samtools index -@ 23 merged_sorted.bam
	ls -lh merged_sorted.bam

		####dedup
	~/1TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/1TB/Tools/picard291.jar MarkDuplicates INPUT=merged_sorted.bam OUTPUT=merged_sorted_dedup.bam METRICS_FILE=metrics.txt REMOVE_DUPLICATES=True ASSUME_SORTED=True
	samtools index -@ 23 merged_sorted_dedup.bam

	####add read group
	~/1TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/1TB/Tools/picard291.jar AddOrReplaceReadGroups \
		INPUT=merged_sorted_dedup.bam \
		OUTPUT=merged_sorted_dedup_RG.bam \
		RGID=${EACHFILE/%.fastq.gz/} \
		RGLB=lib1 \
		RGPL=illumina \
		RGPU=unit1 \
		RGSM=${EACHFILE/%.fastq.gz/}
	samtools index -@ 23 merged_sorted_dedup_RG.bam

	####Re-align reads around known indels
	~/1TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/1TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-R ~/1TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		-I merged_sorted_dedup_RG.bam \
		-known ~/1TB/Databases/hg38/INDEX/Mills_and_1000G_gold_standard.indels.hg38.vcf \
		-nt 23 \
		-o realignment_targets.list

	~/1TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/1TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ~/1TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		-I merged_sorted_dedup_RG.bam \
		-known ~/1TB/Databases/hg38/INDEX/Mills_and_1000G_gold_standard.indels.hg38.vcf \
		-targetIntervals realignment_targets.list \
		-o merged_sorted_dedup_RG_realigned.bam

	samtools index -@ 23 merged_sorted_dedup_RG_realigned.bam

	####Base score recalibration
	~/1TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/1TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-R ~/1TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		-I merged_sorted_dedup_RG_realigned.bam \
		-knownSites ~/1TB/Databases/hg38/INDEX/Mills_and_1000G_gold_standard.indels.hg38.vcf \
		-knownSites ~/1TB/Databases/hg38/INDEX/dbsnp_146.hg38.vcf \
		-o recal_data.table \
		-nct 23

	~/1TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/1TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T PrintReads \
		-R ~/1TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		-I merged_sorted_dedup_RG_realigned.bam \
		-BQSR recal_data.table \
		-o pre_bams/${EACHFILE/%.fastq.gz/}.bam \
		-nct 23
	samtools index pre_bams/${EACHFILE/%.fastq.gz/}.bam

done

