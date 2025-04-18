#!/bin/bash

# Siwei 3 Apr 2018

# Siwei 13 Sept 2018
# Siwei 12 Dec 2019
# Siwei 04 May 2021

jre_8="/home/zhangs3/Data/Tools/jre1.8.0_291/bin/java"

mkdir -p pre_bams

for EACHFILE_1P in *_1.fq.gz
do
	####init
	date
	rm *.bam
	rm *.list
	rm *.table
	rm *.bai

	####variables
        echo $EACHFILE_1P
	echo ${EACHFILE_1P/%_1.fq.gz/_2.fq.gz}  #variable for 2P

	####alignments


	####align for paired reads
	date >> stat.txt
#	echo "Paired"
	echo ${EACHFILE_1P/%_1.fq.gz/} >> stat.txt
	(bowtie2 \
		-p 24 -X 2000 \
		--mm --qc-filter \
		--met 1 -t \
		--sensitive \
		--no-mixed \
		--no-discordant \
		-x "/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa" \
		-1 $EACHFILE_1P \
		-2 ${EACHFILE_1P/%_1.fq.gz/_2.fq.gz} \
		| samtools sort \
			-l 9 -m 10G -@ 20 \
			-o merged_sorted.bam) \
		2>>stat.txt
	ls -lh
	
	cat stat.txt

	####merge reads and sort
#	samtools merge -l 0 -@ 23 -f merged_unsorted.bam Paired.bam Unpaired.bam
#	samtools sort -l 9 -m 2G -@ 23 -o merged_sorted.bam merged_unsorted.bam
	samtools index -@ 20 merged_sorted.bam
	ls -lh

	####dedup
	$jre_8 -Xmx200g -jar \
		/home/zhangs3/Data/Tools/picard291.jar \
		MarkDuplicates \
			INPUT=merged_sorted.bam \
			OUTPUT=merged_sorted_dedup.bam \
			METRICS_FILE=metrics.txt \
			REMOVE_DUPLICATES=True \
			ASSUME_SORTED=True
	samtools index -@ 20 merged_sorted_dedup.bam

	####add read group
	$jre_8 -Xmx200g -jar /home/zhangs3/Data/Tools/picard291.jar \
		AddOrReplaceReadGroups \
			INPUT=merged_sorted_dedup.bam \
			OUTPUT=pre_bams/${EACHFILE_1P/%_1.fq.gz/}.bam \
			RGID=${EACHFILE_1P/%_1.fq.gz/} \
			RGLB=lib1 \
			RGPL=illumina \
			RGPU=unit1 \
			RGSM=${EACHFILE_1P/%_1.fq.gz/}

		samtools index -@ 20 pre_bams/${EACHFILE_1P/%_1.fq.gz/}.bam

done
