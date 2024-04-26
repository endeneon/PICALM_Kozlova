#!/bin/bash

# Siwei 3 Apr 2018

# Siwei 13 Sept 2018
# Siwei 12 Dec 2019
# Siwei 04 May 2021

jre_8="/home/zhangs3/Data/Tools/jre1.8.0_291/bin/java"

# check if output directory exists
if [ ! -d "pre_bams" ]; then
        mkdir pre_bams
fi


for EACHFILE_1P in *_R1_001.fastq.gz
do
	####init
	date
	rm *.bam
	rm *.list
	rm *.table
	rm *.bai

	####variables
        echo $EACHFILE_1P
	echo ${EACHFILE_1P/%_R1_001.fastq.gz/_R2_001.fastq.gz}  #variable for 2P

	####alignments


	####align for paired reads
	date >> stat.txt
#	echo "Paired"
	echo ${EACHFILE_1P/%_R1_001.fastq.gz/} >> stat.txt
	(bowtie2 -p 30 -X 2000 --mm --qc-filter --met 1 -t --sensitive --no-mixed --no-discordant -x /home/zhangs3/Data/Databases/Genomes/hg38/INDEX/Homo_sapiens_assembly38.fasta -1 $EACHFILE_1P -2 ${EACHFILE_1P/%_R1_001.fastq.gz/_R2_001.fastq.gz} | samtools sort -l 9 -m 10G -@ 30 -o merged_sorted.bam) 2>>stat.txt
	ls -lh
	
	cat stat.txt

	####merge reads and sort
#	samtools merge -l 0 -@ 23 -f merged_unsorted.bam Paired.bam Unpaired.bam
#	samtools sort -l 9 -m 2G -@ 23 -o merged_sorted.bam merged_unsorted.bam
	samtools index -@ 30 merged_sorted.bam
	ls -lh

	####dedup
	$jre_8 -Xmx200g -jar /home/zhangs3/Data/Tools/picard291.jar MarkDuplicates INPUT=merged_sorted.bam OUTPUT=merged_sorted_dedup.bam METRICS_FILE=metrics.txt REMOVE_DUPLICATES=True ASSUME_SORTED=True
	samtools index -@ 30 merged_sorted_dedup.bam

	####add read group
	$jre_8 -Xmx200g -jar /home/zhangs3/Data/Tools/picard291.jar AddOrReplaceReadGroups \
		INPUT=merged_sorted_dedup.bam \
		OUTPUT=pre_bams/${EACHFILE_1P/%_R1_001.fastq.gz/}.bam \
		RGID=${EACHFILE_1P/%_R1_001.fastq.gz/} \
		RGLB=lib1 \
		RGPL=illumina \
		RGPU=unit1 \
		RGSM=${EACHFILE_1P/%_R1_001.fastq.gz/}

		samtools index -@ 30 pre_bams/${EACHFILE_1P/%_R1_001.fastq.gz/}.bam

done
