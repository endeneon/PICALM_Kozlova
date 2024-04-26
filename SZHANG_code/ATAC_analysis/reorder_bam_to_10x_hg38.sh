#!/bin/bash

# Siwei 04 Jul 2023
# reorder all previous bams that are mapped using 1000G hg38 fasta
# to the chr order of 10x hg38 fasta

picard_jar="/home/zhangs3/Data/Tools/gatk-4.2.6.1/picard.jar"
GRCh38_10x_ref="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
temp_working_folder="/home/zhangs3/NVME/package_temp/reorder_DN"

mkdir -p temp_working_folder

for eachfile in *.bam
do
	echo $eachfile
	java -jar \
		$picard_jar \
		ReorderSam \
		-I $eachfile \
		-O ${eachfile/%_WASPed.bam/_10x_WASPed.bam} \
		-SD $GRCh38_10x_ref \
		-R $GRCh38_10x_ref \
		-U true \
		-S true \
		--COMPRESSION_LEVEL 9 \
		--CREATE_INDEX true \
		--MAX_RECORDS_IN_RAM 5000000 \
		--TMP_DIR $temp_working_folder
done

rm -r $temp_working_folder
