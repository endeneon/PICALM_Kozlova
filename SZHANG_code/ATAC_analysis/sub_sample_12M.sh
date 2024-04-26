#! /bin/bash

# Siwei 18 May 2021
# subsample all WASPed bam files to 30M reads for MACS2 peak calling

eachfile=$1

	echo $eachfile
	cat <(samtools view -H $eachfile) <(samtools view -@ 20 $eachfile | shuf -n 12000000) \
	       | samtools sort -l 9 -@ 20 -m 10G -o subsampled/${eachfile/%.bam/_12M.bam}


