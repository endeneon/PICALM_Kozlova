#! /bin/bash

# Siwei 18 May 2021
# subsample all WASPed bam files to 30M reads for MACS2 peak calling

for eachfile in *.bam
do
	echo $eachfile
	cat <(samtools view -H $eachfile) <(samtools view -@ 40 $eachfile | shuf -n 30000000) \
	       | samtools sort -l 9 -@ 40 -m 10G -o subsampled/${eachfile/%.bam/_30M.bam}
done


