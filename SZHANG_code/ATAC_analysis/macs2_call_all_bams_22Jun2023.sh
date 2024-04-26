#!/bin/bash

# Siwei 10 May 2022
# use the improved method as in 
# https://www.biorxiv.org/content/10.1101/496521v1

mkdir -p ldsc_peakset

for eachfile in *.bam
do
        echo $eachfile
        samtools index -@ 20 $eachfile
done


macs2 callpeak \
	-t *.bam \
	-f BAMPE \
	-g 2.7e9 \
	-q 0.05 \
	--keep-dup all \
	--nolambda \
	--min-length 100 \
	--max-gap 50 \
	--buffer-size 1000000 \
	-n ldsc_peakset/GABA_20_peak_set_4_LDSC_21Jun2023 \
	--seed 42
