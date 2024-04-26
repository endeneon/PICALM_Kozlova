#!/bin/sh

# Siwei 5 Jun 2018
# Siwei 31 Oct 2018


for EACHFILE in *.bam
do
	echo $EACHFILE
	echo ${EACHFILE/%_WASPed.bam/}
	macs2 callpeak -t $EACHFILE \
		-n ${EACHFILE/%_WASPed.bam/} \
		--outdir MACS2_output/ \
		-f BAMPE \
		-g 3.2e9 \
		--nomodel \
		--nolambda \
		--keep-dup all \
		--call-summits \
		--verbose 3
done

