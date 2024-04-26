#!/bin/bash

for eachfile in R21*.bam
do
	echo $eachfile
	samtools view \
		-H \
		$eachfile \
		| grep "RG"
done


