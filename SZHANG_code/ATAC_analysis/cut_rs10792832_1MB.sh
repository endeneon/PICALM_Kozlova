#!/bin/bash

# Siwei 05 Jul 2023
# cut up/downstream of rs10792832 (chr11:86156833) of 1 MB each

output_dir="rs10792832_1MB"

mkdir -p $output_dir

for eachfile in *.bam
do
	echo $eachfile
	samtools view \
		-h \
		-@ 20 \
		$eachfile \
		chr11:85156833-87156833 \
		| samtools sort \
		-l 9 \
		-m 5G \
		-@ 20 \
		-o $output_dir/${eachfile/%.bam/_rs10792832_1MB.bam}
done


