#!/bin/bash

# Siwei 24 Jun 2022

rm read_sum.txt

for eachfile in *.bam
do
	echo $eachfile
	samtools flagstat -@ 4 $eachfile \
		| grep '0 mapped' \
		| cut -d " " -f 1 \
		>> read_sum.txt
done
