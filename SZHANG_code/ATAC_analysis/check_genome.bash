#!/bin/bash

# Siwei 21 Jan 2019
# Check the genome version of the BAMS 
# by looking at the @PG tag

for eachfile in *.bam
do
	echo ${eachfile/\n/\t}
	samtools view -H $eachfile | grep "^@PG" | grep "38" | grep "truncated"
done

