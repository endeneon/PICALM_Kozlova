#!/bin/bash

for eachfile in *.bam
do
	echo $eachfile

	samtools index -@ 20 $eachfile
	samtools view \
		-@ 10 \
		$eachfile \
		| wc -l
done

