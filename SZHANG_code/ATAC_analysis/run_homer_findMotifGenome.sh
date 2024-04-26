#! /bin/bash
# Siwei 24 Apr 2024

for eachfile in *.bed
do
	echo $eachfile

	findMotifsGenome.pl \
		$eachfile \
		hg38 \
		${eachfile/%.bed/} \
		-size 100 \
		-p 16

done

