#! /bin/bash
# Siwei 26 Mar 2024

# convert the BAM files used to bigwig for GEO submission

mkdir -p bigWig

for eachfile in *.bam
do
	echo $eachfile
	samtools index -@ 20 $eachfile

	bamCoverage \
		--bam $eachfile \
		-o "bigWig/"${eachfile/%.bam/.bw} \
		-of bigwig \
		-bs 50 \
		--normalizeUsing None \
		--effectiveGenomeSize 2913022398 \
		--skipNAs \
		--blackListFileName "/home/zhangs3/Data/Databases/Genomes/hg38/hg38_blacklisted_regions.bed" \
		-p 20 \
		-v
done

