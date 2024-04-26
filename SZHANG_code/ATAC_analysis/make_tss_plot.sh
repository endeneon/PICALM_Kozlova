#!/bin/bash

# Siwei 30 Jun 2021

# use conda env aligners
# make sure R has been installed

# flow
# convert each bam (bamCoverage) bigwig 
# bigwig (ComputeMatrix reference-point) .gz (for plotHeatmap)
# .gz (plotHeatmap) output/pdf

for eachfile in *.bam
do
	echo $eachfile
	## convert bam to bigwig
	bamCoverage \
		-b $eachfile \
		-o tss_plot/bigwig_4_plotting.bigwig \
		-of bigwig \
		-bs 5 \
		-bl ~/Data/Databases/Genomes/hg38/ENCFF356LFX.bed \
		-p 40 \
		--effectiveGenomeSize 2913022398 \
		--normalizeUsing RPGC \
		-ignore chrX chrM

	## compute bigwig to .gz matrix
	computeMatrix reference-point \
		-S tss_plot/bigwig_4_plotting.bigwig \
		-R /home/zhangs3/Data/Databases/Genomes/hg38/gencode.v35.annotation.gtf \
		-a 1000 \
		-b 1000 \
		-o tss_plot/matrix_2_plot.gz \
		--outFileSortedRegions tss_plot/gencode.v35.tss.deeptools.bed \
		--referencePoint TSS \
		-bs 5 \
		--sortRegions keep \
		-bl ~/Data/Databases/Genomes/hg38/ENCFF356LFX.bed \
		-p 40 \
		--metagene
	
	## use the generated .gz matrix to plot heatmap
	plotHeatmap \
		-m tss_plot/matrix_2_plot.gz \
		-o tss_plot/${eachfile/%_new_WASPed\.bam}.pdf \
		--plotFileFormat pdf \
		--plotTitle ${eachfile/%_new_WASPed\.bam} \
		--heatmapHeight 8 \
		--heatmapWidth 4

done

