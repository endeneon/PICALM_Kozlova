#!/bin/bash

# Siwei 06 Jun 2023

# Siwei 07 Jun 2022



mkdir -p annot_file

for eachfile in output_bed4/*.bed
do
	base_bed_name=$(basename -- $eachfile)
	base_bed_name=${base_bed_name/%.bed/}
	echo $base_bed_name

	parallel -j 23 \
		python ..//make_annot.py \
			--bed-file $eachfile \
			--bimfile ../1000G_EUR_Phase3_plink/1000G.EUR.QC.{1}.bim\
			--annot-file annot_file/$base_bed_name.{1}.annot.gz \
			::: {1..22}
done
