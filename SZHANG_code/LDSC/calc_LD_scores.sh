#!/bin/bash

# Siwei 23 Jun 2023
# Siwei 07 Jun 2022


mkdir -p ld_scores

shopt -s nullglob
sumstats_list=(/home/zhangs3/Data/Databases/GWAS/MAGMA_ref_new/ldsc_sumstats/*.sumstats.gz)
shopt -u nullglob

echo "${#sumstats_list[@]}"

for eachfile in output_bed4/*.bed
do
        base_bed_name=$(basename -- $eachfile)
        base_bed_name=${base_bed_name/%.bed/}
        echo $base_bed_name

	parallel -j 23 \
		python ../ldsc.py \
			--l2 \
			--bfile ../1000G_EUR_Phase3_plink/1000G.EUR.QC.{1} \
			--ld-wind-cm 1 \
			--annot annot_file/$base_bed_name.{1}.annot.gz \
			--thin-annot \
			--out annot_file/$base_bed_name.{1} \
			--print-snps ../list.txt \
			::: {1..22}
done

