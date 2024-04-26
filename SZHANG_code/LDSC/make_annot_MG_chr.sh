#!/bin/bash
# Siwei 27 Apr 2021


# ldsc uses numpy, which optimizes threads on its own.
# if you wish to parallelize, set to 1

export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3



# parallel -j22 python ldsc.py \
#	--bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.{1} \
#	--l2 \
#	--ld-wind-cm 1 \
#	--out 1000G_EUR_P3_baseline_27Apr2021/1000G.EUR.QC.{1} \
#	::: {1..22}

 parallel -j20 python make_annot.py \
 	--bed-file bed_hg19_chr/MG_17_lines_to_small_peaks.bed \
	--bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.{1}.bim \
	--annot-file Harry_annot/MG/MG.{1}.annot.gz \
	::: {1..22}


#parallel -j50 python ldsc.py \
#	--l2 \
#	--bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.{1} \
#	--ld-wind-cm 1 \
#	--annot 1000G_EUR_P3_baseline_27Apr2021/1000G.EUR.QC.{1}.annot.gz \
#	--thin-annot \
#	--out 1000G_EUR_P3_baseline_27Apr2021/1000G.EUR.QC.{1} \
#	--print-snps hapmap3_snps/hm.{1}.snp \
#	::: {1..22}


unset OMP_NUM_THREADS
unset MKL_NUM_THREADS

