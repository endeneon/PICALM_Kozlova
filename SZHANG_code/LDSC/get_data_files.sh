#!/bin/bash


## Siwei 14 Apr 2021
## modified from Rob's iMG analysis for Alena

############
# This script assumes you are aware of the tutorial on LDSC for cell specific analyses
# Requires that downloaded files from Broad are unpacked/unzipped in this same folder
# gene lists can be sourced elsewhere
############

# uncomment if you are running in batches
# source activate ldsc || { echo "LDSC Conda environment not activated"; exit; }

# LDSC ON miR137 SET FOR SCZ ### will test different SCZ w versions
cts_name=miR137_SCZ

# GENERATING GENE SPECIFICITY SETS FOR CELL TYPES
#Using same gene sets as 100kb window, modify command_${cts_name}.txt
#mkdir ${cts_name}_1000Gv3_ldscores
#cp Kozlova_top20_1000Gv3_ldscores/*.GeneSet ${cts_name}_1000Gv3_ldscores/
#for k in {1..15} control; do mv Kozlova_top20.${k}.GeneSet ${cts_name}.${k}.GeneSet; done
#Manually modify example Kozlova.ldcts to generate ${cts_name}.ldcts

## GENERATING LD PARTITIONS FOR GENE SETS
# collecting 1000G_phase3 LD baseline and weight files:
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
tar -xvzf weights_hm3_no_hla.tgz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz
tar -xvzf 1000G_Phase3_frq.tgz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
mkdir 1000G_EUR_Phase3_baseline_v2.2
tar -xvzf 1000G_Phase3_baselineLD_v2.2_ldscores.tgz -C 1000G_EUR_Phase3_baseline_v2.2/

# hapmap3 snp set
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 w_hm3.snplist.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz
tar -xvzf hapmap3_snps.tgz

# plink files for LD partitions
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz
tar -xvzf 1000G_Phase3_plinkfiles.tgz

# gene coordinates
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/make_annot_sample_files/ENSG_coord.txt

# generate annot files
# parallel -j100 make_annot.py \
#  --gene-set-file {1}_1000Gv3_ldscores/{1}.{2}.GeneSet \
#  --gene-coord-file ENSG_coord.txt \
#  --windowsize 20000 \
#  --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.{3}.bim \
#  --annot-file {1}_1000Gv3_ldscores/{1}.{2}.{3}.annot.gz \
#  ::: $cts_name ::: {1..15} control ::: {1..22}

#ldsc uses numpy, which optimizes threads on its own.
#if you wish to parallelize, set to 1
#export OMP_NUM_THREADS=1
#LD score computations
#parallel -j100 ldsc.py \
#  --l2 \
#  --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.{3} \
#  --ld-wind-cm 1 \
#  --annot {1}_1000Gv3_ldscores/{1}.{2}.{3}.annot.gz \
#  --thin-annot \
#  --out {1}_1000Gv3_ldscores/{1}.{2}.{3} \
#  --print-snps hapmap3_snps/hm.{3}.snp \
#  ::: $cts_name ::: {1..15} control ::: {1..22}
#unset OMP_NUM_THREADS

##### AD GWAS
## PREP GWAS DATA
# munge_sumstats.py \
#   --sumstats ../R_analysis_redux/AD_sumstats_Jansenetal.txt.gz \
#   --N-col Nsum \
#   --ignore BETA \
#   --merge-alleles w_hm3.snplist \
#   --out Jansen_AD

## RUN THE REGRESSION
#ldsc.py \
#  --h2-cts Jansen_AD.sumstats.gz \
#  --ref-ld-chr 1000G_EUR_Phase3_baseline_v2.2/baselineLD. \
#  --overlap-annot \
#  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
#  --out Jan_AD_${cts_name} \
#  --ref-ld-chr-cts ${cts_name}.ldcts \
#  --w-ld-chr weights_hm3_no_hla/weights.

##### BMI UKBB GWAS
## PREP GWAS DATA
# wget https://data.broadinstitute.org/alkesgroup/UKBB/body_BMIz.sumstats.gz
# munge_sumstats.py \
#   --sumstats body_BMIz.sumstats.gz \
#   --merge-alleles w_hm3.snplist \
#   --out BMI

## RUN THE REGRESSION
#ldsc.py \
#  --h2-cts BMI.sumstats.gz \
#  --ref-ld-chr 1000G_EUR_Phase3_baseline_v2.2/baselineLD. \
#  --overlap-annot \
#  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
#  --out BMI_${cts_name} \
#  --ref-ld-chr-cts ${cts_name}.ldcts \
#  --w-ld-chr weights_hm3_no_hla/weights.

##### PGC SCZ GWAS
## PREP GWAS DATA
# download from UNC
# munge_sumstats.py \
#   --sumstats ckqny.scz2snpres.gz \
#   --N-cas 36989 \
#   --N-con 113075 \
#   --merge-alleles w_hm3.snplist \
#   --out PGC_SCZ

## RUN THE REGRESSION
#ldsc.py \
#  --h2-cts PGC_SCZ.sumstats.gz \
#  --ref-ld-chr 1000G_EUR_Phase3_baseline_v2.2/baselineLD. \
#  --overlap-annot \
#  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
#  --out PGC_SCZ_${cts_name} \
#  --ref-ld-chr-cts ${cts_name}.ldcts \
#  --w-ld-chr weights_hm3_no_hla/weights.

##### Nalls et al 2019 PD GWAS
## PREP GWAS DATA
# download from https://drive.google.com/file/d/1FZ9UL99LAqyWnyNBxxlx6qOUlfAnublN/view?usp=sharing
# significant prep steps (more deatils in 2019-03-12_Magma_aMGL_merged_Nalls_PD.R): 
# extract and rename to ../Magma_analysis folder as PD_sumstats_Nalls.txt
# run MAGMA_Celltyping R function col_headers = format_sumstats_for_magma_linux(gwas_sumstats_path)
# need Nsum column, so make in PD_sumstats_Nalls_Nsum.txt
# awk 'BEGIN{FS=OFS="\t"} {print $0, $10+$11}' PD_sumstats_Nalls.txt > PD_sumstats_Nalls_Nsum.txt
# sed -i '1s/0/N/' PD_sumstats_Nalls_Nsum.txt
# munge_sumstats.py \
#   --sumstats ../Magma_analysis/PD_sumstats_Nalls_Nsum.txt \
#   --merge-alleles w_hm3.snplist \
#   --out Nalls_PD

## RUN THE REGRESSION
#ldsc.py \
#  --h2-cts Nalls_PD.sumstats.gz \
#  --ref-ld-chr 1000G_EUR_Phase3_baseline_v2.2/baselineLD. \
#  --overlap-annot \
#  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
#  --out Nalls_PD_${cts_name} \
#  --ref-ld-chr-cts ${cts_name}.ldcts \
#  --w-ld-chr weights_hm3_no_hla/weights.

