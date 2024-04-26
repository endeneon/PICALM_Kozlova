#! /bin/bash

# Siwei 09 Jul 2021
# added nolambda parameter to make output comparable to previous core calls

# Siwei 18 May 2021
# Call peaks use MACS2
# use conda environment encode-atac-seq-pipeline (python 3.7)
# appearently there are compatibility issues with python 3.8.8


## list all .bam files
shopt -s nullglob
bam_array=(*.bam)
echo "${bam_array[@]}"

macs2 callpeak \
	-t ${bam_array[@]} -f BAMPE -g 2.7e9 -q 0.05 \
	--nomodel \
	--nolambda \
	--keep-dup all -B \
	--shift -75 --extsize 150 \
	-n NGN2_Glut_20_lines_09Jul2021_nolambda --outdir output


