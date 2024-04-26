#!/bin/bash

# Siwei 23 Jun 2023
# Siwei 07 Jun 2022



export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16

mkdir -p results

shopt -s nullglob
sumstats_list=(/home/zhangs3/Data/Databases/GWAS/MAGMA_ref_new/ldsc_sumstats/*.sumstats.gz)
shopt -u nullglob

weights_chr="/home/zhangs3/NVME/python_projects/ldsc/weights_hm3_no_hla/weights."
freq_chr="/home/zhangs3/NVME/python_projects/ldsc/1000G_Phase3_frq/1000G.EUR.QC."

echo "${#sumstats_list[@]}"

for (( i=0; i<${#sumstats_list[@]}; i++ ))
do
	echo ${sumstats_list[$i]}
	base_filename="$(basename -- ${sumstats_list[$i]})"
	echo $base_filename

	python ../ldsc.py \
		--h2 ${sumstats_list[$i]} \
		--ref-ld-chr annot_file/Ast_new_peaks_07Jun2023_peaks_hg19.,annot_file/MG_17_lines_peaks_hg19.,annot_file/GA_20_lines_17Jun2021_peaks_hg19.,annot_file/NGN2_Glut_peak_set_4_LDSC_21Jun2023_peaks_hg19.,annot_file/DN_new_peaks_26Jun2023_peaks_hg19. \
		--w-ld-chr $weights_chr \
		--frqfile-chr $freq_chr \
		--print-coefficients \
		--overlap-annot \
		--out results/$base_filename

done

#		--h2 ${sumstats_list[$i]} \

unset OMP_NUM_THREADS
unset MKL_NUM_THREADS

