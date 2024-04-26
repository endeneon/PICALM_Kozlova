#! /bin/bash
# find all bams matching pattern and make index

mapfile -d $'\0' BAMs_to_index < <(find . -type f -name "Ast*WASPed.bam" -print0)

        macs2 callpeak \
                -t ${BAMs_to_index[@]} \
                -g hs \
                -f BAM \
                -q 0.05 \
                --nomodel \
                --keep-dup all \
                --outdir Ast_peaks \
		-B \
		--verbose 3 \
                -n Ast_peaks/Ast_peaks_27Feb2024 \
                1> Ast_peaks/log.out 2> Ast_peaks/log.err


#for ((i=0; i<${#BAMs_to_index[@]}; i++))
#do
#	echo ${BAMs_to_index[$i]}
#	samtools index -@ 20 \
#		${BAMs_to_index[$i]}
#done
