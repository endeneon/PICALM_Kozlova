#! /bin/bash
# find all bams matching pattern and make index

# mapfile -d $'\0' BAMs_to_index < <(find . -type f -name "*.bam" -print0)

for eachfile in *.bam
do
	echo $eachfile
	samtools index -@ 20 $eachfile

        macs2 callpeak \
                -t $eachfile \
                -g hs \
                -f BAMPE \
                -q 0.05 \
                --nomodel \
                --keep-dup all \
                --outdir peaks_100M \
		-B \
		--verbose 3 \
                -n "peaks_"${eachfile/%.bam/}

done

#                1> Ast_peaks/log.out 2> Ast_peaks/log.err


#for ((i=0; i<${#BAMs_to_index[@]}; i++))
#do
#	echo ${BAMs_to_index[$i]}
#	samtools index -@ 20 \
#		${BAMs_to_index[$i]}
#done
