#! /bin/bash
# collect BAM insert size (1 txt for eachfile)

[[ -d ./insert_sumstat ]] &&
	rm -r insert_sumstats

mkdir -p insert_sumstats

# find all WASPed.bam files
mapfile -d $'\0' BAMs_to_count < <(find . -type f -name "*WASPed.bam" -print0)

for ((i=0; i<${#BAMs_to_count[@]}; i++))
do
	echo ${BAMs_to_count[$i]}

	## test if the index file exists, if not index
	[[ ! -f ${BAMs_to_count[$i]}.bai ]] &&
		samtools index -@ 20 ${BAMs_to_count[$i]}

	## get the basename of BAM file
	bam_base=$(basename -s ".bam" ${BAMs_to_count[$i]})
#	echo $bam_base

	# construct output file name
	output_file_name="insert_sumstats/"$bam_base".txt"
	echo $output_file_name

	samtools stats \
		-@ 30 \
		-d \
		-i 1000 \
		${BAMs_to_count[$i]} \
		| grep "^IS" \
		| cut -f 2- \
		> $output_file_name
done
