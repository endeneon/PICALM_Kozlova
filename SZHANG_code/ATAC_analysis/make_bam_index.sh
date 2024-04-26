#! /bin/bash
# find all bams matching pattern and make index

mapfile -d $'\0' BAMs_to_index < <(find . -type f \( -name "Glut_rapid*WASPed.bam" -o -name "R21*WASPed.bam" \) -print0)

for ((i=0; i<${#BAMs_to_index[@]}; i++))
do
	echo ${BAMs_to_index[$i]}
	samtools index -@ 20 \
		${BAMs_to_index[$i]}
done
