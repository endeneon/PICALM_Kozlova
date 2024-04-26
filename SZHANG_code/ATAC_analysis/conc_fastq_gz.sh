#!/bin/bash

# Siwei 2 Jun 2021
# concacenate fastq files generated from two
# flow cells into one by a list, 
# input file names from the first two columns, 
# output files by the third column

# note use IFS= to preserve white space issue for fail-safe

while IFS=$'\t' read -r part_1 part_2 output_file; do
	echo "$part_1" #"$part_2" "$output_file"
	cat <(cat $part_1 | gzip -d) <(cat $part_2 | gzip -d) | \
		pigz -9 -p 20 > \
		conc_output/$output_file
done < file_list.list


