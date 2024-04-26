#!/bin/bash

# Siwei 2 Jun 2021
# concacenate fastq files generated from two
# flow cells into one by a list, 
# input file names from the first two columns, 
# output files by the third column

# note use IFS= to preserve white space issue for fail-safe

while IFS=$'\t' read -r part_1 part_2 output_file; do
	echo "$part_1" "$part_2" "$output_file"
	cat $part_1 $part_2 | \
		gzip > \
		conc_output/$output_file
done < file_list_no_compress.list


