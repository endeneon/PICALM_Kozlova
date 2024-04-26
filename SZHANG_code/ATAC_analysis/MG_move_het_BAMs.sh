#!/bin/bash

# Siwei 04 Jul 2023
# move all BAM files of Het 10792832 to downsample folder

target_folder="downsample_100M/MG"
mkdir -p $target_folder

het_list_file="MG_rs10792832_het_list.list"

readarray -t het_file_names_array < $het_list_file

for ((i=0; i<${#het_file_names_array[@]}; i++ ))
do
	echo ${het_file_names_array[$i]}
	mv ${het_file_names_array[$i]} $target_folder
	mv ${het_file_names_array[$i]}".bai" $target_folder
done

