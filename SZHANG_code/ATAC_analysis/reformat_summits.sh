#!/bin/bash
# Siwei 22 Apr 2024

mkdir -p summit_cleaned

for eachfile in peaks*summits.bed
do
	echo $eachfile
	cat $eachfile \
		| awk -F '\t' '$1 ~ /chr[0-9]|chr[0-9][0-9]/ {print $0}' \
		> summit_cleaned/$eachfile
done

