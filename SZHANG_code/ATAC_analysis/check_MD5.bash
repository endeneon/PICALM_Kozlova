#!/bin/bash
# Siwei 26 Mar 2019

for eachfile in *.bed
do
	echo $eachfile
	md5sum $eachfile >> BED_md5.txt
done
