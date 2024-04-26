#!/bin/bash

# Siwei 29 Oct 2021

rm MD5_sum.txt
for eachfile in *.bam
do
	echo $eachfile
	printf "$eachfile\t" >> MD5_sum.txt
	md5sum $eachfile >> MD5_sum.txt
done
