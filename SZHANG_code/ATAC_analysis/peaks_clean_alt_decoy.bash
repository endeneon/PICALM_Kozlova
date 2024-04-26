#!/bin/bash

# Siwei 7 Nov 2018

# remove chromosomes of alternative assembly, 
# decoy, chrUn. Anything with chr name includes "_", 
# and sex chromosomes

for EACHFILE in *.bed
do
	echo $EACHFILE
	cat $EACHFILE | grep -v "_alt" | grep -v "_random" | \
		grep -v "decoy" | grep -v "chrUn_" | grep -v "chrX" | \
		grep -v "chrY" | grep -v "chrM" | grep -v "HLA" | grep -v "EBV" \
		> bed_cleanup/cleanup_$EACHFILE
done

