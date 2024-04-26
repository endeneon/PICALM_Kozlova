#!/bin/bash

# Siwei 6 Aug 2021
# Siwei 7 Nov 2018
# Siwei 28 Dec 2018

# remove chromosomes of alternative assembly,
# decoy, chrUn. Anything with chr name includes "_", HLA
# and sex chromosomes

for eachfile in *.bed
do
    echo $eachfile
    cat $eachfile | \
    grep -v "_alt" | grep -v "_random" | grep -v "chrEBV" | \
        grep -v "decoy" | grep -v "chrUn_" | grep -v "chrX" | \
        grep -v "chrY" | grep -v "chrM" | grep -v "HLA" | \
        bedtools slop -b 250 -i - -g /home/zhangs3/Data/Databases/Genomes/hg38/hg38.chrom.sizes > \
        bed_cleanup_250/$eachfile
done
