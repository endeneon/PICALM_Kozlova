#! /bin/bash

# Siwei 14 May 2021
# Run Caper/Cromwell for ATAC-Seq QC

# use conda env encode-atac-seq-pipeline

for eachfile in json_in/*.json
do
	echo $eachfile
	caper run \
		~/Data/Tools/atac-seq-pipeline-2.1.3/atac.wdl \
		-i $eachfile \
		--conda encode-atac-seq-pipeline
done


