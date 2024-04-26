# python code only, will be able to run under command
# Siwei 27 Feb 2024

## imports
import os
import deeptools
import glob

import deeptools.countReadsPerBin as crpb
import pysam

import statistics
import numpy
import pandas as pd

# define files to be counted
# ! note: define peak files here as well !
path = ""
print(path)
bam_files = []
for each_file in glob.glob(path + '**/Ast*WASPed.bam', recursive=True):
    print(each_file)
    bam_files.append(each_file)

peak_file = ["Ast_new_peaks_21Mar2023_peaks.narrowPeak"]

# count FRiP
print(bam_files)
cr = crpb.CountReadsPerBin(bamFilesList=bam_files, # if more than one elements, do not need to relist by []
                           bedFile=peak_file, 
                           numberOfProcessors=40)
reads_at_peaks = cr.run()
reads_at_peaks.shape

# count total reads
total_reads_in_peaks = reads_at_peaks.sum(axis=0)

total_bam_reads = []
for each_file in bam_files:    
    tr = pysam.AlignmentFile(each_file)
    # temp_total_reads = tr.run()
    total_bam_reads.append(tr.mapped)
    # print(tr.mapped)
    # total_bam_reads.appendx)

# count FRiP%
FRiP_count = pd.DataFrame(total_reads_in_peaks)
FRiP_count.index = bam_files
FRiP_count.columns = ["total_reads_peaks"]

FRiP_count['total_bam_reads'] = total_bam_reads

FRiP_count['FRiP'] = [x / y for x, y in zip(total_reads_in_peaks, total_bam_reads)]

# write out
print(FRiP_count)

FRiP_count.to_csv("Ast_Unified_all_lines_FRiP.txt",
                  sep = "\t")