# Siwei 19 Feb 2024
# make peak file contains ASoC SNPs
# ! calculate npglut peaks ! ####

# init ####
{
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(org.Hs.eg.db)



  library(stringr)
  library(future)

  library(readr)

  # library(MASS)

}

# param #####
plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

load("/backup/HiC/hg38_tx_reduced_TSS_region_2k_1k.RData")

# test use 0hr ###
load("~/backuped_space/HiC/ASoC_SNP_list.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_macs2_called_new_peaks_sep_by_timextype.RData")
# Lexi_peak_set <-
#   readRDS("~/NVME/scARC_Duan_018/018-029_combined_analysis/df_Lexi_new_peak_set.RData")

npglut_peaks <-
  vector(mode = "list",
         length = 3L)
names(npglut_peaks) <-
  c("hr_0", "hr_1", "hr_6")

# for (i in 1:length(peaks_uncombined)) {
#   print(paste(i,
#               unique(str_split(string = peaks_uncombined[[i]]@elementMetadata@listData$name,
#                                pattern = "_",
#                                simplify = T)[, 1:2])))
# }
# [1] "1 npglut" "1 1hr"
# [1] "2 GABA" "2 1hr"
# [1] "3 unidentified" "3 6hr"
# [1] "4 unidentified" "4 0hr"
# [1] "5 GABA" "5 0hr"
# [1] "6 GABA" "6 6hr"
# [1] "7 nmglut" "7 1hr"
# [1] "8 npglut" "8 0hr"
# [1] "9 npglut" "9 6hr"
# [1] "10 nmglut" "10 6hr"
# [1] "11 nmglut" "11 0hr"
# [1] "12 unidentified" "12 1hr"

# ! assemble nmglut peaks ! ####
npglut_peaks[[1]] <-
  peaks_uncombined[[8]]
npglut_peaks[[2]] <-
  peaks_uncombined[[1]]
npglut_peaks[[3]] <-
  peaks_uncombined[[9]]
## load the npglut_all peakset
# load("/nvmefs/scARC_Duan_018/018-029_combined_analysis/MACS2_peaks_all_npglut.RData")
# npglut_peaks[[4]] <-
#   peaks_npglut_combined
names(npglut_peaks) <-
  c("hr_0", "hr_1", "hr_6")
# npglut_peaks[[4]]@elementMetadata@listData$ident <-
#   rep_len(x = "npglut_all",
#           length.out = length(npglut_peaks[[4]]@elementMetadata@listData$name))


# saveRDS(npglut_peaks,
#         file = "npglut_peaks_uncombined_016.RData")

## all SNPs #####
# load master SNP set
load("~/NVME/scARC_Duan_018/R_ASoC/Sum_100_lines_21Dec2023.RData")
# npglut_SNPs <-
#   vector(mode = "list",
#          length = 3L)
# names(npglut_SNPs) <-
#   c("hr_0", "hr_1", "hr_6")
npglut_SNPs <-
  list(master_vcf_list[[3]][[1]],
       master_vcf_list[[6]][[1]],
       master_vcf_list[[9]][[1]])
names(npglut_SNPs) <-
  c("hr_0", "hr_1", "hr_6")


## ASoC SNPs only
npglut_ASoC_SNPs <-
  list(master_vcf_list[[3]][[3]],
       master_vcf_list[[6]][[3]],
       master_vcf_list[[9]][[3]])
names(npglut_ASoC_SNPs) <-
  c("hr_0", "hr_1", "hr_6")


## Add a 4th set of merged 0/1/6 ASoC SNPs to test coverage
# npglut_ASoC_SNPs[['all_times']] <-
#   as.data.frame(rbind(npglut_ASoC_SNPs[[1]],
#                       npglut_ASoC_SNPs[[2]],
#                       npglut_ASoC_SNPs[[3]]))
# npglut_ASoC_SNPs[['all_times']] <-
#   npglut_ASoC_SNPs[['all_times']][!duplicated(npglut_ASoC_SNPs[['all_times']]$ID), ]

### enable this snippet if calc non-ASoC SNPs ####
### and recalc from here
npglut_non_ASoC_SNPs <-
  npglut_SNPs
for (i in 1:3) {
  curr_snp_set <- npglut_non_ASoC_SNPs[[i]]
  print(nrow(curr_snp_set))
  curr_snp_set <-
    curr_snp_set[curr_snp_set$pVal > 0.05, ]
  # print(nrow(npglut_ASoC_SNPs[[i]]))
  print(nrow(curr_snp_set))
  # curr_snp_set <-
  #   curr_snp_set[sample(size = nrow(npglut_ASoC_SNPs[[i]]),
  #                       x = nrow(curr_snp_set)), ]
  curr_snp_set$POS.1 <-
    curr_snp_set$POS
  npglut_non_ASoC_SNPs[[i]] <-
    curr_snp_set
}

npglut_ASoC_SNPs <-
  npglut_non_ASoC_SNPs
names(npglut_ASoC_SNPs) <-
  c("hr_0", "hr_1", "hr_6")

###

npglut_ASoC_SNPs_assoc_peaks <-
  vector(mode = "list",
         length = 3L)
names(npglut_ASoC_SNPs_assoc_peaks) <-
  c("hr_0", "hr_1", "hr_6")


## find overlapping peaks with ASoC SNPs #####
for (i in 1:length(npglut_ASoC_SNPs_assoc_peaks)) {
# for (i in 1:1) {
  print(i)
  gRanges_ASoC_SNPs <-
    makeGRangesFromDataFrame(df = npglut_ASoC_SNPs[[i]],
                             keep.extra.columns = T,
                             ignore.strand = T,
                             seqnames.field = "CHROM",
                             start.field = "POS",
                             end.field = "POS.1",
                             starts.in.df.are.0based = F,
                             na.rm = T)
  ASoC_overlap_hits <-
    findOverlaps(query = npglut_peaks[[i]],
                 subject = gRanges_ASoC_SNPs,
                 type = "any",
                 select = "all")

  query_ASoC_overlap_hits <-
    npglut_peaks[[i]][queryHits(ASoC_overlap_hits)]
  subject_gRanges_ASoC_SNPs <-
    gRanges_ASoC_SNPs[subjectHits(ASoC_overlap_hits)]
  # npglut_peaks_overlapped <-
  #   npglut_ASoC_SNPs_assoc_peaks[[i]][subjectHits(overlap_hits[overlap_perc > 0.2])]
  query_ASoC_overlap_hits$rsID <- # this is the peakset
    subject_gRanges_ASoC_SNPs$ID

  # peaks_to_return <-
  #   subsetByOverlaps(x = query_ASoC_overlap_hits,
  #                    ranges = gRanges_ASoC_SNPs,
  #                    type = "any")
  peaks_to_return <-
    query_ASoC_overlap_hits
  peaks_to_return <-
    reduce(peaks_to_return,
           drop.empty.ranges = T)
  peaks_to_return@elementMetadata@listData$name <-
    query_ASoC_overlap_hits$name[!duplicated(query_ASoC_overlap_hits$name)]
  peaks_to_return <-
    peaks_to_return[order(peaks_to_return$name), ]

  test_peaks_metadata <-
    query_ASoC_overlap_hits@elementMetadata@listData
  test_peaks_metadata <-
    as.data.frame(test_peaks_metadata)
  test_peaks_metadata <-
    split(x = test_peaks_metadata,
          f = test_peaks_metadata$name)
  test_peaks_metadata <-
    lapply(1:length(test_peaks_metadata), function(x) {
      return(data.frame(peak_name = unique(unlist(test_peaks_metadata[[x]][, 1])),
                        rsID = paste(unlist(unique(test_peaks_metadata[[x]][, 8])),
                                     sep = "",
                                     collapse = ",")))
    })
  # length(test_peaks_metadata)
  test_peaks_metadata <-
    do.call(what = rbind,
            args = test_peaks_metadata)
  test_peaks_metadata <-
    test_peaks_metadata[order(test_peaks_metadata$peak_name), ]

  peaks_to_return$rsID <-
    test_peaks_metadata$rsID

  npglut_ASoC_SNPs_assoc_peaks[[i]] <-
    peaks_to_return
}


ibed_5kb_all_files <-
  list.files(path = "/home/zhangs3/Data/FASTQ/microC/Capture/chicago",
             pattern = ".*_5000bp.ibed$",
             full.names = T,
             recursive = T)
## try 10000
# ibed_5kb_all_files <-
#   list.files(path = "/home/zhangs3/Data/FASTQ/microC/Capture/chicago",
#              pattern = ".*_5000bp.ibed$",
#              full.names = T,
#              recursive = T)

ibed_5kb_ranges <-
  vector(mode = "list",
         length = 4L)
names(ibed_5kb_ranges) <-
  c("CD_11_0", "CD_11_1", "CD_11_6", "CD_12_0")

npglut_5kb_ASoC_intervals <-
  vector(mode = "list",
         length = 3L)
names(npglut_5kb_ASoC_intervals) <-
  c("hr_0", "hr_1", "hr_6")

for (i in 1:length(npglut_ASoC_SNPs_assoc_peaks)) {
  npglut_5kb_ASoC_intervals[[i]] <-
    vector(mode = "list",
           length = 4L)
  names(npglut_5kb_ASoC_intervals[[i]]) <-
    c("CD_11_0", "CD_11_1", "CD_11_6", "CD_12_0")

  for (j in 1:length(ibed_5kb_ranges)) {
    print(paste(i, j))
    df_2_make_ranges <-
      read_delim(ibed_5kb_all_files[j],
                 delim = "\t",
                 escape_double = FALSE,
                 trim_ws = TRUE)

    ## try score = 7.5? #####
    df_2_make_ranges <-
      df_2_make_ranges[df_2_make_ranges$score > 7.5, ]

    ## remove bait-2-bait interactions
    # df_2_make_ranges <-
    #   df_2_make_ranges[str_detect(string = df_2_make_ranges$otherEnd_name,
    #                               pattern = "^bait.*",
    #                               negate = T), ]
    # print(nrow(df_2_make_ranges))

    ## make 5kb bins
    df_2_make_ranges$bait5k_start <-
      trunc((df_2_make_ranges$bait_start +
               df_2_make_ranges$bait_end) / 2 - 2500)
    df_2_make_ranges$bait5k_end <-
      trunc((df_2_make_ranges$bait_start +
               df_2_make_ranges$bait_end) / 2 + 2500)

    df_2_make_ranges$otherEnd5k_start <-
      trunc((df_2_make_ranges$otherEnd_start +
               df_2_make_ranges$otherEnd_end) / 2 - 2500)
    df_2_make_ranges$otherEnd5k_end <-
      trunc((df_2_make_ranges$otherEnd_start +
               df_2_make_ranges$otherEnd_end) / 2 + 2500)

    ## make target ranges, carry over bait ranges
    target_ranges <-
      makeGRangesFromDataFrame(df_2_make_ranges,
                               keep.extra.columns = T,
                               ignore.strand = T,
                               seqnames.field = "otherEnd_chr",
                               start.field = "otherEnd5k_start",
                               end.field = "otherEnd5k_end")
    overlap_hits <-
      findOverlaps(query = target_ranges,
                   subject = npglut_ASoC_SNPs_assoc_peaks[[i]],
                   type = "any",
                   select = "all")
    overlap_perc <-
      pintersect(x = target_ranges[queryHits(overlap_hits)],
                 y = npglut_ASoC_SNPs_assoc_peaks[[i]][subjectHits(overlap_hits)])
    overlap_perc <-
      width(overlap_perc) / width(npglut_ASoC_SNPs_assoc_peaks[[i]][subjectHits(overlap_hits)])
    final_hits <-
      overlap_perc[overlap_perc > 0.2]
    # target_peaks_overlapped <-
    #   target_ranges[queryHits(overlap_hits[overlap_perc > 0.2])]
    # npglut_peaks_overlapped <-
    #   npglut_ASoC_SNPs_assoc_peaks[[i]][subjectHits(overlap_hits[overlap_perc > 0.2])]

    target_peaks_overlapped <-
      target_ranges[queryHits(overlap_hits)]
    npglut_peaks_overlapped <-
      npglut_ASoC_SNPs_assoc_peaks[[i]][subjectHits(overlap_hits)]

    target_peaks_overlapped$rsID <-
      npglut_peaks_overlapped$rsID
    print(length(target_peaks_overlapped@elementMetadata@listData$bait_chr))

    # find the bait end by metadata ####
    df_bait_end <-
      as.data.frame(target_peaks_overlapped@elementMetadata@listData)
    df_bait_end$otherEnd_chr <-
      unlist(df_bait_end[, 1])
    gRanges_bait_end <-
      makeGRangesFromDataFrame(df = df_bait_end,
                               keep.extra.columns = T,
                               ignore.strand = F,
                               seqnames.field = "bait_chr",
                               start.field = "bait_start",
                               end.field = "bait_end")

    overlap_hits_bait <-
      findOverlaps(query = gRanges_bait_end,
                   subject = hg38_tx_reduced_promoters,
                   type = "any",
                   select = "all")
    bait_end_overlapped <-
      gRanges_bait_end[queryHits(overlap_hits_bait)]
    hg38_promoters_overlapped <-
      hg38_tx_reduced_promoters[subjectHits(overlap_hits_bait)]

    bait_end_overlapped$tx_promoter_annot <-
      hg38_promoters_overlapped$revmap_annot
    # gene_annot <-
    #   select(EnsDb.Hsapiens.v86,
    #          keys = bait_end_overlapped$tx_promoter_annot,
    #          columns = "GENENAME",
    #          keytype = "TXNAME")
    bait_end_overlapped$gene_promoter_annot <-
      hg38_promoters_overlapped$revmap_annot
      # select(EnsDb.Hsapiens.v86,
      #        keys = bait_end_overlapped$tx_promoter_annot,
      #        columns = "GENENAME",
      #        keytype = "TXNAME")

    npglut_5kb_ASoC_intervals[[i]][[j]] <-
      bait_end_overlapped
  }
}

for (i in 1:3) {
  for (j in 1:4) {
    print(paste(names(npglut_5kb_ASoC_intervals)[i],
                names(npglut_5kb_ASoC_intervals[[i]])[j],
                length(unique(c(str_split(string = paste(unlist(npglut_5kb_ASoC_intervals[[i]][[j]]$rsID),
                                                         sep = "",
                                                         collapse = ","),
                                          pattern = ",",
                                          simplify = T)))),
                length(npglut_ASoC_SNPs[[i]]$ID),
                paste0(format(length(unique(c(str_split(string = paste(unlist(npglut_5kb_ASoC_intervals[[i]][[j]]$rsID),
                                                                       sep = "",
                                                                       collapse = ","),
                                                        pattern = ",",
                                                        simplify = T)))) /
                                length(npglut_ASoC_SNPs[[i]]$ID) * 100,
                              digits = 3,
                              nsmall = 1),
                       '%'),
                sep = ","))
  }
}

unique_rsid_collection <-
  vector(mode = "list",
         length = 12L)

collect_counter <- 1
for (i in 1:3) {
  for (j in 1:4) {
    print(paste(names(npglut_5kb_ASoC_intervals)[i],
                names(npglut_5kb_ASoC_intervals[[i]])[j],
                length(unique(c(str_split(string = paste(unlist(npglut_5kb_ASoC_intervals[[i]][[j]]$rsID),
                                                         sep = "",
                                                         collapse = ","),
                                          pattern = ",",
                                          simplify = T)))),
                length(npglut_ASoC_SNPs[[i]]$ID),
                paste0(format(length(unique(c(str_split(string = paste(unlist(npglut_5kb_ASoC_intervals[[i]][[j]]$rsID),
                                                                       sep = "",
                                                                       collapse = ","),
                                                        pattern = ",",
                                                        simplify = T)))) /
                                length(npglut_ASoC_SNPs[[i]]$ID) * 100,
                              digits = 3,
                              nsmall = 1),
                       '%'),
                sep = ","))
    names(unique_rsid_collection)[collect_counter] <-
      paste(names(npglut_5kb_ASoC_intervals)[i],
            names(npglut_5kb_ASoC_intervals[[i]])[j],
            sep = "_")
    unique_rsid_collection[[collect_counter]] <-
      unique(c(str_split(string = paste(unlist(npglut_5kb_ASoC_intervals[[i]][[j]]$rsID),
                                        sep = "",
                                        collapse = ","),
                         pattern = ",",
                         simplify = T)))
    collect_counter <-
      collect_counter + 1
  }
}

{
  print("npglut")


  print(paste("hr_0_rest",
              length(unique(unlist(unique_rsid_collection[(c(1,4)) * 1]))),
              paste0(format(length(unique(unlist(unique_rsid_collection[(c(1,4)) * 1]))) / 11370 * 100,
                            digits = 3,
                            nsmall = 1),
                     '%'),
              sep = ","))
  print(paste("hr_0_stim",
              length(unique(unlist(unique_rsid_collection[(c(2,3) * 1)]))),
              paste0(format(length(unique(unlist(unique_rsid_collection[(c(2,3)) * 1]))) / 11370 * 100,
                            digits = 3,
                            nsmall = 1),
                     '%'),
              sep = ","))

  print(paste("hr_1_rest",
              length(unique(unlist(unique_rsid_collection[(c(1,4)) * 2]))),
              paste0(format(length(unique(unlist(unique_rsid_collection[(c(1,4)) * 2]))) / 22207 * 100,
                            digits = 3,
                            nsmall = 1),
                     '%'),
              sep = ","))
  print(paste("hr_1_stim",
              length(unique(unlist(unique_rsid_collection[(c(2,3) * 2)]))),
              paste0(format(length(unique(unlist(unique_rsid_collection[(c(2,3)) * 2]))) / 22207 * 100,
                            digits = 3,
                            nsmall = 1),
                     '%'),
              sep = ","))

  print(paste("hr_6_rest",
              length(unique(unlist(unique_rsid_collection[(c(1,4)) * 3]))),
              paste0(format(length(unique(unlist(unique_rsid_collection[(c(1,4)) * 3]))) / 21846 * 100,
                            digits = 3,
                            nsmall = 1),
                     '%'),
              sep = ","))
  print(paste("hr_6_stim",
              length(unique(unlist(unique_rsid_collection[(c(2,3) * 3)]))),
              paste0(format(length(unique(unlist(unique_rsid_collection[(c(2,3)) * 3]))) / 21846 * 100,
                            digits = 3,
                            nsmall = 1),
                     '%'),
              sep = ","))

}
saveRDS(unique_rsid_collection,
        file = "npglut_unique_rsid_rest_stim_HiC.RDs")

# "npglut"
# [1] "hr_0_rest,3384,29.8%"
# [1] "hr_0_stim,3380,29.7%"
# [1] "hr_1_rest,5260,23.7%"
# [1] "hr_1_stim,5661,25.5%"
# [1] "hr_6_rest,5368,24.6%"
# [1] "hr_6_stim,7369,33.7%"

## try combined ranges (0 hr x2 vs 1+6 hr) ####

combined_ibed_5kb_ranges <-
  vector(mode = "list",
         length = 2L)
names(combined_ibed_5kb_ranges) <-
  c("rest", "stimulated")

# for (i in 1:length(combined_ibed_5kb_ranges)) {
combined_ibed_5kb_ranges[[1]] <-
  as.data.frame(rbind(read_delim(ibed_5kb_all_files[1],
                                 delim = "\t",
                                 escape_double = FALSE,
                                 trim_ws = TRUE),
                      read_delim(ibed_5kb_all_files[4],
                                 delim = "\t",
                                 escape_double = FALSE,
                                 trim_ws = TRUE)))
combined_ibed_5kb_ranges[[1]]$uuid <-
  str_c(combined_ibed_5kb_ranges[[1]]$bait_chr,
        combined_ibed_5kb_ranges[[1]]$bait_start,
        combined_ibed_5kb_ranges[[1]]$bait_end,
        combined_ibed_5kb_ranges[[1]]$otherEnd_chr,
        combined_ibed_5kb_ranges[[1]]$otherEnd_start,
        combined_ibed_5kb_ranges[[1]]$otherEnd_end,
        sep = "_")
combined_ibed_5kb_ranges[[1]] <-
  combined_ibed_5kb_ranges[[1]][!duplicated(combined_ibed_5kb_ranges[[1]]$uuid), ]

combined_ibed_5kb_ranges[[2]] <-
  as.data.frame(rbind(read_delim(ibed_5kb_all_files[2],
                                 delim = "\t",
                                 escape_double = FALSE,
                                 trim_ws = TRUE),
                      read_delim(ibed_5kb_all_files[3],
                                 delim = "\t",
                                 escape_double = FALSE,
                                 trim_ws = TRUE)))
combined_ibed_5kb_ranges[[2]]$uuid <-
  str_c(combined_ibed_5kb_ranges[[2]]$bait_chr,
        combined_ibed_5kb_ranges[[2]]$bait_start,
        combined_ibed_5kb_ranges[[2]]$bait_end,
        combined_ibed_5kb_ranges[[2]]$otherEnd_chr,
        combined_ibed_5kb_ranges[[2]]$otherEnd_start,
        combined_ibed_5kb_ranges[[2]]$otherEnd_end,
        sep = "_")
combined_ibed_5kb_ranges[[2]] <-
  combined_ibed_5kb_ranges[[2]][!duplicated(combined_ibed_5kb_ranges[[2]]$uuid), ]

# combined_ibed_5kb_ranges[[1]] <-
#   combined_ibed_5kb_ranges[[1]][combined_ibed_5kb_ranges[[1]]$score > 7.5, ]
# combined_ibed_5kb_ranges[[2]] <-
#   combined_ibed_5kb_ranges[[2]][combined_ibed_5kb_ranges[[2]]$score > 7.5, ]

for (i in 1:length(npglut_ASoC_SNPs_assoc_peaks)) {
  npglut_5kb_ASoC_intervals[[i]] <-
    vector(mode = "list",
           length = 2L)
  names(npglut_5kb_ASoC_intervals[[i]]) <-
    c("rest", "stimulated")
    # c("CD_11_0", "CD_11_1", "CD_11_6", "CD_12_0")

  for (j in 1:length(combined_ibed_5kb_ranges)) {
    print(paste(i, j))
    df_2_make_ranges <-
      combined_ibed_5kb_ranges[[j]]
    df_2_make_ranges <-
      df_2_make_ranges[df_2_make_ranges$score > 7.5, ]
    #   read_delim(ibed_5kb_all_files[j],
    #              delim = "\t",
    #              escape_double = FALSE,
    #              trim_ws = TRUE)
    ## remove bait-2-bait interactions
    # df_2_make_ranges <-
    #   df_2_make_ranges[str_detect(string = df_2_make_ranges$otherEnd_name,
    #                               pattern = "^bait.*",
    #                               negate = T), ]
    # print(nrow(df_2_make_ranges))

    ## make 5kb bins
    df_2_make_ranges$bait5k_start <-
      trunc((df_2_make_ranges$bait_start +
               df_2_make_ranges$bait_end) / 2 - 2500)
    df_2_make_ranges$bait5k_end <-
      trunc((df_2_make_ranges$bait_start +
               df_2_make_ranges$bait_end) / 2 + 2500)

    df_2_make_ranges$otherEnd5k_start <-
      trunc((df_2_make_ranges$otherEnd_start +
               df_2_make_ranges$otherEnd_end) / 2 - 2500)
    df_2_make_ranges$otherEnd5k_end <-
      trunc((df_2_make_ranges$otherEnd_start +
               df_2_make_ranges$otherEnd_end) / 2 + 2500)

    # df_2_make_ranges$otherEnd5k_chr <-
    #   df_2_make_ranges$otherEnd_chr
    ## remove inter-chrom links
    df_2_make_ranges <-
      df_2_make_ranges[df_2_make_ranges$bait_chr == df_2_make_ranges$otherEnd_chr, ]

    ## make target ranges, carry over bait ranges
    target_ranges <-
      makeGRangesFromDataFrame(df_2_make_ranges,
                               keep.extra.columns = T,
                               ignore.strand = T,
                               seqnames.field = "otherEnd_chr",
                               start.field = "otherEnd5k_start",
                               end.field = "otherEnd5k_end")
    overlap_hits <-
      findOverlaps(query = target_ranges,
                   subject = npglut_ASoC_SNPs_assoc_peaks[[i]],
                   type = "any",
                   select = "all")
    overlap_perc <-
      pintersect(x = target_ranges[queryHits(overlap_hits)],
                 y = npglut_ASoC_SNPs_assoc_peaks[[i]][subjectHits(overlap_hits)])
    overlap_perc <-
      width(overlap_perc) / width(npglut_ASoC_SNPs_assoc_peaks[[i]][subjectHits(overlap_hits)])
    # final_hits <-
    #   overlap_perc[overlap_perc > 0.2]
    final_hits <-
      overlap_perc > 0.2
    target_peaks_overlapped <-
      target_ranges[queryHits(overlap_hits[final_hits])]
    npglut_peaks_overlapped <-
      npglut_ASoC_SNPs_assoc_peaks[[i]][subjectHits(overlap_hits[final_hits])]
    # target_peaks_overlapped <-
    #   target_ranges[queryHits(overlap_hits[overlap_perc > 0.2])]
    # npglut_peaks_overlapped <-
    #   npglut_ASoC_SNPs_assoc_peaks[[i]][subjectHits(overlap_hits[overlap_perc > 0.2])]

    target_peaks_overlapped <-
      target_ranges[queryHits(overlap_hits)]
    npglut_peaks_overlapped <-
      npglut_ASoC_SNPs_assoc_peaks[[i]][subjectHits(overlap_hits)]

    target_peaks_overlapped$rsID <-
      npglut_peaks_overlapped$rsID
    print(length(target_peaks_overlapped@elementMetadata@listData$bait_chr))

    # find the bait end by metadata ####
    df_bait_end <-
      as.data.frame(target_peaks_overlapped@elementMetadata@listData)
    # df_bait_end$otherEnd_chr <-
    #   unlist(df_bait_end[, 1])
    # df_bait_end <-
    #   df_bait_end[df_bait_end$bait_chr == df_bait_end$otherEnd_chr, ]
    gRanges_bait_end <-
      makeGRangesFromDataFrame(df = df_bait_end,
                               keep.extra.columns = T,
                               ignore.strand = F,
                               seqnames.field = "bait_chr",
                               start.field = "bait_start",
                               end.field = "bait_end")

    overlap_hits_bait <-
      findOverlaps(query = gRanges_bait_end,
                   subject = hg38_tx_reduced_promoters,
                   type = "any",
                   select = "all")
    bait_end_overlapped <-
      gRanges_bait_end[queryHits(overlap_hits_bait)]
    hg38_promoters_overlapped <-
      hg38_tx_reduced_promoters[subjectHits(overlap_hits_bait)]

    bait_end_overlapped$tx_promoter_annot <-
      hg38_promoters_overlapped$revmap_annot
    # gene_annot <-
    #   select(EnsDb.Hsapiens.v86,
    #          keys = bait_end_overlapped$tx_promoter_annot,
    #          columns = "GENENAME",
    #          keytype = "TXNAME")
    bait_end_overlapped$gene_promoter_annot <-
      hg38_promoters_overlapped$revmap_annot
    # select(EnsDb.Hsapiens.v86,
    #        keys = bait_end_overlapped$tx_promoter_annot,
    #        columns = "GENENAME",
    #        keytype = "TXNAME")

    npglut_5kb_ASoC_intervals[[i]][[j]] <-
      bait_end_overlapped
  }
}

# SZ_ASoC_microC_intersect <-
#   readxl::read_excel("SZ_ASoC_microC-intersect.xlsx")
# SZ_ASoC_microC_intersect <-
#   SZ_ASoC_microC_intersect$RSID

for (i in 1:3) {
  for (j in 1:2) {
    print(paste(names(npglut_5kb_ASoC_intervals)[i],
                names(npglut_5kb_ASoC_intervals[[i]])[j],
                length(unique(c(str_split(string = paste(unlist(npglut_5kb_ASoC_intervals[[i]][[j]]$rsID),
                                                         sep = "",
                                                         collapse = ","),
                                          pattern = ",",
                                          simplify = T)))),
                length(npglut_ASoC_SNPs[[i]]$ID),
                paste0(format(length(unique(c(str_split(string = paste(unlist(npglut_5kb_ASoC_intervals[[i]][[j]]$rsID),
                                                                       sep = "",
                                                                       collapse = ","),
                                                        pattern = ",",
                                                        simplify = T)))) /
                                length(npglut_ASoC_SNPs[[i]]$ID) * 100,
                              digits = 3,
                              nsmall = 1),
                       '%'),
                sep = ","))

    # unique_SNP_list <-
    #   unique(c(str_split(string = paste(unlist(npglut_5kb_ASoC_intervals[[i]][[j]]$rsID),
    #                                     sep = "",
    #                                     collapse = ","),
    #                      pattern = ",",
    #                      simplify = T)))
    # print(length(SZ_ASoC_microC_intersect %in% unique_SNP_list))
    # print(SZ_ASoC_microC_intersect[SZ_ASoC_microC_intersect %in% unique_SNP_list])
    # print(SZ_ASoC_microC_intersect[unique_SNP_list %in% SZ_ASoC_microC_intersect])

  }
}

## npglut stats #####
### ASoC SNPs ###
# [1] "hr_0,rest,3321,11370,29.2%"
# [1] "hr_0,stimulated,3226,11370,28.4%"
# [1] "hr_1,rest,5436,22207,24.5%"
# [1] "hr_1,stimulated,5578,22207,25.1%"
# [1] "hr_6,rest,5182,21846,23.7%"
# [1] "hr_6,stimulated,5102,21846,23.4%"
### non-ASoC SNPs ###
# [1] "hr_0,rest,3181,11370,28.0%"
# [1] "hr_0,stimulated,3142,11370,27.6%"
# [1] "hr_1,rest,5710,22207,25.7%"
# [1] "hr_1,stimulated,5682,22207,25.6%"
# [1] "hr_6,rest,5510,21846,25.2%"
# [1] "hr_6,stimulated,5573,21846,25.5%"


## @cutoff=7.5 #####
## ASoC SNPs
# [1] "hr_0,rest,1340,11370,11.8%"
# [1] "hr_0,stimulated,1118,11370,9.83%"
# [1] "hr_1,rest,2146,22207,9.66%"
# [1] "hr_1,stimulated,1864,22207,8.39%"
# [1] "hr_6,rest,1975,21846,9.04%"
# [1] "hr_6,stimulated,1620,21846,7.42%"
## non-ASoC SNPs #
# [1] "hr_0,rest,1303,11370,11.5%"
# [1] "hr_0,stimulated,1197,11370,10.5%"
# [1] "hr_1,rest,2243,22207,10.1%"
# [1] "hr_1,stimulated,2150,22207,9.68%"
# [1] "hr_6,rest,2061,21846,9.43%"
# [1] "hr_6,stimulated,1962,21846,8.98%"


#####
for (i in 1:length(npglut_ASoC_SNPs)) {
  # print(i)
  gRanges_ASoC_SNPs <-
    makeGRangesFromDataFrame(df = npglut_ASoC_SNPs[[i]],
                             keep.extra.columns = T,
                             ignore.strand = T,
                             seqnames.field = "CHROM",
                             start.field = "POS",
                             end.field = "POS.1",
                             starts.in.df.are.0based = F,
                             na.rm = T)
  ASoC_overlap_hits <-
    findOverlaps(query = hg38_tx_reduced_promoters,
                 subject = gRanges_ASoC_SNPs,
                 type = "any",
                 select = "all")

  # query_ASoC_overlap_hits <-
  #   npglut_peaks[[i]][queryHits(ASoC_overlap_hits)]
  subject_gRanges_ASoC_SNPs <-
    gRanges_ASoC_SNPs[subjectHits(ASoC_overlap_hits)]
  print(paste(i,
              length(subject_gRanges_ASoC_SNPs$REF),
              length(gRanges_ASoC_SNPs$REF),
              paste0(length(subject_gRanges_ASoC_SNPs$REF) /
                       length(gRanges_ASoC_SNPs$REF) * 100,
                     '%'),
              sep = ","))
}

# in promoter region, npglut
## ASoC ####
# [1] "1,3145,11370,27.660510114336%"
# [1] "2,4473,22207,20.1422974737695%"
# [1] "3,4296,21846,19.6649272177973%"
# non-ASoC
# [1] "1,39961,141133,28.3144268172575%"
# [1] "2,43168,185151,23.3150239534218%"
# [1] "3,42767,188756,22.6572930132022%"

log2(19348 / (3e3 * 24511 / 2.7e9 * 67380))
