# Siwei 19 Feb 2024
# make peak file contains ASoC SNPs

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

for (i in 1:length(peaks_uncombined)) {
  print(paste(i,
              unique(str_split(string = peaks_uncombined[[i]]@elementMetadata@listData$name,
                               pattern = "_",
                               simplify = T)[, 1:2])))
}
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

npglut_peaks[[1]] <-
  peaks_uncombined[[8]]
npglut_peaks[[2]] <-
  peaks_uncombined[[1]]
npglut_peaks[[3]] <-
  peaks_uncombined[[9]]
## load the npglut_all peakset
load("/nvmefs/scARC_Duan_018/018-029_combined_analysis/MACS2_peaks_all_npglut.RData")
npglut_peaks[[4]] <-
  peaks_npglut_combined
names(npglut_peaks) <-
  c("hr_0", "hr_1", "hr_6", "all_times")
npglut_peaks[[4]]@elementMetadata@listData$ident <-
  rep_len(x = "npglut_all",
          length.out = length(npglut_peaks[[4]]@elementMetadata@listData$name))


saveRDS(npglut_peaks,
        file = "npglut_peaks_uncombined_016.RData")

## all SNPs
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
npglut_ASoC_SNPs[['all_times']] <-
  as.data.frame(rbind(npglut_ASoC_SNPs[[1]],
                      npglut_ASoC_SNPs[[2]],
                      npglut_ASoC_SNPs[[3]]))
npglut_ASoC_SNPs[['all_times']] <-
  npglut_ASoC_SNPs[['all_times']][!duplicated(npglut_ASoC_SNPs[['all_times']]$ID), ]



npglut_ASoC_SNPs_assoc_peaks <-
  vector(mode = "list",
         length = 4L)
names(npglut_ASoC_SNPs_assoc_peaks) <-
  c("hr_0", "hr_1", "hr_6", "all_times")


## find overlapping peaks with ASoC SNPs #####
for (i in 1:length(npglut_ASoC_SNPs_assoc_peaks)) {
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
  query_ASoC_overlap_hits$rsID <-
    subject_gRanges_ASoC_SNPs$ID

  peaks_to_return <-
    subsetByOverlaps(x = query_ASoC_overlap_hits,
                     ranges = gRanges_ASoC_SNPs,
                     type = "any")
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

  peaks_to_return$rsID <-
    test_peaks_metadata$rsID

  npglut_ASoC_SNPs_assoc_peaks[[i]] <-
    peaks_to_return
}

# test_peaks_reduced <-
#   reduce(peaks_to_return,
#          drop.empty.ranges = T)
# test_peaks_reduced@elementMetadata@listData$name <-
#   peaks_to_return$name[!duplicated(peaks_to_return$name)]
# test_peaks_reduced <-
#   test_peaks_reduced[order(test_peaks_reduced$name), ]
#
# test_peaks_metadata <-
#   query_ASoC_overlap_hits@elementMetadata@listData
# test_peaks_metadata <-
#   as.data.frame(test_peaks_metadata)
# # test_peaks_metadata <-
# test_peaks_metadata <-
#   split(x = test_peaks_metadata,
#         f = test_peaks_metadata$name)
# test_peaks_metadata <-
#   lapply(1:length(test_peaks_metadata), function(x) {
#     return(data.frame(peak_name = unique(unlist(test_peaks_metadata[[x]][, 1])),
#                       rsID = paste(unlist(unique(test_peaks_metadata[[x]][, 8])),
#                                   sep = "",
#                                   collapse = ",")))
#   })
# length(test_peaks_metadata)
# test_peaks_metadata <-
#   do.call(what = rbind,
#           args = test_peaks_metadata)
#
#
#
# test_peaks_reduced <-
#   split(x = query_ASoC_overlap_hits,
#         f = query_ASoC_overlap_hits$name,
#         drop = T)


ibed_5kb_all_files <-
  list.files(path = "/home/zhangs3/Data/FASTQ/microC/Capture/chicago",
             pattern = ".*_5000bp.ibed$",
             full.names = T,
             recursive = T)

ibed_5kb_ranges <-
  vector(mode = "list",
         length = 4L)
names(ibed_5kb_ranges) <-
  c("CD_11_0", "CD_11_1", "CD_11_6", "CD_12_0")

npglut_5kb_ASoC_intervals <-
  vector(mode = "list",
         length = 4L)
names(npglut_5kb_ASoC_intervals) <-
  c("hr_0", "hr_1", "hr_6", "all_times")

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

# for (i in 1:3) {
#   for (j in 1:4) {
#     print(paste(names(npglut_5kb_ASoC_intervals)[i],
#                 names(npglut_5kb_ASoC_intervals[[i]])[j],
#                 length(unique(unlist(npglut_5kb_ASoC_intervals[[i]][[j]]$gene_promoter_annot)))))
#   }
# }


for (i in 1:4) {
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


length(unique(c(str_split(string = paste(unlist(gRanges_bait_end$rsID),
                                         sep = "",
                                         collapse = ","),
                          pattern = ",",
                          simplify = T))))
