# Siwei 14 Jul 2023
# Use DiffBind to identify the most accessible peaks of the three
# neuronal types
# for the maximum separation power

# init #####
library(readr)
library(DiffBind)

library(stringr)

library(RColorBrewer)
library(factoextra)

# make the input sample sheet
df_sample_sheet <-
  data.frame(bamReads = c(list.files(path = "/home/zhangs3/Data/FASTQ/Duan_Project_014/GA/macs2", # GA
                                     pattern = ".*WASPed\\.bam$",
                                     full.names = T,
                                     recursive = T),
                          list.files(path = "/home/zhangs3/NVME/Hanwen_ATACseq_test/GA_R21/GA_R21_reprocessed/WASPed_BAMs/GA_BAMs", # GA
                                     pattern = ".*WASPed\\.bam$",
                                     full.names = T,
                                     recursive = T),
                          list.files(path = "/home/zhangs3/NVME/Hanwen_ATACseq_test/GA_R21/new_May2023/WASPed_BAMs/GA", # GA
                                     pattern = ".*WASPed\\.bam$",
                                     full.names = T,
                                     recursive = T),
                          list.files(path = "/home/zhangs3/Data/FASTQ/Duan_Project_014/DN/macs2", # DN
                                     pattern = ".*WASPed\\.bam$",
                                     full.names = T,
                                     recursive = T),
                          list.files(path = "/home/zhangs3/Data/FASTQ/Duan_Project_014/NGN2_neuron_BAMs_20/WASPed_BAMs/BAMs", # NGN2-Glut
                                     pattern = "Glut.*new_WASPed\\.bam$",
                                     full.names = T,
                                     recursive = T),
                          list.files(path = "/home/zhangs3/NVME/Hanwen_ATACseq_test/GA_R21/", # NGN2-Glut
                                     pattern = "R21.*WASPed\\.bam$",
                                     full.names = T,
                                     recursive = T)))

df_sample_sheet$SampleID <-
  unlist(lapply(df_sample_sheet$bamReads, function(x) {
    tail(unlist(str_split(string = x,
                          pattern = "/",
                          simplify = F)),
         1)
  }))

df_sample_sheet$Tissue <-
  str_split(string = str_replace_all(string = df_sample_sheet$SampleID,
                                     pattern = "_",
                                     replacement = "-"),
            pattern = "-",
            simplify = T)[, 1]
df_sample_sheet$Tissue[df_sample_sheet$Tissue %in% "R21"] <- "Glut"
df_sample_sheet$Tissue <-
  factor(df_sample_sheet$Tissue,
         levels = c("Glut", "GA", "DN"))
df_sample_sheet$Condition <- "Neuron"
df_sample_sheet$Replicate <-
  c(1:sum(df_sample_sheet$Tissue %in% "GA"),
    1:sum(df_sample_sheet$Tissue %in% "DN"),
    1:sum(df_sample_sheet$Tissue %in% "Glut"))

df_sample_sheet$Peaks <-
  "/home/zhangs3/backuped_space/R_MG_Ast_GA_NGN2_PCA_16Jun2023/neuron_raw/collapsed_GA_DN_NGN2_neurons_only_501bp_13Jul2023_cleaned.bed.gz"
df_sample_sheet$PeakCaller <- "bed"
# df_sample_sheet$PeakFormat <- "bed"


# make the main DBA object contains the 3 neuron types
dba_GA_DN_NGN2_raw <-
  dba(sampleSheet = df_sample_sheet,
      # mask = NULL,
      config = data.frame(RunParallel = T,
                          # AnalysisMethod = DBA_EDGER_GLM,
                          # bCorPlot = T,
                          # bUsePval = F,
                          bRemoveM = T))

# dba_GA_DN_NGN2_raw <-
#   dba(sampleSheet = df_sample_sheet)

# count reads
dba_GA_DN_NGN2_raw <-
  dba.count(dba_GA_DN_NGN2_raw,
            bRemoveDuplicates = T,
            bParallel = T,
            summits = F,
            minOverlap = 100)
save.image("GA_DN_NGN2_DiffBind_Raw_data_14Jul2023.RData")

## Add contrasts (3) #####
## Run each analysis after established one contrast


### calc Glut vs DN #####
dba_GA_DN_NGN2_raw <-
  dba.contrast(dba_GA_DN_NGN2_raw,
               contrast = c("Tissue",
                            "DN",
                            "Glut"))
# dba_GA_DN_NGN2_results <-
#   dba.analyze(dba_GA_DN_NGN2_raw,
#               method = DBA_DESEQ2,
#               bBlacklist = DBA_BLACKLIST_HG38)

dba_GA_DN_NGN2_raw <-
  dba.contrast(dba_GA_DN_NGN2_raw,
               contrast = c("Tissue",
                            "GA",
                            "Glut"))
# dba_GA_DN_NGN2_results <-
#   dba.analyze(dba_GA_DN_NGN2_raw,
#               method = DBA_DESEQ2,
#               bBlacklist = DBA_BLACKLIST_HG38)

dba_GA_DN_NGN2_raw <-
  dba.contrast(dba_GA_DN_NGN2_raw,
               contrast = c("Tissue",
                            "GA",
                            "DN"))
dba_GA_DN_NGN2_results <-
  dba.analyze(dba_GA_DN_NGN2_raw,
              method = DBA_DESEQ2,
              bBlacklist = DBA_BLACKLIST_HG38)


save.image("GA_DN_NGN2_DiffBind_results_18Jul2023.RData")

## save peak indices into a list (3 indices) #####
## take FDR < 0.01 peaks;
significant_FDR_peak_ids <-
  vector(mode = "list",
         length = 3L)
names(significant_FDR_peak_ids) <-
  c("DNvsGlut",
    "GAvsGlut",
    "GAvsDN")

### calc DN/GA vs Glut, get ~1e6 peaks back each #####
for (i in 1:2) {
  significant_FDR_peak_ids[[i]] <-
    dba_GA_DN_NGN2_results$contrasts[[i]]$DESeq2$de$id[dba_GA_DN_NGN2_results$contrasts[[i]]$DESeq2$de$padj < 1e-12]
}

### DN vs GA peaks need to be calculated separately #####
### for equal distance, also need ~1e6 peaks by controlling
### the FDR cutoff used
significant_FDR_peak_ids[[3]] <-
  dba_GA_DN_NGN2_results$contrasts[[3]]$DESeq2$de$id[dba_GA_DN_NGN2_results$contrasts[[3]]$DESeq2$de$padj < 1e-4]

## check the total number of the union of peaks (148133)
length(unique(c(unlist(significant_FDR_peak_ids))))

## Extract RPKM from each sample for PCA analysis
## Use matrix during binding to save computation time (especially for t() and data.frame())


for (i in 1:length(dba_GA_DN_NGN2_results$peaks)) {
  print(i)
  if (i == 1) {
    df_raw_RPKM <-
      matrix(dba_GA_DN_NGN2_raw$peaks[[i]]$RPKM, 1)
  } else {
    df_raw_RPKM <-
      rbind(df_raw_RPKM,
            matrix(dba_GA_DN_NGN2_raw$peaks[[i]]$RPKM, 1))
  }
}
df_raw_RPKM <- as.data.frame(df_raw_RPKM)

## subset rows (peaks) that are in the significant FDR peak index #####
## note that default colnames start with "V"
df_raw_RPKM_sig_peaks <-
  df_raw_RPKM[, colnames(df_raw_RPKM) %in%
                str_c("V",
                      as.character(unique(c(unlist(significant_FDR_peak_ids)))),
                      sep = "")]


# calculate PCA #####
rpkm_mtable <-
  log1p(df_raw_RPKM_sig_peaks)
rpkm_mtable <-
  scale(rpkm_mtable)

# calc prin comp
prin_comp_rpkm <-
  prcomp(as.matrix(rpkm_mtable),
         center = T,
         scale. = T)

fviz_pca_ind(prin_comp_rpkm,
             geom.ind = "point",
             geom.var = "text",
             habillage = df_sample_sheet$Tissue,
             palette = c(brewer.pal(n = 8, name = "Dark2"),
                         brewer.pal(n = 8, name = "Paired")),
             # label = "ind",
             invisible = (c("quali")))

save.image("GA_DN_NGN2_DiffBind_Raw_data_14Jul2023.RData")
