# count the % of ASoC SNPs that resides within peaks
# that also falls into MicroC bait intervals whose target
# resides in TSS region (-2 kb ~ +1 kb)

# use separate npglut peak sets (0/1/6)
# this will not be necessary since ASoC SNPs are always inside peaks

# Siwei 16 Feb 2024

# init ####
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
# library(GenomeInfoDb)
# library(GenomicFeatures)
# library(AnnotationDbi)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

library(RColorBrewer)

library(stringr)
library(future)

plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load data ####
load("~/NVME/scARC_Duan_018/R_ASoC/Sum_100_lines_21Dec2023.RData")

ASoC_SNP_list <-
  vector(mode = "list",
         length = length(master_vcf_list))
for (i in 1:length(ASoC_SNP_list)) {
  print(i)
  ASoC_SNP_list[[i]] <-
    master_vcf_list[[i]][[3]]
}
save(ASoC_SNP_list,
     file = "ASoC_SNP_list.RData")
