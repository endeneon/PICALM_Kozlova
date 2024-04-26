# 25 May 2022 Siwei

# test vcfR package

# init
library(readr)

library(plyr)
library(dplyr)
library(stringr)

library(Rfast)

library(ggplot2)
library(RColorBrewer)

# init
library(readr)

library(plyr)
library(dplyr)
library(stringr)

library(Rfast)
library(vcfR)

library(ggplot2)
library(RColorBrewer)

##

vcf_file <- "two_batchs_merge/temp_merge/processed/GABA_merged_12_lines_SNP_24May2022.vcf"
dna_file <- "/data/Databases/Genomes/hg38/INDEX/GRCh38.primary_assembly.genome.fa"
gff_file <- "/data/Databases/Genomes/hg38/gencode.v35.annotation.sorted.gtf"

vcf <- read.vcfR(vcf_file)
dna <- ape::read.dna(dna_file,
                     format = "fasta")
gff <- read.table(gff_file,
                  sep = "\t")

chrom_hg38_gencodev35 <- create.chromR(name = 'GRCh38_gencode_v35', 
                                       vcf = vcf, 
                                       seq = dna, 
                                       ann = gff)

AD <- extract.gt(vcf,
                 element = "AD",
                 as.numeric = F, # if as.numeric = T, only the 1st AD count will be extracted
                 extract = T)

AD_test <- as.matrix(AD[1:20, 1:14])
sum(is.na(AD_test))
AD_test[is.na(AD_test)] <- "0,0"

AD_REF_ALT_list <- vector(mode = "list",
                          length = 2L)
AD_REF_ALT_list[[1]] <-
  str_split(AD_test,
            pattern = ",",
            n = 2,
            simplify = T)[, 1]
AD_REF_ALT_list[[1]] <-
  matrix(as.numeric(AD_REF_ALT_list[[1]]),
         ncol = 14)

AD_REF_ALT_list[[2]] <-
  str_split(AD_test,
            pattern = ",",
            n = 2,
            simplify = T)[, 2]
AD_REF_ALT_list[[2]] <-
  matrix(as.numeric(AD_REF_ALT_list[[2]]),
         ncol = 14)

AD_REF_ALT_list[[3]] <- 
  data.frame(rep_len(0, length.out = nrow(AD_REF_ALT_list[[1]])))

AD_REF_ALT_list[[3]]$REF <- rowSums(AD_REF_ALT_list[[1]])
AD_REF_ALT_list[[3]]$ALT <- rowSums(AD_REF_ALT_list[[2]])
AD_REF_ALT_list[[3]]$DP <- sum(AD_REF_ALT_list[[3]]$REF,
                               AD_REF_ALT_list[[3]]$ALT)

