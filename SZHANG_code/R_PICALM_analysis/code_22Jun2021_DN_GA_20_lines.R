# Siwei 22 Jun 2021
# ASoC analysis of new processed GA and DN lines ATAC data

# init
library(readr)

library(plyr)
library(dplyr)
library(stringr)

library(Rfast)

library(ggplot2)
library(RColorBrewer)

# load data
ASoC_df_raw <- 
  read_table2("DP20_data_files/DN_20_lines_merged_SNP_DP_20_10Jun2021.txt_4_R.txt")
ASoC_df_raw$X47 <- NULL



# data cleanup
## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]

# View(ASoC_df[1:10, ])
# View(ASoC_df_raw[ASoC_df_raw$ID == 'rs10792832', ])
# View(ASoC_df_raw[ASoC_df_raw$ID == 'rs61904265', ])
# 
# View(ASoC_df_raw[ASoC_df_raw$ID %in% c("rs10792832", "rs61904265"), ])

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 7:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, 40, 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, 40, 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 20
ASoC_df <- ASoC_df[!(ASoC_df$DP < 20), ]

## remove samples that have either REF_N < 2 or ALT_N < 2
ASoC_df <- ASoC_df[((ASoC_df$REF_N > 1) & (ASoC_df$ALT_N > 1)), ]

## calculate binomial p value (use mapply)
ASoC_binom_output <- 
  mapply(binom.test, ASoC_df$REF_N, ASoC_df$DP,
         MoreArgs = list(alternative = "t"),
         SIMPLIFY = F)

## extract the p values use lapply (p.value is the No.3 element)
ASoC_df$binom.P <- unlist(lapply(ASoC_binom_output, `[[`, 3))
ASoC_df$FDR <- p.adjust(ASoC_df$binom.P, method = "fdr")
ASoC_df$logFDR <- 0 - log10(ASoC_df$FDR)
ASoC_df$ratio <- ASoC_df$REF_N / ASoC_df$DP

## sort ASoC df by logFDR value
ASoC_df <- ASoC_df[order(ASoC_df$logFDR, decreasing = T), ]
ASoC_df_table <- ASoC_df[, -c(7:46)]

## write out the table
write.table(ASoC_df_table, file = "DP20_data_files/DN_20_lines_SNP_174K_22Jun2021.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

## generate ASoC (FDR < 0.05) only files
DN_20_lines_SNP_174K_22Jun2021 <-
  DN_20_lines_SNP_174K_22Jun2021[DN_20_lines_SNP_174K_22Jun2021$FDR < 0.05, ]

write.table(DN_20_lines_SNP_174K_22Jun2021, 
            file = "DN_20_lines_SNP_ASoC_only.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

GA_20_lines_SNP_218K_22Jun2021 <-
  GA_20_lines_SNP_218K_22Jun2021[GA_20_lines_SNP_218K_22Jun2021$FDR < 0.05, ]

write.table(GA_20_lines_SNP_218K_22Jun2021, 
            file = "GA_20_lines_SNP_ASoC_only.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

MG_17_lines_SNP_276K_22Jun2021 <-
  MG_17_lines_SNP_276K_22Jun2021[MG_17_lines_SNP_276K_22Jun2021$FDR < 0.05, ]

write.table(MG_17_lines_SNP_276K_22Jun2021, 
            file = "MG_17_lines_SNP_ASoC_only.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

#########
# load data
ASoC_df_raw <- 
  read_table2("DP20_data_files/GA_20_lines_merged_SNP_01Jun2021.txt_4_R.txt")
ASoC_df_raw$X41 <- NULL



# data cleanup
## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]

View(ASoC_df[1:10, ])
# View(ASoC_df_raw[ASoC_df_raw$ID == 'rs10792832', ])
# View(ASoC_df_raw[ASoC_df_raw$ID == 'rs61904265', ])
# 
# View(ASoC_df_raw[ASoC_df_raw$ID %in% c("rs10792832", "rs61904265"), ])

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 7:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 20
ASoC_df <- ASoC_df[!(ASoC_df$DP < 20), ]

## remove samples that have either REF_N < 2 or ALT_N < 2
ASoC_df <- ASoC_df[((ASoC_df$REF_N > 1) & (ASoC_df$ALT_N > 1)), ]

## calculate binomial p value (use mapply)
ASoC_binom_output <- 
  mapply(binom.test, ASoC_df$REF_N, ASoC_df$DP,
         MoreArgs = list(alternative = "t"),
         SIMPLIFY = F)

## extract the p values use lapply (p.value is the No.3 element)
ASoC_df$binom.P <- unlist(lapply(ASoC_binom_output, `[[`, 3))
ASoC_df$FDR <- p.adjust(ASoC_df$binom.P, method = "fdr")
ASoC_df$logFDR <- 0 - log10(ASoC_df$FDR)
ASoC_df$ratio <- ASoC_df$REF_N / ASoC_df$DP

## sort ASoC df by logFDR value
ASoC_df <- ASoC_df[order(ASoC_df$logFDR, decreasing = T), ]
ASoC_df_table <- ASoC_df[, -c(7:40)]

## write out the table
write.table(ASoC_df_table, file = "DP20_data_files/GA_20_lines_SNP_218K_22Jun2021.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)
