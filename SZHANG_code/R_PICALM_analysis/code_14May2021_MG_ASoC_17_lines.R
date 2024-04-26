# Siwei 14 May 2021
# ASoC analysis of 17 Microglia lines ATAC data

# init
library(readr)

library(plyr)
library(dplyr)
library(stringr)

library(Rfast)

library(ggplot2)
library(RColorBrewer)

# load data
ASoC_df_raw <- read_table2("DP20_data_files/non_500bp_intersected/MG_28_lines_merged_peaks_filtered_03Jun2022_DP_20_4_R.txt")
ASoC_df_raw$X41 <- NULL



# data cleanup
## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]

# View(ASoC_df[1:10, ])
View(ASoC_df_raw[ASoC_df_raw$ID == 'rs10792832', ])
View(ASoC_df_raw[ASoC_df_raw$ID == 'rs61904265', ])

View(ASoC_df_raw[ASoC_df_raw$ID %in% c("rs10792832", "rs61904265"), ])

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
write.table(ASoC_df_table, file = "DP20_data_files/MG_17_lines_SNP_276K_22Jun2021.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)



## make simple plot, remove all p = 1 points
ASoC_df_plot <- ASoC_df[ASoC_df$logFDR > 0, ]

ggplot(ASoC_df_plot, aes(x = ratio,
                         y = logFDR)) +
  geom_point(size = 0.5) +
  theme_classic()

ggplot(ASoC_df, aes(x = ratio,
                    y = logFDR)) +
  geom_point(size = 0.5) +
  theme_classic()

## remove all variants with FDR > 0.05
ASoC_df_table <- 
  ASoC_df_table[ASoC_df_table$FDR < 0.05, ]

ASoC_df_table$CHROM <- 
  str_remove(ASoC_df_table$CHROM, "chr")
ASoC_df_table <- ASoC_df_table[, c(1,2,2,4,5,3,10,11)]

## write out table for annovar
write.table(ASoC_df_table, 
            file = "MG_17_lines_SNP_FDR005_annovar_17May2021.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


binom.test(x = sum(60,81,94,72,142), 
           n = sum(rowSums(ASoC_df_raw[ASoC_df_raw$ID == 'rs10792832', 7:ncol(ASoC_df_raw)]), 
                   -136),
           p = 0.5,
           alternative = "t")

binom.test(x = sum(60,81,94,72,142,56), 
           n = rowSums(ASoC_df_raw[ASoC_df_raw$ID == 'rs10792832', 7:ncol(ASoC_df_raw)]),
           p = 0.5,
           alternative = "t")

binom.test(x = sum(60,81,94,72,142,56, 70 * 4), 
           n = sum(rowSums(ASoC_df_raw[ASoC_df_raw$ID == 'rs10792832', 7:ncol(ASoC_df_raw)]),
                   70 * 4 * 1.5),
           p = 0.5,
           alternative = "t")
binom.test(x = 505,
           n = 939,
           p = 0.5,
           alternative = "t")
