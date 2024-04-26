# 25 May 2022 Siwei

# process new GABA  batch data (May 2022)
# note that this batch of GABA bulk ATAC-seq data, 
# two of the samples has lower quality
# since their original samples had been frozen-thawed.

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

library(ggplot2)
library(RColorBrewer)

# load data
ASoC_df_raw <- 
  read_delim("DP20_data_files/non_500bp_intersected/GABA_34_lines_merged_SNP_in_peaks_03Jun2022_DP20_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

# data cleanup
## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]

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

sum(ASoC_df$binom.P < 0.05)
sum(ASoC_df$FDR < 0.05)

## sort ASoC df by logFDR value
ASoC_df <- ASoC_df[order(ASoC_df$logFDR, decreasing = T), ]
ASoC_df_table <- ASoC_df[, -c(7:68)]


## write out the table
write.table(ASoC_df_table, file = "DP20_data_files/GABA_34_lines_SNP_307K_06Jun2022.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


## make simple plot, remove all p = 1 points
ASoC_df_plot <- ASoC_df[ASoC_df$logFDR > 0, ]

ggplot(ASoC_df_plot, aes(x = ratio,
                         y = logFDR,
                         color = ifelse(FDR < 0.05,
                                        "red",
                                        "black"))) +
  ylim(0, 60) +
  scale_color_manual(values = c("black", "red")) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none")


ggplot(ASoC_df_plot, aes(x = ratio,
                         y = logFDR)) +
  geom_point(size = 0.5) +
  theme_classic()

# remove all variants with FDR > 0.05
ASoC_df_table <- 
  ASoC_df_table[ASoC_df_table$FDR < 0.05, ]

ASoC_df_table$CHROM <- 
  str_remove(ASoC_df_table$CHROM, "chr")
ASoC_df_table <- ASoC_df_table[, c(1,2,2,4,5,3,10,11)]

## write out table for annovar
write.table(ASoC_df_table, 
            file = "GABA_34_lines_SNP_500bp_FDR005_annovar_06Jun2022.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


######### GABA_2_batches_merged


# load data
ASoC_df_raw <- 
  read_delim("two_batchs_merge/GABA_34_lines_merged_SNP_24May2022_DP20_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

# data cleanup
## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]

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
ASoC_df_table <- ASoC_df[, -c(7:68)]


## write out the table
write.table(ASoC_df_table, file = "DP20_data_files/GABA_34_lines_SNP_668K_26May2022.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)


## make simple plot, remove all p = 1 points
ASoC_df_plot <- ASoC_df[ASoC_df$logFDR > 0, ]

ggplot(ASoC_df_plot, aes(x = ratio,
                         y = logFDR)) +
  geom_point(size = 0.5) +
  theme_classic()

# remove all variants with FDR > 0.05
ASoC_df_table <- 
  ASoC_df_table[ASoC_df_table$FDR < 0.05, ]

ASoC_df_table$CHROM <- 
  str_remove(ASoC_df_table$CHROM, "chr")
ASoC_df_table <- ASoC_df_table[, c(1,2,2,4,5,3,10,11)]

## write out table for annovar
write.table(ASoC_df_table, 
            file = "GABA_34_lines_SNP_FDR005_annovar_26May2022.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)
