# 13 Jun 2022 Siwei
# Calc Bulk atac-seq of NGN2-Glut

# init
library(readr)

library(plyr)
library(dplyr)
library(stringr)

library(Rfast)

library(ggplot2)
library(RColorBrewer)


### NGN2_Glut_NEFM_pos
ASoC_df_raw <- 
  read_delim("two_batchs_merge/temp_merge/output/GABA_14_lines_merged_SNP_24May2022_DP20_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 7:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 19
ASoC_df <- ASoC_df[!(ASoC_df$DP < 19), ]

# # subset per-sample REF/ALT counts
# ASoC_per_sample_df <- ASoC_df[, 7:(ncol(ASoC_df) - 3)]

## remove samples that have either REF_N < 2 or ALT_N < 2
ASoC_df <- ASoC_df[((ASoC_df$REF_N > 1) & (ASoC_df$ALT_N > 1)), ]

## calculate binomial p value (use mapply)
ASoC_binom_output <- 
  mapply(binom.test, ASoC_df$REF_N, ASoC_df$DP,
         MoreArgs = list(alternative = "t"),
         SIMPLIFY = F)

## extract the p values use lapply (p.value is the No.3 element)
ASoC_df$binom.P <- unlist(lapply(ASoC_binom_output, `[[`, 3))
ASoC_df$logPval <- 0 - log10(ASoC_df$binom.P)
ASoC_df$FDR <- p.adjust(ASoC_df$binom.P, method = "fdr")
ASoC_df$logFDR <- 0 - log10(ASoC_df$FDR)
ASoC_df$ratio <- ASoC_df$REF_N / ASoC_df$DP

sum(ASoC_df$binom.P < 0.05)
sum(ASoC_df$FDR < 0.05)


## sort ASoC df by logFDR value
ASoC_df <- ASoC_df[order(ASoC_df$logFDR, decreasing = T), ]
ASoC_df_writeout <-
  ASoC_df[, c((1:6), (ncol(ASoC_df) - 6):(ncol(ASoC_df) - 1))]
write.table(ASoC_df_writeout,
            file = paste("/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/SNP_18_liners_output_table",
                         "Glut_NEFMp_scATAC_18_lines_0hr_DP20_full_list.txt",
                         sep = '/'),
            row.names = F, col.names = T, 
            sep = "\t", quote = F)
ASoC_df_writeout <-
  ASoC_df_writeout[ASoC_df_writeout$FDR < 0.05, ]
write.table(ASoC_df_writeout,
            file = paste("/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/SNP_18_liners_output_table",
                         "Glut_NEFMp_scATAC_18_lines_0hr_DP20_FDR_005.txt",
                         sep = '/'),
            row.names = F, col.names = T, 
            sep = "\t", quote = F)
ASoC_df_4_annovar <-
  ASoC_df_writeout[, c(1,2,2,4,5,3,6:ncol(ASoC_df_writeout))]
ASoC_df_4_annovar[[2]] <- ASoC_df_4_annovar[[2]] - 1
ASoC_df_4_annovar$CHROM <- 
  str_remove_all(string = ASoC_df_4_annovar$CHROM,
                 pattern = "^chr")
write.table(ASoC_df_4_annovar,
            file = paste("/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/SNP_18_liners_output_table",
                         "Glut_NEFMp_scATAC_18_lines_0hr_DP20_FDR_005_4_annovar.txt",
                         sep = '/'),
            row.names = F, col.names = T, 
            sep = "\t", quote = F)


ASoC_df_plot <- ASoC_df
ASoC_df_plot$logP <- 0 - log(ASoC_df_plot$binom.P)
ggplot(ASoC_df_plot, aes(x = ratio,
                         y = logP,
                         color = ifelse(FDR < 0.05,
                                        "black",
                                        "red"))) +
  # ylim(0, 50) +
  scale_color_manual(values = c("red", "black")) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0.5,
             linetype = 2) +
  # ylim(0, 20) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle(paste("Bulk ATAC-seq, NGN2-Glut\n", 
                "20 lines, 169531 SNPs, DP >= 20, \nmax DP=6220, median DP=46, FDR<0.05=8124\n"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum(ASoC_df_plot$DP)/nrow(ASoC_df_plot)
sum(ASoC_df_plot$FDR < 0.05)

hist(ASoC_df_plot$DP, breaks = 1000, xlim = c(0, 400))


df_sum <- 
  read_delim("/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_0/WASPed_BAMs/test_call_genotypes_NEFMp/read_sum.txt", 
             delim = "\t", escape_double = FALSE, 
             col_names = F,
             trim_ws = TRUE)
colSums(df_sum) # 105980870


(48198416+20500328+18430019+
    18825446+6341530+11965570+
    69865102+27656471+23128443)/670310156

48198416+18825446+69865102
20500328+6341530+27656471
18430019+11965570+23128443


df_sum <- 
  read_delim("/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_0/read_sum.txt", 
             delim = "\t", escape_double = FALSE, 
             col_names = F,
             trim_ws = TRUE)
colSums(df_sum) # 596354224
