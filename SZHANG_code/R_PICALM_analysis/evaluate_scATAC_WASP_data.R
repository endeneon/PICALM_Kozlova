# 13 Jun 2022 Siwei
# Evaluate WASP calibration results

# init
library(readr)

library(plyr)
library(dplyr)
library(stringr)

library(Rfast)

library(ggplot2)
library(RColorBrewer)


# load data
# ASoC_df_raw <- 
#   read_delim("~/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/CD_50/WASPed_BAMs/vcf_output/het_vcf/output/MG_17_lines_merged_SNP_05May2022_DP15_4_R.txt", 
#              delim = "\t", escape_double = FALSE, 
#              trim_ws = TRUE)
ASoC_df_raw <- 
  read_delim("~/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_0/WASPed_BAMs/vcf_output/hr_0_GABA/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_10x_PBMC_10K/line_51/WASPed_BAMs/vcf_output/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)


# data cleanup

ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_10x_PBMC_10K/line_51/pre_fastq_bams/WASPed_BAMs/vcf_output/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_10x_PBMC_10K/hg38_line_51/WASPed_BAMs/vcf_output/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)


ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_10x_PBMC_10K/hg38_line_51/fastq_out_0hr/pre_bams/WASPed_BAMs/vcf_output/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)


ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_10x_PBMC_10K/hg38_line_51/fastq_out_0hr/pre_bams/WASPed_BAMs/vcf_output/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

# unfiltered
# use both het and alt (0/1 and 1/1 sites) for WASP
# !!! This one works better !!!
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/CD_35_compare_calib/align_results/6hr_unfiltered/vcf_output/output/unfiltered_CD_35_6r_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# use het site for WASP only
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/CD_35_compare_calib/align_results/6hr_cleaned/vcf_output/output/cleaned_CD_35_6hr_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]

# ASoC_df <- ASoC_df_raw
# remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 7:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 7
ASoC_df <- ASoC_df[!(ASoC_df$DP < 10), ]

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





ASoC_df_table <- ASoC_df[, -c(7:ncol(ASoC_df) - 7)]


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
  theme(legend.position = "none") #+
  ggtitle(paste("18 lines, GABA 0hr, 10732 SNPs, DP >= 10, \nmax DP=1879, median DP=17, FDR<0.05=31\n",
                "Nratio < 0.5 = 7549, Nratio >= 0.5 = 3183"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio > 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)

hist(ASoC_df_plot$DP, breaks = 400, xlim = c(0, 400))


### NGN2_Glut_NEFM_pos
ASoC_df_raw <- 
  read_delim("~/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_0/WASPed_BAMs/vcf_output/hr_0_NEFM_pos/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
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
# remove DP < 9
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
ASoC_df$FDR <- p.adjust(ASoC_df$binom.P, method = "fdr")
ASoC_df$logFDR <- 0 - log10(ASoC_df$FDR)
ASoC_df$ratio <- ASoC_df$REF_N / ASoC_df$DP

sum(ASoC_df$binom.P < 0.05)
sum(ASoC_df$FDR < 0.05)


## sort ASoC df by logFDR value
ASoC_df <- ASoC_df[order(ASoC_df$logFDR, decreasing = T), ]
ASoC_df_table <- ASoC_df[, -c(7:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("18 lines, NEFM+ 0hr, 14076 SNPs, DP >= 10, \nmax DP=1550, median DP=18, FDR<0.05=19\n",
                "Nratio < 0.5 = 10576, Nratio >= 0.5 = 3500"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum(ASoC_df_plot$DP)/nrow(ASoC_df_plot)


hist(ASoC_df_plot$DP, breaks = 400, xlim = c(0, 400))

### NGN2_Glut_NEFM_neg
ASoC_df_raw <- 
  read_delim("~/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_0/WASPed_BAMs/vcf_output/hr_0_NEFM_neg/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
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
# remove DP < 7
ASoC_df <- ASoC_df[!(ASoC_df$DP < 7), ]

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
ASoC_df$FDR <- p.adjust(ASoC_df$binom.P, method = "fdr")
ASoC_df$logFDR <- 0 - log10(ASoC_df$FDR)
ASoC_df$ratio <- ASoC_df$REF_N / ASoC_df$DP

sum(ASoC_df$binom.P < 0.05)
sum(ASoC_df$FDR < 0.05)


## sort ASoC df by logFDR value
ASoC_df <- ASoC_df[order(ASoC_df$logFDR, decreasing = T), ]
ASoC_df_table <- ASoC_df[, -c(7:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("18 lines, NEFM- 0hr, 8239 SNPs, DP >= 10, \nmax DP=1506, median DP=18, FDR<0.05=15\n",
                "Nratio < 0.5 = 6127, Nratio >= 0.5 = 2112"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
sum(ASoC_df_plot$DP)/nrow(ASoC_df_plot)
median(ASoC_df_plot$DP)
hist(ASoC_df_plot$DP, breaks = 400, xlim = c(0, 400))


############ 1hr
ASoC_df_raw <- 
  read_delim("~/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_1/WASPed_BAMs/vcf_output/hr_1_GABA/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
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
# remove DP < 7
ASoC_df <- ASoC_df[!(ASoC_df$DP < 7), ]

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
ASoC_df_table <- ASoC_df[, -c(7:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("18 lines, GABA 1hr, 13150 SNPs, DP >= 10, \nmax DP=1952, median DP=18, FDR<0.05=34\n",
                "Nratio < 0.5 = 9298, Nratio >= 0.5 = 3852"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)

hist(ASoC_df_plot$DP, breaks = 400, xlim = c(0, 400))


### NGN2_Glut_NEFM_pos
ASoC_df_raw <- 
  read_delim("~/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_1/WASPed_BAMs/vcf_output/hr_1_NEFM_pos/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
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
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df$FDR <- p.adjust(ASoC_df$binom.P, method = "fdr")
ASoC_df$logFDR <- 0 - log10(ASoC_df$FDR)
ASoC_df$ratio <- ASoC_df$REF_N / ASoC_df$DP

sum(ASoC_df$binom.P < 0.05)
sum(ASoC_df$FDR < 0.05)


## sort ASoC df by logFDR value
ASoC_df <- ASoC_df[order(ASoC_df$logFDR, decreasing = T), ]
ASoC_df_table <- ASoC_df[, -c(7:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("18 lines, NEFM+ 6hr, 20027 SNPs, DP >= 10, \nmax DP=1770, median DP=20, FDR<0.05=24\n",
                "Nratio < 0.5 = 14932, Nratio >= 0.5 = 5095"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.3)
sum(ASoC_df_plot$ratio > 0.7)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum(ASoC_df_plot$DP)/nrow(ASoC_df_plot)


hist(ASoC_df_plot$DP, breaks = 400, xlim = c(0, 400))

### NGN2_Glut_NEFM_neg
ASoC_df_raw <- 
  read_delim("~/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_1/WASPed_BAMs/vcf_output/hr_1_NEFM_neg/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
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
# remove DP < 7
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df$FDR <- p.adjust(ASoC_df$binom.P, method = "fdr")
ASoC_df$logFDR <- 0 - log10(ASoC_df$FDR)
ASoC_df$ratio <- ASoC_df$REF_N / ASoC_df$DP

sum(ASoC_df$binom.P < 0.05)
sum(ASoC_df$FDR < 0.05)


## sort ASoC df by logFDR value
ASoC_df <- ASoC_df[order(ASoC_df$logFDR, decreasing = T), ]
ASoC_df_table <- ASoC_df[, -c(7:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("18 lines, NEFM- 6hr, 13302 SNPs, DP >= 10, \nmax DP=1692, median DP=18, FDR<0.05=17\n",
                "Nratio < 0.5 = 9960, Nratio >= 0.5 = 3342"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
sum(ASoC_df_plot$DP)/nrow(ASoC_df_plot)
median(ASoC_df_plot$DP)

sum(ASoC_df_plot$DP)

hist(ASoC_df_plot$DP, breaks = 400, xlim = c(0, 400))





############ 6hr
ASoC_df_raw <- 
  read_delim("~/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_6/WASPed_BAMs/vcf_output/hr_6_GABA/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
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
# remove DP < 7
ASoC_df <- ASoC_df[!(ASoC_df$DP < 7), ]

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
ASoC_df_table <- ASoC_df[, -c(7:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("18 lines, GABA 6hr, 16163 SNPs, DP >= 10, \nmax DP=1909, median DP=19, FDR<0.05=34\n",
                "Nratio < 0.5 = 11525, Nratio >= 0.5 = 4638"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)

hist(ASoC_df_plot$DP, breaks = 400, xlim = c(0, 400))


### NGN2_Glut_NEFM_pos
ASoC_df_raw <- 
  read_delim("~/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_6/WASPed_BAMs/vcf_output/hr_6_NEFM_pos/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
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
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df$FDR <- p.adjust(ASoC_df$binom.P, method = "fdr")
ASoC_df$logFDR <- 0 - log10(ASoC_df$FDR)
ASoC_df$ratio <- ASoC_df$REF_N / ASoC_df$DP

sum(ASoC_df$binom.P < 0.05)
sum(ASoC_df$FDR < 0.05)


## sort ASoC df by logFDR value
ASoC_df <- ASoC_df[order(ASoC_df$logFDR, decreasing = T), ]
ASoC_df_table <- ASoC_df[, -c(7:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("18 lines, NEFM+ 6hr, 17593 SNPs, DP >= 10, \nmax DP=1476, median DP=19, FDR<0.05=24\n",
                "Nratio < 0.5 = 4694, Nratio >= 0.5 = 1476"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum(ASoC_df_plot$DP)/nrow(ASoC_df_plot)


hist(ASoC_df_plot$DP, breaks = 400, xlim = c(0, 400))

### NGN2_Glut_NEFM_neg
ASoC_df_raw <- 
  read_delim("~/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/hr_6/WASPed_BAMs/vcf_output/hr_6_NEFM_neg/het_vcf/output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt", 
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
# remove DP < 7
ASoC_df <- ASoC_df[!(ASoC_df$DP < 7), ]

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
ASoC_df$FDR <- p.adjust(ASoC_df$binom.P, method = "fdr")
ASoC_df$logFDR <- 0 - log10(ASoC_df$FDR)
ASoC_df$ratio <- ASoC_df$REF_N / ASoC_df$DP

sum(ASoC_df$binom.P < 0.05)
sum(ASoC_df$FDR < 0.05)


## sort ASoC df by logFDR value
ASoC_df <- ASoC_df[order(ASoC_df$logFDR, decreasing = T), ]
ASoC_df_table <- ASoC_df[, -c(7:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("18 lines, NEFM- 6hr, 13377 SNPs, DP >= 10, \nmax DP=1476, median DP=19, FDR<0.05=29\n",
                "Nratio < 0.5 = 9623, Nratio >= 0.5 = 1476"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
sum(ASoC_df_plot$DP)/nrow(ASoC_df_plot)
median(ASoC_df_plot$DP)


hist(ASoC_df_plot$DP, breaks = 400, xlim = c(0, 400))


















#################
## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/GABA/output/GABA__0hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, GABA 0hr, 25046 SNPs, DP >= 10, \nmax DP=466, median DP = 17, FDR<0.05=62\n",
                "Nratio < 0.5 = 12026, Nratio >= 0.5 = 13020"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum((ASoC_df_plot$FDR < 0.05))


## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/NEFMp_CUX2m/output/NEFMp__0hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, NEFM+ 0hr, 57734 SNPs, DP >= 10, \nmax DP=580, median DP=15. FDR<0.05=182\n",
                "Nratio < 0.5 = 28996, Nratio >= 0.5 = 28738"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
sum((ASoC_df_plot$FDR < 0.05))
median(ASoC_df_plot$DP)
sum(ASoC_df_plot$DP)/nrow(ASoC_df_plot)


hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))


## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/NEFMm_CUX2p/output/NEMFm__0hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, NEFM- 0hr, 6704 SNPs, DP >= 10, \nmax DP=162, median DP=15, FDR<0.05=1\n",
                "Nratio < 0.5 = 3149, Nratio >= 0.5 = 3555"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum((ASoC_df_plot$FDR < 0.05))
hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))



## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/GABA/output/GABA__1hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, GABA 1hr, 22447 SNPs, DP >= 10, \nmax DP=337, median DP=16, FDR<0.05=25\n",
                "Nratio < 0.5 = 10637, Nratio >= 0.5 = 11810"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum((ASoC_df_plot$FDR < 0.05))
hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))


## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/NEFMp_CUX2m/output/NEFMp__1hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, NEFM+ 1hr, 45439 SNPs, DP >= 10, \nmax DP=431, FDR<0.05=70\n",
                "Nratio < 0.5 = 28996, Nratio >= 0.5 = 28738"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
sum(ASoC_df$FDR < 0.05)
sum((ASoC_df_plot$FDR < 0.05))
hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))


## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/NEFMm_CUX2p/output/NEMFm__0hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, NEFM- 0hr, 6704 SNPs, DP >= 10, \nmax DP=162, FDR<0.05=1\n",
                "Nratio < 0.5 = 3149, Nratio >= 0.5 = 3555"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
sum((ASoC_df_plot$FDR < 0.05))
hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))

## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/GABA/output/GABA__1hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, GABA 1hr, 22447 SNPs, DP >= 10, \nmax DP=337, median DP=16, FDR<0.05=25\n",
                "Nratio < 0.5 = 10637, Nratio >= 0.5 = 11810"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum((ASoC_df_plot$FDR < 0.05))
hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))


## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/NEFMp_CUX2m/output/NEFMp__1hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, NEFM+ 1hr, 45439 SNPs, DP >= 10, \nmax DP=431, median DP = 15, FDR<0.05=70\n",
                "Nratio < 0.5 = 22183, Nratio >= 0.5 = 23256"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum(ASoC_df$FDR < 0.05)
sum((ASoC_df_plot$FDR < 0.05))
sum(ASoC_df$FDR < 0.05)
hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))


## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/NEFMm_CUX2p/output/NEMFm__1hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, NEFM- 1hr, 7262 SNPs, DP >= 10, \nmax DP=133, FDR<0.05=1\n",
                "Nratio < 0.5 = 3356, Nratio >= 0.5 = 3906"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum((ASoC_df_plot$FDR < 0.05))
hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))

## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/GABA/output/GABA__6hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, GABA 6hr, 18696 SNPs, DP >= 10, \nmax DP=162, median DP=16, FDR<0.05=27\n",
                "Nratio < 0.5 = 8825, Nratio >= 0.5 = 9871"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum((ASoC_df_plot$FDR < 0.05))
hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))


## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/NEFMp_CUX2m/output/NEFMp__6hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


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
  ggtitle(paste("5 lines, NEFM+ 6hr, 43331 SNPs, DP >= 10, \nmax DP=455, median DP = 16, FDR<0.05=72\n",
                "Nratio < 0.5 = 21128, Nratio >= 0.5 = 22203"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum(ASoC_df$FDR < 0.05)
sum((ASoC_df_plot$FDR < 0.05))
sum(ASoC_df$FDR < 0.05)
hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))


## results from the 5 lines
ASoC_df_raw <- 
  read_delim("/home/zhangs3/NVME/scARC_Duan_018/Duan_Project_018-ATAC/VCF_15Feb/het_vcf/NEFMm_CUX2p/output/NEMFm__6hr.table", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
# data cleanup

## remove all SNPs without rs#
ASoC_df <- ASoC_df_raw[!(ASoC_df_raw$ID == '.'), ]
## remove the last col
# ASoC_df <- ASoC_df[, -(ncol(ASoC_df))]

# subset per-sample REF/ALT counts
ASoC_per_sample_df <- ASoC_df[, 6:ncol(ASoC_df)]

## calculate REF_N and ALT_N from per-sample df
ASoC_df$REF_N <- 
  rowSums(ASoC_per_sample_df[, seq(1, ncol(ASoC_per_sample_df), 2)], na.rm = T) # odd columns
ASoC_df$ALT_N <- 
  rowSums(ASoC_per_sample_df[, seq(2, ncol(ASoC_per_sample_df), 2)], na.rm = T) # even columns
ASoC_df$DP <- rowSums(ASoC_per_sample_df, na.rm = T)
# remove DP < 9
ASoC_df <- ASoC_df[!(ASoC_df$DP < 9), ]

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
ASoC_df_table <- ASoC_df[, -c(6:ncol(ASoC_df) - 7)]


ASoC_df_plot <- ASoC_df
ASoC_df_plot$logP <- 0 - log(ASoC_df_plot$binom.P)
ggplot(ASoC_df_plot, aes(x = ratio,
                         y = logP,
                         color = ifelse(FDR < 0.05,
                                        "black",
                                        "red"))) +
  # ylim(0, 50) +
  scale_color_manual(values = c("black")) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0.5,
             linetype = 2) +
  # ylim(0, 20) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle(paste("5 lines, NEFM- 6hr, 5396 SNPs, DP >= 10, \nmax DP=115, median DP=14, FDR<0.05=0\n",
                "Nratio < 0.5 = 2536, Nratio >= 0.5 = 2860"))

sum((ASoC_df_plot$ratio < 0.5) & (ASoC_df_plot$FDR < 0.05))
sum((ASoC_df_plot$ratio > 0.5) & (ASoC_df_plot$FDR < 0.05))
sum(ASoC_df_plot$ratio < 0.5)
sum(ASoC_df_plot$ratio >= 0.5)
max(ASoC_df_plot$DP)
median(ASoC_df_plot$DP)
sum((ASoC_df_plot$FDR < 0.05))
hist(ASoC_df_plot$DP, breaks = 100, xlim = c(0, 400))

####################
