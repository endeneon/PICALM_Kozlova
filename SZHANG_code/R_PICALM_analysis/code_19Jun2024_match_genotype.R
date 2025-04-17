# Siwei 19 Jun 2024
# match RNASeq count matrices to individual and genotype data

# install.packages("LDlinkR")

# init ####
{
  library(readr)

  library(edgeR)
  library(ggplot2)
  library(RColorBrewer)

  library(stringr)
}

# load data ####
df_gene_names <-
  read_csv("freshmicro_counts/RNAseq_GeneInfo.csv")

df_raw_counts <-
  read_csv("freshmicro_counts/RNAseq_count_matrix.csv")
colnames(df_raw_counts)[1] <- "Geneid"
colnames(df_raw_counts) <-
  str_split(string = colnames(df_raw_counts),
            pattern = "__",
            simplify = T)[, 1]

df_gene_info <-
  read_csv("freshmicro_counts/RNAseq_GeneInfo.csv")

df_raw_indiv <-
  read_csv("freshmicroglia_indiv/FreshMicro_individual_metadata.csv")

df_geno_info <-
  readRDS("rs10792832_most_correlated_SNPs_hg38.RDs")
df_geno_meta <-
  as.data.frame(t(df_geno_info[, 10:ncol(df_geno_info)]))
colnames(df_geno_meta) <-
  c("rs561655",
    "rs3851179")
df_geno_meta$individualID <-
  rownames(df_geno_meta)

df_indiv_geno_integrated <-
  merge(x = df_raw_indiv,
        y = df_geno_meta,
        by.x = "individualID",
        by.y = "individualID")

df_indiv_geno_integrated$rs10792832 <- "AG"
df_indiv_geno_integrated$rs10792832[df_indiv_geno_integrated$rs3851179 %in% "1|1"] <- "GG"
df_indiv_geno_integrated$rs10792832[df_indiv_geno_integrated$rs3851179 %in% "0|0"] <- "AA"

raw_gene_list <-
  df_raw_counts$Geneid
df_raw_counts$Geneid <- NULL
rownames(df_raw_counts) <-
  raw_gene_list

df_counts_4_DGE <-
  df_raw_counts[, colnames(df_raw_counts) %in% df_indiv_geno_integrated$individualID]
colnames(df_counts_4_DGE)
df_indiv_geno_integrated$individualID

rownames(df_counts_4_DGE) <-
  raw_gene_list
df_counts_4_DGE <-
  df_counts_4_DGE[, order(colnames(df_counts_4_DGE))]
rownames(df_counts_4_DGE) <-
  raw_gene_list

df_design_matrix <-
  df_indiv_geno_integrated[df_indiv_geno_integrated$individualID %in%
                             colnames(df_counts_4_DGE), ]
df_design_matrix <-
  df_design_matrix[order(df_design_matrix$individualID), ]
df_design_matrix <-
  df_design_matrix[!is.na(df_design_matrix$race), ]
df_design_matrix <-
  df_design_matrix[!is.na(df_design_matrix$Braak), ]
df_design_matrix <-
  df_design_matrix[!(df_design_matrix$diagnosis == 'Not Applicable'), ]
df_design_matrix <-
  df_design_matrix[!is.na(df_design_matrix$brainWeight), ]

unique(df_design_matrix$diagnosis)
df_design_matrix$diagnosis[df_design_matrix$diagnosis == 'Alzheimer Disease'] <-
  'Alzheimer_Disease'
df_design_matrix$diagnosis[df_design_matrix$diagnosis == paste0("Parkinson",
                                                                "'",
                                                                "s",
                                                                " ",
                                                                "disease")] <-
  'Parkinson_Disease'
df_design_matrix$diagnosis[df_design_matrix$diagnosis == 'traumatic brain injury'] <-
  "traumatic_brain_injury"

unique(df_design_matrix$race)
df_design_matrix$race[df_design_matrix$race == 'Black or African American'] <-
  "Black_or_African_American"


df_counts_4_DGE_final <-
  df_counts_4_DGE[, colnames(df_counts_4_DGE) %in%
                    df_design_matrix$individualID]
rownames(df_counts_4_DGE_final) <-
  raw_gene_list


df_DGE <-
  DGEList(counts = as.matrix(df_counts_4_DGE_final),
          samples = colnames(df_counts_4_DGE_final),
          genes = rownames(df_counts_4_DGE_final),
          group = df_design_matrix$rs10792832,
          remove.zeros = T)

df_DGE <-
  calcNormFactors(df_DGE)

df_design_matrix$sex <-
  as.factor(df_design_matrix$sex)
df_design_matrix$diagnosis <-
  as.factor(df_design_matrix$diagnosis)
df_design_matrix$diagnosis <-
  relevel(df_design_matrix$diagnosis,
          ref = "control")

df_design_matrix$race <-
  as.factor(df_design_matrix$race)
df_design_matrix$Braak <-
  as.factor(df_design_matrix$Braak)
df_design_matrix$apoeGenotype <-
  as.factor(df_design_matrix$apoeGenotype)
df_design_matrix$rs10792832 <-
  factor(df_design_matrix$rs10792832,
         levels = c("GG",
                    "AG",
                    "AA"))

matrix_design <-
  model.matrix(~ 0+
                 rs10792832 +
                 # sex +
                 # race +
                 # Braak +
                 # brainWeight +
                 # apoeGenotype +
                 diagnosis,
               data = df_design_matrix)
colnames(matrix_design)


normalised_cpm <-
  cpm(df_DGE,
      normalized.lib.sizes = T,
      log = F)

rowSums(normalised_cpm > 1)

df_DGE <-
  df_DGE[rowSums(normalised_cpm > 1) > 37, ]



df_DGE <-
  estimateDisp(df_DGE,
               design = matrix_design)


df_DGE_QLF <-
  glmQLFit(df_DGE,
           design = matrix_design)

colnames(matrix_design)

contrast_geno <-
  makeContrasts(AGvsGG = rs10792832AG - rs10792832GG,
                AAvsGG = rs10792832AA - rs10792832GG,
                Alzvsctrl = diagnosisAlzheimer_Disease - rs10792832GG,
                levels = matrix_design)

qlf_AGvsGG <-
  glmQLFTest(df_DGE_QLF,
             contrast = contrast_geno[, "AGvsGG"])
qlf_AGvsGG$table[rownames(qlf_AGvsGG$table) == "ENSG00000073921", ]


qlf_AAvsGG <-
  glmQLFTest(df_DGE_QLF,
             contrast = contrast_geno[, "AAvsGG"])
qlf_AAvsGG$table[rownames(qlf_AAvsGG$table) == "ENSG00000073921", ]

qlf_AlzvsCtrl <-
  glmQLFTest(df_DGE_QLF,
             contrast = contrast_geno[, "Alzvsctrl"])
qlf_AlzvsCtrl$table[rownames(qlf_AlzvsCtrl$table) == "ENSG00000073921", ]
View(qlf_AlzvsCtrl$table)

normalised_cpm <-
  as.data.frame(normalised_cpm)

df_2_plot <-
  data.frame(Diag = df_design_matrix$diagnosis,
             Geno = df_design_matrix$rs10792832,
             PICALM = unlist(normalised_cpm[rownames(normalised_cpm) == "ENSG00000073921", ]))

df_2_plot$Geno <-
  factor(df_2_plot$Geno,
         levels = c("GG", "AG", "AA"))
class(df_2_plot$Geno)

# G is risk
df_2_plot$Dosage <- 0
df_2_plot$Dosage[df_2_plot$Geno == "AG"] <- 1
df_2_plot$Dosage[df_2_plot$Geno == "GG"] <- 2


glm_test <-
  summary(glm(PICALM ~ Geno,
                   data = df_2_plot))
coef(glm_test)[, 4]
with(glm_test,
     1 - glm_test$deviance / glm_test$null.deviance)

summary(lm(formula = PICALM ~ Geno,
             data = df_2_plot))

summary(lm(formula = PICALM ~ Dosage,
           data = df_2_plot))
# summary(lm_test)

df_2_plot$metaeQTL_levels <-
  "control"
df_2_plot$metaeQTL_levels[df_2_plot$Diag %in% c("Alzheimer_Disease",
                                                "dementia")] <-
  "Alz"

df_2_lm <-
  df_2_plot[df_2_plot$metaeQTL_levels == "Alz", ]
summary(lm(formula = PICALM ~ Geno,
           data = df_2_lm))
summary(lm(formula = PICALM ~ Dosage,
           data = df_2_lm))
summary(lm(formula = log(PICALM) ~ Dosage,
           data = df_2_lm))

glm_test <-
  summary(glm(PICALM ~ Geno,
              data = df_2_lm))
coef(glm_test)[, 4]
with(glm_test,
     1 - glm_test$deviance / glm_test$null.deviance)


df_2_lm <-
  df_2_plot[df_2_plot$metaeQTL_levels == "control", ]
summary(lm(formula = PICALM ~ Geno,
           data = df_2_lm))
summary(lm(formula = PICALM ~ Dosage,
           data = df_2_lm))
summary(lm(formula = log(PICALM) ~ Dosage,
           data = df_2_lm))
glm_test <-
  summary(glm(PICALM ~ Geno,
              data = df_2_lm))
coef(glm_test)[, 4]
with(glm_test,
     1 - glm_test$deviance / glm_test$null.deviance)


write.table(df_2_plot,
            file = "PICALM_CPM_df_2_plot.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

df_2_plot <-
  read.delim(file = "PICALM_CPM_df_2_plot.tsv",
             header = T,
             sep = "\t")



df_2_plot$metaeQTL_levels <-
  as.factor(df_2_plot$metaeQTL_levels)
df_2_plot$metaeQTL_levels <-
  relevel(df_2_plot$metaeQTL_levels,
          ref = "control")

#
# ggplot(df_2_plot,
#        aes(x = metaeQTL_levels,
#            y = PICALM,
#            fill = metaeQTL_levels)) +
#   geom_boxplot(size = 0.5) +
#   geom_jitter(width = 0.1) +
#   stat_summary(fun.y = mean,
#                geom = "crossbar",
#                colour = "red") +
#   ylim(0, 800) +
#   theme_classic()

ggplot(df_2_plot,
       aes(x = Geno,
           y = PICALM,
           fill = Geno)) +
  geom_boxplot(size = 0.5) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.y = mean,
               geom = "crossbar",
               colour = "red") +
  ylim(0, 800) +
  theme_classic() +
  ggtitle("All diags: AD, PD, SCZ, trauma, control")


df_2_plot <-
  df_2_plot[df_2_plot$Diag %in% c("Alzheimer_Disease",
                                  "control"), ]
ggplot(df_2_plot,
       aes(x = Geno,
           y = PICALM,
           fill = Geno)) +
  geom_boxplot(size = 0.5) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.y = mean,
               geom = "crossbar",
               colour = "red") +
  ylim(0, 800) +
  theme_classic() +
  ggtitle("AD and controls only")

df_2_plot <-
  df_2_plot[df_2_plot$Diag %in% c("control",
                                  "Alzheimer_Disease"), ]
df_2_plot$Diag <-
  factor(df_2_plot$Diag,
         levels = c("control",
                    "Alzheimer_Disease"))


ggplot(df_2_plot,
       aes(x = Diag,
           y = PICALM,
           fill = Diag)) +
  geom_boxplot(size = 0.5,
               width = 0.5,
               outliers = F) +
  geom_jitter(width = 0.1,
              size = 0.5) +
  stat_summary(fun.y = mean,
               geom = "crossbar",
               colour = "red") +
  scale_y_continuous(limits = c(0, 800),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#619CFF",
                               "#F8766D")) +
  xlab("") +
  ylab("PICALM exp. level (CPM)") +
  # ylim(0, 800) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        legend.position = "none")


df_2_plot_regressed <-
  df_2_plot
df_2_plot_regressed <-
  df_2_plot_regressed[df_2_plot_regressed$metaeQTL_levels == "control", ]
# df_2_plot_regressed <-
#   df_2_plot_regressed[df_2_plot_regressed$metaeQTL_levels == "Alz", ]

ggplot(df_2_plot_regressed,
       aes(x = Dosage,
           y = PICALM)) +
  geom_boxplot(aes(group = as.factor(Dosage)),
               size = 0.25,
               linewidth = 0.5,
               outlier.size = 0,
               outlier.colour = "white",
               colour = "black") +
  geom_jitter(aes(colour = as.factor(Dosage)),
              width = 0.1,
              size = 0.5) +
  stat_summary(fun.y = mean,
               geom = "crossbar",
               colour = "red",
               linewidth = 0.2) +
  geom_smooth(method = "lm",
              se = F,
              colour = "darkcyan") +
  ylim(0, 800) +
  scale_colour_manual(values = c("darkblue",
                                 "orange3",
                                 "darkred")) +

  xlab("G dosage") +
  ylab("PICALM exp. level (CPM)") +
  # ggtitle("All individuals\n Adj. sqrt(R) = 0.01549;\np-value = 0.1441.") +
  ggtitle("ctrl individuals\nAdj. sqrt(R) = 0.0030;\np-value = 0.8548.") +
  # ggtitle("Alz individuals\nAdj. sqrt(R) = 0.1673;\np-value = 0.0042.") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        title = element_text(size = 7))



t.test(x = df_2_plot$PICALM[df_2_plot$metaeQTL_levels == "control"],
       y = df_2_plot$PICALM[df_2_plot$metaeQTL_levels == "Alz"],
       paired = F,
       var.equal = T)

results_table <-
  qlf_AlzvsCtrl$table
results_table$FDR <-
  p.adjust(results_table$PValue,
           method = "fdr")
# results_table$Geneid <-
