# Siwei 29 Jan 2025
# plot new lipidomic volcano plot

# init ####
{
  library(readxl)

  library(stringr)
  library(ggplot2)

  library(ggrepel)

  library(scales)
  library(reshape2)

  library(RColorBrewer)
  # library(cm.)
  library(ggpubr)

  library(dplyr)
  library(data.table)

  library(limma)

  library(DescTools)
  library(multcomp)

  library(gridExtra)

  library(gplots)
}

## make columns ####
df_raw <-
  read_excel("Batch_2_of_new_plots/lipid_volcano_plot.xlsx",
             sheet = 1)
df_raw$safe_rownames <-
  make.names(df_raw$LipidID, unique = T)
sum(duplicated(df_raw$safe_rownames))

df_transform <-
  as.data.frame(df_raw)
df_transform$safe_rownames <- NULL
df_transform$LipidID <- NULL
df_transform <-
  as.data.frame(df_transform)
rownames(df_transform) <-
  df_raw$safe_rownames


sum(is.na(rowSums(df_transform)))

df_transform$log2FC <-
  unlist(apply(X = df_transform,
               MARGIN = 1,
               function(x) {
                 log2FC <-
                   log2(mean(x[5:8], na.rm = T) / mean(x[1:4], na.rm = T))
                 return(log2FC)
               }))
sum(is.na(df_transform$log2FC))

# summary(aov(yield ~ block + N*P*K, npk))
# dunnett_output <-
#   as.data.frame(summary(multcomp::glht(aov(yield ~ block + N*P*K,
#                                            npk),
#                                        linfct = mcp(block = "Dunnett"))))
# dunnett_output$test$pvalues

df_transform$dunnett_P_risk_vs_nonrisk <-
  apply(X = df_transform,
        MARGIN = 1,
        function(x) {
          df_4_aov <-
            data.frame(value = x[1:12],
                       genotype = rep(c("non-risk",
                                        "risk",
                                        "risk_OE"),
                                      each = 4))
          df_4_aov$genotype <-
            factor(df_4_aov$genotype,
                   levels = c("non-risk",
                              "risk",
                              "risk_OE"))
          df_4_aov$genotype <-
            relevel(df_4_aov$genotype,
                    ref = "non-risk")
          glht_output <-
            summary(multcomp::glht(aov(value ~ genotype,
                                       data = df_4_aov),
                                   linfct = mcp(genotype = "Dunnett")))
          return(glht_output$test$pvalues[1])
        })

# rownames(df_transform) <-
#            df_raw$safe_rownames


df_transform$FDR <-
  p.adjust(df_transform$dunnett_P_risk_vs_nonrisk,
           method = "fdr")

df_transform <-
  as.data.frame(cbind(df_raw$LipidID,
                      df_transform))
colnames(df_transform)[1] <- "LipidID"
# df_transform$LipidID <-
#   df_raw$LipidID
# df_transform$safe_id <-
#   df_raw$safe_rownames
# df_transform <-
#   df_transform[order(df_transform$significance), ]


write.table(df_transform,
            file = "lipidomics_ANOVA_Dunnett_P_FDR_risk_vs_nonrisk.tsv",
            row.names = F, col.names = T,
            sep = "\t")

df_transform$significance <-
  "FDR non-significant"
df_transform$significance[df_transform$FDR < 0.05] <-
  "FDR significant"

df_transform$significance <-
  factor(df_transform$significance)
df_transform$significance <-
  relevel(df_transform$significance,
          ref = "FDR significant")


df_transform$Lipid <- df_transform$LipidID
df_transform$Lipid[!(df_transform$FDR < 0.05)] <- ""


df_2_plot <-
  df_transform[order(df_transform$dunnett_P_risk_vs_nonrisk), ]
# df_transform$Lipid <-
#   df_2_plot$safe_id

# df_2_plot$Lipid[!(df_transform$FDR < 0.05)] <- ""
df_2_plot$Lipid[21:nrow(df_2_plot)] <- ""

ggplot(df_2_plot,
       aes(x = log2FC,
           y = 0 - log10(dunnett_P_risk_vs_nonrisk),
           colour = significance,
           # fill = significance
           label = Lipid)) +
  geom_point(shape = 20) +
  scale_colour_manual(values = c("darkred",
                                 "darkblue")) +
  ggrepel::geom_text_repel(seed = 42,
                           max.overlaps = 20,
                           force = 100,
                           min.segment.length = 0,
                           show.legend = F,
                           box.padding = 0) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 7)) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-4, 4)) +
  theme_classic() +
  ggtitle("top 20 lipids") +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"))
# x <- unlist(df_transform[1, 1:12])
# df_4_aov <-
#   data.frame(value = x,
#              genotype = rep(c("non-risk",
#                               "risk",
#                               "risk_OE"),
#                             each = 4))
# df_4_aov$genotype <-
#   factor(df_4_aov$genotype,
#          levels = c("non-risk",
#                     "risk",
#                     "risk_OE"))
# df_4_aov$genotype <-
#   relevel(df_4_aov$genotype,
#           ref = "non-risk")
# glht_output <-
#   summary(multcomp::glht(aov(value ~ genotype,
#                              data = df_4_aov),
#                          linfct = mcp(genotype = "Dunnett")))
# glht_output
# glht_output$test$pvalues[1]
#
# data.frame(genotype = rep(c("non-risk",
#                             "risk",
#                             "risk_OE"),
#                           each = 4))

df_z_score <-
  as.data.frame(df_2_plot[(df_transform$FDR < 0.05), 2:13])

saved_lipidID <-
  df_z_score$LipidID

# df_z_score$LipidID <- NULL
# rownames(df_z_score) <-
#   df_2_plot$safe_id[1:20]

df_z_output <-
  apply(X = df_z_score,
        MARGIN = 1,
        FUN = function(x) {
          z_scored_x <-
            (x - mean(x, na.rm = T)) / sd(x,
                                          na.rm = T)
          return(z_scored_x)
        }, simplify = T)
df_z_output_transposed <-
  as.data.frame(t(df_z_output))

# df_z_output <-
#   as.matrix(t(df_z_output))
heatmap_col <-
  rev(colorRampPalette(colors = c("purple2",
                                  "white",
                                  "cyan"))(20))


heatmap.2(as.matrix(df_z_output_transposed))

heatmap.2(as.matrix(df_z_output_transposed),
          Rowv = F,
          labRow = saved_lipidID,
          Colv = F,
          breaks = 21,
          dendrogram = "both",
          # breaks =
          col = rev(colorRampPalette(colors = c("purple2",
          "white",
          "cyan"))(20)),
          cellnote = format(as.matrix(df_z_output_transposed),
                            digits = 1),
          notecol = "black",
          density.info = "none",
          tracecol = "none",
          main = "risk_vs_nonrisk",
          srtCol = (360 - 45),
          trace = "none")


## calc risk-CRISPRa vs risk ####
df_raw <-
  read_excel("Batch_2_of_new_plots/lipid_volcano_plot.xlsx",
             sheet = 1)
df_raw$safe_rownames <-
  make.names(df_raw$LipidID, unique = T)
sum(duplicated(df_raw$safe_rownames))

df_transform <-
  df_raw
df_transform$safe_rownames <- NULL
df_transform$LipidID <- NULL
df_transform <-
  as.data.frame(df_transform)
rownames(df_transform) <-
  df_raw$safe_rownames

sum(is.na(rowSums(df_transform)))

df_transform$log2FC <-
  unlist(apply(X = df_transform,
               MARGIN = 1,
               function(x) {
                 log2FC <-
                   log2(mean(x[9:12], na.rm = T) / mean(x[5:8], na.rm = T))
                 return(log2FC)
               }))
sum(is.na(df_transform$log2FC))

# summary(aov(yield ~ block + N*P*K, npk))
# dunnett_output <-
#   as.data.frame(summary(multcomp::glht(aov(yield ~ block + N*P*K,
#                                            npk),
#                                        linfct = mcp(block = "Dunnett"))))
# dunnett_output$test$pvalues

# df_4_aov <-
#   data.frame(value = unlist(df_transform[1, 1:12]),
#              genotype = rep(c("non-risk",
#                               "risk",
#                               "risk_OE"),
#                             each = 4))
# df_4_aov$genotype <-
#   factor(df_4_aov$genotype,
#          levels = c("risk_OE",
#                     "risk",
#                     "non-risk"))
# df_4_aov$genotype <-
#   relevel(df_4_aov$genotype,
#           ref = "non-risk")
# glht_output <-
#   summary(multcomp::glht(aov(value ~ genotype,
#                              data = df_4_aov),
#                          linfct = mcp(genotype = "Dunnett")))

df_transform$dunnett_P_risk_OE_vs_risk <-
  apply(X = df_transform,
        MARGIN = 1,
        function(x) {
          df_4_aov <-
            data.frame(value = x[1:12],
                       genotype = rep(c("non-risk",
                                        "risk",
                                        "risk_OE"),
                                      each = 4))
          df_4_aov$genotype <-
            factor(df_4_aov$genotype,
                   levels = c("risk_OE",
                              "risk",
                              "non-risk"))
          df_4_aov$genotype <-
            relevel(df_4_aov$genotype,
                    ref = "risk")
          glht_output <-
            summary(multcomp::glht(aov(value ~ genotype,
                                       data = df_4_aov),
                                   linfct = mcp(genotype = "Dunnett")))
          return(glht_output$test$pvalues[1])
        })

df_transform$FDR <-
  p.adjust(df_transform$dunnett_P_risk_OE_vs_risk,
           method = "fdr")

df_transform <-
  as.data.frame(cbind(df_raw$LipidID,
                      df_transform))
colnames(df_transform)[1] <- "LipidID"
# df_transform$LipidID <-
#   df_raw$LipidID
# df_transform$safe_id <-
#   df_raw$safe_rownames
df_transform <-
  df_transform[order(df_transform$dunnett_P_risk_OE_vs_risk), ]


write.table(df_transform,
            file = "lipidomics_ANOVA_Dunnett_P_FDR_risk_OE_vs_risk.tsv",
            row.names = F, col.names = T,
            sep = "\t")
