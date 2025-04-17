# Siwei 10 Feb 2025
# use 425 samples
# make new plots for Alena's revised paper

# init ####
{
  library(ggplot2)
  library(stringr)
  library(lsr)
  library(RColorBrewer)
  library(ggpubr)

  library(readr)
  library(readxl)

  library(stats)
  library(rstatix)

  library(agricolae)
  library(DescTools)

  # library

  library(data.table)

}

# plot 3x3 cell types ####

df_raw <-
  read_excel("Batch_2_of_new_plots/PICALM-400plus_AD-3-type-8other_cell_types.xlsx")

df_2_plot <-
  df_raw

df_2_plot <-
  df_2_plot[, c(9, 11:18)]


df_2_plot <-
  data.table::melt(df_2_plot,
                   na.rm = T)
df_2_plot <-
  df_2_plot[!is.na(df_2_plot$variable), ]
df_2_plot <-
  df_2_plot[!is.na(df_2_plot$value), ]

dev.off()

# i <- 1

for (i in unique(df_2_plot$variable)) {
  print(i)

  df_2_plot_subsetted <-
    df_2_plot[df_2_plot$variable == i, ]
  df_2_plot_subsetted <-
    df_2_plot_subsetted[!is.na(df_2_plot_subsetted$value), ]
  df_2_plot_subsetted <-
    df_2_plot_subsetted[!is.nan(df_2_plot_subsetted$value), ]
  df_2_plot_subsetted <-
    df_2_plot_subsetted[!is.infinite(df_2_plot_subsetted$value), ]
  df_2_plot_subsetted <-
    df_2_plot_subsetted[!(df_2_plot_subsetted$value == 0), ]
#
  # df_2_plot_subsetted <-
  #   df_2_plot[df_2_plot$variable == unique(df_2_plot$variable)[i], ]
  df_2_plot_subsetted$ADdiag3types <-
      factor(df_2_plot_subsetted$ADdiag3types,
             levels = c("earlyAD",
                        "lateAD",
                        "nonAD"))

  pdf(file = paste0("Rev_plots/v3_425_",
                    i,
                    "_Anova_Dunnett.pdf"),
      width = 2,
      height = 4)



# i <-
#   unique(df_2_plot$variable)[1]
# df_2_plot_subsetted <-
#   df_2_plot[df_2_plot$variable == i, ]
#
# aov_results <-
#   anova_test(formula = value ~ ADdiag3types,
#              data = df_2_plot_subsetted,
#              type = 3)
#
# Dunnett_p <-
#   DunnettTest(x = df_2_plot_subsetted$value,
#               g = df_2_plot_subsetted$ADdiag3types,
#               control = "nonAD")
#
# # aov_results <-
#   anova_summary(aov(formula = value ~ ADdiag3types,
#                     data = df_2_plot_subsetted))
# # summary(aov_results)

# pwc <-
  # aov_results %>%
  # print(LSD.test(aov_results,
  #                "ADdiag3types")) #%>%
  # add_xy_position(x = "ADdiag3types")

output_panel <-
  ggplot(df_2_plot_subsetted,
         aes(x = ADdiag3types,
             y = value)) +
  # geom_boxplot(outliers = T,
  #              ) +
  # geom_point(aes(shape = Clone,
  #                colour = Genotype)) +
  geom_jitter(aes(colour = ADdiag3types,
                  fill = ADdiag3types),
              shape = 1,
              width = 0.1) +
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.4,
               linewidth = 0.5) +
  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               geom = "errorbar",
               width = 0.2,
               linewidth = 1) +
  # scale_shape_manual(values = c(1, 2),
  #                    labels = c("Clone 1",
  #                               "Clone 2")) +
  scale_colour_manual(aes(colour = factor(ADdiag3types)),
                      values = c("darkblue",
                                 "darkred",
                                 "darkgreen")) +

  # scale_fill_manual(aes(factor(Genotype,
  #                              levels = c("risk",
  #                                         "non-risk",
  #                                         "risk-CRISPRa"))),
  #                   values = c("darkblue",
  #                              "darkred",
  #                              "darkgreen"),
  #                   breaks = c("risk",
  #                              "non-risk",
  #                              "risk-CRISPRa")) +
  guides(colour = guide_legend(title = "ADdiag3types"),
         fill = "none") +
  # ylim(0, 0.7) +
  xlab("") +
  ylab("PICALM expression in log2CPM") +
  scale_y_continuous(expand = c(0, 0.05)) +
    # stat_compare_means(method = "anova") +
  # geom_linerange()
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12,
                                    face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1),
        axis.title = element_text(size = 12,
                                  face = "bold"))

  if (anova_test(formula = value ~ ADdiag3types,
                 data = df_2_plot_subsetted,
                 type = 3)$p > 0.05) {
    sub_text <- ""
    anova_p_signif = ", not significant"
  } else {
    sub_text <-
      paste0("Dunnett \'s multiple comparison test adjusted P values\n",
             rownames(DunnettTest(x = df_2_plot_subsetted$value,
                                  g = df_2_plot_subsetted$ADdiag3types,
                                  control = "nonAD")[[1]])[1],
             ": ",
             DunnettTest(x = df_2_plot_subsetted$value,
                         g = df_2_plot_subsetted$ADdiag3types,
                         control = "nonAD")[[1]][1, 4],
             "\n",
             rownames(DunnettTest(x = df_2_plot_subsetted$value,
                                  g = df_2_plot_subsetted$ADdiag3types,
                                  control = "nonAD")[[1]])[2],
             ": ",
             DunnettTest(x = df_2_plot_subsetted$value,
                         g = df_2_plot_subsetted$ADdiag3types,
                         control = "nonAD")[[1]][2, 4])
    anova_p_signif = ", significant"

  }

  output_panel <-
    output_panel +
    ggtitle(label = paste0(i,
                           "; Anova \'s P =",
                          anova_test(formula = value ~ ADdiag3types,
                                     data = df_2_plot_subsetted,
                                     type = 3)$p,
                          anova_p_signif),
            subtitle = sub_text)

  print(output_panel)

  dev.off()

}



output_panel +
  # stat_pvalue_manual(Dunnett_p,
  #                    hide.ns = T) +
  # stat_anova_test(method = "one_way",
  #                 group.by = "x.var",
  #                 type = 2,
  #                 correction = "none",
  #                 p.adjust.method = "fdr") +
  ggtitle("MG")


test_var <-
  DunnettTest(x = df_2_plot_subsetted$value,
              g = df_2_plot_subsetted$ADdiag3types,
              control = "nonAD")[[1]]


df_raw <-
  read.table("Immune_cells_P2RY12_pos_cpm_wo_log.tsv")

df_2_plot <-
  df_raw
df_2_plot <-
  df_2_plot[rownames(df_2_plot) == "PICALM", ]

df_2_plot <-
  as.data.frame(t(df_2_plot))
df_2_plot$g_id <-
  rownames(df_2_plot)

df_2_match <-
  read_excel("Batch_2_of_new_plots/PICALM-400plus_AD-3-type-8other_cell_types.xlsx")
df_2_match <-
  df_2_match[, c(2, 9)]

df_2_plot <-
  merge(df_2_plot,
        df_2_match,
        by.x = "g_id",
        by.y = "g_projid")

df_2_plot_subsetted <-
  df_2_plot
df_2_plot_subsetted$PICALM <-
  log2(df_2_plot_subsetted$PICALM)
df_2_plot_subsetted$value <-
  df_2_plot_subsetted$PICALM

df_2_plot_subsetted <-
  df_2_plot_subsetted[!is.na(df_2_plot_subsetted$value), ]
df_2_plot_subsetted <-
  df_2_plot_subsetted[!is.infinite(df_2_plot_subsetted$value), ]

output_panel <-
  ggplot(df_2_plot_subsetted,
         aes(x = ADdiag3types,
             y = value)) +
  # geom_boxplot(outliers = T,
  #              ) +
  # geom_point(aes(shape = Clone,
  #                colour = Genotype)) +
  geom_jitter(aes(colour = ADdiag3types,
                  fill = ADdiag3types),
              shape = 1,
              width = 0.1) +
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.4,
               size = 0.5) +
  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               geom = "errorbar",
               width = 0.2) +
  # scale_shape_manual(values = c(1, 2),
  #                    labels = c("Clone 1",
  #                               "Clone 2")) +
  scale_colour_manual(aes(colour = factor(ADdiag3types)),
                      values = c("darkblue",
                                 "darkred",
                                 "darkgreen")) +

  # scale_fill_manual(aes(factor(Genotype,
  #                              levels = c("risk",
  #                                         "non-risk",
  #                                         "risk-CRISPRa"))),
  #                   values = c("darkblue",
  #                              "darkred",
  #                              "darkgreen"),
  #                   breaks = c("risk",
  #                              "non-risk",
  #                              "risk-CRISPRa")) +
  guides(colour = guide_legend(title = "ADdiag3types"),
         fill = F) +
  # ylim(0, 0.7) +
  xlab("") +
  ylab("PICALM expression in log2CPM") +
  scale_y_continuous(expand = c(0, 0.05),
                     limits = c(8.5, 11.5)) +
  # stat_compare_means(method = "anova") +
  # geom_linerange()
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12,
                                    face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1),
        axis.title = element_text(size = 12,
                                  face = "bold"))

if (rstatix::anova_test(formula = value ~ ADdiag3types,
               data = df_2_plot_subsetted,
               type = 3)$p > 0.05) {
  sub_text <- ""
  anova_p_signif = ", not significant"
} else {
  sub_text <-
    paste0("Dunnett \'s multiple comparison test adjusted P values\n",
           rownames(DunnettTest(x = df_2_plot_subsetted$value,
                                g = df_2_plot_subsetted$ADdiag3types,
                                control = "nonAD")[[1]])[1],
           ": ",
           DunnettTest(x = df_2_plot_subsetted$value,
                       g = df_2_plot_subsetted$ADdiag3types,
                       control = "nonAD")[[1]][1, 4],
           "\n",
           rownames(DunnettTest(x = df_2_plot_subsetted$value,
                                g = df_2_plot_subsetted$ADdiag3types,
                                control = "nonAD")[[1]])[2],
           ": ",
           DunnettTest(x = df_2_plot_subsetted$value,
                       g = df_2_plot_subsetted$ADdiag3types,
                       control = "nonAD")[[1]][2, 4])
  anova_p_signif = ", significant"

}


DunnettTest(x = df_2_plot_subsetted$value,
            g = df_2_plot_subsetted$ADdiag3types,
            control = "nonAD")

mean(df_2_plot_subsetted$value[df_2_plot_subsetted$ADdiag3types == "nonAD"])
mean(df_2_plot_subsetted$value[df_2_plot_subsetted$ADdiag3types == "lateAD"])


output_panel <-
  output_panel +
  ggtitle(label = paste0(i,
                         "; Anova \'s P =",
                         anova_test(formula = value ~ ADdiag3types,
                                    data = df_2_plot_subsetted,
                                    type = 3)$p,
                         anova_p_signif),
          subtitle = "MG")

print(output_panel)


# Krusker-Wallis test


# KW test, plot 3x3 cell types ####

df_raw <-
  read_excel("Batch_2_of_new_plots/PICALM-400plus_AD-3-type-8other_cell_types.xlsx")

df_2_plot <-
  df_raw

df_2_plot <-
  df_2_plot[, c(9, 11:18)]


df_2_plot <-
  data.table::melt(df_2_plot,
                   na.rm = T)
df_2_plot <-
  df_2_plot[!is.na(df_2_plot$variable), ]
df_2_plot <-
  df_2_plot[!is.na(df_2_plot$value), ]

dev.off()

# i <- 1

for (i in unique(df_2_plot$variable)) {
  print(i)

  df_2_plot_subsetted <-
    df_2_plot[df_2_plot$variable == i, ]
  df_2_plot_subsetted <-
    df_2_plot_subsetted[!is.na(df_2_plot_subsetted$value), ]
  df_2_plot_subsetted <-
    df_2_plot_subsetted[!is.nan(df_2_plot_subsetted$value), ]
  df_2_plot_subsetted <-
    df_2_plot_subsetted[!is.infinite(df_2_plot_subsetted$value), ]
  df_2_plot_subsetted <-
    df_2_plot_subsetted[!(df_2_plot_subsetted$value == 0), ]
  #
  # df_2_plot_subsetted <-
  #   df_2_plot[df_2_plot$variable == unique(df_2_plot$variable)[i], ]
  df_2_plot_subsetted$ADdiag3types <-
    factor(df_2_plot_subsetted$ADdiag3types,
           levels = c("earlyAD",
                      "lateAD",
                      "nonAD"))

  pdf(file = paste0("Rev_plots/v3_425_",
                    i,
                    "_KW_Dunnett.pdf"),
      width = 2,
      height = 4)



  # i <-
  #   unique(df_2_plot$variable)[1]
  # df_2_plot_subsetted <-
  #   df_2_plot[df_2_plot$variable == i, ]
  #
  # aov_results <-
  #   anova_test(formula = value ~ ADdiag3types,
  #              data = df_2_plot_subsetted,
  #              type = 3)
  #
  # Dunnett_p <-
  #   DunnettTest(x = df_2_plot_subsetted$value,
  #               g = df_2_plot_subsetted$ADdiag3types,
  #               control = "nonAD")
  #
  # # aov_results <-
  #   anova_summary(aov(formula = value ~ ADdiag3types,
  #                     data = df_2_plot_subsetted))
  # # summary(aov_results)

  # pwc <-
  # aov_results %>%
  # print(LSD.test(aov_results,
  #                "ADdiag3types")) #%>%
  # add_xy_position(x = "ADdiag3types")

  output_panel <-
    ggplot(df_2_plot_subsetted,
           aes(x = ADdiag3types,
               y = value)) +
    # geom_boxplot(outliers = T,
    #              ) +
    # geom_point(aes(shape = Clone,
    #                colour = Genotype)) +
    geom_jitter(aes(colour = ADdiag3types,
                    fill = ADdiag3types),
                shape = 1,
                width = 0.1) +
    stat_summary(fun = "mean",
                 geom = "crossbar",
                 width = 0.4,
                 linewidth = 0.5) +
    stat_summary(fun.data = mean_se,
                 fun.args = list(mult = 1),
                 geom = "errorbar",
                 width = 0.2,
                 linewidth = 1) +
    # scale_shape_manual(values = c(1, 2),
    #                    labels = c("Clone 1",
    #                               "Clone 2")) +
    scale_colour_manual(aes(colour = factor(ADdiag3types)),
                        values = c("darkblue",
                                   "darkred",
                                   "darkgreen")) +

    # scale_fill_manual(aes(factor(Genotype,
    #                              levels = c("risk",
    #                                         "non-risk",
    #                                         "risk-CRISPRa"))),
    #                   values = c("darkblue",
    #                              "darkred",
    #                              "darkgreen"),
    #                   breaks = c("risk",
    #                              "non-risk",
    #                              "risk-CRISPRa")) +
    guides(colour = guide_legend(title = "ADdiag3types"),
           fill = "none") +
    # ylim(0, 0.7) +
    xlab("") +
    ylab("PICALM expression in log2CPM") +
    scale_y_continuous(expand = c(0, 0.05)) +
    # stat_compare_means(method = "anova") +
    # geom_linerange()
    theme_classic() +
    theme(axis.text = element_text(size = 12,
                                   colour = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12,
                                      face = "bold"),
          legend.position = "none",
          axis.text.x = element_text(angle = 315,
                                     hjust = 0,
                                     vjust = 1),
          axis.title = element_text(size = 12,
                                    face = "bold"))

  if (kruskal.test(x = df_2_plot_subsetted$value,
                   g = df_2_plot_subsetted$ADdiag3types)$p.value > 0.05) {
    sub_text <- ""
    KW_p_signif = ", not significant"
  } else {
    sub_text <-
      paste0("Dunnett \'s multiple comparison test adjusted P values\n",
             rownames(DunnettTest(x = df_2_plot_subsetted$value,
                                  g = df_2_plot_subsetted$ADdiag3types,
                                  control = "nonAD")[[1]])[1],
             ": ",
             DunnettTest(x = df_2_plot_subsetted$value,
                         g = df_2_plot_subsetted$ADdiag3types,
                         control = "nonAD")[[1]][1, 4],
             "\n",
             rownames(DunnettTest(x = df_2_plot_subsetted$value,
                                  g = df_2_plot_subsetted$ADdiag3types,
                                  control = "nonAD")[[1]])[2],
             ": ",
             DunnettTest(x = df_2_plot_subsetted$value,
                         g = df_2_plot_subsetted$ADdiag3types,
                         control = "nonAD")[[1]][2, 4])
    KW_p_signif = ", significant"

  }

  output_panel <-
    output_panel +
    ggtitle(label = paste0(i,
                           "; KW \'s P =",
                           kruskal.test(x = df_2_plot_subsetted$value,
                                        g = df_2_plot_subsetted$ADdiag3types)$p.value,
                           KW_p_signif),
            subtitle = sub_text)

  print(output_panel)

  dev.off()

  print(paste0(i,
               "; ",
               kruskal.test(x = df_2_plot_subsetted$value,
                            g = df_2_plot_subsetted$ADdiag3types)$p.value))

}



output_panel +
  # stat_pvalue_manual(Dunnett_p,
  #                    hide.ns = T) +
  # stat_anova_test(method = "one_way",
  #                 group.by = "x.var",
  #                 type = 2,
  #                 correction = "none",
  #                 p.adjust.method = "fdr") +
  ggtitle("MG")


test_var <-
  DunnettTest(x = df_2_plot_subsetted$value,
              g = df_2_plot_subsetted$ADdiag3types,
              control = "nonAD")[[1]]


df_raw <-
  read.table("Immune_cells_P2RY12_pos_cpm_wo_log.tsv")

df_2_plot <-
  df_raw
df_2_plot <-
  df_2_plot[rownames(df_2_plot) == "PICALM", ]

df_2_plot <-
  as.data.frame(t(df_2_plot))
df_2_plot$g_id <-
  rownames(df_2_plot)

df_2_match <-
  read_excel("Batch_2_of_new_plots/PICALM-400plus_AD-3-type-8other_cell_types.xlsx")
df_2_match <-
  df_2_match[, c(2, 9)]

df_2_plot <-
  merge(df_2_plot,
        df_2_match,
        by.x = "g_id",
        by.y = "g_projid")

df_2_plot_subsetted <-
  df_2_plot
df_2_plot_subsetted$PICALM <-
  log2(df_2_plot_subsetted$PICALM)
df_2_plot_subsetted$value <-
  df_2_plot_subsetted$PICALM

df_2_plot_subsetted <-
  df_2_plot_subsetted[!is.na(df_2_plot_subsetted$value), ]
df_2_plot_subsetted <-
  df_2_plot_subsetted[!is.infinite(df_2_plot_subsetted$value), ]

output_panel <-
  ggplot(df_2_plot_subsetted,
         aes(x = ADdiag3types,
             y = value)) +
  # geom_boxplot(outliers = T,
  #              ) +
  # geom_point(aes(shape = Clone,
  #                colour = Genotype)) +
  geom_jitter(aes(colour = ADdiag3types,
                  fill = ADdiag3types),
              shape = 1,
              width = 0.1) +
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.4,
               size = 0.5) +
  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               geom = "errorbar",
               width = 0.2) +
  # scale_shape_manual(values = c(1, 2),
  #                    labels = c("Clone 1",
  #                               "Clone 2")) +
  scale_colour_manual(aes(colour = factor(ADdiag3types)),
                      values = c("darkblue",
                                 "darkred",
                                 "darkgreen")) +

  # scale_fill_manual(aes(factor(Genotype,
  #                              levels = c("risk",
  #                                         "non-risk",
  #                                         "risk-CRISPRa"))),
  #                   values = c("darkblue",
  #                              "darkred",
  #                              "darkgreen"),
  #                   breaks = c("risk",
  #                              "non-risk",
  #                              "risk-CRISPRa")) +
  guides(colour = guide_legend(title = "ADdiag3types"),
         fill = F) +
  # ylim(0, 0.7) +
  xlab("") +
  ylab("PICALM expression in log2CPM") +
  scale_y_continuous(expand = c(0, 0.05),
                     limits = c(8.5, 11.5)) +
  # stat_compare_means(method = "anova") +
  # geom_linerange()
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12,
                                    face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1),
        axis.title = element_text(size = 12,
                                  face = "bold"))

if (rstatix::anova_test(formula = value ~ ADdiag3types,
                        data = df_2_plot_subsetted,
                        type = 3)$p > 0.05) {
  sub_text <- ""
  anova_p_signif = ", not significant"
} else {
  sub_text <-
    paste0("Dunnett \'s multiple comparison test adjusted P values\n",
           rownames(DunnettTest(x = df_2_plot_subsetted$value,
                                g = df_2_plot_subsetted$ADdiag3types,
                                control = "nonAD")[[1]])[1],
           ": ",
           DunnettTest(x = df_2_plot_subsetted$value,
                       g = df_2_plot_subsetted$ADdiag3types,
                       control = "nonAD")[[1]][1, 4],
           "\n",
           rownames(DunnettTest(x = df_2_plot_subsetted$value,
                                g = df_2_plot_subsetted$ADdiag3types,
                                control = "nonAD")[[1]])[2],
           ": ",
           DunnettTest(x = df_2_plot_subsetted$value,
                       g = df_2_plot_subsetted$ADdiag3types,
                       control = "nonAD")[[1]][2, 4])
  anova_p_signif = ", significant"

}


DunnettTest(x = df_2_plot_subsetted$value,
            g = df_2_plot_subsetted$ADdiag3types,
            control = "nonAD")

mean(df_2_plot_subsetted$value[df_2_plot_subsetted$ADdiag3types == "nonAD"])
mean(df_2_plot_subsetted$value[df_2_plot_subsetted$ADdiag3types == "lateAD"])


output_panel <-
  output_panel +
  ggtitle(label = paste0(i,
                         "; Anova \'s P =",
                         anova_test(formula = value ~ ADdiag3types,
                                    data = df_2_plot_subsetted,
                                    type = 3)$p,
                         anova_p_signif),
          subtitle = "MG")

print(output_panel)
#
# [1] "Ast"
# [1] "Ast; 0.841254311208934"
# [1] "OPC"
# [1] "OPC; 0.61173206485834"
# [1] "Oligo"
# [1] "Oligo; 0.109533577280205"
# [1] "Ext-1"
# [1] "Ext-1; 0.0938804963613182"
# [1] "Ext-2"
# [1] "Ext-2; 0.231743511015968"
# [1] "Ext-3"
# [1] "Ext-3; 0.281172037185996"
# [1] "Inh"
# [1] "Inh; 0.141697740780614"
# [1] "Vas"
# [1] "Vas; 0.0395203304357988"
