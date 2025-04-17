# Siwei 24 Jan 2025
# plot new Ext. Fig 5b

# init ####
{
  library(readxl)

  library(stringr)
  library(ggplot2)

  library(scales)
  library(reshape2)

  library(RColorBrewer)
  library(ggpubr)

  library(dplyr)
  library(data.table)

  library(DescTools)
  library(multcomp)

  library(gridExtra)
}

## Fig_3B ####
df_raw <-
  read_excel("Batch_2_of_new_plots/new_Fig_Ex5b.xlsx",
             sheet = 1)

#
# df_2_plot <-
#   df_raw

# df_2_plot <-
#   melt(df_2_plot,
#        value.name = "Value")

## CD04 #####

df_2_plot <-
  df_raw[df_raw$line == "CD04", ]
df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = c("risk",
                    "non-risk",
                    "non-risk-CRISPRoff"))
df_2_plot$genotype <-
  relevel(df_2_plot$genotype,
          ref = "risk")
df_2_plot$Average_Exp <-
  df_2_plot$Average_Exp /
  max(df_2_plot$Average_Exp)

#
# output_results_p1 <-
  summary(multcomp::glht(aov(Average_Exp ~ genotype,
                             data = df_2_plot),
                         linfct = mcp(genotype = "Dunnett")))
p1 <-
ggerrorplot(df_2_plot,
              x = "genotype",
              y = "Average_Exp",
              color = "genotype",
              # shape = "clones",
              # group = "Genotype",
              width = 0.2,
              # facet.by = "Gene",
              # ncol = 4,
              error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue",
                                 "darkorange")) +

  stat_compare_means(label = "p.signif",
                     # label.y.npc = 0.95,
                     # label.x.npc = 0,
                     method = 'wilcox.test',
                     hide.ns = T,
                     ref.group = "non-risk",
                     paired = F) +
  geom_jitter(aes(colour = genotype),
              width = 0.1,
              size = 1,
              shape = 1) +
  # scale_shape_manual(values = c(1)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Normalized PICALM exp.") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                1.1)) +
  # guides(guide_col)
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD04")

## CD09 #####
df_2_plot <-
  df_raw[df_raw$line == "CD09", ]
df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = c("risk",
                    "non-risk",
                    "non-risk-CRISPRoff"))
df_2_plot$genotype <-
  relevel(df_2_plot$genotype,
          ref = "non-risk")
df_2_plot$Average_Exp <-
  df_2_plot$Average_Exp /
  max(df_2_plot$Average_Exp)

#
output_results_p2 <-
  summary(multcomp::glht(aov(Average_Exp ~ genotype,
                             data = df_2_plot),
                         linfct = mcp(genotype = "Dunnett")))
p2 <-
ggerrorplot(df_2_plot,
            x = "genotype",
            y = "Average_Exp",
            color = "genotype",
            # shape = "clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue",
                                 "darkorange")) +

  stat_compare_means(label = "p.signif",
                     # label.y.npc = 0.95,
                     # label.x.npc = 0,
                     method = 'wilcox.test',
                     hide.ns = T,
                     ref.group = "non-risk",
                     paired = F) +
  geom_jitter(aes(colour = genotype),
              width = 0.1,
              size = 1,
              shape = 1) +
  # scale_shape_manual(values = c(1)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                1.1)) +
  # guides(guide_col)
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD09")
#   ggtitle(label = paste0("CD_09, Pa = ",
#                          summary(aov(formula = Average_Exp ~ genotype,
#                                      data = df_2_plot))[[1]][["Pr(>F)"]][1]),
#           subtitle = paste(paste0(rownames(output_results_p1$linfct)[1],
#                                   ", Dunnett\'s P = ",
#                                   output_results_p1$test$pvalues[1]),
#                            paste0(rownames(output_results_p1$linfct)[2],
#                                   ", Dunnett\'s P = ",
#                                   output_results_p1$test$pvalues[2]),
#                            sep = "\n"))

##

#    output_results <-
#     summary(aov(formula = Average_Exp ~ genotype,
#                 data = df_2_plot))
# # output_results <-
#   summary(aov(formula = Average_Exp ~ genotype,
#               data = df_2_plot)[[1]][[1]])#

  # output_results <-
  # summary(multcomp::glht(aov(Average_Exp ~ genotype,
  #                            data = df_2_plot),
  #                        linfct = mcp(genotype = "Dunnett"))$test$pvalues[1])

# 2-panels plot #####
gridExtra::grid.arrange(p1, p2,
                             ncol = 2)

