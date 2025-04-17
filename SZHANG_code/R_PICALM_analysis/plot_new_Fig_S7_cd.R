# Siwei 02 Feb 2025
# plot new S 7c-d

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

# load raw data ####
rm(list = ls())


## Fig. S7c #####
df_raw <-
  read_excel("Batch_2_of_new_plots/new_S_Fig_7c_d.xlsx",
             sheet = 1)

df_2_plot <-
  df_raw

df_2_plot <-
  reshape2::melt(df_2_plot,
                 value.name = "Value")

df_2_plot <-
  df_2_plot[!(is.na(df_2_plot$Value)), ]


df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "risk TrC",
                    "non-risk",
                    "non-risk TrC"))
#
# df_2_plot$Genotype <-
#   factor(df_2_plot$Genotype,
#          levels = c("non-risk",
#                     "non-risk TrC",
#                     "risk",
#                     "risk TrC"))

df_2_plot$variable <-
  factor(df_2_plot$variable,
         levels = c("0 min",
                    "45 min",
                    "90 min",
                    "135 min",
                    "180 min"))



ggplot(df_2_plot,
       aes(x = variable,
           y = Value,
           colour = Genotype,
           fill = Genotype,
           group = Genotype)) +
  # ggerrorplot(df_2_plot,
  #             x = "variable",
  #             y = "Value",
  #             color = "Genotype",
  #             width = 0.2,
  #             # facet.by = "Gene",
  #             # ncol = 4,
  #             error.plot = "errorbar") +
  geom_line(stat = "summary",
            fun = mean,
            linewidth = 1) +
  geom_errorbar(stat = "summary",
                fun.data = mean_se,
                width = 0.3) +
  geom_point(shape = 21,
             stat = "summary",
             fun = mean,
             colour = "black") +


  # geom_point(shape = 21,
  #            colour = "black") +
  scale_fill_manual(values = c("darkred",
                               "green4",
                               "darkblue",
                               "darkorange")) +
  scale_colour_manual(values = c("darkred",
                                 "green4",
                                 "darkblue",
                                 "darkorange")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.5)) +
  labs(x = "",
       y = "pHrodo intensity / Normalized to non-risk (180 min)") +
  theme_classic() +
  # ylim(0, 0.8) +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0))


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "Genotype",
            # shape = "clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "green4",
                                 "darkblue",
                                 "darkorange")) +

  # stat_compare_means(label = "p.signif",
  #                    # label.y.npc = 0.95,
  #                    # label.x.npc = 0,
  #                    method = 'wilcox.test',
  #                    hide.ns = T,
  #                    ref.group = "non-risk",
  #                    paired = F) +
  geom_jitter(aes(colour = Genotype),
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
       y = "pHrodo int./Normalized non-risk (180 min)") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                1.5)) +
  # guides(guide_col)
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1))

## Fig S7d ####
df_raw <-
  read_excel("Batch_2_of_new_plots/new_S_Fig_7c_d.xlsx",
             sheet = 2)

df_2_plot <-
  df_raw

# df_2_plot <-
#   reshape2::melt(df_2_plot,
#                  value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "risk TrC",
                    "non-risk",
                    "non-risk TrC"))

# df_2_plot$Genotype <-
#   factor(df_2_plot$Genotype,
#          levels = c("non-risk",
#                     "non-risk TrC",
#                     "risk",
#                     "risk TrC"))

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Ave_Bodipy",
            color = "Genotype",
            # shape = "clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "green4",
                                 "darkblue",
                                 "darkorange")) +

  # stat_compare_means(label = "p.signif",
  #                    # label.y.npc = 0.95,
  #                    # label.x.npc = 0,
  #                    method = 'wilcox.test',
  #                    hide.ns = T,
  #                    ref.group = "non-risk",
  #                    paired = F) +
  geom_jitter(aes(colour = Genotype),
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
       y = "BODIPY intensity/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                2)) +
  # guides(guide_col)
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1))

df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "non-risk")

summary(aov(Ave_Bodipy ~ Genotype,
            data = df_2_plot))
summary(multcomp::glht(aov(Ave_Bodipy ~ Genotype,
                           data = df_2_plot),
                       linfct = mcp(Genotype = "Dunnett")))

# Multiple Comparisons of Means: Dunnett Contrasts
#
#
# Fit: aov(formula = Ave_Bodipy ~ Genotype, data = df_2_plot)
#
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)
# non-risk TrC - non-risk == 0 -0.13319    0.11050  -1.205    0.515
# risk - non-risk == 0          0.90413    0.09883   9.148   <1e-04 ***
#   risk TrC - non-risk == 0      0.05035    0.11050   0.456    0.945


