# Siwei 26 Jan 2025
# plot Fig. 3e

# init ####
{
  library(readxl)

  library(stringr)
  library(ggplot2)

  library(scales)
  library(reshape2)

  library(RColorBrewer)
  library(ggpubr)
  library(ggridges)

  library(dplyr)
  library(data.table)

  library(DescTools)
  library(multcomp)

  library(gridExtra)

  library(ggallin)
}

set.seed(42)


df_raw <-
  read_excel(path = "Batch_2_of_new_plots/Fig_3e.xlsx")

## CD04 ####
# df_sum_all <-
#   df_raw
# df_sum_all$BR <-
#   str_split(string = df_sum_all$BR_FOV,
#             pattern = "_",
#             simplify = T)[, 1]

# colnames(df_2_plot)[2] <-
#   "Genotype"
df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "risk",
                    "risk-CRISPRa"))
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))



ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "PICALM.actin",
            color = "Genotype",
            # shape = "Clone",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkblue",
                                 "darkred",
                                 "darkgreen"),
                      guide = "none") +
  # stat_compare_means(label = "p.signif",
  #                    label.y.npc = 0.95,
  #                    # label.x.npc = 0,
  #                    method = 't.test',
  #                    hide.ns = F,
  #                    ref.group = "non-risk",
  #                    paired = F) +
  geom_jitter(aes(colour = Genotype),
              width = 0.1,
              shape = 1,
              size = 1) +
  # scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "PICALM/b-actin") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        # legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD04")

df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "risk")

summary(aov(PICALM.actin ~ Genotype,
            data = df_2_plot))
# summary(multcomp::glht(aov(PICALM.actin ~ Genotype,
#                            data = df_2_plot),
#                        linfct = mcp(Genotype = "Dunnett")))
kruskal.test(x = df_2_plot$PICALM.actin,
             g = df_2_plot$Genotype)
DunnTest(PICALM.actin ~ Genotype,
         data = df_2_plot)

# Kruskal-Wallis rank sum test
#
# data:  df_2_plot$PICALM.actin and df_2_plot$Genotype
# Kruskal-Wallis chi-squared = 5.4222, df = 2, p-value = 0.06646

# Dunn's test of multiple comparisons using rank sums : holm
#
#                       mean.rank.diff   pval
# risk-non-risk             -4.6666667 0.1107
# risk-CRISPRa-non-risk     -0.3333333 0.8815
# risk-CRISPRa-risk          4.3333333 0.1107
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# > summary(aov(PICALM.actin ~ Genotype,
#               +             data = df_2_plot))
# Df  Sum Sq Mean Sq F value  Pr(>F)
# Genotype     2 0.21914 0.10957   12.07 0.00789 **
  # Residuals    6 0.05446 0.00908
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)
# non-risk - risk == 0      0.33266    0.07779   4.276  0.00933 **
#   risk-CRISPRa - risk == 0  0.32935    0.07779   4.234  0.00977 **

## CD09 ####
df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]
# colnames(df_2_plot)[2] <-
#   "Genotype"

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "risk",
                    "risk-CRISPRa"))
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))


# summary(aov(PICALM.actin ~ Genotype,
#             data = df_2_plot))
# summary(multcomp::glht(aov(PICALM.actin ~ Genotype,
#                            data = df_2_plot),
#                        linfct = mcp(Genotype = "Dunnett")))
#

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "PICALM.actin",
            color = "Genotype",
            # shape = "Clone",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkblue",
                                 "darkred",
                                 "darkgreen"),
                      guide = "none") +
  # stat_compare_means(label = "p.signif",
  #                    label.y.npc = 0.95,
  #                    # label.x.npc = 0,
  #                    method = 't.test',
  #                    hide.ns = F,
  #                    ref.group = "non-risk",
  #                    paired = F) +
  geom_jitter(aes(colour = Genotype),
              width = 0.1,
              shape = 1,
              size = 1) +
  # scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "PICALM/b-actin") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        # legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD09")

df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "risk")


summary(aov(PICALM.actin ~ Genotype,
            data = df_2_plot))
summary(multcomp::glht(aov(PICALM.actin ~ Genotype,
                           data = df_2_plot),
                       linfct = mcp(Genotype = "Dunnett")))
DunnTest(PICALM.actin ~ Genotype,
         data = df_2_plot)

# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)
# non-risk - risk == 0      0.29853    0.05907   5.054  0.00417 **
#   risk-CRISPRa - risk == 0  0.24679    0.05907   4.178  0.01039 *


# # > summary(aov(PICALM.actin ~ Genotype,
#               +             data = df_2_plot))
# Df Sum Sq Mean Sq F value  Pr(>F)
# Genotype     2 0.1527 0.07635   14.59 0.00496 **
#   Residuals    6 0.0314 0.00523

kruskal.test(x = df_2_plot$PICALM.actin,
             g = df_2_plot$Genotype)
DunnTest(PICALM.actin ~ Genotype,
         data = df_2_plot)


# Kruskal-Wallis rank sum test
#
# data:  df_2_plot$PICALM.actin and df_2_plot$Genotype
# Kruskal-Wallis chi-squared = 5.6, df = 2, p-value = 0.06081

# Dunn's test of multiple comparisons using rank sums : holm
#
#                       mean.rank.diff   pval
# risk-non-risk                     -5 0.0760 .
# risk-CRISPRa-non-risk             -1 0.6547
# risk-CRISPRa-risk                  4 0.1473
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
