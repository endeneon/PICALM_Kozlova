# Siwei 31 Jan 2025
# plot Fig. 6b

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

#
df_raw <-
  read_excel(path = "Batch_2_of_new_plots/new_Fig_6b.xlsx")

# CD04 #####

df_2_plot <-
  df_raw[df_raw$Cell_Line == "CD04", ]

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "risk",
                    "risk-CRISPRa"))
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "value",
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
       y = "Bodipy") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_raw$value) + 0.1)) +
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

summary(aov(value ~ Genotype,
            data = df_2_plot))
summary(multcomp::glht(aov(value ~ Genotype,
                           data = df_2_plot),
                       linfct = mcp(Genotype = "Dunnett")))
#
# Multiple Comparisons of Means: Dunnett Contrasts
#
#
# Fit: aov(formula = value ~ Genotype, data = df_2_plot)
#
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)
# non-risk - risk == 0      -2.6549     0.2775  -9.568 1.74e-07 ***
#   risk-CRISPRa - risk == 0  -2.6403     0.2775  -9.515 1.87e-07 ***

kruskal.test(x = df_2_plot$value,
             g = df_2_plot$Genotype)
DunnTest(value ~ Genotype,
         data = df_2_plot)
# Kruskal-Wallis rank sum test
#
# data:  df_2_plot$value and df_2_plot$Genotype
# Kruskal-Wallis chi-squared = 11.38, df = 2, p-value = 0.003379

# Dunn's test of multiple comparisons using rank sums : holm
#
#                       mean.rank.diff   pval
# non-risk-risk             -9.1666667 0.0088 **
# risk-CRISPRa-risk         -8.8333333 0.0088 **
# risk-CRISPRa-non-risk      0.3333333 0.9139
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# CD09 #####
df_2_plot <-
  df_raw[df_raw$Cell_Line == "CD09", ]

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "risk",
                    "risk-CRISPRa"))
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "value",
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
       y = "Bodipy") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_raw$value) + 0.1)) +
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

summary(aov(value ~ Genotype,
            data = df_2_plot))
summary(multcomp::glht(aov(value ~ Genotype,
                           data = df_2_plot),
                       linfct = mcp(Genotype = "Dunnett")))
#
# Multiple Comparisons of Means: Dunnett Contrasts
#
#
# Fit: aov(formula = value ~ Genotype, data = df_2_plot)
#
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)
# non-risk - risk == 0      -1.4787     0.3397  -4.354  0.00858 **
#   risk-CRISPRa - risk == 0  -1.8061     0.3397  -5.317  0.00324 **
kruskal.test(x = df_2_plot$value,
             g = df_2_plot$Genotype)
DunnTest(value ~ Genotype,
         data = df_2_plot)
# Kruskal-Wallis rank sum test
#
# data:  df_2_plot$value and df_2_plot$Genotype
# Kruskal-Wallis chi-squared = 11.556, df = 2, p-value = 0.003096
#
# Dunn's test of multiple comparisons using rank sums : holm
#
#                       mean.rank.diff   pval
# non-risk-risk              -8.333333 0.0137 *
# risk-CRISPRa-risk          -9.666667 0.0051 **
# risk-CRISPRa-non-risk      -1.333333 0.6653
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
