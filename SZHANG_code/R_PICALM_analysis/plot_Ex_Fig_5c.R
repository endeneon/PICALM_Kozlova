# Siwei 26 Jan 2025
# plot new Ex. Fig 5c

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

rm(list = ls())

df_raw <-
  read_excel("Batch_2_of_new_plots/Alena_TG_quant.xlsx")
colnames(df_raw)[3] <- "genotype"

df_2_plot <-
  df_raw[df_raw$Cell_Line == "CD04", ]
df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = c("non-risk",
                    "risk",
                    "risk-CRISPRa"))
df_2_plot$genotype <-
  relevel(df_2_plot$genotype,
          ref = "risk")

summary(multcomp::glht(aov(value ~ genotype,
                           data = df_2_plot),
                       linfct = mcp(genotype = "Dunnett")))
df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = c("non-risk",
                    "risk",
                    "risk-CRISPRa"))
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)
# non-risk - risk == 0     -1.930e+09  3.790e+08  -5.092  0.00121 **
#   risk-CRISPRa - risk == 0 -3.574e+09  3.790e+08  -9.430  1.1e-05 ***

ggerrorplot(df_2_plot,
            x = "genotype",
            y = "value",
            color = "genotype",
            # shape = "Clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkblue",
                                 "darkred",
                                 "green4"),
                      guide = "none") +
  # stat_compare_means(label = "p.signif",
  #                    label.y.npc = 0.95,
  #                    # label.x.npc = 0,
  #                    method = 'wilcox.test',
  #                    hide.ns = F,
  #                    ref.group = "non-risk",
  #                    paired = F) +
  geom_jitter(aes(colour = genotype),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "TG level") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_raw$value,
                                    na.rm = T) * 1.1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD04")



# CD 09 ####
df_2_plot <-
  df_raw[df_raw$Cell_Line == "CD09", ]
df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = c("non-risk",
                    "risk",
                    "risk-CRISPRa"))
df_2_plot$genotype <-
  relevel(df_2_plot$genotype,
          ref = "risk")

summary(multcomp::glht(aov(value ~ genotype,
                           data = df_2_plot),
                       linfct = mcp(genotype = "Dunnett")))
df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = c("non-risk",
                    "risk",
                    "risk-CRISPRa"))
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)
# non-risk - risk == 0     -3.186e+09  3.642e+08  -8.748 2.03e-05 ***
#   risk-CRISPRa - risk == 0 -3.885e+09  3.642e+08 -10.669 3.94e-06 ***

ggerrorplot(df_2_plot,
            x = "genotype",
            y = "value",
            color = "genotype",
            # shape = "Clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  # stat_compare_means(label = "p.signif",
  #                    label.y.npc = 0.95,
  #                    # label.x.npc = 0,
  #                    method = 'wilcox.test',
  #                    hide.ns = F,
  #                    ref.group = "non-risk",
  #                    paired = F) +
  geom_jitter(aes(colour = genotype),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "TG level") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_raw$value,
                                    na.rm = T) * 1.1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD09")
