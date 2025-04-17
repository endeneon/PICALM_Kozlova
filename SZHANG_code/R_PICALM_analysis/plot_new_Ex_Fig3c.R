# Siwei 24 Jan 2025
# plot new Ex Fig. 3c

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
}

df_raw <-
  read_excel("Batch_2_of_new_plots/new_Fig_Ex3c.xlsx")

df_2_plot <-
  df_raw[df_raw$line == "CD04", ]

df_2_plot$clones <-
  as.numeric(factor(df_2_plot$clone))

df_2_plot$clones[df_2_plot$clones == 3] <- 1
df_2_plot$clones[df_2_plot$clones == 4] <- 2

df_2_plot$clones <-
  str_c("Clone ",
        df_2_plot$clones)

df_2_plot$genotype[df_2_plot$genotype == "nonrisk"] <-
  "non-risk"
df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = c("risk",
                    "non-risk"))

ggerrorplot(df_2_plot,
            x = "genotype",
            y = "Average_Exp",
            color = "genotype",
            shape = "clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +


  stat_compare_means(label = "p.signif",
                     label.y.npc = 0.95,
                     # label.x.npc = 0,
                     method = 'wilcox.test',
                     hide.ns = F,
                     ref.group = "non-risk",
                     paired = F) +
  geom_jitter(aes(colour = genotype,
                  shape = clones),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
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
                                0.033)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD04")

df_2_plot <-
  df_raw[df_raw$line == "CD09", ]

df_2_plot$clones <-
  as.numeric(factor(df_2_plot$clone))

df_2_plot$clones[df_2_plot$clones == 3] <- 1
df_2_plot$clones[df_2_plot$clones == 4] <- 2

df_2_plot$clones <-
  str_c("Clone ",
        df_2_plot$clones)

df_2_plot$genotype[df_2_plot$genotype == "nonrisk"] <-
  "non-risk"
df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = c("risk",
                    "non-risk"))

ggerrorplot(df_2_plot,
            x = "genotype",
            y = "Average_Exp",
            color = "genotype",
            shape = "clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +


  stat_compare_means(label = "p.signif",
                     label.y.npc = 0.95,
                     # label.x.npc = 0,
                     method = 'wilcox.test',
                     hide.ns = F,
                     ref.group = "non-risk",
                     paired = F) +
  geom_jitter(aes(colour = genotype,
                  shape = clones),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
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
                                max(df_2_plot$Average_Exp,
                                    na.rm = T) * 1.1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD09")
