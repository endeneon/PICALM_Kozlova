# Siwei 24 Jan 2025
# plot new Ex Fig. 3b

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

# TREM2, IBA1, PY2 #####

df_raw <-
  read_excel(path = "Batch_2_of_new_plots/new_Fig_Ex3b.xlsx",
             sheet = 1)

df_2_plot <-
  df_raw

# df_2_clones <-

df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = c("risk",
                    "non-risk"))


df_2_plot$clones <-
  str_c("Clone ",
        df_2_plot$clone)


ggerrorplot(df_2_plot,
            x = "genotype",
            y = "value",
            color = "genotype",
            # shape = "clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  # scale_colour_manual(values = brewer.pal(n = 3,
  #                                         name = "Dark2")[c(2, 3)],
  #                     guide = "none") +


  # stat_compare_means(label = "p.signif",
  #                    label.y.npc = 0.95,
  #                    # label.x.npc = 0,
  #                    method = 'wilcox.test',
  #                    hide.ns = F,
  #                    ref.group = "non-risk",
  #                    paired = F) +
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
       y = "TREM2+, IBA1+, PY2+ cells") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                1.1),
                     labels = scales::percent) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("TREM2+, IBA1+, PY2+ cells")


# TREM2, IBA1, TMEM119 #####

df_raw <-
  read_excel(path = "Batch_2_of_new_plots/new_Fig_Ex3b.xlsx",
             sheet = 2)

df_2_plot <-
  df_raw

# df_2_clones <-

df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = c("risk",
                    "non-risk"))


# df_2_plot$Clones <-
#   str_c("Clone ",
#         df_2_plot$clone)


ggerrorplot(df_2_plot,
            x = "genotype",
            y = "value",
            color = "genotype",
            # shape = "clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  # scale_colour_manual(values = brewer.pal(n = 3,
  #                                         name = "Dark2")[c(2, 3)],
  #                     guide = "none") +


  # stat_compare_means(label = "p.signif",
  #                    label.y.npc = 0.95,
  #                    # label.x.npc = 0,
  #                    method = 'wilcox.test',
  #                    hide.ns = F,
  #                    ref.group = "non-risk",
  #                    paired = F) +
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
       y = "TREM2+, IBA1+, TMEM119+ cells") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                1.1),
                     labels = scales::percent) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("TREM2+, IBA1+, TMEM119+ cells")
