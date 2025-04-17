# Siwei 02 Feb 2025
# plot new Fig 7b-c

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

df_raw <-
  read_excel("Batch_2_of_new_plots/new_Fig_Ex_7b.xlsx",
             sheet = 1)

## CD04 #####
## make columns ####

df_sum_all <-
  df_raw[df_raw$Cell_line == "CD04", ]

df_sum_all$BR <-
  str_split(string = df_sum_all$BR_FOV,
            pattern = "_",
            simplify = T)[, 1]

df_sum_all$sum_data_table <-
  str_c(df_sum_all$BR,
        df_sum_all$Genotype,
        # df_sum_all$`Cell Line`,
        df_sum_all$Clone,
        df_sum_all$Batch,
        sep = "_")

dt_sum_all <-
  setDT(df_sum_all[, c(1:2, 9)])
dt_sum_mean <-
  dt_sum_all[, lapply(.SD, mean), keyby = "sum_data_table"]

dt_sum_mean$BR <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 1]
dt_sum_mean$genotype <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 2]
# dt_sum_mean$cell_line <-
#   str_split(dt_sum_mean$sum_data_table,
#             pattern = "_",
#             simplify = T)[, 3]
dt_sum_mean$clone_code <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 3]
dt_sum_mean$batch <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 4]

dt_sum_mean$genotype <-
  factor(dt_sum_mean$genotype,
         levels = c("risk",
                    "non-risk"))

dt_sum_mean$genotype <-
  factor(dt_sum_mean$genotype,
         levels = c("risk",
                    "non-risk"))

colnames(dt_sum_mean) <-
  make.names(colnames(dt_sum_mean))
colnames(dt_sum_mean)

dt_sum_mean <-
  setDF(dt_sum_mean)
df_2_plot_df <-
  dt_sum_mean

### LD area/iMG #####
ggerrorplot(df_2_plot_df,
            x = "genotype",
            y = "LD.area.iMG",
            # color = "genotype",
            shape = "clone_code",
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
                     # var.equal = F,
                     paired = F) +
  geom_jitter(aes(colour = genotype,
                  shape = clone_code),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "LD area/iMG") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        # legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD04 LD area/iMG")

### puncta density/iMG #####
ggerrorplot(df_2_plot_df,
            x = "genotype",
            y = "puncta.iMG",
            # color = "genotype",
            shape = "clone_code",
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
                     # var.equal = F,
                     paired = F) +
  geom_jitter(aes(colour = genotype,
                  shape = clone_code),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Puncta/iMG") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 60)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        # legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD04 Puncta/iMG")


## CD09 #####
## make columns ####

df_sum_all <-
  df_raw[df_raw$Cell_line == "CD09", ]

df_sum_all$BR <-
  str_split(string = df_sum_all$BR_FOV,
            pattern = "_",
            simplify = T)[, 1]

df_sum_all$sum_data_table <-
  str_c(df_sum_all$BR,
        df_sum_all$Genotype,
        # df_sum_all$`Cell Line`,
        df_sum_all$Clone,
        df_sum_all$Batch,
        sep = "_")

dt_sum_all <-
  setDT(df_sum_all[, c(1:2, 9)])
dt_sum_mean <-
  dt_sum_all[, lapply(.SD, mean), keyby = "sum_data_table"]

dt_sum_mean$BR <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 1]
dt_sum_mean$genotype <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 2]
# dt_sum_mean$cell_line <-
#   str_split(dt_sum_mean$sum_data_table,
#             pattern = "_",
#             simplify = T)[, 3]
dt_sum_mean$clone_code <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 3]
dt_sum_mean$batch <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 4]

dt_sum_mean$genotype <-
  factor(dt_sum_mean$genotype,
         levels = c("risk",
                    "non-risk"))

# dt_sum_mean$genotype <-
#   factor(dt_sum_mean$genotype,
#          levels = c("risk",
#                     "non-risk"))

colnames(dt_sum_mean) <-
  make.names(colnames(dt_sum_mean))
colnames(dt_sum_mean)

dt_sum_mean <-
  setDF(dt_sum_mean)
df_2_plot_df <-
  dt_sum_mean

### LD area/iMG #####
ggerrorplot(df_2_plot_df,
            x = "genotype",
            y = "LD.area.iMG",
            # color = "genotype",
            shape = "clone_code",
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
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
                     # var.equal = F,
                     paired = F) +
  geom_jitter(aes(colour = genotype,
                  shape = clone_code),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "LD area/iMG") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        # legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD09 LD area/iMG")

### puncta density/iMG #####
ggerrorplot(df_2_plot_df,
            x = "genotype",
            y = "puncta.iMG",
            # color = "genotype",
            shape = "clone_code",
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
                     # var.equal = F,
                     paired = F) +
  geom_jitter(aes(colour = genotype,
                  shape = clone_code),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Puncta/iMG") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 60)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        # legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD09 Puncta/iMG")
