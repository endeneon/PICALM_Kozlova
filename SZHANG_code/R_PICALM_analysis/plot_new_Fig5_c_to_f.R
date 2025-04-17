# Siwei 26 Jan 2025
# plot new Fig 5c-f

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
## make columns ####
df_raw <-
  read_excel("Batch_2_of_new_plots/Alena_new_LD_iMG_4_plot_16Jan2025.xlsx",
             sheet = 1)

df_sum_all <-
  df_raw
df_sum_all$BR <-
  str_split(string = df_sum_all$BR_FOV,
            pattern = "_",
            simplify = T)[, 1]

df_sum_all$Cell_line <-
  str_split(string = df_sum_all$line_batch_risk,
            pattern = " ",
            simplify = T)[, 1]
df_sum_all$genotype <-
  str_split(string = df_sum_all$line_batch_risk,
            pattern = " ",
            simplify = T)[, 2]
df_sum_all$clone_code <-
  str_split(string = df_sum_all$line_batch_risk,
            pattern = " ",
            simplify = T)[, 3]
df_sum_all$batch <-
  str_split(string = df_sum_all$line_batch_risk,
            pattern = " ",
            simplify = T)[, 4]

unique(df_sum_all$clone_code)

# lme4 ####


df_sum_all$sum_data_table <-
  str_c(df_sum_all$BR,
        df_sum_all$genotype,
        df_sum_all$Cell_line,
        df_sum_all$clone_code,
        df_sum_all$batch,
        sep = "_")

dt_sum_all <-
  setDT(df_sum_all[, c(3:6, 12)])

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
dt_sum_mean$cell_line <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 3]
dt_sum_mean$clone_code <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 4]

dt_sum_mean$Clones <-
  rep_len(rep_len(c("Clone 1",
                    "Clone 2"),
                  length.out = 8),
          length.out = nrow(dt_sum_mean))


dt_sum_mean$genotype <-
  factor(dt_sum_mean$genotype,
         levels = c("risk",
                    "non-risk"))

colnames(dt_sum_mean) <-
  make.names(colnames(dt_sum_mean))
colnames(dt_sum_mean)

dt_sum_mean <-
  setDF(dt_sum_mean)
# df_2_plot <-
#   dt_sum_mean[dt_sum_mean$cell_line == unique(dt_sum_mean$cell_line)[1]]
# colnames(df_2_plot) <-
#   make.names(colnames(df_2_plot))
#
# df_2_plot$genotype <-
#   factor(df_2_plot$genotype,
#          levels = c("risk",
#                     "non-risk"))
# colnames(df_2_plot) <-
#   make.names(colnames(df_2_plot))
# colnames(df_2_plot) <-
#   str_replace_all(string = colnames(df_2_plot),
#                   pattern = "\\.",
#                   replacement = "_")

# try lme4 here
dt_line_mean <-
  dt_sum_mean[dt_sum_mean$cell_line == "CD04", ]
dt_line_mean$Experiments <-
  factor(c(1, 2))

# lm_model <-
#   lmerTest::lmer(puncta.microglia ~
#                    genotype +
#                    (1 | Experiments),
#                  # (2 + Genotype|diff),
#                  data = dt_line_mean)
lm_model <-
  lmerTest::lmer(puncta.microglia ~
                   genotype +
                   (1 | Clones),
                 # (2 + Genotype|diff),
                 data = dt_line_mean)

summary(lm_model)

lm_model <-
  lmerTest::lmer(Ave.size ~
                   genotype +
                   (1 | Clones),
                 # (2 + Genotype|diff),
                 data = dt_line_mean)

summary(lm_model)
# anova(lm_model)

lm_model <-
  lmerTest::lmer(LD.area.iMG ~
                   genotype +
                   (1 | Clones),
                 # (2 + Genotype|diff),
                 data = dt_line_mean)

summary(lm_model)

lm_model <-
  lmerTest::lmer(Fl.int.iMG ~
                   genotype +
                   (1 | Clones),
                 # (2 + Genotype|diff),
                 data = dt_line_mean)

summary(lm_model)

## CD09 #####

dt_line_mean <-
  dt_sum_mean[dt_sum_mean$cell_line == "CD09", ]
dt_line_mean$Experiments <-
  factor(c(1, 2))

# lm_model <-
#   lmerTest::lmer(puncta.microglia ~
#                    genotype +
#                    (1 | Experiments),
#                  # (2 + Genotype|diff),
#                  data = dt_line_mean)
lm_model <-
  lmerTest::lmer(puncta.microglia ~
                   genotype +
                   (1 | Clones),
                 # (2 + Genotype|diff),
                 data = dt_line_mean)



summary(lm_model)
anova(lm_model)

lm_model <-
  lmerTest::lmer(Ave.size ~
                   genotype +
                   (1 | Clones),
                 # (2 + Genotype|diff),
                 data = dt_line_mean)

summary(lm_model)
# anova(lm_model)

lm_model <-
  lmerTest::lmer(LD.area.iMG ~
                   genotype +
                   (1 | Clones),
                 # (2 + Genotype|diff),
                 data = dt_line_mean)

summary(lm_model)

lm_model <-
  lmerTest::lmer(Fl.int.iMG ~
                   genotype +
                   (1 | Clones),
                 # (2 + Genotype|diff),
                 data = dt_line_mean)

summary(lm_model)





df_test <-
  data.frame(clone_number = c(1,1,1,2,2,2,1,1,1,2,2,2),
             value = c(1:12),
             group = c(1,1,1,1,1,1,2,2,2,2,2,2))
t.test(x = df_test$value[1:6],
       y = df_test$value[7:12])

df_test$clone_number <-
  as.factor(df_test$clone_number)
df_test$group <-
  as.factor(df_test$group)
lm_model <-
  lmerTest::lmer(value ~ group +
                   (1 | clone_number),
                 data = df_test)

summary(lm_model)
anova(lm_model)
## make a plotting function ####

### make_sub_plot ####
make_sub_plot <-
  function(i) {
    df_2_plot_sub <-
      dt_sum_mean[dt_sum_mean$cell_line == unique(dt_sum_mean$cell_line)[1], ]
    # colnames(df_2_plot_sub) <-
    #   make.names(colnames(df_2_plot_sub))

    df_2_plot_sub$genotype <-
      factor(df_2_plot_sub$genotype,
             levels = c("risk",
                        "non-risk"))

#
    df_2_plot_df <-
      df_2_plot_sub

    y_label <-
      colnames(df_2_plot_df)[i + 1]



    print(colnames(df_2_plot_df))

    p_CD04 <-
      ggerrorplot(df_2_plot_df,
                  x = "genotype",
                  y = y_label,
                  color = "genotype",
                  shape = "Clones",
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
                      shape = Clones),
                  width = 0.1,
                  size = 1) +
      scale_shape_manual(values = c(1, 2)) +
      stat_summary(geom = "crossbar",
                   fun = "mean",
                   width = 0.5,
                   linewidth = 0.2,
                   colour = "black") +
      labs(x = "",
           y = y_label) +
      scale_y_continuous(expand = c(0, 0),
                         limits = c(0,
                                    max(dt_sum_mean[, (i + 1)],
                                        na.rm = T) * 1.1)) +
      # stat_compare_means(label = "p.signif")
      theme_classic() +
      theme(axis.text = element_text(size = 12,
                                     colour = "black"),
            legend.position = "none",
            axis.text.x = element_text(angle = 315,
                                       hjust = 0,
                                       vjust = 1)) +
      ggtitle(unique(dt_sum_mean$cell_line)[1])


    ### CD09
    df_2_plot_sub <-
      dt_sum_mean[dt_sum_mean$cell_line == unique(dt_sum_mean$cell_line)[2], ]
    # colnames(df_2_plot_sub) <-
    #   make.names(colnames(df_2_plot_sub))

    df_2_plot_sub$genotype <-
      factor(df_2_plot_sub$genotype,
             levels = c("risk",
                        "non-risk"))

    y_label <-
      colnames(df_2_plot_sub)[i + 1]

    df_2_plot_df <-
      df_2_plot_sub
    # df_2_plot_df <-
    #   setDF(df_2_plot_sub)

    p_CD09 <-
      ggerrorplot(df_2_plot_df,
                  x = "genotype",
                  y = y_label,
                  color = "genotype",
                  shape = "Clones",
                  # group = "Genotype",
                  width = 0.2,
                  # facet.by = "Gene",
                  # ncol = 4,
                  error.plot = "errorbar",
                  legend.position = "right") +
      scale_colour_manual(values = c("darkred",
                                     "darkblue")) +
      stat_compare_means(label = "p.signif",
                         label.y.npc = 0.95,
                         # label.x.npc = 0,
                         method = 'wilcox.test',
                         hide.ns = F,
                         ref.group = "non-risk",
                         paired = F) +
      geom_jitter(aes(colour = genotype,
                      shape = Clones),
                  width = 0.1,
                  size = 1) +
      scale_shape_manual(values = c(1, 2)) +
      stat_summary(geom = "crossbar",
                   fun = "mean",
                   width = 0.5,
                   linewidth = 0.2,
                   colour = "black") +
      labs(x = "",
           y = "") +
      scale_y_continuous(expand = c(0, 0),
                         limits = c(0,
                                    max(dt_sum_mean[, (i + 1)],
                                        na.rm = T) * 1.1)) +
      # stat_compare_means(label = "p.signif")
      theme_classic() +
      theme(axis.text = element_text(size = 12,
                                     colour = "black"),
            # legend.position = "none",
            axis.text.x = element_text(angle = 315,
                                       hjust = 0,
                                       vjust = 1)) +
      ggtitle(unique(dt_sum_mean$cell_line)[2])

    p_return <-
      gridExtra::grid.arrange(p_CD04,
                              p_CD09,
                              ncol = 2)

    return(p_return)

    # ggerrorplot(df_2_plot_df,
    #             x = "genotype",
    #             y = "LD_area_iMG")

  }


print(make_sub_plot(i = 1)) # 5F
print(make_sub_plot(i = 2)) # 5D
print(make_sub_plot(i = 3)) # 5G
print(make_sub_plot(i = 4)) # 5E
# plot_2_side_by_side <-
#   function(x) {
#     y_label <-
#       colnames(df_2_plot)[x]
#
#
#   }


## LD area per iMG #####
y_label <-
  colnames(df_2_plot)[2]

ggerrorplot(df_2_plot_df,
            x = "genotype",
            y = y_label,
            color = "genotype",
            shape = "Clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue")) +
  stat_compare_means(label = "p.signif",
                     label.y.npc = 0.95,
                     # label.x.npc = 0,
                     method = 'wilcox.test',
                     hide.ns = F,
                     ref.group = "non-risk",
                     paired = F) +
  geom_jitter(aes(colour = genotype,
                  shape = Clones),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(dt_sum_mean[, (i + 1)],
                                    na.rm = T) * 1.1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle(unique(dt_sum_mean$cell_line)[2])
