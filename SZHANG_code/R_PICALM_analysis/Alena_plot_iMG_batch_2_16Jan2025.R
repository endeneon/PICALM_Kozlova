# Siwei 16 Jan 2025
# reshape Alena's new batch 2 data and plot
{
  library(stringr)
  # library(Seurat)

  library(parallel)
  library(future)


  # library(glmGamPoi)

  # library(edgeR)

  library(data.table)

  library(readr)
  library(readxl)

  library(dplyr)
  library(tidyr)

  library(reshape2)
  library(scales)

  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)

  library(gridExtra)


  plan("multisession", workers = 3)
  # options(mc.cores = 32)
  set.seed(42)
  options(future.globals.maxSize = 429496729600)
}

# load raw data and reshape ####

df_raw_b4_reshape <-
  read_excel(path = "Batch_2_of_new_plots/Alena_new_LD_iMG_4_plot_16Jan2025.xlsx",
             sheet = 1)

df_reshape <-
  df_raw_b4_reshape

df_reshape$meta_space_stripped <-
  str_replace_all(df_reshape$line_batch_risk,
                  pattern = " ",
                  replacement = "_")

df_reshape$uuid_by_diff <-
  paste0(df_reshape$meta_space_stripped,
         "_",
         str_split(df_reshape$BR_FOV,
                   pattern = "_",
                   simplify = T)[, 1])

colnames(df_reshape) <-
  str_replace_all(colnames(df_reshape),
                  pattern = " ",
                  replacement = "_")
colnames(df_reshape) <-
  str_replace_all(colnames(df_reshape),
                  pattern = "\\/",
                  replacement = "_")


df_by_data_type <-
  vector(mode = "list",
         length = 4L)
names(df_by_data_type) <-
  colnames(df_reshape)[3:6]

df_by_data_type[["LD_area_iMG"]] <-
  df_reshape %>%
  group_by(uuid_by_diff) %>%
  dplyr::summarise(Mean_LD_area = mean(LD_area_iMG,
                                       na.rm = T))
df_by_data_type[["puncta_microglia"]] <-
  df_reshape %>%
  group_by(uuid_by_diff) %>%
  dplyr::summarise(puncta_microglia = mean(puncta_microglia,
                                           na.rm = T))
df_by_data_type[["Fl_int_iMG"]] <-
  df_reshape %>%
  group_by(uuid_by_diff) %>%
  dplyr::summarise(Fl_int_iMG = mean(Fl_int_iMG,
                                     na.rm = T))
df_by_data_type[["Ave_size"]] <-
  df_reshape %>%
  group_by(uuid_by_diff) %>%
  dplyr::summarise(Fl_int_iMG = mean(Ave_size,
                                     na.rm = T))

func_split_metadata <-
  function(df) {
    returned_df <- df
    returned_df$cell_line <-
      str_split(string = returned_df$uuid_by_diff,
                pattern = "_",
                simplify = T)[, 1]
    returned_df$pheno <-
      str_split(string = returned_df$uuid_by_diff,
                pattern = "_",
                simplify = T)[, 2]

    returned_df$clone <-
      str_split(string = returned_df$uuid_by_diff,
                pattern = "_",
                simplify = T)[, 3]
    # returned_df$clone[returned_df$clone == unique(returned_df$clone)[1]] <-
    #   "Clone 1"
    # returned_df$clone[returned_df$clone == unique(returned_df$clone)[2]] <-
    #   "Clone 2"

    returned_df$value <-
      unlist(returned_df[, 2])
    return(returned_df)
  }


for (i in 1:length(df_by_data_type)) {
  print(i)
  df_by_data_type[[i]] <-
    func_split_metadata(df_by_data_type[[i]])

}

plot_by_category <-
  function(x) {

    df_2_plot <- x
    df_2_plot$pheno <-
      factor(df_2_plot$pheno,
             levels = c("risk",
                        "non-risk"))

    df_2_plot_CD04 <-
      df_2_plot[df_2_plot$cell_line == "CD04", ]
    df_2_plot_CD04$clone_id <-
      df_2_plot_CD04$clone

    df_2_plot_CD04$clone_id[df_2_plot_CD04$clone == unique(df_2_plot_CD04$clone[df_2_plot_CD04$pheno == "non-risk"])[1]] <-
      "Clone 1"
    df_2_plot_CD04$clone_id[df_2_plot_CD04$clone == unique(df_2_plot_CD04$clone[df_2_plot_CD04$pheno == "non-risk"])[2]] <-
      "Clone 2"
    df_2_plot_CD04$clone_id[df_2_plot_CD04$clone == unique(df_2_plot_CD04$clone[df_2_plot_CD04$pheno == "risk"])[1]] <-
      "Clone 1"
    df_2_plot_CD04$clone_id[df_2_plot_CD04$clone == unique(df_2_plot_CD04$clone[df_2_plot_CD04$pheno == "risk"])[2]] <-
      "Clone 2"

    df_2_plot_CD09 <-
      df_2_plot[df_2_plot$cell_line == "CD09", ]
    df_2_plot_CD09$clone_id <-
      df_2_plot_CD09$clone

    df_2_plot_CD09$clone_id[df_2_plot_CD09$clone == unique(df_2_plot_CD09$clone[df_2_plot_CD09$pheno == "non-risk"])[1]] <-
      "Clone 1"
    df_2_plot_CD09$clone_id[df_2_plot_CD09$clone == unique(df_2_plot_CD09$clone[df_2_plot_CD09$pheno == "non-risk"])[2]] <-
      "Clone 2"
    df_2_plot_CD09$clone_id[df_2_plot_CD09$clone == unique(df_2_plot_CD09$clone[df_2_plot_CD09$pheno == "risk"])[1]] <-
      "Clone 1"
    df_2_plot_CD09$clone_id[df_2_plot_CD09$clone == unique(df_2_plot_CD09$clone[df_2_plot_CD09$pheno == "risk"])[2]] <-
      "Clone 2"

    lm_model <-
      lmerTest::lmer(value ~
                       pheno +
                       (1|Clone/BR_ind),
                     data = df_2_plot_CD04)
    print(summary(lm_model))

    # panel_CD04 <-
    #   ggerrorplot(df_2_plot_CD04,
    #               x = "pheno",
    #               y = "value",
    #               color = "black",
    #               # shape = "clone_id",
    #               # group = "Genotype",
    #               width = 0.2,
    #               # facet.by = "Gene",
    #               # ncol = 4,
    #               error.plot = "errorbar") +
    #   scale_colour_manual(values = c("darkred",
    #                                  "darkblue"),
    #                       guide = "none") +
    #
    #
    #   stat_compare_means(label = "p.signif",
    #                      label.y.npc = 0.95,
    #                      # label.x.npc = 0,
    #                      method = 't.test',
    #                      hide.ns = F,
    #                      ref.group = "non-risk",
    #                      paired = F) +
    #   geom_jitter(aes(colour = pheno,
    #                   shape = clone_id),
    #               width = 0.1,
    #               size = 1) +
    #   scale_shape_manual(values = c(1, 2),
    #                      guide = "none") +
    #   stat_summary(geom = "crossbar",
    #                fun = "mean",
    #                width = 0.5,
    #                linewidth = 0.2,
    #                colour = "black") +
    #   labs(x = "",
    #        y = colnames(df_2_plot_CD04)[2]) +
    #   # ylim(0, 12) +
    #   scale_y_continuous(expand = c(0, 0),
    #                      limits = c(0,
    #                                 max(df_2_plot$value,
    #                                     na.rm = T) * 1.1)) +
    #   # stat_compare_means(label = "p.signif")
    #   theme_classic() +
    #   theme(axis.text = element_text(size = 10,
    #                                  colour = "black"),
    #         axis.text.x = element_text(angle = 315,
    #                                    hjust = 0,
    #                                    vjust = 1)) +
    #   ggtitle("CD04")
    #
    # panel_CD09 <-
    #   ggerrorplot(df_2_plot_CD09,
    #               x = "pheno",
    #               y = "value",
    #               color = "black",
    #               # shape = "clone_id",
    #               # group = "Genotype",
    #               width = 0.2,
    #               # facet.by = "Gene",
    #               # ncol = 4,
    #               error.plot = "errorbar") +
    #   scale_colour_manual(values = c("darkred",
    #                                  "darkblue"),
    #                       guide = "none") +
    #
    #
    #   stat_compare_means(label = "p.signif",
    #                      label.y.npc = 0.95,
    #                      # label.x.npc = 0,
    #                      method = 't.test',
    #                      hide.ns = F,
    #                      ref.group = "non-risk",
    #                      paired = F) +
    #   geom_jitter(aes(colour = pheno,
    #                   shape = clone_id),
    #               width = 0.1,
    #               size = 1) +
    #   scale_shape_manual(values = c(1, 2)) +
    #   stat_summary(geom = "crossbar",
    #                fun = "mean",
    #                width = 0.5,
    #                linewidth = 0.2,
    #                colour = "black") +
    #   labs(x = "",
    #        y = "") +
    #   # ylim(0, 12) +
    #   scale_y_continuous(expand = c(0, 0),
    #                      limits = c(0,
    #                                 max(df_2_plot$value,
    #                                     na.rm = T) * 1.1)) +
    #   # stat_compare_means(label = "p.signif")
    #   theme_classic() +
    #   theme(axis.text = element_text(size = 10,
    #                                  colour = "black"),
    #         axis.text.x = element_text(angle = 315,
    #                                    hjust = 0,
    #                                    vjust = 1)) +
    #   ggtitle("CD09")
    #
    # # panel_assembled <-
    # return(grid.arrange(panel_CD04,
    #                     panel_CD09,
    #                     widths = c(1, 2),
    #                     ncol = 2))
    lm_model <-
      lmerTest::lmer(LD.area.iMG ~
                       Genotype +
                       (1|Clone/BR_ind),
                     data = df_2_plot_df)
    print(summary(lm_model))


  }


dev.off()

for (i in 1:length(df_by_data_type)) {
  print(names(df_by_data_type)[i])
  # pdf(file = paste0("batch_2_new_",
  #                   names(df_by_data_type)[i],
  #                   ".pdf"),
  #     width = 3.5,
  #     height = 3)
  plot_by_category(x = df_by_data_type[[i]])
  dev.off()

}
# panel_assembled
#
# print(panel_assembled)
