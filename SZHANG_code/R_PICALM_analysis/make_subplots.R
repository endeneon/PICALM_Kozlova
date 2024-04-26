# Siwei 11 Mar 2024
# Plot a list of genes from Alena's iMG expression data


# init ####
{
  library(edgeR)
  library(readr)
  library(readxl)

  library(Rfast)
  library(factoextra)
  library(dplyr)
  library(stringr)

  library(ggplot2)
  library(RColorBrewer)
  library(ggpubr)

  library(reshape2)

  library(sva)
}

# load data #####
load("~/backuped_space/Siwei_misc_R_projects/Alena_RNASeq_23Aug2023/Alena_volcano_plot_08Oct2023.RData")

## df_master_4_DGE_combat has been processed, #####
## should be able to use this df to build DGE

df_edgeR_DGE <-
  DGEList(counts = as.matrix(df_master_4_DGE_combat),
          samples = df_metadata_4_DGE$sample_name,
          group = df_metadata_4_DGE$sample_type,
          genes = rownames(df_master_4_DGE_combat),
          remove.zeros = T)

df_edgeR_DGE_sub <-
  df_edgeR_DGE[, str_detect(string = colnames(df_edgeR_DGE),
                            pattern = "off",
                            negate = T)]


sum(df_edgeR_DGE$genes$genes %in% c("HMGCS1",
                                    "MSMO1",
                                    "HMGCR",
                                    "DHCR7",
                                    "FDFT1",
                                    "PICALM"))
unique(df_edgeR_DGE_sub$samples$samples)

cpm_df_edgeR_DGE_sub <-
  as.data.frame(cpm(df_edgeR_DGE_sub,
                    normalized.lib.sizes = T))
rownames(cpm_df_edgeR_DGE_sub)


df_extracted_cpm_df_edgeR_DGE_sub <-
  cpm_gene_count[rownames(cpm_gene_count) %in% c("HMGCS1",
                                                 "MSMO1",
                                                 "HMGCR",
                                                 "DHCR7",
                                                 "FDFT1",
                                                 "PICALM",
                                                 "TREM2",
                                                 "ATP6AP2"), ]
df_extracted_cpm_df_edgeR_DGE_sub <-
  as.data.frame(t(df_extracted_cpm_df_edgeR_DGE_sub))
df_extracted_cpm_df_edgeR_DGE_sub <-
  df_extracted_cpm_df_edgeR_DGE_sub[order(rownames(df_extracted_cpm_df_edgeR_DGE_sub)), ]

df_metadata_4_plot <-
  df_metadata_raw
df_metadata_4_plot$Cell_line <-
  str_split(df_metadata_4_plot$sample_names,
            pattern = "-",
            simplify = T)[, 1]

df_metadata_4_plot <-
  df_metadata_4_plot[order(df_metadata_4_plot$sample_names), ]

df_extracted_cpm_df_edgeR_DGE_sub$sample_type <-
  df_metadata_4_plot$sample_type


melted_df_2_plot <-
  melt(df_extracted_cpm_df_edgeR_DGE_sub,
       value.name = "CPM")

melted_df_2_plot <-
  melted_df_2_plot[!(melted_df_2_plot$sample_type == "non-risk-CRISPRoff"), ]

ggplot(melted_df_2_plot,
       aes(x = factor(variable),
           y = log2(CPM),
           # group = factor(Genotype),
           fill = factor(sample_type))) +
  geom_dotplot(binaxis = "y",
               dotsize = 0.3,
               stackdir = "center",
               binwidth = 0.3,
               # method = "histodot",
               position = position_dodge(width = 0.75),
               stackgroups = T,
               right = F,
               colour = "#00000000") +
  # geom_boxplot(varwidth = F,
  #              position = "dodge",
  #              colour = "black",
  #              outlier.size = 1) +
  # stat_summary(fun = "mean",
  #              shape = 4,
  #              # colour = rep_len("darkorange",
  #              #                  length.out = 6),
  #              position = position_dodge(width = 0.75)) +
  stat_summary(fun.data = ggplot2::median_hilow,
               # shape = 4,
               # fun.args = list(mult = 1),
               geom = "errorbar",
               colour = "red",
               width = 0.2,
               # colour = rep_len("darkorange",
               #                  length.out = 6),
               position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 4,
               size = 3,
               colour = "darkblue",
               position = position_dodge(width = 0.75)) +
  # scale_fill_manual(values = brewer.pal(n = 3,
  #                                       name = "Dark2")) +
# geom_linerange(aes(ymin = after_stat(notchlower),
#                    ymax = after_stat(notchupper)),
#                position = "dodge")
# geom_pointrange(position = "dodge") +

  # = 21,) +
  labs(x = "",
       y = "Counts per million (CPM)") +
  ylim(0, 10) +
  # geom_errorbar(aes(ymin = M))
  # scale_colour_manual(values = brewer.pal(n = 4,
  #                                         name = "Dark2")[c(2,4)],
  #                     guide = guide_legend("Genotype")) +
  # # scale_shape_discrete() +
  scale_fill_manual(values = brewer.pal(n = 4,
                                        name = "Dark2")[c(3,1)],
                    guide = guide_legend("Genotype")) +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))

gene_list <-
  unique(melted_df_2_plot$variable)
# [1] MSMO1  PICALM FDFT1  TREM2  HMGCS1 HMGCR  DHCR7
i <- 1

for (i in 1:length(gene_list)) {
  df_subplot <-
    melted_df_2_plot[melted_df_2_plot$variable == gene_list[i], ]

  df_subplot$sample_type <-
    factor(df_subplot$sample_type,
           levels = c("non-risk",
                      "risk"))
  df_subplot$log2CPM <-
    log2(df_subplot$CPM)

  # ggdotplot(df_subplot,
  #           x = "sample_type",
  #           y = "log2CPM",
  #           color = "sample_type",
  #           fill = "white",
  #           palette = c("darkred",
  #                       "darkblue"),
  #           add = "mean_se")
  #
  pdf(file = paste0(gene_list[i],
                    '.pdf'),
      width = 1.8,
      height = 2.6,
      compress = T)

  output_graph <-
    ggerrorplot(df_subplot,
                x = "sample_type",
                y = "log2CPM",
                color = "black",
                fill = "white",
                width = 0.2,
                # colo
                palette = c("black",
                            "black"),
                # add = "jitter",
                error.plot = "errorbar") +

    scale_colour_manual(values = c("darkred",
                                   "darkblue"),
                        guide = "none") +
    stat_compare_means(label = "p.signif",
                       label.y.npc = .9,
                       hide.ns = T) +
    geom_jitter(aes(colour = sample_type),
                width = 0.05,
                size = 2,
                shape = 1) +
    stat_summary(geom = "crossbar",
                 fun = "mean",
                 width = 0.5,
                 linewidth = 0.2,
                 colour = "black") +
    labs(x = "",
         y = "log2(CPM)") +
    # ylim(0, 12) +
    scale_y_continuous(limits = c(0, 12),
                       expand = c(0, 0)) +
    theme_classic() +
    theme(axis.text = element_text(size = 10,
                                   colour = "black")) +
    ggtitle(gene_list[i])

  #
  # output_graph <-
  #   ggplot(df_subplot,
  #          aes(x = factor(sample_type),
  #              y = log2(CPM),
  #              # group = factor(Genotype),
  #              fill = factor(sample_type),
  #              colour = factor(sample_type))) +
  #   stat_summary(geom = "crossbar",
  #                fun = "mean",
  #                width = 0.5,
  #                linewidth = 0.5,
  #                colour = "black") +
  #   geom_jitter(width = 0.1,
  #               size = 2,
  #               shape = 1) +
  #   scale_fill_manual(values = c("darkred",
  #                                "darkblue"),
  #                     guide = "none") +
  #   scale_colour_manual(values = c("darkred",
  #                                  "darkblue"),
  #                       guide = "none") +
  #   labs(x = "",
  #        y = "log2(CPM)") +
  #   # ylim(0, 12) +
  #   scale_y_continuous(limits = c(0, 12),
  #                      expand = c(0, 0)) +
  #   theme_classic() +
  #   theme(axis.text = element_text(size = 10,
  #                                  colour = "black")) +
  #   ggtitle(gene_list[i])
  print(output_graph)
  dev.off()

}
