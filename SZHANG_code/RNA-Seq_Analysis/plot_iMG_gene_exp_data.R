# Siwei 08 Mar 2024
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
  cpm_df_edgeR_DGE_sub[rownames(cpm_df_edgeR_DGE_sub) %in% c("HMGCS1",
                                                             "MSMO1",
                                                             "HMGCR",
                                                             "DHCR7",
                                                             "FDFT1",
                                                             "PICALM"), ]
df_extracted_cpm_df_edgeR_DGE_sub <-
  as.data.frame(t(df_extracted_cpm_df_edgeR_DGE_sub))
df_extracted_cpm_df_edgeR_DGE_sub$Genotype <- "G"
df_extracted_cpm_df_edgeR_DGE_sub$Genotype[str_detect(string = rownames(df_extracted_cpm_df_edgeR_DGE_sub),
                                                      pattern = "\\-A")] <- "A"

df_extracted_cpm_df_edgeR_DGE_sub$Cell_line <- "Line_09"
df_extracted_cpm_df_edgeR_DGE_sub$Cell_line[str_detect(string = rownames(df_extracted_cpm_df_edgeR_DGE_sub),
                                                      pattern = "^CD04")] <- "Line_04"


melted_df_2_plot <-
  melt(df_extracted_cpm_df_edgeR_DGE_sub,
       value.name = "CPM")

ggplot(melted_df_2_plot,
        aes(x = factor(variable),
            y = CPM,
            # group = factor(Genotype),
            fill = factor(Genotype),
            colour = factor(Cell_line))) +
  # geom_boxplot(varwidth = F,
  #              position = "dodge",
  #              colour = "black",
  #              outlier.size = 1) +
  # stat_summary_bin(fun = "mean",
  #              shape = 4,
  #              # colour = rep_len("darkorange",
  #              #                  length.out = 6),
  #              position = position_dodge(width = 0.75)) +
  # scale_fill_manual(values = brewer.pal(n = 3,
  #                                       name = "Dark2")) +
  # geom_linerange(aes(ymin = after_stat(notchlower),
  #                    ymax = after_stat(notchupper)),
  #                position = "dodge")
  # geom_pointrange(position = "dodge") +
  geom_dotplot(binaxis = "y",
               dotsize = 0.5,
               stackdir = "center",
               # method = "histodot",
               position = position_dodge(width = 0.75),
               stackgroups = T,
               right = F) +
                # = 21,) +
  labs(x = "",
       y = "Counts per million (CPM)") +
  scale_colour_manual(values = brewer.pal(n = 4,
                                          name = "Dark2")[c(2,4)],
                      guide = guide_legend("Genotype")) +
  # scale_shape_discrete() +
  scale_fill_manual(values = brewer.pal(n = 4,
                                        name = "Dark2")[c(1,3)],
                    guide = guide_legend("Genotype")) +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))

