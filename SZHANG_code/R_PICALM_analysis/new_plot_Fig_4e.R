# Siwei 31 Jan 2025
# plot Fig. 4e

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

  library(ggrepel)
}

set.seed(42)

#
df_raw <-
  read_excel(path = "Batch_2_of_new_plots/Fig_4e.xlsx")

df_2_plot <-
  df_raw

df_2_plot$label <- ""
df_2_plot$label[df_2_plot$gene...1 %in% c("SULF1",
                                          "CH25H",
                                          "TCIRG1",
                                          "VAMP1",
                                          "OAS3",
                                          "CD74",
                                          "MKI67")] <-
  df_2_plot$gene...1[df_2_plot$gene...1 %in% c("SULF1",
                                               "CH25H",
                                               "TCIRG1",
                                               "VAMP1",
                                               "OAS3",
                                               "CD74",
                                               "MKI67")]

ggplot(df_2_plot,
       aes(x = logFC,
           y = mouse_log2_fc,
           label = label)) +
  geom_point(colour = "darkred") +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(-3, 4)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(-2.5, 3)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = lm,
              se = F,
              fullrange = T,
              linetype = 2) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text = element_text(size = 12)) +
  geom_text_repel(min.segment.length = 0)

  # scale_x_continuous()
