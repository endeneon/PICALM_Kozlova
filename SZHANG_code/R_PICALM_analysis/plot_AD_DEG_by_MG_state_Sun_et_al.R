# Siwei 10 Feb 2025
# plot use Sun et al. 2023

# init ####
{
  library(ggplot2)
  library(stringr)
  library(lsr)
  library(RColorBrewer)
  library(ggpubr)

  library(readr)
  library(readxl)

  library(stats)
  library(rstatix)

  library(agricolae)
  library(DescTools)

  # library

  library(data.table)

}


df_raw <-
  read_excel("Batch_2_of_new_plots/Sun_et al_Cell 2023_PICALM-AD-DEG-by_MG_state.xlsx")
df_raw$row.names <- NULL

# df_2_plot <-
#   reshape2::melt(df_raw)

df_2_plot <-
  df_raw
df_2_plot$log2FDR <-
  0 - log10(df_2_plot$fdr)

ggplot(df_2_plot,
       aes(x = groupID,
           y = coef)) +


  scale_size(limits = c(0, 0.1),
                  breaks = c(0, 0.3, 0.6, 0.9) / 10,
                  labels = c("0", "10", "20", "30"),
                  name = "-log10FDR") +
  geom_linerange(aes(ymax = ci.hi,
                     ymin = ci.lo),
                 colour = "black",
                 size = 1) +
  geom_point(aes(size = log2FDR / 500),
             shape = 21,
             colour = "black",
             fill = "darkred") +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(-0.25, 0.75)) +
  coord_flip() +
  geom_hline(yintercept = 0) +
  theme_classic()

max(df_2_plot$log2FDR)
