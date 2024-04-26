# Siwei 07 Mar 2024
# plot FRiP sumstats of 5 cell types

# init ####
library(readr)
library(ggplot2)

library(RColorBrewer)

library(stringr)
library(reshape2)
library(plyr)

# load raw data ####
df_raw <-
  read_delim("FRiP/FRiP_sumstat.sumstat",
             delim = "\t", escape_double = FALSE,
             col_names = T, trim_ws = TRUE)

df_2_plot <-
  df_raw[!(df_raw$total_reads_peaks == "total_reads_peaks"), ]
colnames(df_2_plot)[1] <- "Sample_name"
df_2_plot$FRiP <-
  as.numeric(df_2_plot$total_reads_peaks) / as.numeric(df_2_plot$total_bam_reads)


df_2_plot$Samples <-
  str_split(string = df_2_plot$Sample_name,
            pattern = "\\/",
            simplify = T)[, 3]

df_2_plot$Samples <-
  sapply(X = 1:nrow(df_2_plot),
        # MARGIN = 1,
        FUN = function(x) {
          vector_2_split <-
            df_2_plot[[1]][x]
          vector_2_split <-
            unlist(str_split(string = vector_2_split,
                             pattern = "\\/",
                             simplify = T))
          return_value <-
            vector_2_split[length(vector_2_split)]
          return(return_value)
        })
df_2_plot$Samples <-
  str_split(string = df_2_plot$Samples,
            pattern = "_S|\\-S",
            simplify = T)[, 1]
df_2_plot$Samples <-
  str_split(string = df_2_plot$Samples,
            pattern = "_10x",
            simplify = T)[, 1]
df_2_plot <-
  df_2_plot[order(df_2_plot$Samples), ]
df_2_plot <-
  df_2_plot[!duplicated(df_2_plot$Samples), ]
df_2_plot$cell_type <-
  str_split(string = df_2_plot$Samples,
            pattern = "\\-|\\_",
            simplify = T)[, 1]

df_2_plot$cell_type[df_2_plot$cell_type == "R21"] <- "Glut"

df_2_plot$FRiP <-
  as.numeric(df_2_plot$FRiP)

unique(df_2_plot$cell_type)

df_2_plot$cell_type <-
  plyr::mapvalues(x = df_2_plot$cell_type,
                   from = unique(df_2_plot$cell_type),
                   to = c("iAst",
                          "iDN", "iGA", "iGlut",
                          "iMG"))

df_2_plot$cell_type <-
  factor(df_2_plot$cell_type,
         levels = c("iAst", "iMG",
                    "iDN", "iGA", "iGlut"))

ggplot(df_2_plot,
       aes(x = cell_type,
           y = FRiP * 100,
           group = cell_type,
           fill = cell_type,
           colour = cell_type)) +
  geom_boxplot(alpha = 0.4,
               linewidth = 0.5,
               outlier.colour = "transparent") +
  # stat_summary(geom = "crossbar",
  #              fun = "median",
  #              width = 0.5,
  #              linewidth = 0.5,
  #              colour = "black") +
  # geom_linerange() +
  # geom_jitter(stat = "identity",
  #             width = 0.1,
  #             size = 0.5,
  #             colour = "black",
  #             fill = "black") +
  stat_summary(geom = "linerange",
               fun.min = "min",
               fun.max = "max",
               # colour = "darkred",
               size = 1) +
  stat_summary(geom = "crossbar",
               fun = "min",
               # colour = "darkred",
               size = 0.5,
               width = 0.5) +
  stat_summary(geom = "crossbar",
               fun = "max",
               # colour = "darkred",
               size = 0.5,
               width = 0.5) +
  # geom_dotplot(binaxis = "y",
  #              stackdir = "center",
  #              dotsize = 0.2,
  #              binwidth = 5,
  #              # shape = 20,
  #              colour = "black",
  #              fill = "black") +
  stat_summary(geom = "point",
               fun = "mean",
               shape = 23,
               size = 2,
               colour = "darkred",
               fill = "orange") +

  scale_fill_manual(values = brewer.pal(n = 5,
                                        name = "Dark2"),
                    guide = "none") +
  scale_colour_manual(values = brewer.pal(n = 5,
                                        name = "Dark2"),
                    guide = "none") +
  # ylim(0, 100) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100)) +
  labs(x = "Cell Type",
       y = "FRiP Percentage") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black")) #+
  ggtitle("FRiP %")
