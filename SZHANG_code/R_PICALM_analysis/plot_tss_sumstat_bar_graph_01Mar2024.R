# Siwei 01 Mar 2024
# plot TSS sumstats of 5 cell types

# init ####
library(readr)
library(ggplot2)

library(RColorBrewer)

library(stringr)

# load raw data ####
df_raw <-
  read_delim("tss_sumstats/tss_sumstat.sumstat",
             delim = "\t", escape_double = FALSE,
             col_names = FALSE, trim_ws = TRUE)
colnames(df_raw) <-
  c("Sample", "TSS_enrich")

df_2_plot <-
  df_raw
df_2_plot <-
  df_2_plot[str_detect(string = df_2_plot$Sample,
                       pattern = "^CN|^NSC",
                       negate = T), ]
df_2_plot$Cell_type <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\-|\\_",
            simplify = T)[, 1]
# table(df_2_plot$Cell_type)

df_2_plot$Cell_type[df_2_plot$Cell_type == "R21"] <- "Glut"
unique(df_2_plot$Cell_type)
# table(df_2_plot$Cell_type)

df_2_plot$Cell_type <-
  plyr::mapvalues(x = df_2_plot$Cell_type,
                  from = unique(df_2_plot$Cell_type),
                  to = c("iAst",
                         "iDN", "iGA", "iGlut",
                         "iMG"))

df_2_plot$Cell_type <-
  factor(df_2_plot$Cell_type,
         levels = c("iAst", "iMG",
                    "iDN", "iGA", "iGlut"))

# df_2_plot$Cell_type <-
#   factor(df_2_plot$Cell_type,
#          levels = sort(unique(df_2_plot$Cell_type)))


ggplot(df_2_plot,
       aes(x = Cell_type,
           y = TSS_enrich,
           group = Cell_type,
           fill = Cell_type,
           colour = Cell_type)) +
  geom_boxplot(alpha = 0.4,
               linewidth = 0.5,
               outlier.colour = "transparent") +
  # stat_summary(geom = "crossbar",
  #              fun = "median",
  #              width = 0.5,
  #              linewidth = 0.5,
  #              colour = "black") +
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
  # geom_linerange() +
  # stat_summary(geom = "linerange",
  #              fun.min = "min",
  #              fun.max = "max",
  #              # colour = "darkred",
  #              size = 1) +
  stat_summary(geom = "point",
               fun = "mean",
               shape = 23,
               size = 2,
               colour = "darkred",
               fill = "orange") +
  # geom_jitter(stat = "identity",
  #             width = 0.1,
  #             size = 0.5,
  #             colour = "black",
  #             fill = "black") +

  # geom_dotplot(binaxis = "y",
  #              stackdir = "center",
  #              dotsize = 0.2,
  #              binwidth = 5,
  #              # shape = 20,
  #              colour = "black",
  #              fill = "black") +


  scale_fill_manual(values = brewer.pal(n = 5,
                                        name = "Dark2"),
                    guide = "none") +
  scale_colour_manual(values = brewer.pal(n = 5,
                                          name = "Dark2"),
                      guide = "none") +
  # ylim(0, 100) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 10)) +
  labs(x = "Cell Type",
       y = "TSS enrichment fold") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"))

ggplot(df_2_plot,
       aes(x = Cell_type,
           y = TSS_enrich,
           group = Cell_type,
           fill = Cell_type)) +
  geom_boxplot(alpha = 0.7) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.2,
               # shape = 20,
               colour = "black",
               fill = "black") +
  stat_summary(geom = "point",
               fun = "mean",
               shape = 4,
               size = 2,
               colour = "darkred") +
  scale_fill_manual(values = brewer.pal(n = 5,
                                        name = "Dark2")) +
  ylim(0, 10) +
  labs(x = "Cell Type",
       y = "TSS enrichment fold") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"))
