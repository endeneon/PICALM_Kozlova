# make figures for Alena PICALM paper
# Siwei 22 May 2024

# init ####
{
  library(readxl)
  
  library(stringr)
  library(ggplot2)
  
  library(scales)
  library(reshape2)
  
  library(RColorBrewer)
  library(ggpubr)
}

## Western blot quantify ####
df_raw <-
  read_excel("table_western.xlsx",
             sheet = 1)

df_melted <-
  melt(df_raw, 
       variable.name = "Genotype")

df_2_plot <-
  df_melted[str_detect(string = df_melted$Genotype,
                       pattern = "^CD04.*"), ]
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c(unique(df_2_plot$Genotype)[2],
                    unique(df_2_plot$Genotype)[1]))


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     # ref.group = "non-risk",
                     paired = F) +
  geom_jitter(aes(colour = Genotype),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Relative abundance") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.2),
                     labels = waiver()) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 0)) +
  ggtitle("CD04")

#####
df_melted <-
  melt(df_raw, 
       variable.name = "Genotype")

df_2_plot <-
  df_melted[str_detect(string = df_melted$Genotype,
                       pattern = "^CD09.*"), ]
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c(unique(df_2_plot$Genotype)[2],
                    unique(df_2_plot$Genotype)[1]))


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     # ref.group = "non-risk",
                     paired = F) +
  geom_jitter(aes(colour = Genotype),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Relative abundance") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.2),
                     labels = waiver()) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 0)) +
  ggtitle("CD09")


# plot together
df_raw <-
  read_excel("table_western.xlsx",
             sheet = 2,
             col_names = T)

df_2_plot <- df_raw

# df_2_plot$Genotype <-
#   factor(df_2_plot$Genotype,
#          levels = c(unique(df_2_plot$Genotype)[2],
#                     unique(df_2_plot$Genotype)[1]))
# df_2_plot$Cell_Line <-
#   factor(df_2_plot$Cell_Line)

# df_2_plot$Genotype[df_2_plot$Genotype %in% "risk"] <- "darkred"
# df_2_plot$Genotype[df_2_plot$Genotype %in% "non-risk"] <- "darkblue"
# 
# df_2_plot$Cell_Line[df_2_plot$Cell_Line == "CD04"] <- 1
# df_2_plot$Cell_Line[df_2_plot$Cell_Line == "CD09"] <- 2

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))
df_2_plot$Cell_Line <-
  factor(df_2_plot$Cell_Line)

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "value",
            # shape = "Cell_Line",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            add = "jitter",
            # facet.by = "Gene",
            # ncol = 4,
            add.params = list(color = "Genotype",
                              shape = "Cell_Line"),
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  scale_shape_manual(values = c(1, 2)) +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     # ref.group = "non-risk",
                     paired = F) +
  # geom_jitter(aes(colour = Genotype),
  #             width = 0.05,
  #             size = 2,
  #             shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Relative abundance") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.2),
                     labels = waiver()) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 0)) 
