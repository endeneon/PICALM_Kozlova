# make figures for Alena PICALM paper
# Siwei 18 Apr 2024

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

## Fig_5F ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 1)

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))

df_2_plot$Cell_line <-
  factor(df_2_plot$Cell_line,
         levels = c("CD04",
                    "CD09"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD04", ]
# df_2_plot <-
#   df_2_plot[!(df_2_plot$Genotype == "risk"), ]

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
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD puncta/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 200),
                     labels = waiver()) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD04")

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD09", ]

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
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD area/iMG (um2)") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 200)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD09")

## Fig_5D ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 2)

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))
# 
# df_2_plot$Cell_line <-
#   factor(df_2_plot$Cell_line,
#          levels = c("CD04",
#                     "CD09"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD04", ]

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
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD puncta/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 200),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD04")

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD09", ]

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
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD puncta/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 200)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD09")

## Fig_5E ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 3)

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
colnames(df_2_plot) <- 
  c("Cell_line",
    "Genotype",
    "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))
# 
# df_2_plot$Cell_line <-
#   factor(df_2_plot$Cell_line,
#          levels = c("CD04",
#                     "CD09"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD04", ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD size/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 2),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD04")

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
# df_2_plot <-
#   melt(df_2_plot)
# colnames(df_2_plot) <- 
#   c("Cell_line",
#     "Genotype",
#     "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD09", ]

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
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD size/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 5)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD09")

## Fig_5G ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 4)

df_2_plot <- df_raw
df_2_plot <-
  melt(df_2_plot)
colnames(df_2_plot) <- 
  c("Cell_line",
    "Genotype",
    "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))
# 
# df_2_plot$Cell_line <-
#   factor(df_2_plot$Cell_line,
#          levels = c("CD04",
#                     "CD09"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD04", ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  geom_jitter(aes(colour = Genotype),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
                     paired = F) +
  labs(x = "",
       y = "Fluor. int/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 550),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD04")

df_2_plot <- df_raw
df_2_plot <-
  melt(df_2_plot)
colnames(df_2_plot) <- 
  c("Cell_line",
    "Genotype",
    "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD09", ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  geom_jitter(aes(colour = Genotype),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
                     paired = F) +
  labs(x = "",
       y = "Fluor. int/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1250),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD09")

## Fig_S11B ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 5)

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))
# 
# df_2_plot$Cell_line <-
#   factor(df_2_plot$Cell_line,
#          levels = c("CD04",
#                     "CD09"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD04", ]

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
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD puncta/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 150),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD04")

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD09", ]

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
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD puncta/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 150)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD09")


## Fig_S11C ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 6)

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
colnames(df_2_plot) <- 
  c("Cell_line",
    "Genotype",
    "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))
# 
# df_2_plot$Cell_line <-
#   factor(df_2_plot$Cell_line,
#          levels = c("CD04",
#                     "CD09"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD04", ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD size/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 2.6),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD04")

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
# df_2_plot <-
#   melt(df_2_plot)
# colnames(df_2_plot) <- 
#   c("Cell_line",
#     "Genotype",
#     "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD09", ]

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
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD size/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 2.6)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD09")


## Fig_S11D ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 7)

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))

df_2_plot$Cell_line <-
  factor(df_2_plot$Cell_line,
         levels = c("CD04",
                    "CD09"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD04", ]
# df_2_plot <-
#   df_2_plot[!(df_2_plot$Genotype == "risk"), ]

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
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD area/iMG (um2)") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 200),
                     labels = waiver()) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD04")

df_2_plot <- df_raw

df_2_plot <-
  melt(df_2_plot, 
       variable.name = "Genotype")
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD09", ]

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
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
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
       y = "LD area/iMG (um2)") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 200)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD09")

## Fig_S11E ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 8)

df_2_plot <- df_raw
df_2_plot <-
  melt(df_2_plot)
colnames(df_2_plot) <- 
  c("Cell_line",
    "Genotype",
    "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))
# 
# df_2_plot$Cell_line <-
#   factor(df_2_plot$Cell_line,
#          levels = c("CD04",
#                     "CD09"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD04", ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  geom_jitter(aes(colour = Genotype),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
                     paired = F) +
  labs(x = "",
       y = "Fluor. int/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1000),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD04")

df_2_plot <- df_raw
df_2_plot <-
  melt(df_2_plot)
colnames(df_2_plot) <- 
  c("Cell_line",
    "Genotype",
    "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "CRISPRoff"))

df_2_plot <-
  df_2_plot[df_2_plot$Cell_line == "CD09", ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue",
                                 "green4"),
                      guide = "none") +
  geom_jitter(aes(colour = Genotype),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F,
                     ref.group = "non-risk",
                     paired = F) +
  labs(x = "",
       y = "Fluor. int/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1000),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD09")

## Fig_5I ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 9)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk_TrC",
                    '',
                    "non-risk", "non-risk_TrC"))


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "green4",
                                 "darkblue",
                                 "orange4"),
                      guide = "none") +
  geom_jitter(aes(colour = Genotype),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  stat_compare_means(label = "p.signif",
                     # label.y = 20,
                     # label.y.npc = c(0.5, 0.6),
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("risk", "risk_TrC"),
                                        c("non-risk", "non-risk_TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "LD puncta/iMG") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 15),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

## Fig_5J ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 10)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk_TrC",
                    '',
                    "non-risk", "non-risk_TrC"))


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "green4",
                                 "darkblue",
                                 "orange4"),
                      guide = "none") +
  geom_jitter(aes(colour = Genotype),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  stat_compare_means(label = "p.signif",
                     # label.y = 20,
                     # label.y.npc = c(0.5, 0.6),
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("risk", "risk_TrC"),
                                        c("non-risk", "non-risk_TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "LD size/iMG (um2)") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 35),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

## Fig_5K ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 11)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk_TrC",
                    '',
                    "non-risk", "non-risk_TrC"))


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "green4",
                                 "darkblue",
                                 "orange4"),
                      guide = "none") +
  geom_jitter(aes(colour = Genotype),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  stat_compare_means(label = "p.signif",
                     # label.y = 20,
                     # label.y.npc = c(0.5, 0.6),
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("risk", "risk_TrC"),
                                        c("non-risk", "non-risk_TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "LD size/iMG (um2)") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 10),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

## Fig_5L ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 12)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk_TrC",
                    '',
                    "non-risk", "non-risk_TrC"))


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "green4",
                                 "darkblue",
                                 "orange4"),
                      guide = "none") +
  geom_jitter(aes(colour = Genotype),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  stat_compare_means(label = "p.signif",
                     # label.y = 20,
                     # label.y.npc = c(0.5, 0.6),
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("risk", "risk_TrC"),
                                        c("non-risk", "non-risk_TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "CellRox LD area/iMG (um2)") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 35),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

## Fig_S14B ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 13)


df_2_plot <-
  df_raw
# 
# df_2_plot <-
#   melt(df_2_plot,
#        value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "risk_TrC"))


df_2_plot$variable <-
  factor(df_2_plot$variable,
         levels = c("0 min",
                    "45 min",
                    "90 min",
                    "135 min",
                    "180 min"))

ggplot(df_2_plot,
       aes(x = variable,
           y = Value,
           colour = Genotype,
           fill = Genotype,
           group = Genotype)) +
  # ggerrorplot(df_2_plot,
  #             x = "variable",
  #             y = "Value",
  #             color = "Genotype",
  #             width = 0.2,
  #             # facet.by = "Gene",
  #             # ncol = 4,
  #             error.plot = "errorbar") +
  geom_line(stat = "summary",
            fun = mean,
            linewidth = 1) +
  geom_errorbar(stat = "summary",
                fun.data = mean_cl_boot,
                width = 0.3,
                colour = "black") +
  geom_point(shape = 21,
             stat = "summary",
             fun = mean,
             colour = "black") +
  
  
  # geom_point(shape = 21,
  #            colour = "black") +
  scale_fill_manual(values = c("darkred", 
                               "darkblue",
                               "green4")) +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue",
                                 "green4")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1),
                     na.value = NA) +
  labs(x = "",
       y = "Myelin pHrodo intensity / cell\nNormalized to non-risk (180 min)") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) #+
  ggtitle(unique(df_2_plot$Line))

  ## Fig_S14C ####
df_raw <-
  read_excel("tables_4_plot_v3.xlsx",
             sheet = 14)


df_2_plot <-
  df_raw

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value")

df_2_plot$variable <-
  factor(df_2_plot$variable,
         levels = c("risk",
                    "risk_TrC",
                    # "",
                    "non-risk"))

ggerrorplot(df_2_plot,
            x = "variable",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "green4",
                                 "darkblue",
                                 "orange4"),
                      guide = "none") +
  geom_jitter(aes(colour = variable),
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  stat_compare_means(label = "p.signif",
                     # label.y = 20,
                     label.y.npc = 0.5,
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("risk", "risk_TrC"),
                                        # c("non-risk", "non-risk_TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "Myelin fluor. intensity/iMG") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 2),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 
