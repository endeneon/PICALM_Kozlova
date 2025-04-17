# Siwei 04 Feb 2025
# plot new Fig Ex 7bc

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

df_raw <-
  read_excel("Batch_2_of_new_plots/Fig_Ex_7bc.xlsx",
             sheet = 1)

df_2_plot <-
  df_raw
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))

df_2_plot$BR_ind <-
  str_split(string = df_2_plot$BR,
            pattern = "_",
            simplify = T)[, 1]

# CD04 ####
df_2_plot_df <-
  df_2_plot[df_2_plot$Cell.Line == "CD04", ]
# df_2_plot_df$uuid <-
#   str_c(df_2_plot_df$)

df_2_plot_df$Genotype <-
  factor(df_2_plot_df$Genotype,
         levels = c("risk",
                    "non-risk"))
df_2_plot_df$Clone <-
  tolower(df_2_plot_df$Clone)
unique(df_2_plot_df$Clone)

lm_model <-
  lmerTest::lmer(LD.area.iMG ~
                   Genotype +
                   (1|Clone/BR_ind),
                 data = df_2_plot_df)
summary(lm_model)

lm_model <-
  lmerTest::lmer(Fl.iMG ~
                   Genotype +
                   (1|Clone/BR_ind),
                 data = df_2_plot_df)
summary(lm_model)

## PLIN2 LD area ####
ggerrorplot(df_2_plot_df,
            x = "Genotype",
            y = "LD.area.iMG",
            color = "Genotype",
            # shape = "Clone",
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
  geom_jitter(aes(colour = Genotype,
                  shape = Clone),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "LD area/iMG") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_2_plot$LD.area.iMG * 1.1))) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0)) +
  ggtitle("CD04 PLIN2 LD area")

## PLIN2 Fluor intensity area ####
ggerrorplot(df_2_plot_df,
            x = "Genotype",
            y = "Fl.iMG",
            color = "Genotype",
            # shape = "Clone",
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
  geom_jitter(aes(colour = Genotype,
                  shape = Clone),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "LD area/iMG") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_2_plot$Fl.iMG * 1.1))) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0)) +
  ggtitle("CD04 PLIN2 Fluor ")


# CD09 ####
df_2_plot_df <-
  df_2_plot[df_2_plot$Cell.Line == "CD09", ]
# df_2_plot_df$uuid <-
#   str_c(df_2_plot_df$)

df_2_plot_df$Genotype <-
  factor(df_2_plot_df$Genotype,
         levels = c("risk",
                    "non-risk"))
df_2_plot_df$Clone <-
  tolower(df_2_plot_df$Clone)
unique(df_2_plot_df$Clone)

lm_model <-
  lmerTest::lmer(LD.area.iMG ~
                   Genotype +
                   (1|Clone/BR_ind),
                 data = df_2_plot_df)
summary(lm_model)

lm_model <-
  lmerTest::lmer(Fl.iMG ~
                   Genotype +
                   (1|Clone/BR_ind),
                 data = df_2_plot_df)
summary(lm_model)

## PLIN2 LD area ####
ggerrorplot(df_2_plot_df,
            x = "Genotype",
            y = "LD.area.iMG",
            color = "Genotype",
            # shape = "Clone",
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
  geom_jitter(aes(colour = Genotype,
                  shape = Clone),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "LD area/iMG") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_2_plot$LD.area.iMG * 1.1))) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0)) +
  ggtitle("CD09 PLIN2 LD area")

## PLIN2 Fluor intensity area ####
ggerrorplot(df_2_plot_df,
            x = "Genotype",
            y = "Fl.iMG",
            color = "Genotype",
            # shape = "Clone",
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
  geom_jitter(aes(colour = Genotype,
                  shape = Clone),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Fluor. intensity") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_2_plot$Fl.iMG * 1.1))) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0)) +
  ggtitle("CD09 PLIN2 Fluor ")


## shape_guide ####
ggerrorplot(df_2_plot_df,
            x = "Genotype",
            y = "Fl.iMG",
            color = "Genotype",
            # shape = "Clone",
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
  geom_jitter(aes(colour = Genotype,
                  shape = Clone),
              width = 0.1,
              size = 1) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Fluor. intensity") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_2_plot$Fl.iMG * 1.1))) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        # legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0)) +
  ggtitle("CD09 PLIN2 Fluor ")
