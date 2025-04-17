# Siwei 26 Jan 2025
# plot Fig. 2f

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
}

set.seed(42)


df_raw <-
  read_excel(path = "Batch_2_of_new_plots/Fig_2F.xlsx")

## CD04 ####
df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]
colnames(df_2_plot)[2] <-
  "Genotype"

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))
df_2_plot$batch <-
  c(1,2,3,4,5,6,7,
    1,2,3,4,5,6,7)

lm_model <-
  lmerTest::lmer(PICALM.actin ~
                   Genotype +
                   (1 | batch),
                 # (2 + Genotype|diff),
                 data = df_2_plot)



summary(lm_model)
anova(lm_model)


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "PICALM.actin",
            color = "Genotype",
            shape = "Clone",
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
       y = "PICALM/b-actin") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        # legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD04")


## CD09 ####
df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]
colnames(df_2_plot)[2] <-
  "Genotype"

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))

df_2_plot$batch <-
  c(1,2,3,4,5,6,7,
    1,2,3,4,5,6,7)

lm_model <-
  lmerTest::lmer(PICALM.actin ~
                   Genotype +
                   (1 | batch),
                 # (2 + Genotype|diff),
                 data = df_2_plot)



summary(lm_model)
anova(lm_model)

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "PICALM.actin",
            color = "Genotype",
            shape = "Clone",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = 0.8,
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
       y = "PICALM/b-actin") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        # legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD09")
