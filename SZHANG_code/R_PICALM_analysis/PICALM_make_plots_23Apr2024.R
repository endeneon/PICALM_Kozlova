# make figures for Alena PICALM paper, additional data
# Siwei 23 Apr 2024

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

## Fig_5I ####
df_raw <-
  read_excel("tables_4_plot_v4.xlsx",
             sheet = 1)

df_2_plot <- df_raw

df_2_plot <-
  reshape2::melt(df_2_plot,
                 value.name = "Value")
colnames(df_2_plot) <-
  c("Genotype", "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk_TrC",
                    # '',
                    "non-risk", "non-risk_TrC"))
df_2_plot$BR <-
  rep_len(x = c(1,1,1,1,
                2,2,2,2),
          length.out = 32)

lm_model <-
  lmerTest::lmer(Value ~ 
                   Genotype + 
                   (1|BR),
                 data = df_2_plot)
summary(lm_model)





ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkblue",
                                 "green4",
                                 "darkred",
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
                     # method = 't.test',
                     method = "wilcox.test",
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
  # scale_y_continuous(expand = c(0, 0),
  #                    limits = c(0, 120),
  #                    na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

## Fig_5J ####
df_raw <-
  read_excel("tables_4_plot_v4.xlsx",
             sheet = 2)

df_2_plot <- df_raw

df_2_plot <-
  reshape2::melt(df_2_plot,
                 value.name = "Value")
colnames(df_2_plot) <-
  c("Genotype", "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk_TrC",
                    # '',
                    "non-risk", "non-risk_TrC"))

df_2_plot$BR <-
  rep_len(x = c(1,1,1,1,
                2,2,2,2),
          length.out = 32)

lm_model <-
  lmerTest::lmer(Value ~ 
                   Genotype + 
                   (1|BR),
                 data = df_2_plot)
summary(lm_model)
anova(lm_model)

lm_model <-
  lmerTest::lmer(Value ~ 
                   Genotype,
                 data = df_2_plot)
summary(lm_model)
anova(lm_model)

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk", "non-risk_TrC",
                    "risk", "risk_TrC"))

lm_model <-
  lmerTest::lmer(Value ~ 
                   Genotype + 
                   (1|BR),
                 data = df_2_plot)
summary(lm_model)
anova(lm_model)

lm_model <-
  glm(Value ~ 
                   Genotype,
                 data = df_2_plot)
summary(lm_model)
anova(lm_model)


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkblue",
                                 "green4",
                                 "darkred",
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
       y = "LD area/iMG") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 120),
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
  read_excel("tables_4_plot_v4.xlsx",
             sheet = 3)

df_2_plot <- df_raw

df_2_plot <-
  reshape2::melt(df_2_plot,
                 value.name = "Value")
colnames(df_2_plot) <-
  c("Genotype", "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk_TrC",
                    # '',
                    "non-risk", "non-risk_TrC"))


df_2_plot$BR <-
  rep_len(x = c(1,1,1,1,
                2,2,2,2),
          length.out = 32)

lm_model <-
  lmerTest::lmer(Value ~ 
                   Genotype + 
                   (1|BR),
                 data = df_2_plot)

summary(lm_model)
anova(lm_model)

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk", "non-risk_TrC",
                    "risk", "risk_TrC"))

lm_model <-
  lmerTest::lmer(Value ~ 
                   Genotype + 
                   (1|BR),
                 data = df_2_plot)
summary(lm_model)
anova(lm_model)

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkblue",
                                 "green4",
                                 "darkred",
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
       y = "CellRox LD puncta/iMG") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100),
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
  read_excel("tables_4_plot_v4.xlsx",
             sheet = 4)

df_2_plot <- df_raw

df_2_plot <-
  reshape2::melt(df_2_plot,
                 value.name = "Value")
colnames(df_2_plot) <-
  c("Genotype", "Value")

# df_2_plot$Genotype <-
#   factor(df_2_plot$Genotype,
#          levels = c("risk", "risk_TrC",
#                     '',
#                     "non-risk", "non-risk_TrC"))

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk_TrC",
                    # '',
                    "non-risk", "non-risk_TrC"))


df_2_plot$BR <-
  rep_len(x = c(1,1,1,1,
                2,2,2,2),
          length.out = 32)

lm_model <-
  lmerTest::lmer(Value ~ 
                   Genotype + 
                   (1|BR),
                 data = df_2_plot)

summary(lm_model)
anova(lm_model)

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk", "non-risk_TrC",
                    "risk", "risk_TrC"))

lm_model <-
  lmerTest::lmer(Value ~ 
                   Genotype + 
                   (1|BR),
                 data = df_2_plot)
summary(lm_model)
anova(lm_model)

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkblue",
                                 "green4",
                                 "darkred",
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
       y = "CellRox LD area/iMG") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 
