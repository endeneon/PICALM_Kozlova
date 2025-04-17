# Siwei 11 Mar 2025
# plot for Nature appeal letter

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
}

rm(list = ls())

set.seed(42)

# PY2 #####
df_raw <-
  read_excel("plot_PY2_TMEM119.xlsx",
             sheet = 1)

df_2_plot <-
  df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

# lm_model <-
lm_model <-
  lmerTest::lmer(Per_cover_slip ~
                   Genotype +
                   (1 | Cell_line) +
                   (1 | Clone),
                 # (2 + Genotype|diff),
                 data = df_2_plot)

summary(lm_model)

ggplot(df_2_plot,
       aes(x = Genotype,
           y = Per_cover_slip,
           colour = Genotype,
           shape = Cell_line)) +
  geom_jitter(width = 0.1,
              height = 0) +
  scale_shape_manual(values = c(1, 2)) +
  scale_colour_manual(values = c("darkred",
                                 "darkblue")) +
  stat_summary(aes(group = Genotype),
               geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  geom_errorbar(stat = "summary",
                aes(group = Genotype),
                # group = "Genotypeaes(group = Genotype),",
                fun.data = mean_cl_boot,
                width = 0.1,
                colour = "black",
                size = 0.5) +
  guides(colour = F) +
  # stat_summary(mappaes(y = all_pos_cells / nuclei_count),
  #              geom = "segment",
  #              fun = "mean",
  #              mapping = aes(xend = after_stat(x) - 0.1,
  #                            yend = after_stat(y)),
  #              colour = "black") +
  # stat_summary(geom = "segment",
  #              fun = "mean",
  #              mapping = aes(xend = after_stat(x) + 0.1,
  #                            yend = after_stat(y)),
  #              colour = "black") +
  scale_y_continuous(limits = c(0, 1.2),
                     labels = scales::percent,
                     expand = c(0, 0),

                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                     # breaks = c(0, 0.5, 0.9, 1)) +
  # geom_hline(yintercept = 0.9) +
  ylab("TREM+, CD45+, PY2+ cells per coverslip") +
  xlab("") +

  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"))



#
# TMEM119 #####
df_raw <-
  read_excel("plot_PY2_TMEM119.xlsx",
             sheet = 2)

df_2_plot <-
  df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

ggplot(df_2_plot,
       aes(x = Genotype,
           y = Per_cover_slip,
           colour = Genotype,
           shape = Cell_line)) +
  geom_jitter(width = 0.1,
              height = 0) +
  scale_shape_manual(values = c(1, 2)) +
  scale_colour_manual(values = c("darkred",
                                 "darkblue")) +
  stat_summary(aes(group = Genotype),
               geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  geom_errorbar(aes(group = Genotype),
                stat = "summary",
                fun.data = mean_cl_boot,
                width = 0.1,
                colour = "black",
                size = 0.5) +
  guides(colour = F) +
  # stat_summary(mappaes(y = all_pos_cells / nuclei_count),
  #              geom = "segment",
  #              fun = "mean",
  #              mapping = aes(xend = after_stat(x) - 0.1,
  #                            yend = after_stat(y)),
  #              colour = "black") +
  # stat_summary(geom = "segment",
  #              fun = "mean",
  #              mapping = aes(xend = after_stat(x) + 0.1,
  #                            yend = after_stat(y)),
  #              colour = "black") +
  scale_y_continuous(limits = c(0, 1.2),
                     labels = scales::percent,
                     expand = c(0, 0),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  # geom_hline(yintercept = 0.8) +
  ylab("TREM+, CD45+, TMEM119+ cells per coverslip") +
  xlab("") +

  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"))

t.test(x = (df_2_plot$all_pos_cells / df_2_plot$nuclei_count)[df_2_plot$Genotype == "risk"],
       y = (df_2_plot$all_pos_cells / df_2_plot$nuclei_count)[df_2_plot$Genotype == "non-risk"],
       paired = F, var.equal = F)

# PICALM_AD #####
rm(list = ls())
df_raw <-
  read_excel("PICALM_reduced_in_AD_proteomics.xlsx",
             sheet = 2)
df_2_plot <-
  df_raw

df_2_plot$Phenotype[df_2_plot$Phenotype == "Control"] <- "nonAD"
df_2_plot$Phenotype <-
  factor(df_2_plot$Phenotype,
         levels = c("nonAD",
                    "AsymAD",
                    "AD"))


ggplot(df_2_plot,
       aes(x = Phenotype,
           y = Value,
           colour = Phenotype)) +
  geom_point(aes(colour = Phenotype),
             shape = 20,
             position = position_dodge2(width = 0.1)) +
  # geom_jitter(width = 0.1,
  #             size = 0.1) +
  # scale_shape_manual(values = c(1, 2)) +
  # scale_colour_manual(values = c("darkred",
  #                                "darkblue")) +
  stat_summary(aes(group = Phenotype),
               geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black",
               position = position_dodge2(width = 0.2)) +
  geom_errorbar(stat = "summary",
                fun.data = mean_cl_boot,
                width = 0.1,
                colour = "black",
                size = 0.5,
                position = position_dodge2(width = 0.5)) +

  # geom_
  # stat_summary(mappaes(y = all_pos_cells / nuclei_count),
  #              geom = "segment",
  #              fun = "mean",
  #              mapping = aes(xend = after_stat(x) - 0.1,
  #                            yend = after_stat(y)),
  #              colour = "black") +
  # stat_summary(geom = "segment",
  #              fun = "mean",
  #              mapping = aes(xend = after_stat(x) + 0.1,
  #                            yend = after_stat(y)),
  #              colour = "black") +
  scale_y_continuous(limits = c(-0.4, 0.4),
                     # labels = scales::percent,
                     expand = c(0, 0)) +
  scale_colour_manual(values = brewer.pal(n = 3,
                                          name = "Dark2")) +
  guides(colour = F) +

  ylab("Regressed clean log2(-centered abundance)") +
  xlab("") +

  theme_classic() +
  theme(axis.text = element_text(size = 10)) +
  ggtitle("Anova's P value = 4.32e-6")


summary(aov(glm(Value ~ Phenotype,
                data = df_2_plot)))
anova(lm(Value ~ Phenotype,
         data = df_2_plot))
DunnettTest(x = df_2_plot$Value,
            g = df_2_plot$Phenotype)
summary(anova(lm(Value ~ Phenotype,
                  data = df_2_plot)))
