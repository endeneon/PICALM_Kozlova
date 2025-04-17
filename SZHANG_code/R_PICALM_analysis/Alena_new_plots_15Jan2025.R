# Siwei 15 Jan 2025
# make new plots for Alena's revised paper

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

  library(lme4)
  library(lmerTest)

}

# qPCR ####
df_all_sources_qPCR <-
  vector(mode = "list",
         length = 4L)

source_xlsx <-
  "Rev_plots/box_whiskers.xlsx"

for (i in 1:4) {
  df_all_sources_qPCR[[i]] <-
    read_excel(path = source_xlsx,
               sheet = i)
}
names(df_all_sources_qPCR) <-
  excel_sheets(source_xlsx)

## panel 1 of 4 (qPCR) ####

df_2_plot <-
  df_all_sources_qPCR[[1]]
t.test(df_2_plot$Average_Exp[df_2_plot$Genotype == "non-risk"],
       df_2_plot$Average_Exp[df_2_plot$Genotype == "risk"])

df_2_plot$Genotype <-
  as.factor(df_2_plot$Genotype)
df_2_plot$clone_new_old <-
  c(rep_len("new", length.out = 6),
    rep_len("old", length.out = 9),
    rep_len("old", length.out = 9),
    rep_len("new", length.out = 6))
# df_2_plot$round <-
#   c(1,1,1,2,2,2,1,1,1,2,2,2,3,3,3,
#     1,1,1,2,2,2,3,3,3,1,1,1,2,2,2)
# df_2_plot$round <-
#   c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,
#     3,3,3,4,4,4,5,5,5,1,1,1,2,2,2)
# df_2_plot$round <-
#   as.factor(df_2_plot$round)

df_2_plot$batch <-
  c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,
    3,3,3,1,1,1,2,2,2,1,1,1,2,2,2)

lm_model <-
  lmerTest::lmer(Average_Exp ~
                   Genotype +
                   (1 | batch/clone_new_old),
                 # (2 + Genotype|diff),
                 data = df_2_plot)
summary(lm_model)

lm_model <-
  lme4::glmer(Average_Exp ~
                   Genotype +
                   (1 | batch/clone_new_old),
                 # (2 + Genotype|diff),
                 data = df_2_plot)
summary(lm_model)

lm_model <-
  lmerTest::lmer(Average_Exp ~
                   Genotype +
                   (1 | clone_new_old) +
                   (1 | batch),
                 # (2 + Genotype|diff),
                 data = df_2_plot)
View(summary(lm_model))
str(summary(lm_model))


lm_model <-
  lmerTest::lmer(Average_Exp ~
                   Genotype +
                   (1 | batch),
                 # (2 + Genotype|diff),
                 data = df_2_plot)

summary(lm_model)
anova(lm_model)

df_2_plot$Clone <-
  as.factor(df_2_plot$Clone)
df_2_plot$diff <-
  as.factor(c(1,1,2,2,3,3,1,1,1,2,2,2,3,3,3,
              1,1,1,2,2,2,3,3,3,1,1,2,2,3,3))

# lm_model <-
#   lmerTest::lmer(Average_Exp ~
#                    Genotype +
#                    (1 + Genotype|Clone/diff),
#                    # (2 + Genotype|diff),
#                  data = df_2_plot)
lm_model <-
  lmerTest::lmer(Average_Exp ~
                   Genotype +
                   (1 | Clone/diff),
                 # (2 + Genotype|diff),
                 data = df_2_plot)


summary(lm_model)
anova(lm_model)

# Random effects:
#   Groups     Name        Variance Std.Dev.
# diff:Clone (Intercept) 0.002398 0.04897
# Clone      (Intercept) 0.000713 0.02670
# Residual               0.003250 0.05701
# Number of obs: 30, groups:  diff:Clone, 12; Clone, 4
#
# Fixed effects:
#   Estimate Std. Error       df t value Pr(>|t|)
# (Intercept)   0.40762    0.03130  1.98316  13.023  0.00604 **
#   Genotyperisk -0.19756    0.04427  1.98316  -4.463  0.04744 *

lm_model <-
  lmerTest::lmer(Average_Exp ~
                   Genotype +
                   (1 | diff),
                 # (2 + Genotype|diff),
                 data = df_2_plot)
summary(lm_model)
anova(lm_model)

# Random effects:
#   Groups   Name        Variance  Std.Dev.
# diff     (Intercept) 7.463e-05 0.008639
# Residual             5.800e-03 0.076158
# Number of obs: 30, groups:  diff, 3
#
# Fixed effects:
#   Estimate Std. Error       df t value Pr(>|t|)
# (Intercept)   0.41394    0.02029  6.70883  20.404 2.74e-07 ***
#   Genotyperisk -0.20129    0.02781 26.00000  -7.238 1.09e-07 ***

df_2_plot$Clone[df_2_plot$Clone %in% c("B4", "A5")] <- 1
df_2_plot$Clone[df_2_plot$Clone %in% c("G3", "A6")] <- 2

ggplot(df_2_plot,
       aes(x = factor(Genotype,
                      levels = c("risk",
                                 "non-risk")),
           y = Average_Exp)) +
  # geom_boxplot(outliers = T,
  #              ) +
  # geom_point(aes(shape = Clone,
  #                colour = Genotype)) +
  geom_jitter(aes(shape = Clone,
                  colour = Genotype,
                  fill = Genotype),
              width = 0.1) +
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.4,
               size = 0.5) +
  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               geom = "errorbar",
               width = 0.2) +
  scale_shape_manual(values = c(1, 2),
                     labels = c("Clone 1",
                                "Clone 2")) +
  scale_colour_manual(aes(factor(Genotype,
                                 levels = c("risk",
                                            "non-risk"))),
                      values = c("darkblue",
                                 "darkred"),
                      breaks = c("risk",
                                 "non-risk")) +

  scale_fill_manual(aes(factor(Genotype,
                               levels = c("risk",
                                          "non-risk"))),
                    values = c("darkblue",
                               "darkred"),
                    breaks = c("risk",
                               "non-risk")) +
  guides(colour = guide_legend(title = "Genotype"),
         fill = guide_legend(title = "Genotype"),
         shape = guide_legend(title = "Clone")) +
  ylim(0, 0.7) +
  xlab("") +
  ylab("PICALM expression") +
  # geom_linerange()
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12,
                                    face = "bold"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1),
        axis.title = element_text(size = 12,
                                  face = "bold")) +
  ggtitle(label = unique(df_2_plot$line)[1],
          subtitle = paste0("Wilcox's test; P = ",
                            wilcox.test(x = df_2_plot$Average_Exp[df_2_plot$Genotype == unique(df_2_plot$Genotype)[1]],
                                        y = df_2_plot$Average_Exp[df_2_plot$Genotype == unique(df_2_plot$Genotype)[2]],
                                        alternative = "two.sided", paired = F)$p.value))


## panel 2 of 4 (qPCR) ####

df_2_plot <-
  df_all_sources_qPCR[[2]]

df_2_plot$Genotype <-
  as.factor(df_2_plot$Genotype)
df_2_plot$Clone <-
  as.factor(df_2_plot$Clone)

df_2_plot$batch <-
  c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,
    1,1,1,2,2,2,3,3,3,4,4,4,5,5,5)
# df_2_plot$diff <-
#   as.factor(c(1,1,1,2,2,2,3,3,3,1,1,2,2,3,3,
#               1,1,1,2,2,2,3,3,3,1,1,2,2,3,3))
# df_2_plot$diff <-
#   as.factor(c(1,1,1,2,2,2,3,3,3,1,1,2,2,3,3,
#               4,4,4,5,5,5,6,6,6,4,4,5,5,6,6))

lm_model <-
  lmerTest::lmer(Average_Exp ~
                   Genotype +
                   (1 | Clone/batch),
                 # (2 + Genotype|diff),
                 data = df_2_plot)

summary(lm_model)
anova(lm_model)


lm_model <-
  lmerTest::lmer(Average_Exp ~
                   Genotype +
                   (1 | batch),
                 # (2 + Genotype|diff),
                 data = df_2_plot)



summary(lm_model)
anova(lm_model)




df_2_plot$Clone[df_2_plot$Clone %in% c("A1", "A3")] <- 1
df_2_plot$Clone[df_2_plot$Clone %in% c("F6", "A9")] <- 2

ggplot(df_2_plot,
       aes(x = factor(Genotype,
                      levels = c("risk",
                                 "non-risk")),
           y = Average_Exp)) +
  # geom_boxplot(outliers = T,
  #              ) +
  # geom_point(aes(shape = Clone,
  #                colour = Genotype)) +
  geom_jitter(aes(shape = Clone,
                  colour = Genotype,
                  fill = Genotype),
              width = 0.1) +
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.4,
               size = 0.5) +
  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               geom = "errorbar",
               width = 0.2) +
  scale_shape_manual(values = c(1, 2),
                     labels = c("Clone 1",
                                "Clone 2")) +
  scale_colour_manual(aes(factor(Genotype,
                                 levels = c("risk",
                                            "non-risk"))),
                      values = c("darkblue",
                                 "darkred"),
                      breaks = c("risk",
                                 "non-risk")) +

  scale_fill_manual(aes(factor(Genotype,
                               levels = c("risk",
                                          "non-risk"))),
                    values = c("darkblue",
                               "darkred"),
                    breaks = c("risk",
                               "non-risk")) +
  guides(colour = guide_legend(title = "Genotype"),
         fill = guide_legend(title = "Genotype"),
         shape = guide_legend(title = "Clone")) +

  xlab("") +
  ylab("PICALM expression") +
  ylim(0, 0.7) +
  # geom_linerange()
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12,
                                    face = "bold"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1),
        axis.title = element_text(size = 12,
                                  face = "bold")) +
  ggtitle(label = unique(df_2_plot$line)[1],
          subtitle = paste0("Wilcox's test; P = ",
                            wilcox.test(x = df_2_plot$Average_Exp[df_2_plot$Genotype == unique(df_2_plot$Genotype)[1]],
                                        y = df_2_plot$Average_Exp[df_2_plot$Genotype == unique(df_2_plot$Genotype)[2]],
                                        alternative = "two.sided", paired = F)$p.value))

## panel 3 of 4 (CRISPR) ####

df_2_plot <-
  df_all_sources_qPCR[[3]]

df_2_plot$Genotype[df_2_plot$Genotype == "risk_CRISPRa"] <-
  "risk-CRISPRa"

df_2_plot$diff <-
  factor(c(1,1,1,2,2,2,
           1,1,1,2,2,2,
           1,1,1,2,2,2))

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "risk-CRISPRa"))
df_2_plot$Clone <-
  as.factor(df_2_plot$Clone)


# lm_model <-
#   lmerTest::lmer(Average_Exp ~
#                    Genotype +
#                    (1 | Clone/diff),
#                  # (2 + Genotype|diff),
#                  data = df_2_plot)
lm_model <-
  lmerTest::lmer(Average_Exp ~
                   Genotype +
                   (1 | diff),
                 # (2 + Genotype|diff),
                 data = df_2_plot)



summary(lm_model)
anova(lm_model)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Average_Exp ~ Genotype + (1 | diff)
#    Data: df_2_plot
#
# REML criterion at convergence: -26.1
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -1.9753 -0.4579  0.2023  0.6454  1.7242
#
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  diff     (Intercept) 0.004846 0.06961
#  Residual             0.006261 0.07912
# Number of obs: 18, groups:  diff, 2
#
# Fixed effects:
#                      Estimate Std. Error       df t value Pr(>|t|)
# (Intercept)           0.21169    0.05888  1.55811   3.596  0.09903 .
# Genotypenon-risk      0.23329    0.04568 14.00000   5.107  0.00016 ***
# Genotyperisk-CRISPRa  0.16313    0.04568 14.00000   3.571  0.00307 **

output_panel <-
  ggplot(df_2_plot,
       aes(x = Genotype,
           y = Average_Exp)) +
  # geom_boxplot(outliers = T,
  #              ) +
  # geom_point(aes(shape = Clone,
  #                colour = Genotype)) +
  geom_jitter(aes(colour = Genotype,
                  fill = Genotype),
              shape = 1,
              width = 0.1) +
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.4,
               size = 0.5) +
  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               geom = "errorbar",
               width = 0.2) +
  # scale_shape_manual(values = c(1, 2),
  #                    labels = c("Clone 1",
  #                               "Clone 2")) +
  scale_colour_manual(aes(Genotype),
                      values = c("darkblue",
                                 "darkred",
                                 "darkgreen")) +

  scale_fill_manual(aes(Genotype),
                    values = c("darkblue",
                               "darkred",
                               "darkgreen")) +
  guides(colour = guide_legend(title = "Genotype"),
         fill = F) +
  # ylim(0, 0.6) +
  scale_y_continuous(limits = c(0, 0.7),
                     expand = c(0, 0)) +
  xlab("") +
  ylab("PICALM expression") +
  # geom_linerange()
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12,
                                    face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1),
        axis.title = element_text(size = 12,
                                  face = "bold"))

kruskal.test(x = df_2_plot$Average_Exp,
             g = df_2_plot$Genotype)


if (anova_test(formula = Average_Exp ~ Genotype,
               data = df_2_plot,
               type = 3)$p > 0.05) {
  sub_text <- ""
  anova_p_signif = ", not significant"
} else {
  sub_text <-
    paste0("Dunnett \'s multiple comparison test adjusted P values\n",
           rownames(DunnettTest(x = df_2_plot$Average_Exp,
                                g = df_2_plot$Genotype,
                                control = "non-risk")[[1]])[1],
           ": ",
           DunnettTest(x = df_2_plot$Average_Exp,
                       g = df_2_plot$Genotype,
                       control = "non-risk")[[1]][1, 4],
           "\n",
           rownames(DunnettTest(x = df_2_plot$Average_Exp,
                                g = df_2_plot$Genotype,
                                control = "non-risk")[[1]])[2],
           ": ",
           DunnettTest(x = df_2_plot$Average_Exp,
                       g = df_2_plot$Genotype,
                       control = "non-risk")[[1]][2, 4])
  anova_p_signif = ", significant"

}

output_panel <-
  output_panel +
  ggtitle(label = paste0(unique(df_2_plot$line)[1],
                         "; Anova \'s P =",
                         anova_test(formula = Average_Exp ~ Genotype,
                                    data = df_2_plot,
                                    type = 3)$p,
                         anova_p_signif),
          subtitle = sub_text)

print(output_panel)

#
# df_2_plot$Genotype <-
#   relevel(df_2_plot$Genotype,
#           ref = "risk")
#
#
# summary(aov(Average_Exp ~ Genotype,
#             data = df_2_plot))
# summary(multcomp::glht(aov(Average_Exp ~ Genotype,
#                            data = df_2_plot),
#                        linfct = mcp(Genotype = "Dunnett")))



## panel 4 of 4 (CRISPR) ####

df_2_plot <-
  df_all_sources_qPCR[[4]]
#
df_2_plot$Genotype[df_2_plot$Genotype == "risk_CRISPRa"] <-
  "risk-CRISPRa"

# df_2_plot$Genotype <-
#   factor(df_2_plot$Genotype,
#          levels = c("non-risk",
#                     "risk",
#                     "risk-CRISPRa"))


df_2_plot$diff <-
  factor(c(1,1,1,2,2,2,
           1,1,1,2,2,2,
           1,1,1,2,2,2))

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "risk-CRISPRa"))
df_2_plot$Clone <-
  as.factor(df_2_plot$Clone)


lm_model <-
  lmerTest::lmer(Average_Exp ~
                   Genotype +
                   (1 | Clone/diff),
                 # (2 + Genotype|diff),
                 data = df_2_plot)
lm_model <-
  lmerTest::lmer(Average_Exp ~
                   Genotype +
                   (1 | diff),
                 # (2 + Genotype|diff),
                 data = df_2_plot)



summary(lm_model)
anova(lm_model)



output_panel <-
  ggplot(df_2_plot,
         aes(x = Genotype,
             y = Average_Exp)) +
  # geom_boxplot(outliers = T,
  #              ) +
  # geom_point(aes(shape = Clone,
  #                colour = Genotype)) +
  geom_jitter(aes(colour = Genotype,
                  fill = Genotype),
              shape = 1,
              width = 0.1) +
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.4,
               size = 0.5) +
  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               geom = "errorbar",
               width = 0.2) +
  # scale_shape_manual(values = c(1, 2),
  #                    labels = c("Clone 1",
  #                               "Clone 2")) +
  scale_colour_manual(aes(colour = Genotype),
                      values = c("darkblue",
                                 "darkred",
                                 "darkgreen")) +

  # scale_fill_manual(aes(factor(Genotype,
  #                              levels = c("risk",
  #                                         "non-risk",
  #                                         "risk-CRISPRa"))),
  #                   values = c("darkblue",
  #                              "darkred",
  #                              "darkgreen"),
  #                   breaks = c("risk",
  #                              "non-risk",
  #                              "risk-CRISPRa")) +
  guides(colour = guide_legend(title = "Genotype"),
         fill = F) +
  # ylim(0, 0.7) +
  xlab("") +
  ylab("PICALM expression") +
  scale_y_continuous(limits = c(0, 0.7),
                     expand = c(0, 0)) +
  # geom_linerange()
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12,
                                    face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1),
        axis.title = element_text(size = 12,
                                  face = "bold"))

if (anova_test(formula = Average_Exp ~ Genotype,
               data = df_2_plot,
               type = 3)$p > 0.05) {
  sub_text <- ""
  anova_p_signif = ", not significant"
} else {
  sub_text <-
    paste0("Dunnett \'s multiple comparison test adjusted P values\n",
           rownames(DunnettTest(x = df_2_plot$Average_Exp,
                                g = df_2_plot$Genotype,
                                control = "non-risk")[[1]])[1],
           ": ",
           DunnettTest(x = df_2_plot$Average_Exp,
                       g = df_2_plot$Genotype,
                       control = "non-risk")[[1]][1, 4],
           "\n",
           rownames(DunnettTest(x = df_2_plot$Average_Exp,
                                g = df_2_plot$Genotype,
                                control = "non-risk")[[1]])[2],
           ": ",
           DunnettTest(x = df_2_plot$Average_Exp,
                       g = df_2_plot$Genotype,
                       control = "non-risk")[[1]][2, 4])
  anova_p_signif = ", significant"

}

output_panel <-
  output_panel +
  ggtitle(label = paste0(unique(df_2_plot$line)[1],
                         "; Anova \'s P =",
                         anova_test(formula = Average_Exp ~ Genotype,
                                    data = df_2_plot,
                                    type = 3)$p,
                         anova_p_signif),
          subtitle = sub_text)

print(output_panel)


df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "risk")

#
# summary(aov(Average_Exp ~ Genotype,
#             data = df_2_plot))
# summary(multcomp::glht(aov(Average_Exp ~ Genotype,
#                            data = df_2_plot),
#                        linfct = mcp(Genotype = "Dunnett")))




# plot 3x3 cell types ####

df_raw <-
  read_excel("Rev_plots/Sun_Cell_2023_MG_PICALM_cpm_indiv_w_metadata_28Jun2024_JD2025_plot.xlsx")

df_2_plot <-
  df_raw

df_2_plot <-
  data.table::melt(df_2_plot,
                   na.rm = T)

dev.off()

# i <- 1

for (i in unique(df_2_plot$variable)) {
  print(i)

  df_2_plot_subsetted <-
    df_2_plot[df_2_plot$variable == i, ]
#
  df_2_plot_subsetted <-
    df_2_plot[df_2_plot$variable == unique(df_2_plot$variable)[i], ]

  pdf(file = paste0("Rev_plots/v2_",
                    i,
                    "_Anova_Dunnett.pdf"),
      width = 2,
      height = 4)



# i <-
#   unique(df_2_plot$variable)[1]
# df_2_plot_subsetted <-
#   df_2_plot[df_2_plot$variable == i, ]
#
# aov_results <-
#   anova_test(formula = value ~ ADdiag3types,
#              data = df_2_plot_subsetted,
#              type = 3)
#
# Dunnett_p <-
#   DunnettTest(x = df_2_plot_subsetted$value,
#               g = df_2_plot_subsetted$ADdiag3types,
#               control = "nonAD")
#
# # aov_results <-
#   anova_summary(aov(formula = value ~ ADdiag3types,
#                     data = df_2_plot_subsetted))
# # summary(aov_results)

# pwc <-
  # aov_results %>%
  # print(LSD.test(aov_results,
  #                "ADdiag3types")) #%>%
  # add_xy_position(x = "ADdiag3types")

output_panel <-
  ggplot(df_2_plot_subsetted,
         aes(x = ADdiag3types,
             y = value)) +
  # geom_boxplot(outliers = T,
  #              ) +
  # geom_point(aes(shape = Clone,
  #                colour = Genotype)) +
  geom_jitter(aes(colour = ADdiag3types,
                  fill = ADdiag3types),
              shape = 1,
              width = 0.1) +
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.4,
               size = 0.5) +
  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               geom = "errorbar",
               width = 0.2) +
  # scale_shape_manual(values = c(1, 2),
  #                    labels = c("Clone 1",
  #                               "Clone 2")) +
  scale_colour_manual(aes(colour = factor(ADdiag3types)),
                      values = c("darkblue",
                                 "darkred",
                                 "darkgreen")) +

  # scale_fill_manual(aes(factor(Genotype,
  #                              levels = c("risk",
  #                                         "non-risk",
  #                                         "risk-CRISPRa"))),
  #                   values = c("darkblue",
  #                              "darkred",
  #                              "darkgreen"),
  #                   breaks = c("risk",
  #                              "non-risk",
  #                              "risk-CRISPRa")) +
  guides(colour = guide_legend(title = "ADdiag3types"),
         fill = F) +
  # ylim(0, 0.7) +
  xlab("") +
  ylab("PICALM expression in log2CPM") +
  scale_y_continuous(expand = c(0, 0.05)) +
    # stat_compare_means(method = "anova") +
  # geom_linerange()
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12,
                                    face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1),
        axis.title = element_text(size = 12,
                                  face = "bold"))

  if (anova_test(formula = value ~ ADdiag3types,
                 data = df_2_plot_subsetted,
                 type = 3)$p > 0.05) {
    sub_text <- ""
    anova_p_signif = ", not significant"
  } else {
    sub_text <-
      paste0("Dunnett \'s multiple comparison test adjusted P values\n",
             rownames(DunnettTest(x = df_2_plot_subsetted$value,
                                  g = df_2_plot_subsetted$ADdiag3types,
                                  control = "nonAD")[[1]])[1],
             ": ",
             DunnettTest(x = df_2_plot_subsetted$value,
                         g = df_2_plot_subsetted$ADdiag3types,
                         control = "nonAD")[[1]][1, 4],
             "\n",
             rownames(DunnettTest(x = df_2_plot_subsetted$value,
                                  g = df_2_plot_subsetted$ADdiag3types,
                                  control = "nonAD")[[1]])[2],
             ": ",
             DunnettTest(x = df_2_plot_subsetted$value,
                         g = df_2_plot_subsetted$ADdiag3types,
                         control = "nonAD")[[1]][2, 4])
    anova_p_signif = ", significant"

  }

  output_panel <-
    output_panel +
    ggtitle(label = paste0(i,
                           "; Anova \'s P =",
                          anova_test(formula = value ~ ADdiag3types,
                                     data = df_2_plot_subsetted,
                                     type = 3)$p,
                          anova_p_signif),
            subtitle = sub_text)

  print(output_panel)

  dev.off()

}



output_panel +
  # stat_pvalue_manual(Dunnett_p,
  #                    hide.ns = T) +
  # stat_anova_test(method = "one_way",
  #                 group.by = "x.var",
  #                 type = 2,
  #                 correction = "none",
  #                 p.adjust.method = "fdr") +
  ggtitle(i)


test_var <-
  DunnettTest(x = df_2_plot_subsetted$value,
              g = df_2_plot_subsetted$ADdiag3types,
              control = "nonAD")[[1]]


