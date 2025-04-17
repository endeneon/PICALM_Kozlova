# make figures for Alena PICALM paper
# Siwei 19 Mar 2024

# init ####
{
  library(readxl)
  
  library(stringr)
  library(ggplot2)
  
  library(scales)
  library(reshape2)
  
  library(RColorBrewer)
  library(ggpubr)
  
  library(lme4)
  library(lmerTest)
}
# library(seqLogo)

## Fig_1e ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 1)

df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]

ggplot(df_2_plot,
       aes(x = Genotype,
           y = Avg_exp,
           group = Genotype,
           colour = Genotype,
           fill = Genotype)) +
  stat_summary(geom = "crossbar",
               fun = "median",
               width = 0.5,
               linewidth = 0.5,
               colour = "black") +
  geom_jitter(width = 0.1,
              size = 2,
              shape = 1) +
  scale_fill_manual(values = c("darkred", 
                               "darkblue"),
                    guide = "none") +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue"),
                      guide = "none") +
  # ylim(0, 100) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.6)) +
  labs(x = "",
       y = "Normalized exp. value") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black")) +
  ggtitle(label = unique(df_2_plot$Line))

df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]

ggplot(df_2_plot,
       aes(x = Genotype,
           y = Avg_exp,
           group = Genotype,
           colour = Genotype,
           fill = Genotype)) +
  stat_summary(geom = "crossbar",
               fun = "median",
               width = 0.5,
               linewidth = 0.5,
               colour = "black") +
  geom_jitter(width = 0.1,
              size = 2,
              shape = 1) +
  scale_fill_manual(values = c("darkred", 
                               "darkblue"),
                    guide = "none") +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue"),
                      guide = "none") +
  # ylim(0, 100) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.6)) +
  labs(x = "",
       y = "Normalized exp. value") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black")) +
  ggtitle(label = unique(df_2_plot$Line))

## Fig_2f ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 2)

df_2_plot <- df_raw
  # df_raw[df_raw$Line == "CD04", ]
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Avg_exp",
            # group = "Line",
            color = "black",
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
                     hide.ns = F) +
  geom_jitter(aes(colour = Genotype,
                  shape = Line),
              width = 0.05,
              size = 2) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Relative quantity") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.25)) +
  scale_shape_manual(values = c(1, 2)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0))

ggplot(df_2_plot,
       aes(x = Genotype,
           y = Avg_exp,
           group = Genotype,
           colour = Genotype,
           fill = Genotype,
           shape = Line)) +
  stat_summary(geom = "crossbar",
               fun = "median",
               width = 0.5,
               linewidth = 0.5,
               colour = "black") +
  geom_jitter(width = 0.1,
              size = 2) +
  scale_fill_manual(values = c("darkred", 
                               "darkblue"),
                    guide = "none") +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue"),
                      guide = "none") +
  scale_shape_manual(values = c(1,2)) +
  # ylim(0, 100) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  labs(x = "",
       y = "Relative quantity") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black")) #+
  ggtitle(label = unique(df_2_plot$Line))

df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]

ggplot(df_2_plot,
       aes(x = Genotype,
           y = Avg_exp,
           group = Genotype,
           colour = Genotype,
           fill = Genotype)) +
  stat_summary(geom = "crossbar",
               fun = "median",
               width = 0.5,
               linewidth = 0.5,
               colour = "black") +
  geom_jitter(width = 0.1,
              size = 2,
              shape = 1) +
  scale_fill_manual(values = c("darkred", 
                               "darkblue"),
                    guide = "none") +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue"),
                      guide = "none") +
  # ylim(0, 100) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  labs(x = "",
       y = "Relative quantity") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black")) +
  ggtitle(label = unique(df_2_plot$Line))

df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]

ggplot(df_2_plot,
       aes(x = Genotype,
           y = Avg_exp,
           group = Genotype,
           colour = Genotype,
           fill = Genotype)) +
  stat_summary(geom = "crossbar",
               fun = "median",
               width = 0.5,
               linewidth = 0.5,
               colour = "black") +
  geom_jitter(width = 0.1,
              size = 2,
              shape = 1) +
  scale_fill_manual(values = c("darkred", 
                               "darkblue"),
                    guide = "none") +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue"),
                      guide = "none") +
  # ylim(0, 100) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  labs(x = "",
       y = "Relative quantity") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black")) +
  ggtitle(label = unique(df_2_plot$Line))

# make PU1 motif ####
df_raw <-
  read.table(file = "~/Data/Tools/HOMER/motifs/pu1.motif",
             header = F,
             skip = 1,
             sep = "\t")

df_4_seqLogo <-
  as.data.frame(t(df_raw))
df_seqLogo_pwm <-
  makePWM(df_4_seqLogo)
seqLogo(df_seqLogo_pwm,
        ic.scale = F)


## Fig.2E ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 9)

df_2_plot <-
  df_raw

df_2_plot <-
  df_2_plot[df_2_plot$Line == "CD04", ]

df_2_plot <-
  df_2_plot[!(df_2_plot$Genotype == "CRISPRoff"), ]

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

df_2_plot$Avg_exp_normalised <-
  df_2_plot$Avg_exp / mean(df_2_plot$Avg_exp[df_2_plot$Genotype == "risk"])

# df_2_plot$Genotype <-
#   factor(df_2_plot$Genotype,
#          levels = c("risk",
#                     "non-risk"))

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Avg_exp_normalised",
            color = "black",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 'wilcox.test',
                     hide.ns = F) +
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
       y = "Normalized Expression Level") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 5),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) +
  ggtitle(unique(df_2_plot$Line))


df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]
# df_2_plot$Genotype <-
#   factor(df_2_plot$Genotype,
#          levels = unique(df_2_plot$Genotype))
df_2_plot <-
  df_2_plot[!(df_2_plot$Genotype == "CRISPRoff"), ]

df_2_plot$Avg_exp_normalised <-
  df_2_plot$Avg_exp / mean(df_2_plot$Avg_exp[df_2_plot$Genotype == "risk"])

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Avg_exp_normalised",
            color = "black",
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
                     hide.ns = F) +
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
       y = "Normalized Expression Level") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 5)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) +
  ggtitle(unique(df_2_plot$Line))

{
  library(lme4)
  library(lmerTest)
  library(rstatix)
  }


## Fig.2J ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 3)
df_2_plot <-
  df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

df_2_plot <-
  df_2_plot[df_2_plot$Line == "CD04", ]

df_2_plot$batch <-
  as.factor(c(1, 2, 1, 2))

lm_model <-
  lmerTest::lmer(Abundance ~ Genotype + (1|batch),
             data = df_2_plot)
summary(lm_model)
anova(lm_model)

# > summary(lm_model)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Abundance ~ Genotype + (1 | batch)
#    Data: df_2_plot
# 
# REML criterion at convergence: -20.2
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -0.5073 -0.4963  0.0000  0.4963  0.5073 
# 
# Random effects:
#  Groups   Name        Variance  Std.Dev. 
#  batch    (Intercept) 4.048e-05 0.0063621
#  Residual             1.747e-08 0.0001322
# Number of obs: 4, groups:  batch, 2
# 
# Fixed effects:
#                   Estimate Std. Error        df t value Pr(>|t|)   
# (Intercept)      0.0878964  0.0044997 1.0004315   19.53  0.03252 * 
# Genotypenon-risk 0.0542748  0.0001322 1.0000000  410.67  0.00155 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr)
# Gntypnn-rsk -0.015

# diag(length(coef(lm_model)))[-1, ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Abundance",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
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
                     # method = 't.test',
                     method = 'wilcox.test',
                     hide.ns = F,
                     # ref.group = "risk",
                     paired = F) +
  labs(x = "",
       y = "Normalized abundance") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.2),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle(unique(df_2_plot$Line))



df_2_plot <-
  df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

df_2_plot <-
  df_2_plot[df_2_plot$Line == "CD09", ]

df_2_plot$batch <-
  as.factor(c(1, 2, 1, 2))

lm_model <-
  lmerTest::lmer(Abundance ~ Genotype + (1|batch),
                 data = df_2_plot)
# lm_model <-
#   lmerTest::lmer(Abundance ~ Genotype + (1 + Genotype|batch),
#                  data = df_2_plot)
# > lm_model <-
# +   lmerTest::lmer(Abundance ~ Genotype + (1 + Genotype|batch),
# +                  data = df_2_plot)
# Error: number of observations (=4) <= number of random effects (=4) for term (1 + Genotype | batch); 
# the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
summary(lm_model)
anova(lm_model)

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Abundance",
            color = "black",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
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
                     # method = 't.test',
                     method = 'wilcox.test',
                     hide.ns = F,
                     # ref.group = "risk",
                     paired = F) +
  labs(x = "",
       y = "Normalized abundance") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.2),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle(unique(df_2_plot$Line))


## Fig_S5b ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 4)

df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]

ggplot(df_2_plot,
       aes(x = Genotype,
           y = Avg_exp,
           group = Genotype,
           colour = Genotype,
           fill = Genotype)) +
  stat_summary(geom = "crossbar",
               fun = "median",
               width = 0.5,
               linewidth = 0.5,
               colour = "black") +
  geom_jitter(width = 0.1,
              size = 2,
              shape = 1) +
  scale_fill_manual(values = c("darkred", 
                               "darkblue"),
                    guide = "none") +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue"),
                      guide = "none") +
  # ylim(0, 100) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.04)) +
  labs(x = "",
       y = "Normalized exp. value") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        axis.text.x = element_text(size = 12,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle(label = unique(df_2_plot$Line))

df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]

ggplot(df_2_plot,
       aes(x = Genotype,
           y = Avg_exp,
           group = Genotype,
           colour = Genotype,
           fill = Genotype)) +
  stat_summary(geom = "crossbar",
               fun = "median",
               width = 0.5,
               linewidth = 0.5,
               colour = "black") +
  geom_jitter(width = 0.1,
              size = 2,
              shape = 1) +
  scale_fill_manual(values = c("darkred", 
                               "darkblue"),
                    guide = "none") +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue"),
                      guide = "none") +
  # ylim(0, 100) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.04)) +
  labs(x = "",
       y = "Normalized exp. value") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        axis.text.x = element_text(size = 12,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle(label = unique(df_2_plot$Line))


## Fig_S5d ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 5)

df_2_plot <-
  df_raw

ggplot(df_2_plot,
       aes(x = Seq_type,
           y = Value,
           # group = Seq_type,
           # colour = Genotype,
           fill = Genotype)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("darkred", 
                               "darkblue")) +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue"),
                      guide = "none") +
  # ylim(0, 100) +
  scale_y_continuous(expand = c(0, 0),
                     # limits = c(0, 1),
                     labels = scales::percent) +
  labs(x = "Sample type",
       y = "Area % of the signal for each base") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black")) 


# Fig 3B ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 7)

df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("Non-risk",
                    "Risk",
                    "CRISPRoff"))
df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "risk")


df_2_plot$variable <-
  factor(df_2_plot$variable,
         levels = c("0 min",
                    "45 min",
                    "90 min",
                    "135 min",
                    "180 min"))
df_2_plot$batch <-
  as.factor(rep_len(x = c(1, 2, 3, 4),
                    length.out = 60))

for (cell_time in unique(df_2_plot$variable)[2:5]) {
  df_2_plot_sub <-
    df_2_plot[df_2_plot$variable == cell_time, ]
  df_2_plot_sub$Genotype <-
    factor(df_2_plot_sub$Genotype,
           levels = c("Non-risk",
                      "Risk",
                      "CRISPRoff"))
  df_2_plot_sub$Genotype <-
    relevel(df_2_plot_sub$Genotype,
            ref = "Risk")
  # prin
  lm_model <-
    lmerTest::lmer(Value ~ Genotype + (1|batch),
                   data = df_2_plot_sub)
  # summary(lm_model)
  lm_output <-
    summary(lm_model)
  #   anova(lm_model)
  print(cell_time)
  print(lm_output$coefficients)
}
# 
# boundary (singular) fit: see help('isSingular')
# [1] "45 min"
# Estimate  Std. Error df   t value     Pr(>|t|)
# (Intercept)        0.34704374 0.006020324  9 57.645363 7.164501e-13
# GenotypeNon-risk   0.13780837 0.008514023  9 16.186045 5.815428e-08
# GenotypeCRISPRoff -0.04772096 0.008514023  9 -5.604984 3.321474e-04
# boundary (singular) fit: see help('isSingular')
# [1] "90 min"
# Estimate Std. Error df   t value     Pr(>|t|)
# (Intercept)       0.56687695 0.01018560  9 55.654742 9.821920e-13
# GenotypeNon-risk  0.17990745 0.01440461  9 12.489571 5.473251e-07
# GenotypeCRISPRoff 0.03461152 0.01440461  9  2.402808 3.971462e-02
# boundary (singular) fit: see help('isSingular')
# [1] "135 min"
# Estimate  Std. Error df   t value     Pr(>|t|)
# (Intercept)       0.65944885 0.009781905  9 67.415176 1.756050e-13
# GenotypeNon-risk  0.17540786 0.013833703  9 12.679748 4.809140e-07
# GenotypeCRISPRoff 0.07176904 0.013833703  9  5.187985 5.731909e-04
# [1] "180 min"
# Estimate Std. Error       df   t value     Pr(>|t|)
# (Intercept)       0.79988168 0.01104544 7.669769 72.417366 3.730130e-12
# GenotypeNon-risk  0.20011832 0.01312056 6.000000 15.252270 5.014839e-06
# GenotypeCRISPRoff 0.07110677 0.01312056 6.000000  5.419493 1.633575e-03

summary(lm_model)
# summary(lm_output)
lmerTest::ranova(lm_model, reduce.terms = F)
lmerTest::difflsmeans(lm_model,
                      test.effs)

ggplot(df_2_plot,
       aes(x = variable,
           y = Value,
           colour = Genotype,
           fill = Genotype,
           group = Genotype)) +
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
                     limits = c(0, 1.1)) +
  labs(x = "",
       y = "pHrodo intensity / cell\nNormalized to non-risk (180 min)") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle(unique(df_2_plot$Line))

df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("Non-risk",
                    "Risk",
                    "CRISPRoff"))


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
  scale_fill_manual(values = c("darkred", 
                               "darkblue",
                               "green4")) +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue",
                                 "green4")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  labs(x = "",
       y = "pHrodo intensity / cell\nNormalized to non-risk (180 min)") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle(unique(df_2_plot$Line))

df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "risk",
                    "CRISPRoff"))
df_2_plot$Time <-
  factor(df_2_plot$Time,
         levels = c("0 min",
                    "45 min",
                    "90 min",
                    "135 min",
                    "180 min"))

ggplot(df_2_plot,
       aes(x = Time,
           y = Value,
           colour = Genotype,
           fill = Genotype,
           group = Genotype)) +
  geom_line(linewidth = 1) +
  geom_point(shape = 21,
             colour = "black") +
  scale_fill_manual(values = c("darkred", 
                               "darkblue",
                               "green4")) +
  scale_colour_manual(values = c("darkred", 
                                 "darkblue",
                                 "green4")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  labs(x = "",
       y = "pHrodo intensity / cell\nNormalized to non-risk (180 min)") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle(unique(df_2_plot$Line))

# Fig S7B ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 8)

df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("Non-risk",
                    "Risk",
                    "CRISPRoff"))


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
                     limits = c(0, 1),
                     labels = scales::percent) +
  labs(x = "",
       y = "% of positive Ab-pHrodo cells") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle(unique(df_2_plot$Line))


df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("Non-risk",
                    "Risk",
                    "CRISPRoff"))


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
                     limits = c(0, 1),
                     labels = scales::percent) +
  labs(x = "",
       y = "% of positive Ab-pHrodo cells") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle(unique(df_2_plot$Line))


## Fig_3D ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 9)

df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = unique(df_2_plot$Genotype))
df_2_plot <-
  df_2_plot[!(df_2_plot$Genotype == "risk"), ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Avg_exp",
            color = "black",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F) +
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
       y = "Normalized Expression Level") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) +
  ggtitle(unique(df_2_plot$Line))


df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = unique(df_2_plot$Genotype))
df_2_plot <-
  df_2_plot[!(df_2_plot$Genotype == "risk"), ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Avg_exp",
            color = "black",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "green4"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F) +
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
       y = "Normalized Expression Level") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) +
  ggtitle(unique(df_2_plot$Line))


## new Fig E6A-B ####
df_raw <-
  read_excel("new_Fig_E6ab_by_clone.xlsx",
             sheet = 1)

df_2_plot <-
  df_raw
df_2_plot <-
  melt(df_2_plot)

colnames(df_2_plot) <-
  c("Cell_line", 
    "Genotype", 
    "Gene", 
    "Value")


df_2_plot <-
  df_2_plot[(df_2_plot$Gene == "ATP6AP2"), ]

df_2_plot <-
  df_2_plot[!(df_2_plot$Gene == "ATP6AP2"), ]
df_2_plot$Gene <-
  factor(df_2_plot$Gene,
         levels = c("NEAT1",
                    "LDLR",
                    "HMGCS1",
                    "DHCR7"))

df_2_plot$Cell_line <-
  as.factor(df_2_plot$Cell_line)
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))
# df_2_plot <-
#   df_2_plot[!(df_2_plot$Gene == "ATP6AP2"), ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            # shape = "Cell_line",
            color = "black",
            width = 0.2,
            facet.by = "Gene",
            ncol = 4,
            # add = "jitter",
            # add.params = list(shape = c(1, 2)),
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +

  stat_compare_means(label = "p.signif",
                     label.y.npc = .8,
                     method = 't.test',
                     method.args = list(var.equal = T),
                     hide.ns = F) +
  geom_jitter(aes(colour = Genotype,
                  shape = Cell_line),
              width = 0.05,
              size = 2) +
  scale_shape_manual(values = c(21, 24)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Normalized Exp. Level") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.5)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) 

df_test_lmm <-
  df_2_plot[df_2_plot$Gene == "LDLR", ]
lm_model <-
  lmerTest::lmer(Value ~ Genotype + (1|Cell_line),
                 data = df_test_lmm)

lm_results <-
  summary(lm_model)
lm_results$coefficients[, 1]

for (i in unique(df_2_plot$Gene)) {
  print(i)
  df_test_lmm <-
    df_2_plot[df_2_plot$Gene == i, ]
  print(summary(lmerTest::lmer(Value ~ Genotype + (1|Cell_line),
                               data = df_test_lmm))$coefficients[, 5])
}





## Fig_S10A ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 10)

# test if stat_compare_means can be used with ggplot?
df_2_plot <-
  df_raw

# df_2_plot$rep <-
  





df_2_plot <-
  melt(df_2_plot)
colnames(df_2_plot) <-
  c("Genotype", "Gene", "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

# df_2_plot <-
#   df_2_plot[!(df_2_plot$variable == '0 min'), ]

for (gene_name in unique(df_2_plot$Gene)) {
  df_2_plot_sub <-
    df_2_plot[df_2_plot$Gene == gene_name, ]
  print(gene_name)
  # df_2_plot_sub$Genotype <-
  #   factor(df_2_plot_sub$Genotype,
  #          levels = c("Non-risk",
  #                     "Risk",
  #                     "CRISPRoff"))
  df_2_plot_sub$Genotype <-
    relevel(df_2_plot_sub$Genotype,
            ref = "risk")
  df_2_plot_sub$Batch <-
    rep_len(x = c(1,1,1,2,2,2),
            length.out = nrow(df_2_plot_sub))
  # prin
  lm_model <-
    lmerTest::lmer(Value ~ Genotype + (1|Batch),
                   data = df_2_plot_sub)
  # summary(lm_model)
  lm_output <-
    summary(lm_model)
  print(lm_output)
  #   anova(lm_model)
  # print(cell_time)
  # print(lm_output$coefficients)
}




# df_2_plot <-
#   df_2_plot[(df_2_plot$Gene == "ATP6AP2"), ]

df_2_plot <-
  df_2_plot[!(df_2_plot$Gene == "ATP6AP2"), ]
df_2_plot$Gene <-
  factor(df_2_plot$Gene,
         levels = c("NEAT1",
                    "LDLR",
                    "HMGCS1",
                    "DHCR7"))

# df_2_plot <-
#   df_2_plot[!(df_2_plot$Gene == "ATP6AP2"), ]

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            width = 0.2,
            facet.by = "Gene",
            ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F) +
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
       y = "Normalized Expression Level") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.5)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) 

## Fig_S10D ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 11)

# test if stat_compare_means can be used with ggplot?
df_2_plot <-
  df_raw
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

df_2_plot$Gene <-
  factor(df_raw$Gene,
         levels = c("ATP6AP2",
                    "VAMP1",
                    "HMGCR",
                    "CD74"))


for (gene_name in unique(df_2_plot$Gene)) {
  df_2_plot_sub <-
    df_2_plot[df_2_plot$Gene == gene_name, ]
  print(gene_name)
  # df_2_plot_sub$Genotype <-
  #   factor(df_2_plot_sub$Genotype,
  #          levels = c("Non-risk",
  #                     "Risk",
  #                     "CRISPRoff"))
  df_2_plot_sub$Genotype <-
    relevel(df_2_plot_sub$Genotype,
            ref = "risk")
  df_2_plot_sub$Batch <-
    rep_len(x = c(1,1,1,2,2,2),
            length.out = nrow(df_2_plot_sub))
  # prin
  lm_model <-
    lmerTest::lmer(Value ~ Genotype + (1|Batch),
                   data = df_2_plot_sub)
  # summary(lm_model)
  lm_output <-
    summary(lm_model)
  print(lm_output)
  #   anova(lm_model)
  # print(cell_time)
  # print(lm_output$coefficients)
}


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            width = 0.2,
            facet.by = "Gene",
            ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .5,
                     method = 't.test',
                     hide.ns = F) +
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
       y = "Normalized Fluor. Intensity") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 400)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) 


## new_Fig_E6D ####
df_raw <-
  read_excel("Alena_DEG_staining_Fig_E6D_16Apr2025.xlsx",
             sheet = 1)

# test if stat_compare_means can be used with ggplot?
df_2_plot <-
  df_raw

df_2_plot$BR <-
  str_split(df_2_plot$Replicates,
            pattern = '_',
            simplify = T)[, 1]
df_2_plot$UUID <-
  str_c(df_2_plot$Genotype,
        df_2_plot$Gene,
        df_2_plot$Cell_line,
        df_2_plot$BR,
        sep = '.')
df_2_plot$plot_value <-
  df_2_plot$MEAN

# df_2_plot_backup <-
#   df_2_plot
df_2_plot <-
  df_2_plot_backup

df_2_plot <-
  df_2_plot[, c(1,2,9,10,11,12)]

df_2_plot_BR <-
  df_2_plot %>%
  group_by(UUID) %>%
  dplyr::summarise(plot_value_BR = mean(plot_value,
                                        na.rm = T))

df_2_plot_BR$Genotype <-
  str_split(df_2_plot_BR$UUID,
            pattern = '\\.',
            simplify = T)[, 1]
df_2_plot_BR$Gene <-
  str_split(df_2_plot_BR$UUID,
            pattern = '\\.',
            simplify = T)[, 2]
df_2_plot_BR$Cell_line <-
  str_split(df_2_plot_BR$UUID,
            pattern = '\\.',
            simplify = T)[, 3]
df_2_plot_BR$BR <-
  # df_2_plot_BR$Cell_line <-
  str_split(df_2_plot_BR$UUID,
            pattern = '\\.',
            simplify = T)[, 4]

df_2_plot_BR$UUID <- NULL

df_2_plot_BR$Gene[df_2_plot_BR$Gene == "ATP"] <- "ATP6AP2"





df_2_plot_BR$Genotype <-
  factor(df_2_plot_BR$Genotype,
         levels = c("risk",
                    "non-risk"))

df_2_plot_BR$Gene <-
  factor(df_2_plot_BR$Gene,
         levels = c("ATP6AP2",
                    "VAMP1",
                    "HMGCR",
                    "CD74"))

colnames(df_2_plot_BR)[1] <- "Value"

for (gene_name in unique(df_2_plot_BR$Gene)) {
  # df_2_plot_sub <-
  #   df_2_plot_BR[df_2_plot_BR$Gene == gene_name, ]
  print(gene_name)

  # summary(lm_model)
  
  df_test_lmm <-
    df_2_plot_BR[df_2_plot_BR$Gene == gene_name, ]
  print(summary(lmerTest::lmer(Value ~ Genotype + (1|Cell_line),
                               data = df_test_lmm))$coefficients[, 5])
  
  # 
  # lm_output <-
  #   summary(lm_model)
  # print(lm_output)
  #   anova(lm_model)
  # print(cell_time)
  # print(lm_output$coefficients)
}


ggerrorplot(df_2_plot_BR,
            x = "Genotype",
            y = "Value",
            color = "black",
            width = 0.2,
            facet.by = "Gene",
            ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = .75,
                     method = 't.test',
                     hide.ns = F) +
  geom_jitter(aes(colour = Genotype,
                  shape = Cell_line),
              width = 0.05,
              size = 2) +
  scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Normalized Fluor. Intensity") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 400)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) 


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            # shape = "Cell_line",
            color = "black",
            width = 0.2,
            facet.by = "Gene",
            ncol = 4,
            # add = "jitter",
            # add.params = list(shape = c(1, 2)),
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  
  stat_compare_means(label = "p.signif",
                     label.y.npc = .8,
                     method = 't.test',
                     method.args = list(var.equal = T),
                     hide.ns = F) +
  geom_jitter(aes(colour = Genotype,
                  shape = Cell_line),
              width = 0.05,
              size = 2) +
  scale_shape_manual(values = c(21, 24)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Normalized Exp. Level") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.5)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) 


## Fig_S10D ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 12)

df_2_plot <- df_raw
df_2_plot <-
  melt(df_2_plot)
# iMG

df_2_plot <-
  df_2_plot[df_2_plot$variable == "iMG_purity", ]
df_2_plot$

ggerrorplot(df_2_plot,
            x = "variable",
            y = "value",
            color = "black",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkblue"),
                      guide = "none") +
  # stat_compare_means(label = "p.signif",
  #                    label.y.npc = .5,
  #                    method = 't.test',
  #                    hide.ns = F) +
  geom_jitter(colour = "darkblue",
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "% of TREM+ CD45+ PU1+") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     labels = scales::percent) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) 

df_2_plot <- df_raw
df_2_plot <-
  melt(df_2_plot)
# iMG

df_2_plot <-
  df_2_plot[df_2_plot$variable %in% c("iAst_Vim_S100b",
                                      "iAst_GFAP"), ]

ggerrorplot(df_2_plot,
            x = "variable",
            y = "value",
            color = "black",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkblue"),
                      guide = "none") +
  # stat_compare_means(label = "p.signif",
  #                    label.y.npc = .5,
  #                    method = 't.test',
  #                    hide.ns = F) +
  geom_jitter(colour = "darkblue",
              width = 0.05,
              size = 2,
              shape = 1) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "% of positive cells") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     labels = scales::percent) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) 

## Fig_5B ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 13)

df_2_plot <- df_raw

df_2_plot$Cell_line <-
  c(rep_len(x = "CD04",
            length.out = 6),
    rep_len(x = "CD09",
            length.out = 6),
    rep_len(x = "CD04",
            length.out = 4),
    rep_len(x = "CD09",
            length.out = 4))

df_2_plot$BR <-
  c(1,1,1,2,2,2,
    1,1,1,2,2,2,
    1,1,2,2,
    1,1,2,2)

lm_model <-
  lmerTest::lmer(`Relative Fluorescence` ~ 
                   Genotype + 
                   (1|Cell_line) +
                   (1|BR),
                 data = df_2_plot)
summary(lm_model)

# iMG

# df_2_plot <-
#   df_2_plot[df_2_plot$variable == "iMG_purity", ]




ggerrorplot(df_2_plot,
            x = "Genotype",
            y = 'Relative Fluorescence',
            color = "black",
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
                     hide.ns = F) +
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
       y = "Relative Fluorescence") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.5),
                     labels = scales::percent) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 0)) 

## Fig_5D ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 14)

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
                     # method = 't.test',
                     method = 'wilcox.test',
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
                     limits = c(0, 15)) +
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
                     # method = 't.test',
                     method = 'wilcox.test',
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
       y = "Puncta area (um2)") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 15)) +
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
  read_excel("tables_4_plot.xlsx",
             sheet = 15)

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
       y = "LD puncta size/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 5),
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
       y = "Puncta area (um2)") +
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

## Fig_5F ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 16)

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
       y = "LD area/iMG") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 5),
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
       y = "LD size (um2)") +
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
  read_excel("tables_4_plot.xlsx",
             sheet = 17)

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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 30),
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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 30),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD09")


## Fig_Ex_7e-h ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 18)

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
  df_2_plot[df_2_plot$Cell_line == "CD04", ]


df_2_plot$BR <-
  as.factor(c(1,1,1,
              2,2,2,
              3,3,3,3))

lm_model <-
  lmerTest::lmer(Value ~ Genotype + (1|Batch),
                 data = df_2_plot)
# summary(lm_model)
lm_output <-
  summary(lm_model)
print(lm_output)

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
       y = "Value") +
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
       y = "Value") +
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
                                   vjust = 0)) +
  ggtitle("CD09")


## Fig_S11C ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 19)

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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 5),
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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 5),
                     na.value = NA) +
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
  read_excel("tables_4_plot.xlsx",
             sheet = 20)

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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 250),
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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 250),
                     na.value = NA) +
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
  read_excel("tables_4_plot.xlsx",
             sheet = 21)

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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1600),
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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1600),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD09")

## Fig_S12C ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 22)

df_2_plot <- df_raw
# df_2_plot <-
#   melt(df_2_plot)
# colnames(df_2_plot) <- 
#   c("Cell_line",
#     "Genotype",
#     "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))


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
       y = "Value") +
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

df_2_plot$Batch <-
  rep_len(x = c(1,1,1,1,1,
                2,2,2,2,2,
                3,3,3,3,3),
          length.out = nrow(df_2_plot))

lm_model <-
  lmerTest::lmer(Value ~
                   Genotype +
                   (1|Batch),
                 data = df_2_plot)
summary(lm_model)





## Fig_S12D ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 23)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))


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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 450),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

## Fig_S12F ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 24)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "CRISPRoff"))


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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 25),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 


df_2_plot$Batch <-
  rep_len(x = c(1,1,2,2,3,3,
                1,1,2,2,3,3,3),
          length.out = nrow(df_2_plot))

lm_model <-
  lmerTest::lmer(Value ~
                   Genotype +
                   (1|Batch),
                 data = df_2_plot)
summary(lm_model)








## Fig_S12G ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 25)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "CRISPRoff"))


p_output <-
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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 75),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

print(p_output)
df_2_plot$Batch <-
  rep_len(x = c(1,1,2,2,3,3,
                1,1,2,2,3,3,3),
          length.out = nrow(df_2_plot))

lm_model <-
  lmerTest::lmer(Value ~
                   Genotype +
                   (1|Batch),
                 data = df_2_plot)
summary(lm_model)

dev.off()
pdf(file = "updated_Fig_E8J_13Mar2025.pdf",
    width = 1.15,
    height = 2.54)
print(p_output)
dev.off()


## Fig_S12H ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 26)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "CRISPRoff"))


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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 25),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

df_2_plot$Batch <-
  rep_len(x = c(1,1,2,2,3,3,
                1,1,2,2,3,3,3),
          length.out = nrow(df_2_plot))

lm_model <-
  lmerTest::lmer(Value ~
                   Genotype +
                   (1|Batch),
                 data = df_2_plot)
summary(lm_model)


## Fig_S12I ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 27)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "CRISPRoff"))

pdf(file = "updated_Fig_E8L_13Mar2025.pdf",
    width = 1.15,
    height = 2.54)


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
       y = "Value") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 75),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

# print(p_output)
dev.off()


df_2_plot$Batch <-
  rep_len(x = c(1,1,2,2,3,3,
                1,1,2,2,3,3,3),
          length.out = nrow(df_2_plot))

lm_model <-
  lmerTest::lmer(Value ~
                   Genotype +
                   (1|Batch),
                 data = df_2_plot)
summary(lm_model)



## Fig_5I ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 28)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk-TrC",
                    '',
                    "non-risk", "non-risk-TrC"))


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
                     comparisons = list(c("risk", "risk-TrC"),
                                        c("non-risk", "non-risk-TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "Value") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 20),
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
  read_excel("tables_4_plot.xlsx",
             sheet = 29)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk-TrC",
                    '',
                    "non-risk", "non-risk-TrC"))


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
                     label.y.npc = 0.5,
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("risk", "risk-TrC"),
                                        c("non-risk", "non-risk-TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "Value") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 30),
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
  read_excel("tables_4_plot.xlsx",
             sheet = 30)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk-TrC",
                    '',
                    "non-risk", "non-risk-TrC"))


pdf("Fig_5K.pdf",
    width = 2.40,
    height = 2.19,
    compress = T)

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
                     label.y.npc = 0.5,
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("risk", "risk-TrC"),
                                        c("non-risk", "non-risk-TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "Value") +
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
dev.off()

## Fig_5L ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 31)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk-TrC",
                    '',
                    "non-risk", "non-risk-TrC"))


pdf("Fig_5L.pdf",
    width = 2.40,
    height = 2.19,
    compress = T)

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
                     label.y.npc = 0.5,
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("risk", "risk-TrC"),
                                        c("non-risk", "non-risk-TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "Value") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 30),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 
dev.off()

## Fig_S13C ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 32)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk-TrC",
                    '',
                    "non-risk", "non-risk-TrC"))


# pdf("Fig_5L.pdf",
#     width = 2.40,
#     height = 2.19,
#     compress = T)

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
                     label.y.npc = 0.5,
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("risk", "risk-TrC"),
                                        c("non-risk", "non-risk-TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "Value") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 20),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 
# dev.off()


df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk-TrC",
                    '',
                    "non-risk", "non-risk-TrC"))
df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "non-risk")
df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "risk")

df_2_plot$Batch <-
  rep_len(x = c(1,1,1,2,2,2,3,3,
                1,1,1,2,2,2,3,3,3),
          length.out = nrow(df_2_plot))

lm_model <-
  lmerTest::lmer(Value ~
                   Genotype +
                   (1|Batch),
                 data = df_2_plot)
summary(lm_model)






## Fig_S13D ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 33)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk", "risk-TrC",
                    '',
                    "non-risk", "non-risk-TrC"))

# 
# pdf("Fig_5L.pdf",
#     width = 2.40,
#     height = 2.19,
#     compress = T)

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
                     label.y.npc = 0.5,
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("risk", "risk-TrC"),
                                        c("non-risk", "non-risk-TrC"),
                                        c("risk", "non-risk"))) +
  labs(x = "",
       y = "Value") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 25),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

# df_2_plot$Genotype <-
#   factor(df_2_plot$Genotype,
#          levels = c("risk", "risk-TrC",
#                     '',
#                     "non-risk", "non-risk-TrC"))
df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "non-risk")
df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "risk")

df_2_plot$Batch <-
  rep_len(x = c(1,1,1,2,2,2,3,3,
                1,1,1,2,2,2,3,3,3),
          length.out = nrow(df_2_plot))

lm_model <-
  lmerTest::lmer(Value ~
                   Genotype +
                   (1|Batch),
                 data = df_2_plot)
summary(lm_model)


# dev.off()

## Fig_5H ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 34)

df_2_plot <- df_raw

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = unique(df_2_plot$Genotype))

unique(df_2_plot$Genotype)
df_2_plot$Category <-
  c(rep_len("untreated", 
            length.out = 8),
    rep_len("CytoD", 
            length.out = 6))

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            shape = "Category",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkblue",
                                 "orange4",
                                 "darkblue",
                                 "orange4"),
                      guide = "none") +
  scale_shape_manual(values = c(1, 2)) +
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
                     label.y.npc = 0.5,
                     # method = 't.test',
                     method = "wilcox.test",
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("control", "KO"),
                                        c("control", "control_CytoD"),
                                        c("KO", "KO_CytoD"))) +
  labs(x = "",
       y = "Fluor. intensity/cell") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  # scale_y_continuous(expand = c(0, 0),
  #                    limits = c(0, 15000),
  #                    na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

summary(aov(Value ~ Genotype,
            data = df_2_plot))
summary(multcomp::glht(aov(Value ~ Genotype,
                           data = df_2_plot),
                       linfct = mcp(Genotype = "Dunnett")))

kruskal.test(x = df_2_plot$Value,
             g = df_2_plot$Genotype)
DescTools::DunnTest(Value ~ Genotype,
                    data = df_2_plot)

# Dunn's test of multiple comparisons using rank sums : holm  
# 
#                        mean.rank.diff   pval    
# KO-control                       -4.0 0.4770    
# control_CytoD-control            -8.5 0.0390 *  
# KO_CytoD-control                 -9.5 0.0177 *  
# control_CytoD-KO                 -4.5 0.4770    
# KO_CytoD-KO                      -5.5 0.3407    
# KO_CytoD-control_CytoD           -1.0 0.7697    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


## Fig_S12F ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 36)

df_2_plot <- df_raw

# df_2_plot$Genotype <-
#   factor(df_2_plot$Genotype,
#          levels = c("non-risk",
#                     "CRISPRoff"))


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
                                 "orange"),
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
                     ref.group = "control",
                     paired = F) +
  labs(x = "",
       y = "Median droplets/cell ") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 50),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

## Fig.2D ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 37)

df_2_plot <-
  df_raw

# df_2_plot <-
#   df_2_plot[df_2_plot$Line == "CD04", ]
# 
# df_2_plot <-
#   df_2_plot[!(df_2_plot$Genotype == "CRISPRoff"), ]

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

# df_2_plot$Avg_exp_normalised <-
#   df_2_plot$Avg_exp / mean(df_2_plot$Avg_exp[df_2_plot$Genotype == "risk"])

# df_2_plot$Genotype <-
#   factor(df_2_plot$Genotype,
#          levels = c("risk",
#                     "non-risk"))

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "black",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  # stat_compare_means(label = "p.signif",
  #                    label.y.npc = .5,
  #                    method = 't.test',
  #                    hide.ns = F) +
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
       y = '% of TREM2+,CD45+,PU1+ cells') +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0))

## Fig.2K ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 38)

df_2_plot <-
  df_raw

ggplot(df_2_plot,
       aes(x = Sample,
           y = Value,
           fill = Genotype)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  scale_fill_manual(values = c("darkred",
                               "darkblue")) +
  scale_y_continuous(expand = c(0,0),
                     labels = scales::percent,
                     na.value = NA) +
  labs(x = "",
       y = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 0,
                                   vjust = 0,
                                   hjust = 0))
