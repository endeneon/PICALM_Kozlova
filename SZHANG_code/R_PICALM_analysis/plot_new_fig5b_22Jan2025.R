# Siwei 22 Jan 2025
# plot new Fig. 5b

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

# CD04 ####
df_raw <-
  read_excel("Batch_2_of_new_plots/Fig_5b.xlsx",
             sheet = 1)

df_2_plot <-
  df_raw
df_2_plot$`0 min` <- 0

df_2_plot$Geno_batch <-
  str_c(df_2_plot$Genotype,
        df_raw$Batch,
        str_split(string = df_2_plot$BR_FOV,
                  pattern = "_",
                  simplify = T)[, 1],
        sep = "_")

## if each BR_FOV as one data point #####
df_2_plot$Batch <- NULL
df_2_plot$BR_FOV <- NULL
# df_2_plot$Geno_batch <- NULL

df_2_plot <-
  reshape2::melt(df_2_plot,
                 value.name = "Value")

df_2_plot$Genotype[df_2_plot$Genotype == "CRISPRa"] <-
  "risk-CRISPRa"

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "risk-CRISPRa"))

df_2_plot$variable <-
  factor(df_2_plot$variable,
         levels = c('0 min',
                    '45 min',
                    '90 min',
                    '135 min',
                    '180 min'))

# df_2_plot$batch <-
#   as.factor(rep_len(x = c(1, 2, 3, 4),
#                     length.out = 60))

# for (cell_time in unique(df_2_plot$variable)[2:5]) {
#   df_2_plot_sub <-
#     df_2_plot[df_2_plot$variable == cell_time, ]
#   # df_2_plot_sub$Genotype <-
#   #   factor(df_2_plot_sub$Genotype,
#   #          levels = c("Non-risk",
#   #                     "Risk",
#   #                     "CRISPRoff"))
#   df_2_plot_sub$Genotype <-
#     relevel(df_2_plot_sub$Genotype,
#             ref = "Risk")
#   # prin
#   lm_model <-
#     lmerTest::lmer(Value ~ Genotype + (1|batch),
#                    data = df_2_plot_sub)
#   # summary(lm_model)
#   lm_output <-
#     summary(lm_model)
#   #   anova(lm_model)
#   print(cell_time)
#   print(lm_output$coefficients)
# }

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
                fun.data = mean_se,
                width = 0.3) +
  geom_point(shape = 21,
             stat = "summary",
             fun = mean,
             colour = "black") +


  # geom_point(shape = 21,
  #            colour = "black") +
  scale_fill_manual(values = c("darkblue",
                               "darkred",
                               "green4")) +
  scale_colour_manual(values = c("darkblue",
                                 "darkred",
                                 "green4")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  labs(x = "",
       y = "pHrodo intensity / cell\nNormalized to non-risk (180 min)") +
  theme_classic() +
  # ylim(0, 0.8) +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  ggtitle("CD04")

ggplot(df_2_plot,
       aes(x = variable,
           y = Value,
           # colour = Genotype,
           # fill = Genotype,
           )) +
  geom_boxplot(aes(fill = Genotype),
               outliers = F,
               alpha = 0.5) +
  geom_point(aes(colour = Genotype),
             position = position_jitterdodge()) +
  scale_fill_manual(values = c("darkblue",
                               "darkred",
                               "green4")) +
  scale_colour_manual(values = c("black",
                                 "black",
                                 "black")) +
  theme_classic() +
  ggtitle("CD04")

df_2_plot <-
  df_2_plot[!(df_2_plot$variable == '0 min'), ]

### 2-way anova #####
summary(aov(Value ~ Genotype + variable + Genotype:variable,
            data = df_2_plot))
summary(aov(Value ~ Genotype * variable,
            data = df_2_plot))

View(TukeyHSD(aov(Value ~ Genotype + variable + Genotype:variable,
                  data = df_2_plot))[[3]])


### Dunnett's test (against single reference group) #####
df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "risk")

summary(multcomp::glht(aov(Value ~ Genotype * variable,
                           data = df_2_plot),
                       linfct = mcp(Genotype = "Dunnett")))
# write.table(x = TukeyHSD(aov(Value ~ Genotype + variable + Genotype:variable,
#                   data = df_2_plot))[[3]],
#             file = "ANOVA_post_test_CD09.tsv",
#             quote = F, sep = "\t")


## take mean of each two replicates #####
df_2_plot_test <-
  df_2_plot %>%
  group_by(Geno_batch) %>%
  dplyr::reframe(across(c('45 min',
                          '90 min',
                          '135 min',
                          '180 min',
                          '0 min')),
                 mean,
                 .unpack = T,
                 .by = Geno_batch)
# df_2_plot_mean <-
#   data.table::setDT(df_2_plot)


df_2_plot <-
  data.table::setDT(df_2_plot)[ , lapply(X = .SD,
                                        FUN = base::mean),
                               keyby = Geno_batch]


df_2_plot$BR_FOV <- NULL
df_2_plot$Genotype <-
  str_split(string = df_2_plot$Geno_batch,
            pattern = "_",
            simplify = T)[, 1]
df_2_plot$Batch <-
  str_split(string = df_2_plot$Geno_batch,
            pattern = "_",
            simplify = T)[, 2]
df_2_plot$Geno_batch <- NULL
df_2_plot$variable <- NULL

df_2_plot <-
  reshape2::melt(df_2_plot,
                 value.name = "Value")

df_2_plot <-
  df_2_plot[!(df_2_plot$variable == '0 min'), ]
df_2_plot$Genotype[df_2_plot$Genotype == "CRISPRa"] <-
  "risk-CRISPRa"

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "risk-CRISPRa"))

df_2_plot$variable <-
  factor(df_2_plot$variable,
         levels = c('45 min',
                    '90 min',
                    '135 min',
                    '180 min'))

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
                fun.data = mean_se,
                width = 0.3) +
  geom_point(shape = 21,
             stat = "summary",
             fun = mean,
             colour = "black") +


  # geom_point(shape = 21,
  #            colour = "black") +
  scale_fill_manual(values = c("darkblue",
                               "darkred",
                               "green4")) +
  scale_colour_manual(values = c("darkblue",
                                 "darkred",
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
  ggtitle("CD04")

summary(aov(Value ~ Genotype + variable + Genotype:variable,
            data = df_2_plot))

TukeyHSD(aov(Value ~ Genotype + variable + Genotype:variable,
             data = df_2_plot))

for (cell_time in unique(df_2_plot$variable)) {
  df_2_plot_sub <-
    df_2_plot[df_2_plot$variable == cell_time, ]
  # df_2_plot_sub$Genotype <-
  #   factor(df_2_plot_sub$Genotype,
  #          levels = c("Non-risk",
  #                     "Risk",
  #                     "CRISPRoff"))
  df_2_plot_sub$Genotype <-
    relevel(df_2_plot_sub$Genotype,
            ref = "risk")
  # prin
  lm_model <-
    lmerTest::lmer(Value ~ Genotype + (1|Batch),
                   data = df_2_plot_sub)
  # summary(lm_model)
  lm_output <-
    summary(lm_model)
  #   anova(lm_model)
  print(cell_time)
  print(lm_output$coefficients)
}

# [1] "45 min"
# Estimate Std. Error        df   t value    Pr(>|t|)
# (Intercept)          0.009023137 0.02661816  1.603541 0.3389843 0.773733101
# Genotypenon-risk     0.073792350 0.02124519 14.000000 3.4733675 0.003727388
# Genotyperisk-CRISPRa 0.024785643 0.02124519 14.000000 1.1666473 0.262841671
# [1] "90 min"
# Estimate Std. Error        df  t value    Pr(>|t|)
# (Intercept)          0.04431799 0.03662750  1.324666 1.209965 0.399215430
# Genotypenon-risk     0.07724542 0.02303702 14.000000 3.353099 0.004734678
# Genotyperisk-CRISPRa 0.06187567 0.02303702 14.000000 2.685924 0.017738767
# [1] "135 min"
# Estimate Std. Error        df   t value     Pr(>|t|)
# (Intercept)          0.08421299 0.12756941  1.121402 0.6601346 0.6185955768
# Genotypenon-risk     0.22257990 0.05219316 14.000000 4.2645416 0.0007856786
# Genotyperisk-CRISPRa 0.15473381 0.05219316 14.000000 2.9646377 0.0102446670
# [1] "180 min"
# Estimate Std. Error        df   t value     Pr(>|t|)
# (Intercept)          0.1885682 0.23175515  1.085369 0.8136528 0.5563319668
# Genotypenon-risk     0.3418142 0.08047691 14.000001 4.2473571 0.0008122977
# Genotyperisk-CRISPRa 0.2087193 0.08047691 14.000001 2.5935298 0.0212431040

# write.table(df_2_plot,
#             file = "df_2_test_")

# TukeyHSD(Value ~ Genotype + variable + Genotype::variable,
#          data = df_2_plot)


# CD09 ####
df_raw <-
  read_excel("Batch_2_of_new_plots/Fig_3b.xlsx",
             sheet = 2)

df_2_plot <-
  df_raw
df_2_plot$`0 min` <- 0

df_2_plot$Geno_batch <-
  str_c(df_2_plot$Genotype,
        df_raw$Batch,
        str_split(string = df_2_plot$BR_FOV,
                  pattern = "_",
                  simplify = T)[, 1],
        sep = "_")

## if each BR_FOV as one data point #####
df_2_plot$Batch <- NULL
df_2_plot$BR_FOV <- NULL
df_2_plot$Geno_batch <- NULL

df_2_plot <-
  reshape2::melt(df_2_plot,
                 value.name = "Value")

df_2_plot$Genotype[df_2_plot$Genotype == "CRISPRa"] <-
  "risk-CRISPRa"

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "risk-CRISPRa"))

df_2_plot$variable <-
  factor(df_2_plot$variable,
         levels = c('0 min',
                    '45 min',
                    '90 min',
                    '135 min',
                    '180 min'))


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
                fun.data = mean_se,
                width = 0.3) +
  geom_point(shape = 21,
             stat = "summary",
             fun = mean,
             colour = "black") +


  # geom_point(shape = 21,
  #            colour = "black") +
  scale_fill_manual(values = c("darkblue",
                               "darkred",
                               "green4")) +
  scale_colour_manual(values = c("darkblue",
                                 "darkred",
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
  ggtitle("CD09")

ggplot(df_2_plot,
       aes(x = variable,
           y = Value,
           # colour = Genotype,
           # fill = Genotype,
       )) +
  geom_boxplot(aes(fill = Genotype),
               outliers = F,
               alpha = 0.5) +
  geom_point(aes(colour = Genotype),
             position = position_jitterdodge()) +
  scale_fill_manual(values = c("darkblue",
                               "darkred",
                               "green4")) +
  scale_colour_manual(values = c("black",
                                 "black",
                                 "black")) +
  theme_classic() +
  ggtitle("CD09")

df_2_plot <-
  df_2_plot[!(df_2_plot$variable == '0 min'), ]

summary(aov(Value ~ Genotype + variable + Genotype:variable,
            data = df_2_plot))

TukeyHSD(aov(Value ~ Genotype + variable + Genotype:variable,
             data = df_2_plot))

## take mean of each two replicates #####
# df_2_plot <-
#   df_2_plot %>%
#   # group_by(batch_BR_FOV) %>%
#   dplyr::reframe(across(c('45 min',
#                           '90 min',
#                           '135 min',
#                           '180 min',
#                           '0 min')),
#                  mean,
#                  .unpack = T,
#                  .by = batch_BR_FOV)
# df_2_plot_mean <-
#   data.table::setDT(df_2_plot)


df_2_plot <-
  data.table::setDT(df_2_plot)[ , lapply(X = .SD,
                                         FUN = base::mean),
                                keyby = Geno_batch]

df_2_plot$BR_FOV <- NULL
df_2_plot$Genotype <-
  str_split(string = df_2_plot$Geno_batch,
            pattern = "_",
            simplify = T)[, 1]
df_2_plot$Batch <-
  str_split(string = df_2_plot$Geno_batch,
            pattern = "_",
            simplify = T)[, 2]
df_2_plot$Geno_batch <- NULL


df_2_plot <-
  reshape2::melt(df_2_plot,
                 value.name = "Value")

df_2_plot$Genotype[df_2_plot$Genotype == "CRISPRa"] <-
  "risk-CRISPRa"

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "risk-CRISPRa"))

df_2_plot$variable <-
  factor(df_2_plot$variable,
         levels = c('0 min',
                    '45 min',
                    '90 min',
                    '135 min',
                    '180 min'))

df_2_plot$Batch <-
  as.factor(df_2_plot$Batch)

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
                fun.data = mean_se,
                width = 0.3) +
  geom_point(shape = 21,
             stat = "summary",
             fun = mean,
             colour = "black") +


  # geom_point(shape = 21,
  #            colour = "black") +
  scale_fill_manual(values = c("darkblue",
                               "darkred",
                               "green4")) +
  scale_colour_manual(values = c("darkblue",
                                 "darkred",
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
  ggtitle("CD09")

summary(aov(Value ~ Genotype + variable + Genotype:variable,
            data = df_2_plot))

TukeyHSD(aov(Value ~ Genotype + variable + Genotype:variable,
             data = df_2_plot))


for (cell_time in unique(df_2_plot$variable)) {
  df_2_plot_sub <-
    df_2_plot[df_2_plot$variable == cell_time, ]
  # df_2_plot_sub$Genotype <-
  #   factor(df_2_plot_sub$Genotype,
  #          levels = c("Non-risk",
  #                     "Risk",
  #                     "CRISPRoff"))
  df_2_plot_sub$Genotype <-
    relevel(df_2_plot_sub$Genotype,
            ref = "risk")
  # prin
  lm_model <-
    lmerTest::lmer(Value ~ Genotype + (1|Batch),
                   data = df_2_plot_sub)
  # summary(lm_model)
  lm_output <-
    summary(lm_model)
  #   anova(lm_model)
  print(cell_time)
  print(lm_output$coefficients)
}

# [1] "45 min"
# Estimate Std. Error       df   t value    Pr(>|t|)
# (Intercept)          0.01059811 0.01724882  1.68017 0.6144252 0.611669410
# Genotypenon-risk     0.03958965 0.01435894 14.00000 2.7571423 0.015427071
# Genotyperisk-CRISPRa 0.05404332 0.01435894 14.00000 3.7637389 0.002095685
# [1] "90 min"
# Estimate Std. Error        df  t value     Pr(>|t|)
# (Intercept)          0.04886575 0.03994521  1.122744 1.223319 0.4188345062
# Genotypenon-risk     0.08492588 0.01642586 14.000000 5.170254 0.0001421169
# Genotyperisk-CRISPRa 0.08571552 0.01642586 14.000000 5.218327 0.0001301531
# [1] "135 min"
# Estimate Std. Error        df  t value   Pr(>|t|)
# (Intercept)          0.1078011 0.08089823  1.229728 1.332553 0.37759811
# Genotypenon-risk     0.1528737 0.04400205 14.000000 3.474241 0.00372092
# Genotyperisk-CRISPRa 0.1181202 0.04400205 14.000000 2.684425 0.01779087
# [1] "180 min"
# Estimate Std. Error        df  t value    Pr(>|t|)
# (Intercept)          0.2253050 0.14697945  1.351573 1.532901 0.318430293
# Genotypenon-risk     0.3730340 0.09547784 14.000000 3.907022 0.001579774
# Genotyperisk-CRISPRa 0.2472478 0.09547784 14.000000 2.589583 0.021406825
