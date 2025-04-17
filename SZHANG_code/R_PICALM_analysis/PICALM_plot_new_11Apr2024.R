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
  library(lmtest)
}
# library(seqLogo)

## Fig_Ex_5F ####
df_raw <-
  read_excel("tables_4_plot_new.xlsx",
             sheet = 1)


df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
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
  df_2_plot[!(df_2_plot$variable == '0 min'), ]

for (cell_time in unique(df_2_plot$variable)) {
  df_2_plot_sub <-
    df_2_plot[df_2_plot$variable == cell_time, ]
  print(cell_time)
  # df_2_plot_sub$Genotype <-
  #   factor(df_2_plot_sub$Genotype,
  #          levels = c("Non-risk",
  #                     "Risk",
  #                     "CRISPRoff"))
  df_2_plot_sub$Genotype <-
    relevel(df_2_plot_sub$Genotype,
            ref = "risk")
  df_2_plot_sub$Batch <-
    rep_len(x = c(1,1,2,2,3,3),
            length.out = nrow(df_2_plot_sub))
  # prin
  lm_model <-
    lmerTest::lmer(Value ~ Genotype + 
                     (1|Batch),
                   data = df_2_plot_sub)
  # summary(lm_model)
  lm_output <-
    summary(lm_model)
  print(lm_output)
  #   anova(lm_model)
  # print(cell_time)
  # print(lm_output$coefficients)
}


## CD09 #####

df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
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
  df_2_plot[!(df_2_plot$variable == '0 min'), ]

for (cell_time in unique(df_2_plot$variable)) {
  df_2_plot_sub <-
    df_2_plot[df_2_plot$variable == cell_time, ]
  print(cell_time)
  # df_2_plot_sub$Genotype <-
  #   factor(df_2_plot_sub$Genotype,
  #          levels = c("Non-risk",
  #                     "Risk",
  #                     "CRISPRoff"))
  df_2_plot_sub$Genotype <-
    relevel(df_2_plot_sub$Genotype,
            ref = "risk")
  df_2_plot_sub$Batch <-
    rep_len(x = c(1,1,2,2,3,3),
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


## Fig_E5h ####
df_raw <-
  read_excel("tables_4_plot_new.xlsx",
             sheet = 2)


df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
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
  df_2_plot[!(df_2_plot$variable == '0 min'), ]

for (cell_time in unique(df_2_plot$variable)) {
  df_2_plot_sub <-
    df_2_plot[df_2_plot$variable == cell_time, ]
  print(cell_time)
  # df_2_plot_sub$Genotype <-
  #   factor(df_2_plot_sub$Genotype,
  #          levels = c("Non-risk",
  #                     "Risk",
  #                     "CRISPRoff"))
  df_2_plot_sub$Genotype <-
    relevel(df_2_plot_sub$Genotype,
            ref = "risk")
  df_2_plot_sub$Batch <-
    rep_len(x = c(1,1,2,2,3,3),
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


df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
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
  df_2_plot[!(df_2_plot$variable == '0 min'), ]

for (cell_time in unique(df_2_plot$variable)) {
  df_2_plot_sub <-
    df_2_plot[df_2_plot$variable == cell_time, ]
  print(cell_time)
  # df_2_plot_sub$Genotype <-
  #   factor(df_2_plot_sub$Genotype,
  #          levels = c("Non-risk",
  #                     "Risk",
  #                     "CRISPRoff"))
  df_2_plot_sub$Genotype <-
    relevel(df_2_plot_sub$Genotype,
            ref = "risk")
  df_2_plot_sub$Batch <-
    rep_len(x = c(1,1,2,2,3,3),
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

## Fig.6C ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 9)

df_2_plot <-
  df_raw

df_2_plot <-
  df_2_plot[df_2_plot$Line == "CD04", ]

df_2_plot <-
  df_2_plot[!(df_2_plot$Genotype == "risk"), ]

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "CRISPRoff"))

df_2_plot$Avg_exp_normalised <-
  df_2_plot$Avg_exp / mean(df_2_plot$Avg_exp[df_2_plot$Genotype == "non-risk"])

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
                     limits = c(0, 2),
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
  df_2_plot[!(df_2_plot$Genotype == "risk"), ]

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("non-risk",
                    "CRISPRoff"))

# df_2_plot <-
#   df_2_plot[!(df_2_plot$Genotype == "CRISPRoff"), ]

df_2_plot$Avg_exp_normalised <-
  df_2_plot$Avg_exp / mean(df_2_plot$Avg_exp[df_2_plot$Genotype == "non-risk"])

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

## Fig_6B ####
df_raw <-
  read_excel("tables_4_plot_new.xlsx",
             sheet = 3)


df_2_plot <-
  df_raw

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk",
                    "risk-TrC"))


df_2_plot$variable <-
  factor(df_2_plot$variable,
         levels = c("0 min",
                    "45 min",
                    "90 min",
                    "135 min",
                    "180 min"))

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
         levels = c("risk",
                    "non-risk",
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

## Fig_6C ####
df_raw <-
  read_excel("tables_4_plot_new.xlsx",
             sheet = 4)


df_2_plot <-
  df_raw

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value", 
       na.rm = T)

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
       y = "BODIPY fluor. intensity/iMG") +
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


## Fig_S14C ####
df_raw <-
  read_excel("tables_4_plot_new.xlsx",
             sheet = 6)


df_2_plot <-
  df_raw

df_2_plot <-
  melt(df_2_plot,
       value.name = "Value", 
       na.rm = T)

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
       y = "BODIPY fluor. intensity/iMG") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 2.5),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

df_2_plot$BR <-
  c(1,1,2,2,2,2,
     1,1,2,2,2,2,
     1,1,2,2)
lm_model <-
  lmerTest::lmer(Value ~ 
                   variable + 
                   (1|BR),
                 data = df_2_plot)

summary(lm_model)
anova(lm_model)
