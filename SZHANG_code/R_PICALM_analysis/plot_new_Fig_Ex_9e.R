# Siwei 02 Feb 2025
# plot new Ex 9e

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

# load raw data ####
set.seed(42)
rm(list = ls())


## Fig. Ex 9e #####
df_raw <-
  read_excel("Batch_2_of_new_plots/new_Fig_Ex_9e.xlsx",
             sheet = 1)

df_sum_all <-
  df_raw
# df_sum_all$BR <-
#   str_split(string = df_sum_all$Replicates,
#             pattern = "_",
#             simplify = T)[, 1]
df_sum_all$Replicates <- NULL

dt_sum_all <-
  setDT(df_sum_all)
# dt_sum_mean <-
#   dt_sum_all[, lapply(.SD, mean), keyby = "BR"]

df_2_plot <-
  reshape2::melt(setDF(dt_sum_all),
                 value.name = "Value")

colnames(df_2_plot)[1] <- "Genotype"

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "risk TrC",
                    "non-risk",
                    "non-risk TrC"
         ))


ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Value",
            color = "Genotype",
            # shape = "clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue",
                                 "green4",
                                 "darkorange")) +

  # stat_compare_means(label = "p.signif",
  #                    # label.y.npc = 0.95,
  #                    # label.x.npc = 0,
  #                    method = 'wilcox.test',
  #                    hide.ns = T,
  #                    ref.group = "non-risk",
  #                    paired = F) +
  geom_jitter(aes(colour = Genotype),
              width = 0.1,
              size = 1,
              shape = 1) +
  # scale_shape_manual(values = c(1)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "BD.11 Ox/BD.11 non-ox") +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                7)) +
  # guides(guide_col)
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1))


df_2_plot$Batch <-
  as.factor(rep_len(x = c(1,1,1,2,2,2,3,3,3),
                    length.out = nrow(df_2_plot)))

df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "non-risk")
df_2_plot$Genotype <-
  relevel(df_2_plot$Genotype,
          ref = "risk")

lm_model <-
  lmerTest::lmer(Value ~
                   Genotype +
                   (1|Batch),
                 data = df_2_plot)
summary(lm_model)


summary(aov(Value ~ Genotype,
            data = df_2_plot))
summary(multcomp::glht(aov(Value ~ Genotype,
                           data = df_2_plot),
                       linfct = mcp(Genotype = "Dunnett")))
