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
                     method = 't.test',
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
  df_2_plot[df_2_plot$Line == "CD09", ]

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
                     method = 't.test',
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

## Fig_S10A ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 10)

# test if stat_compare_means can be used with ggplot?
df_2_plot <-
  df_raw

df_2_plot <-
  melt(df_2_plot)
colnames(df_2_plot) <-
  c("Genotype", "Gene", "Value")

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

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


## Fig_S11B ####
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
                     limits = c(0, 20),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

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
                     limits = c(0, 120),
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

## Fig_S12G ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 25)

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
                     limits = c(0, 75),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

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

## Fig_S12I ####
df_raw <-
  read_excel("tables_4_plot.xlsx",
             sheet = 27)

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
                     limits = c(0, 75),
                     na.value = NA) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0)) 

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
                     method = 't.test',
                     paired = F,
                     hide.ns = F,
                     # ref.group = "non-risk",
                     comparisons = list(c("control", "KO"),
                                        c("control", "control_CytoD"),
                                        c("KO", "KO_cytoD"))) +
  labs(x = "",
       y = "Fluor. intensity/cell") +
  scale_x_discrete(drop = F) +
  # ylim(0, 12) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 15000),
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
