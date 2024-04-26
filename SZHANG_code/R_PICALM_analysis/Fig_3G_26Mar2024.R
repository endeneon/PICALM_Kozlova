# Siwei 26 Mar 2024
# plot Fig2G

# init ####
library(readxl)
library(stringr)
library(RColorBrewer)
library(ggplot2)

library(scales)

# load raw data ####
df_raw <-
  read_excel("Fig_2G.xlsx")

df_2_plot <-
  df_raw
df_2_plot$Age <-
  trunc(df_2_plot$Age_raw / 10) * 10
df_2_plot$Diagnosis <-
  factor(df_2_plot$Diagnosis,
         levels = c("Control", "AD"))

df_2_plot <-
  df_2_plot[, c(2, 4:8)]
df_2_plot$Age <-
  as.factor(df_2_plot$Age)
df_2_plot$group_info <-
  "group"

ggplot(df_2_plot, 
        aes(x = Diagnosis,
            y = ddCt_PICALM_exon1,
            size = Age,
            shape = Repository,
            colour = Sex,
            fill = Sex)) +
  geom_jitter(width = 0.1,
              # shape = c(21, 24),
              colour = "black") +
  stat_summary(aes(group = Diagnosis),
               geom = "crossbar",
               fun = "mean",
               width = 0.3,
               linewidth = 0.5,
               colour = "black") +
  geom_errorbar(aes(group = Diagnosis),
                stat = "summary",
                fun.data = mean_cl_boot,
                width = 0.1,
                colour = "black",
                size = 1) +

  scale_shape_manual(values = c(21, 24)) +
  scale_size_manual(values = c(2, 2.5 , 3, 3.5)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 2.5)) +
  scale_colour_manual(values = brewer.pal(n = 3,
                                          name = "Set1")) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                          name = "Set1")) +
  labs(x = "",
       y = "Normalized exp. level") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10))


# plot Fig S6C ####
df_raw <-
  read_excel("Fig_S6C.xlsx")

df_2_plot <-
  df_raw

ggplot(df_2_plot, 
       aes(x = cell.type,
           y = `Log2(level)`,
           # group = genotype,
           colour = genotype)) +
  geom_point(# shape = c(21, 24),
              # colour = "black",
              position = position_dodge2(width = 1),
              inherit.aes = T) +
  # geom_dotplot(binaxis = "y",
  #              stackdir = "center",
  #              binpositions = "all",
  #              position = position_dodge2(width = 1),
  #              # colour = "black",
  #              # position = position_dodge2(width = 0.5),
  #              inherit.aes = T) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.3,
               linewidth = 0.5,
               colour = "black",
               position = position_dodge2(width = 0.5)) +
  geom_errorbar(stat = "summary",
                fun.data = mean_cl_boot,
                width = 0.1,
                colour = "black",
                size = 1,
                position = position_dodge2(width = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = brewer.pal(n = 3,
                                          name = "Set1")) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Set1")) +
  labs(x = "",
       y = "Normalized exp. level") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10))

df_2_plot$genotype <-
  factor(df_2_plot$genotype,
         levels = unique(df_2_plot$genotype))

ggplot(df_2_plot,
       aes(x = genotype,
           y = `Log2(level)`,
           group = genotype,
           colour = genotype)) +
  geom_point(position = position_dodge2(width = 0.1),
             inherit.aes = T) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.3,
               linewidth = 0.2,
               colour = "black",
               position = position_dodge2(width = 0.2)) +
  geom_errorbar(stat = "summary",
                fun.data = mean_cl_boot,
                width = 0.1,
                colour = "black",
                size = 0.5,
                position = position_dodge2(width = 0.5)) +
  scale_colour_manual(values = brewer.pal(n = 3,
                                          name = "Set1")) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Set1")) +
  labs(x = "",
       y = "Normalized exp. level") +
  ylim(12, 15) +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  facet_grid(~ cell.type)

