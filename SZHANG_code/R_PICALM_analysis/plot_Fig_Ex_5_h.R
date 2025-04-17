# Siwei 27 Jan 2025
# plot Fig. Ex 5h

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

# CD04 ####
df_raw <-
  read_excel("Batch_2_of_new_plots/Fig_Ex_5h.xlsx",
             sheet = 1)

df_sum_all <-
  df_raw
df_sum_all$BR <-
  str_split(string = df_sum_all$Replicates,
            pattern = "_",
            simplify = T)[, 1]

df_sum_all$sum_data_table <-
  str_c(df_sum_all$Genotype,
        df_sum_all$Time,
        df_sum_all$Batch,
        df_sum_all$BR,
        sep = "_")

dt_sum_all <-
  setDT(df_sum_all[, -c(1:5, 8:10)])

dt_sum_all$F_cells <-
  dt_sum_all$`Mean_Abeta-phrodo` / dt_sum_all$`Cell count`
dt_sum_all$`Mean_Abeta-phrodo` <- NULL
dt_sum_all$`Cell count` <- NULL


dt_sum_mean <-
  dt_sum_all[, lapply(.SD, mean), keyby = "sum_data_table"]
dt_sum_mean$normalised <-
  dt_sum_mean$F_cells / max(dt_sum_all$F_cells)

dt_sum_mean$Genotype <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 1]
dt_sum_mean$Time <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 2]


df_2_plot <-
  setDF(dt_sum_mean)
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))

df_2_plot <-
  rbind(df_2_plot,
        c("non-risk",
          0,0,
          "non-risk",
          '0 min'))
df_2_plot <-
  rbind(df_2_plot,
        c("risk",
          0,0,
          "risk",
          '0 min'))

df_2_plot$normalised <-
  as.numeric(df_2_plot$normalised)
sum(is.na(df_2_plot$normalised))

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

df_2_plot$Time <-
  factor(df_2_plot$Time,
         levels = c('0 min',
                    '45 min',
                    '90 min',
                    '135 min',
                    '180 min'))

ggplot(df_2_plot,
       aes(x = Time,
           y = normalised,
           colour = Genotype,
           fill = Genotype,
           group = Genotype)) +

  geom_line(stat = "summary",
            fun = mean,
            linewidth = 1) +

  stat_summary(fun.data = mean_se,
                fun.args = list(mult = 1),
                width = 0.3,
               geom = "errorbar") +
  geom_point(shape = 21,
             stat = "summary",
             fun = mean,
             colour = "black") +


  # geom_point(shape = 21,
  #            colour = "black") +
  scale_fill_manual(values = c("darkred",
                               "darkblue")) +
  scale_colour_manual(values = c("darkred",
                                 "darkblue")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  labs(x = "",
       y = "Ab intensity / cell\nNormalized to non-risk (180 min)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        axis.text.x = element_text(size = 12,
                                   angle = 315,
                                   hjust = 0)) +
  ggtitle("CD04, means used")

for (i in c('45 min',
            '90 min',
            '135 min',
            '180 min')) {
  print(i)
  print(t.test(x = df_2_plot$normalised[(df_2_plot$Time == i) &
                                               (df_2_plot$Genotype == "risk")],
               y = df_2_plot$normalised[(df_2_plot$Time == i) &
                                          (df_2_plot$Genotype == "non-risk")],
               alternative = "t",
               paired = F,
               var.equal = F)$p.value)
}

df_2_plot$batch <-
  str_split(df_2_plot$sum_data_table,
            pattern = "_",
            simplify = T)[, 3]

for (i in c('45 min',
            '90 min',
            '135 min',
            '180 min')) {
  print(i)
  print(t.test(x = df_2_plot$normalised[(df_2_plot$Time == i) &
                                          (df_2_plot$Genotype == "risk")],
               y = df_2_plot$normalised[(df_2_plot$Time == i) &
                                          (df_2_plot$Genotype == "non-risk")],
               alternative = "t",
               paired = F,
               var.equal = F)$p.value)
}

for (i in c('45 min',
            '90 min',
            '135 min',
            '180 min')) {
  print(i)
  df_2_plot_sub <-
    df_2_plot[df_2_plot$Time == i, ]
  # df_2
  print(summary(lmerTest::lmer(normalised ~
                                 Genotype +
                                 (1 | batch),
                               data = df_2_plot_sub, REML = T)))
  # print(t.test(x = df_2_plot$normalised[(df_2_plot$Time == i) &
  #                                         (df_2_plot$Genotype == "risk")],
  #              y = df_2_plot$normalised[(df_2_plot$Time == i) &
  #                                         (df_2_plot$Genotype == "non-risk")],
  #              alternative = "t",
  #              paired = F,
  #              var.equal = F)$p.value)
}


# CD09 ####
df_raw <-
  read_excel("Batch_2_of_new_plots/Fig_Ex_5h.xlsx",
             sheet = 2)

df_sum_all <-
  df_raw
df_sum_all$BR <-
  str_split(string = df_sum_all$Replicates,
            pattern = "_",
            simplify = T)[, 1]

df_sum_all$sum_data_table <-
  str_c(df_sum_all$Genotype,
        df_sum_all$Time,
        df_sum_all$Batch,
        df_sum_all$BR,
        sep = "_")

dt_sum_all <-
  setDT(df_sum_all[, -c(1:5, 8:10)])

dt_sum_all$F_cells <-
  dt_sum_all$`Mean_Abeta-phrodo` / dt_sum_all$`Cell count`
dt_sum_all$`Mean_Abeta-phrodo` <- NULL
dt_sum_all$`Cell count` <- NULL


dt_sum_mean <-
  dt_sum_all[, lapply(.SD, mean), keyby = "sum_data_table"]
dt_sum_mean$normalised <-
  dt_sum_mean$F_cells / max(dt_sum_all$F_cells)

dt_sum_mean$Genotype <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 1]
dt_sum_mean$Time <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 2]


df_2_plot <-
  setDF(dt_sum_mean)
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))

df_2_plot <-
  rbind(df_2_plot,
        c("non-risk",
          0,0,
          "non-risk",
          '0 min'))
df_2_plot <-
  rbind(df_2_plot,
        c("risk",
          0,0,
          "risk",
          '0 min'))

df_2_plot$normalised <-
  as.numeric(df_2_plot$normalised)
sum(is.na(df_2_plot$normalised))

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

df_2_plot$Time <-
  factor(df_2_plot$Time,
         levels = c('0 min',
                    '45 min',
                    '90 min',
                    '135 min',
                    '180 min'))

ggplot(df_2_plot,
       aes(x = Time,
           y = normalised,
           colour = Genotype,
           fill = Genotype,
           group = Genotype)) +

  geom_line(stat = "summary",
            fun = mean,
            linewidth = 1) +

  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               width = 0.3,
               geom = "errorbar") +
  geom_point(shape = 21,
             stat = "summary",
             fun = mean,
             colour = "black") +


  # geom_point(shape = 21,
  #            colour = "black") +
  scale_fill_manual(values = c("darkred",
                               "darkblue")) +
  scale_colour_manual(values = c("darkred",
                                 "darkblue")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  labs(x = "",
       y = "Ab intensity / cell\nNormalized to non-risk (180 min)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        axis.text.x = element_text(size = 12,
                                   angle = 315,
                                   hjust = 0)) +
  ggtitle("CD09, means used")

for (i in c('45 min',
            '90 min',
            '135 min',
            '180 min')) {
  print(i)
  print(t.test(x = df_2_plot$normalised[(df_2_plot$Time == i) &
                                          (df_2_plot$Genotype == "risk")],
               y = df_2_plot$normalised[(df_2_plot$Time == i) &
                                          (df_2_plot$Genotype == "non-risk")],
               alternative = "t",
               paired = F,
               var.equal = F)$p.value)
}


# without taking means #####
## CD04 ####
df_raw <-
  read_excel("Batch_2_of_new_plots/Fig_Ex_5h.xlsx",
             sheet = 1)

df_sum_all <-
  df_raw
df_sum_all$BR <-
  str_split(string = df_sum_all$Replicates,
            pattern = "_",
            simplify = T)[, 1]

df_sum_all$sum_data_table <-
  str_c(df_sum_all$Genotype,
        df_sum_all$Time,
        df_sum_all$Batch,
        df_sum_all$BR,
        sep = "_")

dt_sum_all <-
  setDT(df_sum_all[, -c(1:5, 8:10)])

dt_sum_all$F_cells <-
  dt_sum_all$`Mean_Abeta-phrodo` / dt_sum_all$`Cell count`
dt_sum_all$`Mean_Abeta-phrodo` <- NULL
dt_sum_all$`Cell count` <- NULL


# dt_sum_mean <-
#   dt_sum_all[, lapply(.SD, mean), keyby = "sum_data_table"]
dt_sum_mean <-
  dt_sum_all
dt_sum_mean$normalised <-
  dt_sum_mean$F_cells / max(dt_sum_all$F_cells)
# dt_sum_mean$F_cells <- NULL

dt_sum_mean$Genotype <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 1]
dt_sum_mean$Time <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 2]


df_2_plot <-
  setDF(dt_sum_mean)
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))

df_2_plot <-
  rbind(df_2_plot,
        c("non-risk",
          0,0,
          "non-risk",
          '0 min'))
df_2_plot <-
  rbind(df_2_plot,
        c("risk",
          0,0,
          "risk",
          '0 min'))

df_2_plot$normalised <-
  as.numeric(df_2_plot$normalised)
sum(is.na(df_2_plot$normalised))

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

df_2_plot$Time <-
  factor(df_2_plot$Time,
         levels = c('0 min',
                    '45 min',
                    '90 min',
                    '135 min',
                    '180 min'))

ggplot(df_2_plot,
       aes(x = Time,
           y = normalised,
           colour = Genotype,
           fill = Genotype,
           group = Genotype)) +

  geom_line(stat = "summary",
            fun = mean,
            linewidth = 1) +

  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               width = 0.3,
               geom = "errorbar") +
  geom_point(shape = 21,
             stat = "summary",
             fun = mean,
             colour = "black") +


  # geom_point(shape = 21,
  #            colour = "black") +
  scale_fill_manual(values = c("darkred",
                               "darkblue")) +
  scale_colour_manual(values = c("darkred",
                                 "darkblue")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  labs(x = "",
       y = "Ab intensity / cell\nNormalized to non-risk (180 min)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        axis.text.x = element_text(size = 12,
                                   angle = 315,
                                   hjust = 0)) +
  ggtitle("CD04, means not used")

for (i in c('45 min',
            '90 min',
            '135 min',
            '180 min')) {
  print(i)
  print(t.test(x = df_2_plot$normalised[(df_2_plot$Time == i) &
                                          (df_2_plot$Genotype == "risk")],
               y = df_2_plot$normalised[(df_2_plot$Time == i) &
                                          (df_2_plot$Genotype == "non-risk")],
               alternative = "t",
               paired = F,
               var.equal = F)$p.value)
}
#
# [1] "45 min"
# [1] 0.02266595
# [1] "90 min"
# [1] 0.0156134
# [1] "135 min"
# [1] 0.0231373
# [1] "180 min"
# [1] 0.004247324

## CD09 ####
df_raw <-
  read_excel("Batch_2_of_new_plots/Fig_Ex_5h.xlsx",
             sheet = 2)

df_sum_all <-
  df_raw
df_sum_all$BR <-
  str_split(string = df_sum_all$Replicates,
            pattern = "_",
            simplify = T)[, 1]

df_sum_all$sum_data_table <-
  str_c(df_sum_all$Genotype,
        df_sum_all$Time,
        df_sum_all$Batch,
        df_sum_all$BR,
        sep = "_")

dt_sum_all <-
  setDT(df_sum_all[, -c(1:5, 8:10)])

dt_sum_all$F_cells <-
  dt_sum_all$`Mean_Abeta-phrodo` / dt_sum_all$`Cell count`
dt_sum_all$`Mean_Abeta-phrodo` <- NULL
dt_sum_all$`Cell count` <- NULL


# dt_sum_mean <-
#   dt_sum_all[, lapply(.SD, mean), keyby = "sum_data_table"]
dt_sum_mean <-
  dt_sum_all
dt_sum_mean$normalised <-
  dt_sum_mean$F_cells / max(dt_sum_all$F_cells)

dt_sum_mean$Genotype <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 1]
dt_sum_mean$Time <-
  str_split(dt_sum_mean$sum_data_table,
            pattern = "_",
            simplify = T)[, 2]


df_2_plot <-
  setDF(dt_sum_mean)
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))

df_2_plot <-
  rbind(df_2_plot,
        c("non-risk",
          0,0,
          "non-risk",
          '0 min'))
df_2_plot <-
  rbind(df_2_plot,
        c("risk",
          0,0,
          "risk",
          '0 min'))

df_2_plot$normalised <-
  as.numeric(df_2_plot$normalised)
sum(is.na(df_2_plot$normalised))

df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

df_2_plot$Time <-
  factor(df_2_plot$Time,
         levels = c('0 min',
                    '45 min',
                    '90 min',
                    '135 min',
                    '180 min'))

ggplot(df_2_plot,
       aes(x = Time,
           y = normalised,
           colour = Genotype,
           fill = Genotype,
           group = Genotype)) +

  geom_line(stat = "summary",
            fun = mean,
            linewidth = 1) +

  stat_summary(fun.data = mean_se,
               fun.args = list(mult = 1),
               width = 0.3,
               geom = "errorbar") +
  geom_point(shape = 21,
             stat = "summary",
             fun = mean,
             colour = "black") +


  # geom_point(shape = 21,
  #            colour = "black") +
  scale_fill_manual(values = c("darkred",
                               "darkblue")) +
  scale_colour_manual(values = c("darkred",
                                 "darkblue")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.1)) +
  labs(x = "",
       y = "Ab intensity / cell\nNormalized to non-risk (180 min)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        axis.text.x = element_text(size = 12,
                                   angle = 315,
                                   hjust = 0)) +
  ggtitle("CD09, means not used")


df_2_plot$batch <-
  str_split(df_2_plot$sum_data_table,
            pattern = "_",
            simplify = T)[, 3]


for (i in c('45 min',
            '90 min',
            '135 min',
            '180 min')) {
  print(i)
  df_2_plot_sub <-
    df_2_plot[df_2_plot$Time == i, ]
  # df_2
  print(summary(lmerTest::lmer(normalised ~
                                 Genotype +
                                 (1 | batch),
                               data = df_2_plot_sub, REML = T)))
  # print(t.test(x = df_2_plot$normalised[(df_2_plot$Time == i) &
  #                                         (df_2_plot$Genotype == "risk")],
  #              y = df_2_plot$normalised[(df_2_plot$Time == i) &
  #                                         (df_2_plot$Genotype == "non-risk")],
  #              alternative = "t",
  #              paired = F,
  #              var.equal = F)$p.value)
}

for (i in c('45 min',
            '90 min',
            '135 min',
            '180 min')) {
  print(i)
  print(t.test(x = df_2_plot$normalised[(df_2_plot$Time == i) &
                                          (df_2_plot$Genotype == "risk")],
               y = df_2_plot$normalised[(df_2_plot$Time == i) &
                                          (df_2_plot$Genotype == "non-risk")],
               alternative = "t",
               paired = F,
               var.equal = F)$p.value)
}
#
# [1] "45 min"
# [1] 0.000132693
# [1] "90 min"
# [1] 0.002017698
# [1] "135 min"
# [1] 0.001094981
# [1] "180 min"
# [1] 0.0001623474
