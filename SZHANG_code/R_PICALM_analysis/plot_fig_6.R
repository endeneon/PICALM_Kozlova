# Siwei 16 Jan 2024
# plot fig. 6b

# init ####
library(ggplot2)
library(RColorBrewer)
library(readxl)

library(stringr)

# function to plot excel-type quantile #####
xl.type.whisker <-
  function(d) {
    xl_median <-
      stats::quantile(d,
                      probs = c(0.5),
                      type = 6,
                      na.rm = T)
    xl_quartile <-
      stats::quantile(d,
                      probs = c(0.25, 0.75),
                      type = 6,
                      na.rm = T)
    xl_IQR <-
      stats::IQR(d,
                 na.rm = T,
                 type = 6)
    return_data <-
      data.frame(ymin = min(d[d > (xl_quartile[1] - 1.5 * xl_IQR)]),
                 ymax = max(d[d < (xl_quartile[2] + 1.5 * xl_IQR)]),
                 lower = xl_quartile[1],
                 upper = xl_quartile[2],
                 middle = xl_median)
    return(return_data)

  }


# load data ####
df_raw <-
  read_excel("Fig6b_12.xlsx")


# Fig 6b #1 #####
df_to_plot <-
  df_raw[, c(1, 2)]

ggplot(df_to_plot,
       aes(x = factor(Gene,
                      levels = c("NTC",
                                 "CUL1",
                                 "SHANK3")),
           y = Mean_Neurite_Outgrowth,
           fill = factor(Gene,
                         levels = c("NTC",
                                    "CUL1",
                                    "SHANK3")))) +
  stat_summary(fun.data = xl.type.whisker,
               geom = "errorbar",
               size = 0.5,
               width = 0.5) +
  # geom_errorbar(size = 1) +
  stat_summary(fun.data = xl.type.whisker,
               geom = "boxplot",
               width = 0.8) +
  geom_jitter(size = 0.75,
              width = 0.1) +
  # geom_dotplot(position = "identity",
  #              method = "histodot") +
  # stat_summary(fun = mean,
  #              geom = "point",
  #              shape = 4,
  #              colour = "black",
  #              size = 2) +
  # stat_summary(fun.data = mean_se,
  #              geom = "errorbar") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Set3")) +
  ylab("Mean neurite outgrowth/cell") +
  ylim(0, 2000) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.title.y = element_text(size = 10))


# Fig 6b #2 #####
df_to_plot <-
  df_raw[, c(1, 3)]



ggplot(df_to_plot,
       aes(x = factor(Gene,
                      levels = c("NTC",
                                 "CUL1",
                                 "SHANK3")),
           y = Branches_Per_Cell,
           fill = factor(Gene,
                         levels = c("NTC",
                                    "CUL1",
                                    "SHANK3")))) +
  stat_summary(fun.data = xl.type.whisker,
               geom = "errorbar",
               size = 0.5,
               width = 0.5) +
  # geom_errorbar(size = 1) +
  stat_summary(fun.data = xl.type.whisker,
               geom = "boxplot",
               width = 0.8) +
  # geom_boxplot(outlier.shape = NA,
  #              notch = F,
  #              # notchwidth = 0.25,
  #              width = 0.8) +
  # stat_boxplot(geom = "errorbar",
  #              width = 0.25,
  #              type = 6) +
  geom_jitter(size = 0.75,
              width = 0.1) +
  # geom_dotplot(position = "identity",
  #              method = "histodot") +
  stat_summary(fun = mean,
               geom = "point",
               shape = 4,
               colour = "black",
               size = 2) +
  # stat_summary(fun.data = mean_se,
  #              geom = "errorbar") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Set3")) +
  ylab("Mean neurite branches/cell") +
  ylim(0, 250) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.title.y = element_text(size = 10))

median_hilow(df_to_plot$Branches_Per_Cell[df_to_plot$Gene %in% "CUL1"])
boxplot.stats(df_to_plot$Branches_Per_Cell[df_to_plot$Gene %in% "CUL1"])
fivenum(df_to_plot$Branches_Per_Cell[df_to_plot$Gene %in% "CUL1"])

quantile(df_to_plot$Branches_Per_Cell[df_to_plot$Gene %in% "CUL1"],
         # probs = c()
         type = 6)
IQR(df_to_plot$Branches_Per_Cell[df_to_plot$Gene %in% "CUL1"],
    type = 6)

df_raw <-
  read_excel("Fig_6b3.xlsx")
df_to_plot <- df_raw

ggplot(df_to_plot,
       aes(x = factor(Gene,
                      levels = c("NTC",
                                 "CUL1",
                                 "SHANK3")),
           y = Puncta_Density,
           fill = factor(Gene,
                         levels = c("NTC",
                                    "CUL1",
                                    "SHANK3")))) +

  stat_summary(fun.data = xl.type.whisker,
               geom = "errorbar",
               size = 0.5,
               width = 0.5) +
  # geom_errorbar(size = 1) +
  stat_summary(fun.data = xl.type.whisker,
               geom = "boxplot",
               width = 0.8) +
  geom_jitter(size = 0.75,
              width = 0.1) +
  # geom_dotplot(position = "identity",
  #              method = "histodot") +
  # stat_summary(fun = mean,
  #              geom = "point",
  #              shape = 4,
  #              colour = "black",
  #              size = 2) +
  # stat_summary(fun.data = mean_se,
  #              geom = "errorbar") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Set3")) +
  ylab("Puncta density") +
  ylim(0, 0.08) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.title.y = element_text(size = 10))
