# Siwei 16 Jan 2024
# plot fig. 6b use plotly

# init ####
library(ggplot2)
library(RColorBrewer)
library(readxl)

library(stringr)
library(plotly)

# load data ####
df_raw <-
  read_excel("Fig6b_12.xlsx")

# Fig 6b #1 #####
df_to_plot <-
  df_raw[, c(1, 2)]


df_to_plot$Gene <-
  factor(df_to_plot$Gene,
         levels = c("NTC",
                    "CUL1",
                    "SHANK3"))

# "#8DD3C7" "#FFFFB3" "#BEBADA"
plot_ly(data = df_to_plot,
        y = ~ Mean_Neurite_Outgrowth,
        marker = list(colors = c("#8DD3C7", "#FFFFB3", "#BEBADA")),
        fillcolor = ~ Gene,
        type = "box",
        boxpoints = "all",
        jitter = 0.3,
        pointpos = 0,
        size = 1)

df_to_plot <-
  df_raw[, c(1, 3)]
plot_ly(data = df_to_plot,
        y = ~ Branches_Per_Cell,
        marker = list(colors = c("#8DD3C7", "#FFFFB3", "#BEBADA")),
        fillcolor = ~ Gene,
        type = "box",
        boxpoints = "all",
        jitter = 0.3,
        pointpos = 0,
        size = 1)
