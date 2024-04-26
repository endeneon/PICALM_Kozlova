# Siwei 28 Mar 2024
# plot Fig 4E

# init ####
library(readxl)
library(scales)

library(ggplot2)
library(RColorBrewer)

# load data ####
df_raw <-
  read_excel("lysosomal_genes.xlsx",
             sheet = 3)

df_2_plot <-
  df_raw

ggplot(df_2_plot,
       aes(x = risk_vs_nonrisk,
           y = mouse_aging_LDAM)) +
  geom_point(colour = "magenta4") +
  stat_smooth(method = "lm",
              se = F) +
  # stat_smooth(method = "gam",
  #             colour = "darkred") +
  theme_bw()
