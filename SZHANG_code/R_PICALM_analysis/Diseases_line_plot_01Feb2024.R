# Siwei 01 Feb 2024
# Plot Jubao's Diseases line graph

# init ####
library(readxl)

library(reshape2)

library(ggplot2)
library(scales)

library(plyr)

library(RColorBrewer)

# load data #####
df_raw <-
  read_excel("Jubao_Diseases_plot_line_graph.xlsx")

df_melt <-
  reshape2::melt(data = df_raw,
                 id = c(1, 11))
colnames(df_melt)[3] <- "Dynamic_Pattern"
colnames(df_melt)[4] <- "Fold_Change"

unique(df_melt$Dynamic_Pattern)
df_melt$Dynamic_Pattern <-
  factor(df_melt$Dynamic_Pattern,
         levels = c("up-up",
                    "up-flat",
                    "up-down",
                    "flat-up",
                    "down-up",
                    "flat-down",
                    "down-down",
                    "down-flat",
                    "unchanged"))

df_melt$Diseases <-
  factor(df_melt$Diseases,
         levels = unique(df_melt$Diseases))

ggplot(df_melt,
       aes(x = factor(Dynamic_Pattern),
           y = Fold_Change,
           # group = Diseases,
           colour = factor(Diseases))) +
  # geom_line() +
  geom_path(aes(group = Diseases)) +
  scale_colour_manual(values = brewer.pal(n = 8,
                                        name = "Dark2"),
                    name = "Diseases") +
  facet_wrap(~ Cell_type,
             nrow = 1) +
  labs(y = "Fold Change",
       x = "Dynamic Patterns") +
  # guides(fill = "Dynamic patterns") +
  # scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))

### make cumulative plot
df_cum <-
  df_raw
df_cum <-
  split(x = df_cum,
        f = df_cum$Cell_type)

for (i in 1:length(df_cum)) {
  print(i)
  df_to_operate <-
    df_cum[[i]]
  print(df_to_operate$Diseases)
  df_to_operate <-
    df_to_operate[rev(1:nrow(df_to_operate)), ] # reverse rows
  Diseases_vector <-
    df_to_operate$Diseases
  cell_type <-
    unique(df_to_operate$Cell_type)

  df_to_operate$Diseases <- NULL
  df_to_operate$Cell_type <- NULL

  # Then add rows in a rolling pattern
  for (j in 1:(nrow(df_to_operate) - 1)) {
    df_to_operate[j + 1, ] <-
      df_to_operate[j, ] +
      df_to_operate[j + 1, ]
  }

  # add back the two columns
  df_to_operate$Diseases <-
    Diseases_vector
    # rownames(df_to_operate)
  df_to_operate$Cell_type <-
    cell_type

  df_cum[[i]] <-
    df_to_operate
}

# ! use l_ply from plyr !
df_cum_2_plot <-
  ldply(df_cum,
        data.frame)
colnames(df_cum_2_plot)

df_cum_2_plot$.id <- NULL
# df_cum_2_plot$unchanged <- NULL
# df_cum_2_plot$unchanged <- 1

df_cum_2_plot_melt <-
  reshape2::melt(data = df_cum_2_plot,
                 id = c(10, 11))
colnames(df_cum_2_plot_melt)[3] <- "Dynamic_Pattern"
colnames(df_cum_2_plot_melt)[4] <- "Fold_Change"


df_cum_2_plot_melt$Dynamic_Pattern <-
  factor(df_cum_2_plot_melt$Dynamic_Pattern,
         levels = c("up.up",
                    "up.flat",
                    "up.down",
                    "flat.up",
                    "down.up",
                    "flat.down",
                    "down.down",
                    "down.flat",
                    "unchanged"))

df_cum_2_plot_melt$Diseases <-
  factor(df_cum_2_plot_melt$Diseases,
         levels = rev(unique(df_cum_2_plot_melt$Diseases)))

ggplot(df_cum_2_plot_melt,
       aes(x = factor(Dynamic_Pattern),
           y = Fold_Change,
           # group = Diseases,
           colour = factor(Diseases))) +
  # geom_line() +
  geom_path(aes(group = Diseases)) +
  scale_colour_manual(values = brewer.pal(n = 8,
                                          name = "Dark2"),
                      name = "Diseases") +
  facet_wrap(~ Cell_type,
             nrow = 1) +
  labs(y = "Cumulative Fold Change",
       x = "Dynamic Patterns") +
  # guides(fill = "Dynamic patterns") +
  # scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))

