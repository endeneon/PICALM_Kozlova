# Siwei 01 Feb 2024
# Plot Jubao's Diseases bar plot

# init ####
library(readxl)

library(reshape2)

library(ggplot2)
library(scales)

library(RColorBrewer)

# load data #####
df_raw <-
  read_excel("Jubao_Diseases_bar_plot_v2.xlsx")

df_melt <-
  reshape2::melt(data = df_raw,
                 id = c(1, 10))
colnames(df_melt)[3] <- "Diseases"
# colnames(df_melt)[4] <- ""

df_melt$Dynamic_Pattern <-
  factor(df_melt$Dynamic_Pattern,
         levels = df_melt$Dynamic_Pattern[1:9])
df_melt$Diseases <-
  factor(df_melt$Diseases,
         levels = unique(df_melt$Diseases))

ggplot(df_melt,
        aes(x = Diseases,
            y = value,
            # group = Diseases,
            fill = Dynamic_Pattern)) +
  geom_bar(position = "fill",
           stat = "identity") +
  scale_fill_manual(values = brewer.pal(n = 9,
                                        name = "Set1"),
                    name = "Dynamic patterns") +
  facet_wrap(~ Cell_Type,
             nrow = 1) +
  labs(y = "Percentage") +
  # guides(fill = "Dynamic patterns") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))
