# Siwei 09 Oct 2023
# Make bubble plot of Alena's PICALM MAGMA result

# init #####
{
  library(readr)
  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
  library(scales)
}

# load data #####

df_raw <-
  read_table("assembled_MAGMA_results.txt")

df_to_plot <-
  df_raw
df_to_plot$`-logP` <-
  0 - log10(df_to_plot$P)
df_to_plot$GWAS_SET <-
  str_split(string = df_to_plot$GWAS_SET,
            pattern = "\\.",
            simplify = T)[, 1]

df_to_plot <-
  df_to_plot[!(df_to_plot$GWAS_SET %in% "IBD_2015"), ]

ggplot(data = df_to_plot,
       aes(x = factor(GWAS_SET),
           y = factor(Category,
                      levels = c("up_regulated",
                                 "down_regulated")),
           size = abs(BETA),
           # alpha = 0 - log10(FDR),
           fill = `-logP`)) +
  geom_point(shape = 21) +
  scale_size_area() +
  scale_fill_gradientn(colours = c("darkblue",
                                   "white",
                                   "darkred"),
                       values = rescale(c(min(df_to_plot$`-logP`),
                                          0,
                                          max(df_to_plot$`-logP`))),
                       limits = c(min(df_to_plot$`-logP`),
                                  max(df_to_plot$`-logP`)
                       )) +
labs(fill = "-logP",
     size = "Beta",
     # fill = "Direction",
     x = "GWAS set",
     y = "Direction") +
  guides() +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")

ggplot(data = df_to_plot,
       aes(x = factor(GWAS_SET),
           y = factor(Category,
                      levels = c("up_regulated",
                                 "down_regulated")),
           size = abs(BETA),
           # alpha = 0 - log10(FDR),
           fill = `-logP`)) +
  geom_point(shape = 21) +
  scale_size_area() +
  scale_fill_gradientn(colours = c("white",
                                   "darkred"),
                       # values = rescale(c(min(df_to_plot$`-logP`),
                       #                    0,
                       #                    max(df_to_plot$`-logP`))),
                       limits = c(0,
                                  max(df_to_plot$`-logP`)
                       )) +
  labs(fill = "-logP",
       size = "Beta",
       # fill = "Direction",
       x = "GWAS set",
       y = "Direction") +
  guides() +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")



df_to_plot_alt <-
  df_to_plot
df_to_plot_alt$`-logP`[df_to_plot_alt$BETA < 0] <-
  0 - df_to_plot_alt$`-logP`[df_to_plot_alt$BETA < 0]

ggplot(data = df_to_plot_alt,
       aes(x = factor(GWAS_SET,
                      levels = c("Alz_Jansenetal_2019", "Alz_Jansenetal_2021",
                                 "SCZw3_2022", "UKBMDD_2019", "Bipolar_2019",
                                 "ASD_2019", "Neurotics_2018", "PGCAlcohol_2018",
                                 "Crohn_2018")),
           y = factor(Category,
                      levels = c("up_regulated",
                                 "down_regulated")),
           size = abs(BETA),
           # alpha = 0 - log10(FDR),
           fill = `-logP`)) +
  geom_point(shape = 21) +
  scale_size_area() +
  scale_fill_gradientn(colours = c("darkblue",
                                   "white",
                                   "darkred"),
                       values = rescale(c(min(df_to_plot_alt$`-logP`),
                                          0,
                                          max(df_to_plot_alt$`-logP`))),
                       limits = c(min(df_to_plot_alt$`-logP`),
                                  max(df_to_plot_alt$`-logP`)
                       )) +
  labs(fill = "-logP",
       size = "Beta",
       # fill = "Direction",
       x = "GWAS set",
       y = "Direction") +
  guides(fill = "none") +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right") +
  ggtitle("PICALM risk vs non-risk,\nDE genes with FDR < 0.05")


ggplot(data = df_to_plot_alt,
       aes(x = factor(GWAS_SET,
                      levels = c("Alz_Jansenetal_2019", "Alz_Jansenetal_2021",
                                 "SCZw3_2022", "UKBMDD_2019", "Bipolar_2019",
                                 "ASD_2019", "Neurotics_2018", "PGCAlcohol_2018",
                                 "Crohn_2018")),
           y = factor(Category,
                      levels = c("up_regulated",
                                 "down_regulated")),
           size = abs(BETA),
           # alpha = 0 - log10(FDR),
           fill = `-logP`)) +
  geom_point(shape = 21) +
  scale_size_area() +
  scale_fill_gradientn(colours = c("white",
                                   "darkred"),
                       # values = rescale(c(min(df_to_plot_alt$`-logP`),
                       #                    0,
                       #                    max(df_to_plot_alt$`-logP`))),
                       limits = c(0,
                                  max(df_to_plot_alt$`-logP`)
                       )) +
  labs(fill = "-logP",
       size = "Beta",
       # fill = "Direction",
       x = "GWAS set",
       y = "Direction") +
  # guides(fill = "none") +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")


ggplot(data = df_to_plot_alt,
       aes(x = factor(GWAS_SET),
           y = factor(Category,
                      levels = c("up_regulated",
                                 "down_regulated")),
           size = abs(BETA),
           # alpha = 0 - log10(FDR),
           fill = `-logP`)) +
  geom_point(shape = 21) +
  scale_size_area() +
  scale_fill_gradientn(colours = c("white", "darkblue"),
                       # values = rescale(c(min(df_to_plot_alt$`-logP`),
                       #                    0,
                       #                    max(df_to_plot_alt$`-logP`))),
                       limits = c(0,
                                  0 - min(df_to_plot_alt$`-logP`)
                       )) +
  labs(fill = "-logP",
       size = "Beta",
       # fill = "Direction",
       x = "GWAS set",
       y = "Direction") +
  guides() +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")


## Use size as -logP #####
df_to_plot_alt <-
  df_to_plot
df_to_plot_alt <-
  df_to_plot_alt[!df_to_plot_alt$GWAS_SET %in% "Alz_Jansenetal_2021", ]

ggplot(data = df_to_plot_alt,
       aes(x = factor(GWAS_SET,
                      levels = c("Alz_Jansenetal_2019",
                                 # "Alz_Jansenetal_2021",
                                 "SCZw3_2022", "UKBMDD_2019", "Bipolar_2019",
                                 "ASD_2019", "Neurotics_2018", "PGCAlcohol_2018",
                                 "Crohn_2018")),
           y = factor(Category,
                      levels = c("up_regulated",
                                 "down_regulated")),
           fill = BETA,
           # alpha = 0 - log10(FDR),
           size = `-logP`)) +
  geom_point(shape = 21) +
  scale_radius() +
  scale_fill_gradientn(colours = c("darkblue",
                                   "white",
                                   "darkred"),
                       values = rescale(c(min(df_to_plot_alt$BETA),
                                          0,
                                          max(df_to_plot_alt$BETA)))
                       # ,
                       # limits = c(min(df_to_plot_alt$`-logP`),
                       #            max(df_to_plot_alt$`-logP`)
                       ) +
  labs(fill = "Beta",
       size = "-log10P",
       # fill = "Direction",
       x = "GWAS set",
       y = "Direction") +
  # guides(fill = "none") +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.2,1,0,0),
                           "cm"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "left") +
  ggtitle("PICALM risk vs non-risk,\nDE genes with FDR < 0.05")


## make 3-group plot #####
df_raw <-
  read_table("Alena_PICALM_groups_3x_4_R.txt")

df_to_plot <-
  df_raw
df_to_plot$`-logP` <-
  0 - log10(df_to_plot$P)
df_to_plot$GWAS_SET <-
  str_split(string = df_to_plot$GWAS_SET,
            pattern = "\\.",
            simplify = T)[, 1]

## Also use size as -logP #####
df_to_plot_alt <-
  df_to_plot

df_to_plot_alt <-
  df_to_plot_alt[df_to_plot_alt$GWAS_SET %in% "Alz_Jansenetal_2021", ]

ggplot(data = df_to_plot_alt,
       aes(x = factor(GWAS_SET,
                      levels = c("Alz_Jansenetal_2019", "Alz_Jansenetal_2021",
                                 "SCZw3_2022", "UKBMDD_2019", "Bipolar_2019",
                                 "ASD_2019", "Neurotics_2018", "PGCAlcohol_2018",
                                 "Crohn_2018")),
           y = factor(Category,
                      levels = c("up-regulated",
                                 "non-significant",
                                 "down-regulated")),
           fill = BETA,
           # alpha = 0 - log10(FDR),
           size = `-logP`)) +
  geom_point(shape = 21) +
  scale_radius() +
  scale_fill_gradientn(colours = c("darkblue",
                                   "white",
                                   "darkred"),
                       values = rescale(c(min(df_to_plot_alt$BETA),
                                          0,
                                          max(df_to_plot_alt$BETA)))
                       # ,
                       # limits = c(min(df_to_plot_alt$`-logP`),
                       #            max(df_to_plot_alt$`-logP`)
  ) +
  labs(fill = "Beta",
       size = "-log10P",
       # fill = "Direction",
       x = "GWAS set",
       y = "Direction") +
  # guides(fill = "none") +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.2,1,0,0),
                           "cm"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "left") +
  ggtitle("PICALM risk vs non-risk,\nDE genes with FDR < 0.05")

