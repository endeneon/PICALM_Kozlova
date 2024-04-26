# Siwei 03 Jul 2023
# plot Xiaotong's data in dot plot format

# init
library(ggplot2)
library(ggnewscale)
library(scales)

library(readxl)
library(RColorBrewer)
library(stringr)

# plot the 3 cell types #####
Xiaotong_3_cell_types <-
  read_excel("plot_Xiaotong/Xiaotong_all_cell_types.xlsx")

df_to_plot <- Xiaotong_3_cell_types
df_to_plot$Traits <-
  factor(df_to_plot$Traits,
         levels = c("Schizophrenia",
                    "Bipolar",
                    "Autism",
                    "Depression",
                    "ADHD",
                    "Neuroticism",
                    "Intelligence",
                    "Insomnia",
                    "Alzheimer",
                    "Alzheimer_newer",
                    "Parkinson's disease",
                    "BMI",
                    "T2D"))
df_to_plot$SNP <-
  factor(df_to_plot$SNP,
         levels = c("NGN2",
                    "GABA",
                    "DOPA",
                    "iAstro",
                    "iMG",
                    "hMGeQTL",
                    "hMGcaQTL"))

ggplot(df_to_plot,
       aes(x = Traits,
           y = SNP,
           fill = estimate,
           size = `-log10p`)) +
  geom_point(shape = 21) +
  # scale_fill_gradient(low = "white",
  #                     high = "darkred") +
  # scale_fill_gradient2(low = "darkblue",
  #                      mid = "white",
  #                      high = "darkred",
  #                      # mid = 0,
  #                      breaks = c(min(df_to_plot$estimate),
  #                                 # 0,
  #                                 max(df_to_plot$estimate))) +
  scale_fill_gradientn(colours = c("darkblue",
                                   "white",
                                   "darkred"),
                       values = rescale(c(min(df_to_plot$estimate),
                                          0,
                                          max(df_to_plot$estimate))),
                       limits = c(min(df_to_plot$estimate),
                                  # 0,
                                  max(df_to_plot$estimate))) +
  # scale_fill_distiller(palette = "Spectral",
  #                   type = "seq") +
  # scale_radius() +
  xlab("Traits") +
  # ylab("Cell Type") +
  labs(fill = "Odds ratio",
       size = "-log10p") +
  # scale_y_discrete(limits = rev) +
  # scale_alpha_continuous(df_to_plot$`-logP`,
  #                        range = c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0),
        axis.text = element_text(size = 12)) #+
  ggtitle("Enrichment")

# use new_scale_fill() #####

df_to_plot_est_more_0 <-
  df_to_plot
df_to_plot_est_more_0$estimate[df_to_plot_est_more_0$estimate < 0] <- NA

df_to_plot_est_less_0 <-
  df_to_plot
df_to_plot_est_less_0$estimate[df_to_plot_est_less_0$estimate > 0] <- NA

ggplot() +
  geom_point(data = df_to_plot_est_more_0,
             aes(x = Traits,
                 y = SNP,
                 size = abs(estimate),
                 fill = `-log10p`),
             shape = 21) +
  scale_fill_gradient(low = "white",
                      high = "darkred") +
  scale_radius() +
  new_scale_fill() +
  geom_point(data = df_to_plot_est_less_0,
             aes(x = Traits,
                 y = SNP,
                 size = abs(estimate),
                 fill = `-log10p`),
             shape = 21) +
  scale_fill_gradient(low = "white",
                      high = "steelblue",
                      limits = c(0,
                                 max(df_to_plot_est_less_0$`-log10p`[!is.na(df_to_plot_est_less_0$estimate)]))) +
  scale_radius() +
  xlab("Traits") +
  ylab("") +
  labs(size = "estimate",
       fill = "-log10p") +
  # scale_y_discrete(limits = rev) +
  # scale_alpha_continuous(df_to_plot$`-logP`,
  #                        range = c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0),
        axis.text = element_text(size = 12)) +
  ggtitle("Enrichment")
