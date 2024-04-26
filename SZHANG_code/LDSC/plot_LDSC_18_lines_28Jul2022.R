# make disease/gene group plots for LDSC enrichment on scARC
# peaks enriched from psedobulk exp
# Siwei 28 Jul 2022

# init
library(ggplot2)
library(readr)
library(RColorBrewer)
library(stringr)

# import data

df_raw <-
  read_delim("LDSC_output_4_R.tsv",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

df_raw <-
  read_delim("Lexi_0.01_LDSC_output_4_R.tsv",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

df_raw <-
  read_delim("Lexi_0.02_LDSC_output_4_R.tsv",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

df_to_plot <- df_raw

df_to_plot <-
  df_to_plot[!(df_to_plot$Category %in% "baselineLD"), ]
df_to_plot <-
  df_to_plot[!(df_to_plot$Disease %in% "AD_jansen_from_txt"), ]
df_to_plot <-
  df_to_plot[!(df_to_plot$Disease %in% "Alz_Jansen_etal_2019"), ]
df_to_plot$Enrichment[df_to_plot$Enrichment < 0] <- NA

df_to_plot$logP <- 0 - log10(df_to_plot$Enrichment_p)
sum(is.na(df_to_plot$logP))

df_to_plot$Category <-
  factor(df_to_plot$Category,
         levels = unique(df_to_plot$Category)[c(1, 13, 7, 2, 14, 8,
                                               3, 15, 9, 4, 16, 10,
                                               5, 17, 11, 6, 18, 12)])
df_to_plot$Disease <-
  factor(df_to_plot$Disease,
         levels = unique(df_to_plot$Disease)[c(3, 1, 2, 5, 6,
                                               4, 7, 8)])

ggplot(df_to_plot,
       aes(x = Disease,
           y = Category,
           size = Enrichment,
           fill = logP)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low = "white",
                      high = "darkred") +
  # scale_size_continuous() +
  scale_radius(range = c(0, 6)) +
  xlab("Diseases") +
  labs(fill = "-log10P") +
  scale_y_discrete(limits = rev) +
  # scale_alpha_continuous(df_to_plot$`-logP`,
  #                        range = c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) +
  ggtitle("LDSC peaks enrichment")

