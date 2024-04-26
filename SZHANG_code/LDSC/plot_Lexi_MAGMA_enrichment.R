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
#   read_delim("results_4_R.txt",
#              delim = "\t", escape_double = FALSE,
#              trim_ws = TRUE)
# df_raw <-
#   read_delim("results_4_R_08Aug2022.txt",
#              delim = "\t", escape_double = FALSE,
#              trim_ws = TRUE)

df_raw <-
  read_table2("results_4_R_08Aug2022.txt",
              col_names = FALSE)

df_to_plot <- df_raw

colnames(df_to_plot) <- c("GENE_SET", "TYPE", "NGENES", "BETA",
                          "BETA_STD", "SE", "P", "FULL_NAME", "DISEASES")



# df_to_plot <-
#   df_to_plot[!(df_to_plot$Category %in% "baselineLD"), ]
# df_to_plot <-
#   df_to_plot[!(df_to_plot$Disease %in% "AD_jansen_from_txt"), ]
# df_to_plot <-
#   df_to_plot[!(df_to_plot$Disease %in% "Alz_Jansen_etal_2019"), ]
df_to_plot$BETA[df_to_plot$BETA < 0] <- NA

# df_to_plot$logP <- 0
df_to_plot$logP <- 0 - log10(df_to_plot$P)
sum(is.na(df_to_plot$logP))

df_to_plot$DISEASES <-
  str_remove_all(string = df_to_plot$DISEASES,
                 pattern = "_20k_5k\\.$")

# df_to_plot$DISEASES <-
#   factor(df_to_plot$DISEASES,
#          levels = unique(df_to_plot$DISEASES)[c(1, 13, 7, 2, 14, 8,
#                                                 3, 15, 9, 4, 16, 10,
#                                                 5, 17, 11, 6, 18, 12)])
# df_to_plot$Disease <-
#   factor(df_to_plot$Disease,
#          levels = unique(df_to_plot$Disease)[c(3, 1, 2, 5, 6,
#                                                4, 7, 8)])

df_to_plot <-
  df_to_plot[!(df_to_plot$DISEASES %in% "Ulcerative_colitis"), ]

df_to_plot$DISEASES <-
  factor(df_to_plot$DISEASES,
         levels = unique(df_to_plot$DISEASES)[c(13, 15, 3, 4, 9,
                                                10, 12, 1,
                                                2, 5, 6, 7, 8,
                                                11, 14)])

df_to_plot$FULL_NAME <-
  factor(df_to_plot$FULL_NAME,
         levels = unique(df_to_plot$FULL_NAME)[c(1, 3, 4, 6,
                                                 7, 9, 10, 12,
                                                 13, 15, 16, 18,
                                                 2, 5, 8, 11, 14, 17)])


ggplot(df_to_plot,
       aes(x = DISEASES,
           y = FULL_NAME,
           size = BETA,
           fill = logP)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low = "white",
                      high = "darkred") +
  # scale_size_continuous() +
  scale_radius() +
  xlab("Diseases") +
  ylab("Gene category") +
  labs(fill = "-log10P") +
  scale_y_discrete(limits = rev) +
  # scale_alpha_continuous(df_to_plot$`-logP`,
  #                        range = c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) +
ggtitle("MAGMA enrichment of Lexi's gene list, 20K, 5K")
