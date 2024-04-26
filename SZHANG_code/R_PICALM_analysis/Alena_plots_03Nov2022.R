# Siwei 03 Nov 2022
# make plots for Alena

library(ggplot2)
library(readxl)
library(stringr)

Alena_table <- read_excel("Alena_table.xlsx")

# add a small value to BinomFDRQ for -log10 transformation
Alena_table$BinomFdrQ <-
  Alena_table$BinomFdrQ + 1e-301

df_to_plot <-
  Alena_table[1:20, ]
df_to_plot$GO_disc <-
  str_c(df_to_plot$ID, 
        df_to_plot$Desc,
        sep = ":")

df_to_plot$`-logFDR` <-
  0 - log10(df_to_plot$BinomFdrQ)

df_to_plot$GO_disc <-
  factor(df_to_plot$GO_disc,
         levels = df_to_plot$GO_disc[order(df_to_plot$`-logFDR`)])

ggplot(df_to_plot,
       aes(fill = RegionFoldEnrich,
           colour = RegionFoldEnrich,
           x = GO_disc,
           ymin = 0,
           ymax = `-logFDR`)) +
  geom_linerange(colour = "grey10") +
  geom_point(aes(x = GO_disc,
                 y = `-logFDR`),
             shape = 21,
             colour = "grey50",
             size = 4) +
  # scale_colour_gradient(low = "grey90",
  #                       high = "darkred") +
  scale_fill_gradient(low = "grey90",
                      high = "darkred",
                      name = "Regional Enrichment Fold") +
  labs(x = "Top 20 GO terms ranked by enrichment FDR value", 
       y = "-log10FDR enrichment") +
  theme_classic() +
  coord_flip()


### volcano plot of the 360k MG ASoC SNPs
library(readr)

df_raw <-
  read_delim("DP20_data_files/MG_28_lines_SNP_361KK_06Jun2022.txt", 
                delim = "\t", escape_double = FALSE, 
                trim_ws = TRUE)

df_raw$FDR <- df_raw$FDR + 1e-295
df_raw$logFDR <- 0 - log10(df_raw$FDR)

df_raw$sigif <- 'Non-significant (FDR > 0.05)'
df_raw$sigif[df_raw$FDR < 0.05] <- 'Significant (FDR < 0.05)'

ggplot(df_raw,
       aes(x = ratio,
           y = logFDR,
           colour = sigif)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("black", "red"),
                     name = "SNP ASoC significance") +
  ylim(0, 100) +
  labs(x = "reference allele ratio",
       y = '-log10 FDR') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14))


Jeff_set <-
  1:14423
Other_set <-
  (14423 - 12363):(14423 - 12363 + 16966)

for (i in 1:10000) {
  print(i)
  Jeff_subset <-
    Jeff_set[sample(size = 2584,
                    x = length(Jeff_set))]
  Other_subset <-
    Other_set[sample(size = 3771,
                     x = length(Other_set))]
  shared_count <-
    sum(Jeff_subset %in% Other_subset)
  if (i == 1) {
    hist_plot <- shared_count
  } else {
    hist_plot <-
      c(hist_plot, 
        shared_count)
  }
}

hist(hist_sum,
     breaks = 100)
