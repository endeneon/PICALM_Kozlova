# Siwei 15 Mar 2024
# Make Ast plots
# make plots for Alena

library(readr)
library(ggplot2)
library(readxl)
library(stringr)

# Alena_table <- read_excel("Alena_table.xlsx")
Alena_table <-
  read_delim("~/backuped_space/Siwei_misc_R_projects/R_MG_17/Ast_GREAT_enrich_100.tsv", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE, skip = 1)

# add a small value to BinomFDRQ for -log10 transformation
Alena_table$`Binom FDR Q-Val` <-
  Alena_table$`Binom FDR Q-Val` + 1e-250

df_to_plot <-
  Alena_table[1:20, ]
df_to_plot$GO_disc <-
  str_c(df_to_plot$`Term ID`,
        df_to_plot$`# Term Name`, 
        sep = ":")

df_to_plot$`-logFDR` <-
  0 - log10(df_to_plot$`Binom FDR Q-Val`)

df_to_plot$GO_disc <-
  factor(df_to_plot$GO_disc,
         levels = df_to_plot$GO_disc[order(df_to_plot$`-logFDR`)])

ggplot(df_to_plot,
       aes(fill = `Binom Fold Enrichment`,
           colour = `Binom Fold Enrichment`,
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
  # ylim(0, max(df_to_plot$`-logFDR`)) +
  theme_classic() +
  coord_flip()

max(df_to_plot$`-logFDR`)

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
