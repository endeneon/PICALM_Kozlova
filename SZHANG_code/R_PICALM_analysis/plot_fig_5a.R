# Siwei 15 Jan 2024
# plot fig. 5a

# init ####
library(gplots)
library(RColorBrewer)
library(readxl)

library(stringr)

# load data
df_raw <-
  read_excel("Fig_5a.xlsx")

genotype_info <-
  df_raw$Genotype
gene_names <-
  df_raw$Gene_symbol

df_to_plot <-
  df_raw
df_to_plot$Gene_symbol <- NULL
df_to_plot$Genotype <- NULL
df_to_plot$Control <- NULL
rownames(df_to_plot) <-
  gene_names

genotype_colour <-
  genotype_info
genotype_colour <-
  str_replace_all(string = genotype_colour,
                  pattern = "homo",
                  replacement = '#B3E2CD')
genotype_colour <-
  str_replace_all(string = genotype_colour,
                  pattern = "Hetero",
                  replacement = '#FDCDAC')

heatmap.2(as.matrix(df_to_plot),
          Rowv = F,
          Colv = F,
          dendrogram = "none",
          scale = "row",
          col = colorRampPalette(colors = c("white",
                                            "purple4"))(200),
          sepcolor = "none",
          cellnote = format(as.matrix(df_to_plot),
                            digits = 2),
          notecex = 1,
          notecol = "black",
          trace = "none",
          labRow = gene_names,
          RowSideColors = genotype_colour)
