# Siwei 09 Jun 2022
# prepare data and plot for MAGMA analysis

# init

# load data
library(readr)

library(ggplot2)
library(RColorBrewer)

Ast_MG_RGlut_magma_Gene_list <- 
  read_delim("Ast_MG_RGlut_magma_Gene_list.txt", 
             delim = "\t", escape_double = FALSE, 
             col_names = FALSE, trim_ws = TRUE)
Ast_MG_RGlut_magma_Gene_list$category <-
  c(rep_len(x = "Ast", length.out = 5196),
    rep_len(x = "MG", length.out = 5727),
    rep_len(x = "RGlut", length.out = 5096))

write.table(Ast_MG_RGlut_magma_Gene_list,
            file = "Ast_MG_RGlut_magma.set.annot",
            quote = F, sep = "\t",
            row.names = F, col.names = F)
