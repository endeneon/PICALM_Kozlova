# Siwei 20 Jun 2023 #####
# plot a large PCA include MG, Ast, GA, and possibly NGN2
# ATAC-Seq data using the count matrix of Kosoy et al. (syn26207321)
# (microglia regulome)

# init #####
{
  library(readr)

  library(edgeR)
  library(Rfast)
  library(factoextra)

  library(Rtsne)
  library(irlba)

  library(stringr)

  library(dplyr)
  library(tidyverse)

  library(ggplot2)
  library(RColorBrewer)
}

# load data #####
df_raw <-
  read_delim("raw_data/collapsed_Ast_MG_GA_DN_NGN2_501bp_count_matrix_featureCounts_29Jun2023.txt",
             delim = "\t", escape_double = FALSE, skip = 1,
             trim_ws = TRUE)
temp_rownames <-
  unlist(df_raw[, 1])
df_raw <-
  df_raw[, -c(1:6)]
## Ast samples were counted twice (1-18), need to be removed #####
df_raw <-
  df_raw[, -c(1:18)]
rownames(df_raw) <- temp_rownames

# isolate colnames (sample names) and make them look better #####
raw_colnames <-
  colnames(df_raw)
temp_colnames <-
  str_split(string = raw_colnames,
            pattern = "/",
            simplify = F)
temp_colnames <-
  lapply(X = temp_colnames,
         FUN = function(x) {
           return(tail(x = x,
                       n = 1L,
                       keepnums = F))
         })
temp_colnames <-
  unlist(temp_colnames)
temp_colnames <-
  str_split(string = temp_colnames,
            pattern = "_S[0-9]_trimmed",
            simplify = T)[, 1]
temp_colnames <-
  str_split(string = temp_colnames,
            pattern = "_S[0-9][0-9]_trimmed",
            simplify = T)[, 1]
temp_colnames <-
  str_split(string = temp_colnames,
            pattern = "_S[0-9]_WASPed",
            simplify = T)[, 1]
temp_colnames <-
  str_split(string = temp_colnames,
            pattern = "_S[0-9][0-9]_WASPed",
            simplify = T)[, 1]
temp_colnames <-
  str_replace_all(string = temp_colnames,
                  pattern = "Glut_rapid_neuron20",
                  replacement = "NGN2-Glut")
temp_colnames <-
  str_replace_all(string = temp_colnames,
                  pattern = "R21-",
                  replacement = "NGN2-Glut-")
temp_colnames <-
  str_split(string = temp_colnames,
            pattern = "_S[0-9]",
            simplify = T)[, 1]
temp_colnames <-
  str_split(string = temp_colnames,
            pattern = "_WASPed",
            simplify = T)[, 1]
temp_colnames <-
  str_replace_all(string = temp_colnames,
                  pattern = "_",
                  replacement = "-")

## assign cleaned colnames back, store as a new df #####
df_4_DGE <- df_raw
colnames(df_4_DGE) <- temp_colnames
## reorder samples by name #####
df_4_DGE <-
  df_4_DGE[, order(colnames(df_4_DGE))]

# make the master DGE table

DGE_master_table <-
  DGEList(counts = as.matrix(df_4_DGE),
          samples = colnames(df_4_DGE),
          genes = rownames(df_4_DGE),
          remove.zeros = T)

cpm_mtable <-
  as.data.frame(cpm(DGE_master_table))
# hist(unlist(cpm_mtable))

## construct a metadata table
master_metadata_table <-
  data.frame(sample_name = colnames(cpm_mtable))
master_metadata_table$cell_type <-
  str_replace_all(string = master_metadata_table$sample_name,
                  pattern = "GA36",
                  replacement = "GA")
master_metadata_table$cell_type <-
  str_replace_all(string = master_metadata_table$cell_type,
                  pattern = "NGN2-",
                  replacement = "NGN2_")
master_metadata_table$cell_type <-
  str_split(string = master_metadata_table$cell_type,
            pattern = "-",
            simplify = T)[, 1]
master_metadata_table$cell_type <-
  str_replace_all(string = master_metadata_table$cell_type,
                  pattern = "NGN2_",
                  replacement = "NGN2-")

# test different normalisation methods
cpm_mtable_transformed <-
  log1p(cpm_mtable)
cpm_mtable_transformed <-
  scale(cpm_mtable_transformed)

# calc prin comp
prin_comp_cpm <-
  prcomp(as.matrix(t(cpm_mtable_transformed)),
         center = T,
         scale. = T)

fviz_pca_ind(prin_comp_cpm,
             geom.ind = "point", # "point" or "text"
             geom.var = "text",
             habillage = master_metadata_table$cell_type,
             palette = c(brewer.pal(n = 8, name = "Dark2"),
                         "darkblue", "gold3"),
             label = "ind",
             invisible = (c("quali")))

save.image("Ast_MG_GA_DN_NGN2_501bp_intervals_30Jun2023.RData")
# phs001373_metadata$cell_label <-
#   recode(.x = phs001373_metadata$histological_type,
#          Neuron = "hNeuron",
#          Microglia = "hMG",
#          Oligodendrocyte = "hOlig",
#          Astrocyte = "hAst")

# try voom normalisation -- not good in Ast group #####
# voom accepts count as DGEgroup

DGE_master_table <-
  master_raw_table[rowSums(master_raw_table, na.rm = F) != 0, ]
DGE_master_table <-
  DGE_master_table + 1

DGE_master_table <-
  DGEList(counts = as.matrix(DGE_master_table),
          samples = master_metadata_table$samples,
          genes = rownames(DGE_master_table),
          group = master_metadata_table$cell_label,
          remove.zeros = T)
DGE_voom <-
  voom(DGE_master_table)


prin_comp_cpm <-
  prcomp(as.matrix(t(DGE_voom$E)),
         center = F,
         scale. = T)

fviz_pca_ind(prin_comp_cpm,
             geom = "point",
             habillage = master_metadata_table$cell_label,
             invisible = (c("quali")))
