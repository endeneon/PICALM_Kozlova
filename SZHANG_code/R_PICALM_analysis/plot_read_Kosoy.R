# Siwei 20 Jun 2023 #####
# plot a large PCA include MG, Ast, GA, and possibly NGN2
# ATAC-Seq data using the count matrix of Kosoy et al. (syn26207321)
# (microglia regulome)

# init #####
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

# load data #####

# load data of GSE143666 (Hauberg ME et al., Nat Comm 2020) #####
# ! GSE143666 and microglia_regulome use two different peak sets!
# ! DO NOT MIX!
GSE143666_raw <-
  read_delim("OCR_microglia_regulome/GSE143666_microglia_regulome_featureCounts_20Jun2023.txt",
             delim = "\t", escape_double = FALSE, skip = 1,
             trim_ws = TRUE)
temp_rownames <-
  unlist(GSE143666_raw[, 1])
GSE143666_raw <-
  GSE143666_raw[, -c(1:6)]
rownames(GSE143666_raw) <- temp_rownames
colnames(GSE143666_raw) <-
  str_split(string = colnames(GSE143666_raw),
            pattern = "\\.",
            simplify = T)[, 1]


GSE143666_metadata <-
  read_delim("OCR_microglia_regulome/GSE143666_MG_Ast_GABA_metadata.txt",
             delim = ",", escape_double = FALSE,
             trim_ws = TRUE)
GSE143666_metadata <-
  GSE143666_metadata[, c(1, 11, 30)]
GSE143666_metadata$assembled_names <-
  str_c(GSE143666_metadata$Cell_type,
        GSE143666_metadata$source_name,
        sep = "_")
GSE143666_metadata$assembled_names <-
  str_replace_all(string = GSE143666_metadata$assembled_names,
                  pattern = " ",
                  replacement = "_")
GSE143666_metadata$assembled_names <-
  str_c(GSE143666_metadata$Run,
        GSE143666_metadata$assembled_names,
        sep = "_")


colnames(GSE143666_raw) <-
  GSE143666_metadata$assembled_names[match(colnames(GSE143666_raw),
                                           GSE143666_metadata$Run)]


GSE143666_metadata$cell_label <-
  c(rep_len("GABA", length.out = 22),
    rep_len("MGAS", length.out = 24))

# load iMG+iAst
iMG_Ast_raw <-
  read_delim("OCR_microglia_regulome/all_MGs_Asts_OCR_microglia_regulome_count_matrix_featureCounts_01Jun2023.txt",
             delim = "\t", escape_double = FALSE, skip = 1,
             trim_ws = TRUE)
temp_rownames <-
  unlist(iMG_Ast_raw[, 1])
iMG_Ast_raw <-
  iMG_Ast_raw[, -c(1:6)]
rownames(iMG_Ast_raw) <- temp_rownames
colnames(iMG_Ast_raw) <-
  str_split(string = colnames(iMG_Ast_raw),
            pattern = "_S",
            simplify = T)[, 1]
colnames(iMG_Ast_raw) <-
  str_replace(string = colnames(iMG_Ast_raw),
              pattern = "_",
              replacement = "-")

iMG_Ast_metadata <-
  data.frame(samples = colnames(iMG_Ast_raw))
iMG_Ast_metadata$samples <-
  str_replace(string = iMG_Ast_metadata$samples,
              pattern = "_",
              replacement = "-")
iMG_Ast_metadata$cell_label <-
  str_split(string = iMG_Ast_metadata$samples,
            pattern = "-",
            simplify = T)[, 1]
iMG_Ast_metadata$cell_label <-
  str_c("i",
        iMG_Ast_metadata$cell_label,
        sep = "")


# load Novikova et al. 2021 (all iMGs)
iMG_Novikova <-
  read_delim("OCR_microglia_regulome/iPSC_microglia_Novikova_etal_2021_OCR_microglia_regulome_featureCounts_19Jun2023.txt",
             delim = "\t", escape_double = FALSE, skip = 1,
             trim_ws = TRUE)
temp_rownames <-
  unlist(iMG_Novikova[, 1])
iMG_Novikova <-
  iMG_Novikova[, -c(1:6)]
rownames(iMG_Novikova) <- temp_rownames
colnames(iMG_Novikova) <-
  c("iMG_Novikova_1",
    "iMG_Novikova_2",
    "iMG_Novikova_3")

iMG_Novikova_metadata <-
  data.frame(samples = colnames(iMG_Novikova),
             cell_label = "iMG_Novikova")

# load phs001373 (real human MG, Ast, etc.), the original microglia_regulome
# interval set (Kosoy et al.)
phs001373_raw <-
  read_delim("OCR_microglia_regulome/microglia_astrocyte_phs001373_featureCounts_19Jun2023.txt",
             delim = "\t", escape_double = FALSE, skip = 1,
             trim_ws = TRUE)
temp_rownames <-
  unlist(phs001373_raw[, 1])
phs001373_raw <-
  phs001373_raw[, -c(1:6)]
rownames(phs001373_raw) <- temp_rownames
colnames(phs001373_raw) <-
  str_split(string = colnames(phs001373_raw),
            pattern = "\\.",
            simplify = T)[, 1]
colnames(phs001373_raw) <-
  str_split(string = colnames(phs001373_raw),
            pattern = "\\/",
            simplify = T)[, 2]

phs001373_metadata <-
  read_delim("OCR_microglia_regulome/microglia_astrocyte_metadata.txt",
             delim = ",", escape_double = FALSE,
             trim_ws = TRUE)
phs001373_metadata <-
  phs001373_metadata[, c(1, 20)]
phs001373_metadata$assembled_names <-
  str_c(phs001373_metadata$Run,
        phs001373_metadata$histological_type,
        sep = "_")

colnames(phs001373_raw) <-
  phs001373_metadata$assembled_names[match(colnames(phs001373_raw),
                                           phs001373_metadata$Run)]
phs001373_metadata$cell_label <-
  recode(.x = phs001373_metadata$histological_type,
         Neuron = "hNeuron",
         Microglia = "hMG",
         Oligodendrocyte = "hOlig",
         Astrocyte = "hAst")

master_metadata_table <-
  data.frame(samples = c(colnames(GSE143666_raw),
                         colnames(iMG_Ast_raw),
                         colnames(iMG_Novikova),
                         colnames(phs001373_raw)))
master_metadata_table$cell_label <- ""

master_metadata_table$cell_label <-
  c(GSE143666_metadata$cell_label[match(colnames(GSE143666_raw),
                                        GSE143666_metadata$assembled_names)],
    iMG_Ast_metadata$cell_label[match(colnames(iMG_Ast_raw),
                                      iMG_Ast_metadata$samples)],
    iMG_Novikova_metadata$cell_label,
    phs001373_metadata$cell_label[match(colnames(phs001373_raw),
                                        phs001373_metadata$assembled_names)])
master_metadata_table$cell_label[master_metadata_table$cell_label %in% "hMG"] <-
  "hMG-2"
master_metadata_table$cell_label[120:135][master_metadata_table$cell_label[120:135] %in% "hMG-2"] <- "hMG-1"

# make master data table and normalise #####
master_raw_table <-
  as.data.frame(cbind(GSE143666_raw,
                      iMG_Ast_raw,
                      iMG_Novikova,
                      phs001373_raw))
rownames(master_raw_table) <-
  temp_rownames

save.image("Kosoy_interval_20Jun2023.RData")

# use conventional normalisation (log1p) #####
DGE_master_table <-
  DGEList(counts = as.matrix(master_raw_table),
          samples = master_metadata_table$samples,
          genes = rownames(master_raw_table),
          remove.zeros = T)

cpm_mtable <-
  as.data.frame(cpm(DGE_master_table))
# hist(unlist(cpm_mtable))

# test different normalisation methods
# cpm_mtable <-
#   log1p(cpm_mtable)
# cpm_mtable <-
#   scale(cpm_mtable)

# calc prin comp
prin_comp_cpm <-
  prcomp(as.matrix(t(cpm_mtable)),
         center = T,
         scale. = T)

fviz_pca_ind(prin_comp_cpm,
             geom.ind = "point",
             geom.var = "text",
             habillage = master_metadata_table$cell_label,
             palette = c(brewer.pal(n = 8, name = "Dark2"),
                         "darkblue", "gold3"),
             # label = "ind",
             invisible = (c("quali")))


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
