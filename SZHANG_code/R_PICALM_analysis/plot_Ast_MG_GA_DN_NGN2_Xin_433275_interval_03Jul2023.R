# Siwei 03 Jul 2023 #####
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

# -46 columns from GSE143666

# load iMG+iAst
iMG_Ast_raw <-
  read_delim("raw_data/OCRs_Jun2022_peaks_Xin_433275_featureCounts_03Jul2023.txt",
             delim = "\t", escape_double = FALSE, skip = 1,
             trim_ws = TRUE)
temp_rownames <-
  unlist(iMG_Ast_raw[, 1])
iMG_Ast_raw <-
  iMG_Ast_raw[, -c(1:6)]
rownames(iMG_Ast_raw) <- temp_rownames

# isolate colnames (sample names) and make them look better #####
raw_colnames <-
  colnames(iMG_Ast_raw)
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
colnames(iMG_Ast_raw) <-
  temp_colnames

iMG_Ast_metadata <-
  data.frame(samples = colnames(iMG_Ast_raw))
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
  read_delim("raw_data/iPSC_microglia_Novikova_etal_Xin_433275_featureCounts_03Jul2023.txt",
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
  read_delim("raw_data/phs001373_featureCounts_Xin_433275_peakset_03Jul2023.txt",
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
  data.frame(samples = c(colnames(iMG_Ast_raw),
                         colnames(iMG_Novikova),
                         colnames(phs001373_raw)))
master_metadata_table$cell_label <- ""

master_metadata_table$cell_label <-
  c(iMG_Ast_metadata$cell_label[match(colnames(iMG_Ast_raw),
                                      iMG_Ast_metadata$samples)],
    iMG_Novikova_metadata$cell_label,
    phs001373_metadata$cell_label[match(colnames(phs001373_raw),
                                        phs001373_metadata$assembled_names)])

master_metadata_table$cell_label[master_metadata_table$cell_label %in% "hMG"] <-
  "hMG-2"
master_metadata_table$cell_label[165:179][master_metadata_table$cell_label[165:179] %in% "hMG-2"] <- "hMG-1"
master_metadata_table$cell_label[master_metadata_table$cell_label %in% "iGA36"] <- "iGA"

# make master data table and normalise #####
master_raw_table <-
  as.data.frame(cbind(iMG_Ast_raw,
                      iMG_Novikova,
                      phs001373_raw))
rownames(master_raw_table) <-
  temp_rownames

save.image("Xin_433275_interval_5_celltypes_03Jul2023.RData")

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
cpm_mtable <-
  log1p(cpm_mtable)
cpm_mtable <-
  scale(cpm_mtable)

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
                         brewer.pal(n = 8, name = "Paired")),
             # label = "ind",
             invisible = (c("quali")))

save.image("Ast_MG_GA_DN_NGN2_Xin_433275_intervals_30Jun2023.RData")

load("~/Data/FASTQ/Duan_Project_014/Xin_s_original_peakset_sub/all_5_cell_types_with_microglia_666614.RData")
# try voom normalisation -- not good in Ast group #####
# voom accepts count as DGEgroup

original_samples <-
  colnames(summit_binding_matrix)
original_samples <-
  original_samples[str_detect(string = original_samples,
                              pattern = "iPS*",
                              negate = T)]
original_samples <-
  original_samples[str_detect(string = original_samples,
                              pattern = "MG*",
                              negate = T)]
original_samples <-
  str_replace_all(string = original_samples,
                  pattern = "_",
                  replacement = "-")
original_samples <-
  original_samples[str_detect(string = original_samples,
                              pattern = "NPC*",
                              negate = T)]
original_samples <-
  str_replace_all(string = original_samples,
                  pattern = "Glut",
                  replacement = "NGN2")
original_samples <-
  str_replace_all(string = original_samples,
                  pattern = "NGN2",
                  replacement = "NGN2-Glut")

master_metadata_table$sample_original <- "follow-up"
master_metadata_table$sample_original[master_metadata_table$samples %in% original_samples] <- "core"
master_metadata_table$shape_num <-
  as.numeric(as.factor(master_metadata_table$cell_label))

fviz_pca_ind(prin_comp_cpm,
             geom.ind = "point",
             geom.var = "text",
             # habillage = master_metadata_table$sample_original,
             palette = c(brewer.pal(n = 8, name = "Dark2"),
                         brewer.pal(n = 8, name = "Paired")),
             # label = "ind",
             invisible = (c("quali"))) +
  geom_point(aes(shape = as.factor(master_metadata_table$cell_label),
                 colour = as.factor(master_metadata_table$sample_original))) +
  scale_shape_manual(values = 1:13) +
  guides(shape = guide_legend(title = "Cell Type"),
         colour = guide_legend(title = "Samplei Batch"))


######
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
