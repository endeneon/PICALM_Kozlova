# Siwei 16 Jun 2023 #####
# plot a large PCA include MG, Ast, GA, and possibly NGN2
# ATAC-Seq data using the count matrix of Hauberg ME et al. (GSE143666)
# need to add data from Kosoy R et al. (syn26207321)

# init #####
library(readr)

library(edgeR)
library(Rfast)
library(factoextra)

library(Rtsne)
library(irlba)

library(stringr)

# load data of GSE143666 (Hauberg ME et al., Nat Comm 2020) #####
# ! GSE143666 and microglia_regulome use two different peak sets!

# ! DO NOT MIX!
GSE143666_raw <-
  read_delim("raw_data/GSE143666_Count_matrix.tsv",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

MG_Ast_raw <-
  read_delim("raw_data/all_Ast_MG_use_GSE143666_interval_09Jun2023.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

MG_Ast <-
  MG_Ast_raw[, -c(2:6)]
colnames(MG_Ast) <-
  str_split(string = colnames(MG_Ast),
            pattern = "_S",
            simplify = T)[, 1]

master_table <-
  merge(x = GSE143666_raw,
        y = MG_Ast,
        by = "Peak_ID")
Peak_ID_list <-
  master_table$Peak_ID
master_table <-
  master_table[, -1]
rownames(master_table) <- Peak_ID_list

DGE_master_table <-
  DGEList(counts = as.matrix(master_table),
          samples = colnames(master_table),
          genes = rownames(master_table),
          remove.zeros = T)

cpm_mtable <-
  as.data.frame(cpm(DGE_master_table))
# hist(unlist(cpm_mtable))

# sum(rowMeans(cpm_mtable) > 3)
# cpm_mtable <-
#   cpm_mtable[rowMeans(cpm_mtable) > 1, ]

cpm_mtable <-
  log1p(cpm_mtable)

# calc prin comp

prin_comp_cpm <-
  prcomp(as.matrix(t(cpm_mtable)),
         center = T,
         scale. = T)

colour_seq <-
  data.frame(raw_name = colnames(cpm_mtable))
colour_seq$seq <-
  str_remove_all(string = colour_seq$raw_name,
                 pattern = "^S[0-9]_")
colour_seq$seq_final <-
  str_split(string = colour_seq$seq,
            pattern = "\\-",
            simplify = T)[, 1]
colour_seq$seq_final <-
  str_replace_all(string = colour_seq$seq_final,
                  pattern = "MG_28",
                  replacement = "MG")

fviz_pca_ind(prin_comp_cpm,
             geom = "point",
             habillage = colour_seq$seq_final,
             invisible = (c("quali")))
# View(USArrests)


### plot a tsne as in Hauberg ME #####
tsne_comp_cpm <-
  Rtsne(as.matrix(t(cpm_mtable)),
        pca = F,
        partial_pca = T,
        pca_scale = T,
        verbose = T,
        num_threads = 6)

df_to_plot <-
  data.frame(x = tsne_comp_cpm$Y[, 1],
             y = tsne_comp_cpm$Y[, 2],
             cell_type = colour_seq$seq_final)

df_to_plot <-
  df_to_plot[!df_to_plot$cell_type %in% c("MG", "Ast"), ]

ggplot(df_to_plot,
       aes(x = x,
           y = y,
           colour = cell_type)) +
  geom_point() +
  theme_classic()


# try voom normalisation #####
# voom accepts count as DGEgroup

DGE_master_table <-
  master_table[rowSums(master_table, na.rm = F) != 0, ]
DGE_master_table <-
  DGE_master_table + 1

DGE_master_table <-
  DGEList(counts = as.matrix(DGE_master_table),
          samples = colnames(DGE_master_table),
          genes = rownames(DGE_master_table),
          group = colour_seq$seq_final,
          remove.zeros = T)
DGE_voom <-
  voom(DGE_master_table)


prin_comp_cpm <-
  prcomp(as.matrix(t(DGE_voom$E)),
         center = T,
         scale. = T)

# colour_seq <-
#   data.frame(raw_name = colnames(cpm_mtable))
# colour_seq$seq <-
#   str_remove_all(string = colour_seq$raw_name,
#                  pattern = "^S[0-9]_")
# colour_seq$seq_final <-
#   str_split(string = colour_seq$seq,
#             pattern = "\\-",
#             simplify = T)[, 1]
# colour_seq$seq_final <-
#   str_replace_all(string = colour_seq$seq_final,
#                   pattern = "MG_28",
#                   replacement = "MG")

fviz_pca_ind(prin_comp_cpm,
             geom = "point",
             habillage = colour_seq$seq_final,
             invisible = (c("quali")))
