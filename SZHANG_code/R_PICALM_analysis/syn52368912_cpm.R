# Extract Microglia , use syn52368912



{
  library(stringr)
  library(Seurat)

  library(parallel)
  library(future)


  library(glmGamPoi)

  library(edgeR)

  library(data.table)

  library(readr)

  library(readxl)

  library(stringr)
  library(ggplot2)

  library(scales)
  library(reshape2)

  library(RColorBrewer)
  library(ggpubr)

  library(dplyr)
  library(data.table)

  library(DescTools)
  library(multcomp)

  library(gridExtra)

  library(harmony)

  plan("multisession", workers = 3)
  # options(mc.cores = 32)
  set.seed(42)
  options(future.globals.maxSize = 429496729600)

  setwd("~/backuped_space/Siwei_misc_R_projects/R_Alena_PICALM")
}

all_RDS_files <-
  list.files("~/Data/FASTQ/sage_synapse/syn52368912",
             pattern = ".*.rds",
             full.names = T,
             recursive = F)
names(all_RDS_files) <-
  str_split(string = all_RDS_files,
            pattern = '/',
            simplify = T)[, 8]
# names(all_RDS_files)


PICALM_results <-
  vector(mode = "list",
         length = length(all_RDS_files))
names(PICALM_results) <-
  names(all_RDS_files)

Immune_cells <-
  readRDS(all_RDS_files[5])
unique(Immune_cells$cell_type_high_resolution)
Immune_cells <-
  Immune_cells[, Immune_cells$cell_type_high_resolution == 'Mic P2RY12']

pseudo_projid_counts <-
  AggregateExpression(Immune_cells,
                      group.by = "projid",
                      return.seurat = T)

as.matrix(pseudo_projid_counts[rownames(pseudo_projid_counts) == "PICALM", ])

pseudo_counts_4_edgeR <-
  pseudo_projid_counts@assays$RNA@layers$counts
rownames(pseudo_counts_4_edgeR) <-
  rownames(pseudo_projid_counts)
colnames(pseudo_counts_4_edgeR) <-
  colnames(pseudo_projid_counts)

pseudo_DGE <-
  DGEList(counts = pseudo_counts_4_edgeR)
rownames(pseudo_DGE)

pseudo_DGE_cpm <-
  as.data.frame(cpm(pseudo_DGE))
pseudo_DGE_cpm[rownames(pseudo_DGE_cpm) == "PICALM", ]

write.table(pseudo_DGE_cpm,
      file = "Immune_cells_P2RY12_pos_cpm_wo_log.tsv",
      sep = "\t", quote = F)

write.table(Immune_cells_P2RY12_meta_diag,
            file = "all_immune_cells_w_diag.tsv",
            row.names = F, col.names = T,
            sep = "\t", quote = F)
