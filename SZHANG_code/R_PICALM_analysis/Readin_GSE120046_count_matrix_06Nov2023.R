# Siwei 06 Nov 2023

# init ####
{
  library(Seurat)
  library(Signac)
  library(readr)
  library(future)
  library(parallel)
  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
}

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
options(future.seed = T)
set.seed(42)

# Load data ####
# load Raw UMI
raw_UMI <-
  CreateAssayObject(counts = as.matrix(read.table(gzfile("GSE120046/GSE120046_brain_all_UMIcounts.txt.gz",
                                                         'rt'),
                                                  header = T,
                                                  sep = "\t")),
                    check.matrix = T)

normalised_UMI <-
  CreateAssayObject(data = as.matrix(read.table(gzfile("GSE120046/GSE120046_brain_normlized.txt.gz",
                                                         'rt'),
                                                  header = T,
                                                  sep = "\t")),
                    check.matrix = F)

seurat_GSE120046_normalised_UMI <-
  CreateSeuratObject(normalised_UMI,
                     project = "GSE120046")
seurat_GSE120046_normalised_UMI$weeks <-
  str_replace(string = seurat_GSE120046_normalised_UMI$orig.ident,
              pattern = "W",
              replacement = "W_")
seurat_GSE120046_normalised_UMI$weeks <-
  str_split(string = seurat_GSE120046_normalised_UMI$weeks,
            pattern = "_",
            simplify = T)[, 1]

## run sctransform
Assays(seurat_GSE120046_normalised_UMI)
seurat_GSE120046_normalised_UMI <-
  FindVariableFeatures(seurat_GSE120046_normalised_UMI,
                       verbose = T)
seurat_GSE120046_normalised_UMI <-
  ScaleData(seurat_GSE120046_normalised_UMI,
            verbose = T)
seurat_GSE120046_normalised_UMI <-
  RunPCA(seurat_GSE120046_normalised_UMI,
         verbose = T)
seurat_GSE120046_normalised_UMI <-
  RunUMAP(seurat_GSE120046_normalised_UMI,
          dims = 1:30,
          verbose = T)

seurat_GSE120046_normalised_UMI <-
  FindNeighbors(seurat_GSE120046_normalised_UMI,
                dims = 1:30,
                verbose = T)
seurat_GSE120046_normalised_UMI <-
  FindClusters(seurat_GSE120046_normalised_UMI,
               verbose = T)

DimPlot(seurat_GSE120046_normalised_UMI,
        group.by = "weeks")

cell_metadata <-
  read.table(gzfile("GSE120046/GSE120046_metadata.txt.gz",
                    'rt'),
             header = T,
             sep = "\t")

seurat_GSE120046_normalised_UMI$Pos <-
  cell_metadata$Pos
seurat_GSE120046_normalised_UMI$type <-
  cell_metadata$type
seurat_GSE120046_normalised_UMI$subtype <-
  cell_metadata$subtype
seurat_GSE120046_normalised_UMI$tSNE_1 <-
  cell_metadata$tSNE_1
seurat_GSE120046_normalised_UMI$tSNE_2 <-
  cell_metadata$tSNE_2

DimPlot(seurat_GSE120046_normalised_UMI,
        group.by = "type",
        cols = brewer.pal(n = 11,
                          name = "Set3"),
        label = T,
        repel = T)
DimPlot(seurat_GSE120046_normalised_UMI,
        group.by = "subtype",
        # cols = brewer.pal(n = 11,
        #                   name = "Set3"),
        label = T,
        repel = T)
# DimPlot(seurat_GSE120046_normalised_UMI,
#         dims = c(cell_metadata$tSNE_1,
#                  cell_metadata$tSNE_2),
#         group.by = "weeks")

raw_UMI_filtered <-
  raw_UMI[, colnames(raw_UMI) %in% colnames(seurat_GSE120046_normalised_UMI)]
raw_UMI_filtered <-
  raw_UMI_filtered[rownames(raw_UMI_filtered) %in%
                     rownames(normalised_UMI), ]

pre_merge_GSE120046 <-
  CreateSeuratObject(raw_UMI_filtered,
                     project = "GSE120046")

pre_merge_GSE120046$Pos <-
  cell_metadata$Pos[colnames(normalised_UMI) %in%
                      colnames(pre_merge_GSE120046)]
pre_merge_GSE120046$type <-
  cell_metadata$type[colnames(normalised_UMI) %in%
                       colnames(pre_merge_GSE120046)]
pre_merge_GSE120046$subtype <-
  cell_metadata$subtype[colnames(normalised_UMI) %in%
                          colnames(pre_merge_GSE120046)]
pre_merge_GSE120046$tSNE_1 <-
  cell_metadata$tSNE_1[colnames(normalised_UMI) %in%
                         colnames(pre_merge_GSE120046)]
pre_merge_GSE120046$tSNE_2 <-
  cell_metadata$tSNE_2[colnames(normalised_UMI) %in%
                         colnames(pre_merge_GSE120046)]


saveRDS(pre_merge_GSE120046,
        file = "pre_merge_GSE120046_cleaned.RDs")
