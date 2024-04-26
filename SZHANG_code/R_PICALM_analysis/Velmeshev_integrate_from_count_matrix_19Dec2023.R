# 19 Dec 2023 Siwei
# Reconstruct the Velmeshev raw object from scratch (mtx files)
# Will integrate as instructed in their github code
# https://github.com/velmeshevlab/dev_hum_cortex/blob/main/snRNAseq_integration

# init ####
{
  library(Seurat)
  library(Signac)
  library(readr)
  library(future)
  library(parallel)
  library(harmony)
  library(data.table)

  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
  library(gridExtra)
  library(viridis)
}

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
options(Seurat.object.assay.version = "v3")
options(future.seed = T)
set.seed(42)

#
mat <- Read10X("raw_data_sets",
               gene.column = 1)
meta <-
  read.table("raw_data_sets/meta.tsv",
             header = T,
             sep = "\t",
             as.is = T,
             row.names = 1)
Velmeshev_raw <-
  CreateSeuratObject(counts = mat,
                     project = "Velmeshev_709K_from_count_mtx",
                     meta.data = meta)
# saveRDS(Velmeshev_raw,
#         file = "Velmeshev_709K_raw_constructed_from_count_mtx_19Dec2023.RDs")
# rm(mat)

Velmeshev_normalised <-
  NormalizeData(Velmeshev_raw)
Velmeshev_normalised <-
  FindVariableFeatures(Velmeshev_normalised,
                       selection.method = "mean.var.plot")
Velmeshev_normalised <-
  ScaleData(Velmeshev_normalised,
            features = VariableFeatures(Velmeshev_normalised))
Velmeshev_normalised <-
  RunPCA(Velmeshev_normalised,
         features = VariableFeatures(Velmeshev_normalised),
         verbose = T)
ElbowPlot(object = Velmeshev_normalised,
          ndims = 50)

Velmeshev_normalised@reductions$pca2 <-
  Velmeshev_normalised@reductions$pca
Velmeshev_normalised@reductions$pca2@cell.embeddings <-
  Velmeshev_normalised@reductions$pca2@cell.embeddings[, 1:30]

Velmeshev_normalised <-
  RunHarmony(Velmeshev_normalised,
             group.by.vars = "chemistry",
             theta = 2,
             max.iter = 20,
             reduction.use = 'pca2',
             ncores = 32,
             dims.use = 1:30)
Velmeshev_normalised <-
  RunUMAP(Velmeshev_normalised,
          dims = 1:30,
          reduction = 'harmony',
          return.model = T)
Velmeshev_normalised <-
  FindNeighbors(Velmeshev_normalised,
                reduction = "harmony",
                dims = 1:30) %>%
  FindClusters()

# "10-20 years"   "Adult"         "4-10 years"    "2-4 years"     "0-1 years"
# "3rd trimester" "1-2 years"     "2nd trimester"
Velmeshev_normalised$age_range <-
  factor(Velmeshev_normalised$age_range,
         levels = c("2nd trimester",
                    "3rd trimester",
                    "0-1 years",
                    "1-2 years",
                    "2-4 years",
                    "4-10 years",
                    "10-20 years",
                    "Adult"))
# saveRDS(Velmeshev_normalised,
#         file = "Velmeshev_709K_processed_constructed_from_count_mtx_19Dec2023.RDs")



DimPlot(Velmeshev_normalised,
        reduction = "umap")
FeaturePlot(Velmeshev_normalised,
            features = c("SLC17A7", "GAD1"))
DimPlot(Velmeshev_normalised,
        reduction = "umap",
        group.by = "age_range",
        cols = plasma(n = length(unique(Velmeshev_normalised$age_range))))
DimPlot(Velmeshev_normalised,
        reduction = "umap",
        group.by = "lineage")

unique(Velmeshev_normalised$age)
unique(Velmeshev_normalised$age_range)
unique(Velmeshev_normalised$age.days.)
unique(Velmeshev_normalised$lineage)

unique(Velmeshev_normalised$age[Velmeshev_normalised$age_range == '0-1 years'])


### load in cells from 76 lines
# process cells_3_types_2K
cell_3_types_2K_each <-
  readRDS("cell_3_types_from_76_lines_2K_each_4_Velmeshev.RDs")

## SCT on each of the list items ####
Integrate_raw_df <-
  vector(mode = "list",
         length = 3L)
transfer_anchor_list <-
  vector(mode = "list",
         length = 3L)

for (i in 1:length(cell_3_types_2K_each)) {
  print(i)
  new_raw_x <- cell_3_types_2K_each[[i]]

  new_raw_x <-
    NormalizeData(new_raw_x)
  new_raw_x <-
    FindVariableFeatures(new_raw_x,
                         selection.method = "mean.var.plot")
  new_raw_x <-
    ScaleData(new_raw_x,
              features = VariableFeatures(new_raw_x))
  Integrate_raw_df[[i]] <- new_raw_x

  transfer_anchor_list[[i]] <-
    FindTransferAnchors(reference = Velmeshev_normalised,
                        query = new_raw_x,
                        dims = 1:30,
                        # normalization.method = "LogNormalize",
                        reference.reduction = "harmony")
}

names(Integrate_raw_df) <-
  names(cell_3_types_2K_each)
names(transfer_anchor_list) <-
  names(cell_3_types_2K_each)

#####
plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[1]],
           reference = Velmeshev_normalised,
           query = Integrate_raw_df[[1]],
           refdata = list(celltype = "age_range"),
           reference.reduction = "harmony",
           reduction.model = "umap")
unique(plot_query$predicted.celltype)
plot_query$predicted.celltype <-
  factor(plot_query$predicted.celltype,
         levels = c("2nd trimester",
                    "3rd trimester",
                    "0-1 years",
                    "1-2 years",
                    "2-4 years",
                    "4-10 years",
                    "10-20 years"))
DimPlot(plot_query,
        reduction = "ref.umap",
        group.by = "predicted.celltype",
        cols = plasma(n = length(unique(Velmeshev_normalised$age_range)))[c(1:7)],
        alpha = 0.5,
        pt.size = 0.4) +
  ggtitle(names(cell_3_types_2K_each)[1])
FeaturePlot(plot_query,
            features = c("GAD1"),
            alpha = 0.5,
            pt.size = 0.1,
            max.cutoff = 2,
            cols = c("lightgrey", "darkblue"))
FeaturePlot(plot_query,
            features = c("SLC17A7"),
            alpha = 0.5,
            pt.size = 0.1,
            max.cutoff = 5,
            cols = c("lightgrey", "darkblue"))

plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[2]],
           reference = Velmeshev_normalised,
           query = Integrate_raw_df[[2]],
           refdata = list(celltype = "age_range"),
           reference.reduction = "harmony",
           reduction.model = "umap")
unique(plot_query$predicted.celltype)
plot_query$predicted.celltype <-
  factor(plot_query$predicted.celltype,
         levels = c("2nd trimester",
                    "3rd trimester",
                    "0-1 years",
                    "1-2 years",
                    "2-4 years",
                    "4-10 years"))
DimPlot(plot_query,
        reduction = "ref.umap",
        group.by = "predicted.celltype",
        cols = plasma(n = length(unique(Velmeshev_normalised$age_range)))[c(1:6)],
        alpha = 0.5,
        pt.size = 0.4,
        cols = c("lightgrey", "darkblue")) +
  ggtitle(names(cell_3_types_2K_each)[2])




plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[3]],
           reference = Velmeshev_normalised,
           query = Integrate_raw_df[[3]],
           refdata = list(celltype = "age_range"),
           reference.reduction = "harmony",
           reduction.model = "umap")
unique(plot_query$predicted.celltype)
plot_query$predicted.celltype <-
  factor(plot_query$predicted.celltype,
         levels = c("2nd trimester",
                    "3rd trimester",
                    "0-1 years",
                    "1-2 years"))
DimPlot(plot_query,
        reduction = "ref.umap",
        group.by = "predicted.celltype",
        cols = plasma(n = length(unique(Velmeshev_normalised$age_range)))[c(1:6)],
        alpha = 0.5,
        pt.size = 0.4) +
  ggtitle(names(cell_3_types_2K_each)[3])

FeaturePlot(plot_query,
            features = c("SLC17A7", "GAD1"))

DimPlot(Velmeshev_normalised,
        reduction = "umap",
        cols = c("#FFFFFF00",
                 "#FFFFFF00",
                 "#FFFFFF00",
                 "#003285FF",
                 "#FFFFFF00",
                 "#FFFFFF00",
                 "#FFFFFF00",
                 "#FFFFFF00",
                 "#FFFFFF00"),
        group.by = "lineage")
DimPlot(Velmeshev_normalised,
        reduction = "umap",
        cols = c("#FFEEEE05",
                 "#FFDDCC01",
                 "#FFEEEE05",
                 "#003285FF",
                 "#FFEEEE05",
                 "#FFEEEE05",
                 "#FFEEEE05",
                 "#FFEEEE05",
                 "#FFEEEE05"),
        group.by = "lineage")
FeaturePlot(Velmeshev_normalised,
            features = "GAD1",
            max.cutoff = 2,
            cols = c("lightgrey", "darkblue"))
FeaturePlot(Velmeshev_normalised,
            features = "SLC17A7",
            # max.cutoff = 2,
            cols = c("lightgrey", "darkblue"))

DimPlot(Velmeshev_normalised,
        reduction = "umap",
        group.by = "lineage")
