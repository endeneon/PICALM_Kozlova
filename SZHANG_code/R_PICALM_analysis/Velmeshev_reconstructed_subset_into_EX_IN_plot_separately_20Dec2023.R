# 20 Dec 2023 Siwei

# Will integrate as instructed in their github code
# https://github.com/velmeshevlab/dev_hum_cortex/blob/main/snRNAseq_integration

# Subset the reconstructed Velmeshev object into Ex and IN subsets
# Need to make the plots separately

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
  library(scales)
  library(viridis)
}

plan("multisession", workers = 12)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
options(Seurat.object.assay.version = "v3")
options(future.seed = T)
set.seed(42)

Velmeshev_normalised <-
  readRDS(file = "Velmeshev_709K_processed_constructed_from_count_mtx_19Dec2023.RDs")
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

Idents(Velmeshev_normalised) <- "lineage"

# ExNeu_Velmeshev <-
#   subset(Velmeshev_normalised,
#          idents = "ExNeu")
# IN_Velmeshev <-
#   subset(Velmeshev_normalised,
#          idents = "IN")

# saveRDS(ExNeu_Velmeshev,
#         file = "Velmeshev_processed_ExNeu_20Dec2023.RDs")
# saveRDS(IN_Velmeshev,
#         file = "Velmeshev_processed_IN_20Dec2023.RDs")

##### make from raw
Velmeshev_raw <-
  readRDS(file = "Velmeshev_709K_raw_constructed_from_count_mtx_19Dec2023.RDs")
Velmeshev_raw$age_range <-
  factor(Velmeshev_raw$age_range,
         levels = c("2nd trimester",
                    "3rd trimester",
                    "0-1 years",
                    "1-2 years",
                    "2-4 years",
                    "4-10 years",
                    "10-20 years",
                    "Adult"))

Idents(Velmeshev_raw) <- "lineage"
ExNeu_Velmeshev <-
  subset(Velmeshev_raw,
         idents = "ExNeu")
IN_Velmeshev <-
  subset(Velmeshev_raw,
         idents = "IN")

## Process ExNeu ####
ExNeu_Velmeshev <-
  NormalizeData(ExNeu_Velmeshev)
ExNeu_Velmeshev <-
  FindVariableFeatures(ExNeu_Velmeshev,
                       selection.method = "mean.var.plot")
ExNeu_Velmeshev <-
  ScaleData(ExNeu_Velmeshev,
            features = VariableFeatures(ExNeu_Velmeshev))
ExNeu_Velmeshev <-
  RunPCA(ExNeu_Velmeshev,
         features = VariableFeatures(ExNeu_Velmeshev),
         verbose = T)
ElbowPlot(object = ExNeu_Velmeshev,
          ndims = 50)

ExNeu_Velmeshev@reductions$pca2 <-
  ExNeu_Velmeshev@reductions$pca
ExNeu_Velmeshev@reductions$pca2@cell.embeddings <-
  ExNeu_Velmeshev@reductions$pca2@cell.embeddings[, 1:30]

ExNeu_Velmeshev <-
  RunHarmony(ExNeu_Velmeshev,
             group.by.vars = "chemistry",
             theta = 2,
             max.iter = 20,
             reduction.use = 'pca2',
             ncores = 32,
             dims.use = 1:30)
ExNeu_Velmeshev <-
  RunUMAP(ExNeu_Velmeshev,
          dims = 1:30,
          reduction = 'harmony',
          return.model = T)
ExNeu_Velmeshev <-
  FindNeighbors(ExNeu_Velmeshev,
                reduction = "harmony",
                dims = 1:30) %>%
  FindClusters()
DimPlot(ExNeu_Velmeshev,
        reduction = "umap",
        group.by = "age_range",
        cols = plasma(n = length(unique(ExNeu_Velmeshev$age_range))))
# saveRDS(ExNeu_Velmeshev,
#         file = "Velmeshev_processed_ExNeu_20Dec2023.RDs")

## Process IN ####
plan("multisession", workers = 32)
IN_Velmeshev <-
  NormalizeData(IN_Velmeshev,
                verbose = T)
IN_Velmeshev <-
  FindVariableFeatures(IN_Velmeshev,
                       selection.method = "mean.var.plot",
                       verbose = T)
IN_Velmeshev <-
  ScaleData(IN_Velmeshev,
            features = VariableFeatures(IN_Velmeshev),
            verbose = T)
IN_Velmeshev <-
  RunPCA(IN_Velmeshev,
         features = VariableFeatures(IN_Velmeshev),
         verbose = T)
ElbowPlot(object = IN_Velmeshev,
          ndims = 50)

IN_Velmeshev@reductions$pca2 <-
  IN_Velmeshev@reductions$pca
IN_Velmeshev@reductions$pca2@cell.embeddings <-
  IN_Velmeshev@reductions$pca2@cell.embeddings[, 1:30]

IN_Velmeshev <-
  RunHarmony(IN_Velmeshev,
             group.by.vars = "chemistry",
             theta = 2,
             max.iter = 20,
             reduction.use = 'pca2',
             ncores = 32,
             dims.use = 1:30)
IN_Velmeshev <-
  RunUMAP(IN_Velmeshev,
          dims = 1:30,
          reduction = 'harmony',
          return.model = T)
IN_Velmeshev <-
  FindNeighbors(IN_Velmeshev,
                reduction = "harmony",
                dims = 1:30) %>%
  FindClusters()
DimPlot(IN_Velmeshev,
        reduction = "umap",
        group.by = "age_range",
        cols = plasma(n = length(unique(IN_Velmeshev$age_range))),
        alpha = 0.5) +
  ggtitle("Interneurons")

# saveRDS(IN_Velmeshev,
#         file = "Velmeshev_processed_IN_20Dec2023.RDs")


# load in cells from 76 lines ####
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


### Process Interneuron (GABA, #1) #####
for (i in 1:1) {
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
  new_raw_x <-
    RunPCA(new_raw_x,
           features = VariableFeatures(new_raw_x))
  Integrate_raw_df[[i]] <- new_raw_x

  transfer_anchor_list[[i]] <-
    FindTransferAnchors(reference = IN_Velmeshev,
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
           reference = IN_Velmeshev,
           query = Integrate_raw_df[[1]],
           refdata = list(celltype = "age_range"),
           reference.reduction = "harmony",
           reduction.model = "umap")
length(unique(plot_query$predicted.celltype))
unique(plot_query$predicted.celltype)
plot_query$predicted.celltype <-
  factor(plot_query$predicted.celltype,
         levels = c("2nd trimester",
                    "3rd trimester",
                    "0-1 years",
                    "1-2 years",
                    "2-4 years",
                    "4-10 years",
                    "10-20 years",
                    "Adult"))
DimPlot(plot_query,
        reduction = "ref.umap",
        group.by = "predicted.celltype",
        cols = plasma(n = length(unique(IN_Velmeshev$age_range))),
        alpha = 0.5,
        pt.size = 0.1) +
  ggtitle(names(cell_3_types_2K_each)[1])
table(plot_query$predicted.celltype)
# table(plot_query$predicted.celltype)

cell_prediction_count_list[[1]] <-
  as.data.frame(table(plot_query$predicted.celltype))

#####
### Process ExNeu (Glut, #2-3) #####
DimPlot(ExNeu_Velmeshev,
        reduction = "umap",
        group.by = "age_range",
        cols = plasma(n = length(unique(ExNeu_Velmeshev$age_range))),
        alpha = 0.5) +
  ggtitle("Excitatory Neurons")


for (i in 2:3) {
  print(i)
  new_raw_x <- cell_3_types_2K_each[[i]]

  new_raw_x <-
    NormalizeData(new_raw_x,
                  verbose = T)
  new_raw_x <-
    FindVariableFeatures(new_raw_x,
                         selection.method = "mean.var.plot",
                         verbose = T)
  new_raw_x <-
    ScaleData(new_raw_x,
              features = VariableFeatures(new_raw_x),
              verbose = T)
  new_raw_x <-
    RunPCA(new_raw_x,
           features = VariableFeatures(new_raw_x),
           verbose = T)
  Integrate_raw_df[[i]] <- new_raw_x

  transfer_anchor_list[[i]] <-
    FindTransferAnchors(reference = ExNeu_Velmeshev,
                        query = new_raw_x,
                        dims = 1:30,
                        # normalization.method = "LogNormalize",
                        reference.reduction = "harmony")
}

for (i in 2:3) {
  transfer_anchor_list[[i]] <-
    FindTransferAnchors(reference = ExNeu_Velmeshev,
                        query = Integrate_raw_df[[i]],
                        dims = 1:30,
                        # normalization.method = "LogNormalize",
                        reference.reduction = "harmony")
}

#####
plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[2]],
           reference = ExNeu_Velmeshev,
           query = Integrate_raw_df[[2]],
           refdata = list(celltype = "age_range"),
           reference.reduction = "harmony",
           reduction.model = "umap")
length(unique(plot_query$predicted.celltype))
unique(plot_query$predicted.celltype)
plot_query$predicted.celltype <-
  factor(plot_query$predicted.celltype,
         levels = c("2nd trimester",
                    "3rd trimester",
                    "0-1 years",
                    "1-2 years",
                    "2-4 years",
                    "4-10 years",
                    "10-20 years",
                    "Adult"))
DimPlot(plot_query,
        reduction = "ref.umap",
        group.by = "predicted.celltype",
        cols = plasma(n = length(unique(ExNeu_Velmeshev$age_range))),
        alpha = 0.5,
        pt.size = 0.1) +
  ggtitle(names(cell_3_types_2K_each)[2])
table(plot_query$predicted.celltype)

cell_prediction_count_list <-
  vector(mode = "list",
         length = 3L)
names(cell_prediction_count_list) <-
  names(cell_3_types_2K_each)
cell_prediction_count_list[[2]] <-
  as.data.frame(table(plot_query$predicted.celltype))


#####
plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[3]],
           reference = ExNeu_Velmeshev,
           query = Integrate_raw_df[[3]],
           refdata = list(celltype = "age_range"),
           reference.reduction = "harmony",
           reduction.model = "umap")
length(unique(plot_query$predicted.celltype))
unique(plot_query$predicted.celltype)
plot_query$predicted.celltype <-
  factor(plot_query$predicted.celltype,
         levels = c("2nd trimester",
                    "3rd trimester",
                    "0-1 years",
                    "1-2 years",
                    "2-4 years",
                    "4-10 years",
                    "10-20 years",
                    "Adult"))
DimPlot(plot_query,
        reduction = "ref.umap",
        group.by = "predicted.celltype",
        cols = plasma(n = length(unique(ExNeu_Velmeshev$age_range))),
        alpha = 0.5,
        pt.size = 0.1) +
  ggtitle(names(cell_3_types_2K_each)[3])
table(plot_query$predicted.celltype)

cell_prediction_count_list[[3]] <-
  as.data.frame(table(plot_query$predicted.celltype))


### Prepare bar plot #####
cell_prediction_count_list[[1]]$Cell_type <- "GABA"
cell_prediction_count_list[[2]]$Cell_type <- "nmGlut"
cell_prediction_count_list[[3]]$Cell_type <- "npGlut"

df_bar_plot <-
  as.data.frame(rbind(cell_prediction_count_list[[1]],
                      cell_prediction_count_list[[2]],
                      cell_prediction_count_list[[3]]))
colnames(df_bar_plot)[1] <-
  "Stage"

df_bar_plot$Stage <-
  factor(df_bar_plot$Stage)

ggplot(data = df_bar_plot,
       aes(x = Cell_type,
           y = Freq,
           fill = Stage)) +
  geom_bar(position = "fill",
           stat = "identity",
           width = 0.5) +
  scale_fill_manual(values = brewer.pal(n = 8,
                                        name = "Dark2")) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(fill = "Predicted neurons at stage",
       x = "Cell type",
       y = "Percentage") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"))

