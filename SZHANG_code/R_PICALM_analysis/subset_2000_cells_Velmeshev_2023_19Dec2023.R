# 10 Nov 2023 Siwei
# sample GABA, nmglut, and npglut neurons
# Use 2000 cells per type each
# project 2000 cells, do not integrate

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
  library(gridExtra)
  library(viridis)
}

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
options(Seurat.object.assay.version = "v3")
options(future.seed = T)
set.seed(42)

## subset scARC_Velmeshev
scARC_Velmeshev_2023 <-
  readRDS("~/Data/Databases/scARC_prenatal_postnatal_human_cort_dev_Velmeshev_2023/local.rds")

load("~/backuped_space/Siwei_misc_R_projects/Alena_RNASeq_23Aug2023/ENSG_gene_index.RData")
Idents(scARC_Velmeshev_2023) <- 'cell_type'

DimPlot(scARC_Velmeshev_2023,
        reduction = "umap",
        group.by = "cell_type")

neuron_subsetted_scARC_Velmeshev_2023 <-
  subset(x = scARC_Velmeshev_2023,
         # idents = 'neural cell',
         downsample = 20000)

DimPlot(neuron_subsetted_scARC_Velmeshev_2023,
        reduction = "umap",
        group.by = "cell_type")

### assign dev stages
unique(neuron_subsetted_scARC_Velmeshev_2023$development_stage)


# assemble developmental stages
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified <-
  "prenatal"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "[123]\\-month\\-old.*")] <-
  "First_trimester"
unique(neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified)
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "[456]\\-month\\-old.*")] <-
  "Second_trimester"
unique(neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified)
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'fourth LMP')] <-
  "Second_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'fifth LMP')] <-
  "Second_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'sixth LMP')] <-
  "Second_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'seventh LMP')] <-
  "Third_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'eighth LMP')] <-
  "Third_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'ninth LMP')] <-
  "Third_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'under-1-year-old')] <-
  "1-3_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "newborn.*")] <-
  "1-3_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^1-year-old.*")] <-
  "1-3_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^2-year-old.*")] <-
  "1-3_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^3-year-old.*")] <-
  "1-3_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^[4-9]-year-old.*")] <-
  "4-9_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^1[0-9]-year-old.*")] <-
  "10-19_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^2[0-9]-year-old.*")] <-
  "20-29_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^3[0-9]-year-old.*")] <-
  "more_than_30_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^4[0-9]-year-old.*")] <-
  "more_than_30_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^5[0-9]-year-old.*")] <-
  "more_than_30_year_old"

neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified <-
  factor(neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified,
         levels = c("First_trimester",
                    "Second_trimester",
                    "Third_trimester",
                    "1-3_year_old",
                    "4-9_year_old",
                    "10-19_year_old",
                    "20-29_year_old",
                    "more_than_30_year_old"))

DimPlot(neuron_subsetted_scARC_Velmeshev_2023,
        reduction = "umap",
        group.by = "dev_stages_simplified",
        cols = plasma(n = length(unique(neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified))))

FeaturePlot(neuron_subsetted_scARC_Velmeshev_2023,
            features = c("ENSG00000091664", "ENSG00000128683"))

neuron_subsetted_scARC_Velmeshev_2023_backup <-
  neuron_subsetted_scARC_Velmeshev_2023

neuron_subsetted_scARC_Velmeshev_2023 <-
  NormalizeData(neuron_subsetted_scARC_Velmeshev_2023)
neuron_subsetted_scARC_Velmeshev_2023 <-
  FindVariableFeatures(neuron_subsetted_scARC_Velmeshev_2023,
                       selection.method = "mean.var.plot")
neuron_subsetted_scARC_Velmeshev_2023 <-
  ScaleData(neuron_subsetted_scARC_Velmeshev_2023,
            features = VariableFeatures(neuron_subsetted_scARC_Velmeshev_2023))
neuron_subsetted_scARC_Velmeshev_2023 <-
  RunPCA(neuron_subsetted_scARC_Velmeshev_2023,
         features = VariableFeatures(neuron_subsetted_scARC_Velmeshev_2023),
         verbose = T)
ElbowPlot(object = neuron_subsetted_scARC_Velmeshev_2023,
          ndims = 50)
neuron_subsetted_scARC_Velmeshev_2023 <-
  RunUMAP(neuron_subsetted_scARC_Velmeshev_2023,
          dims = 1:40,
          reduction = "pca",
          return.model = T)
neuron_subsetted_scARC_Velmeshev_2023 <-
  FindNeighbors(neuron_subsetted_scARC_Velmeshev_2023,
                reduction = "pca",
                dims = 1:40) %>%
  FindClusters()

DimPlot(neuron_subsetted_scARC_Velmeshev_2023,
        reduction = "umap",
        group.by = "dev_stages_simplified",
        cols = plasma(n = length(unique(neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified))))
FeaturePlot(neuron_subsetted_scARC_Velmeshev_2023,
            features = c("ENSG00000091664", "ENSG00000128683"))

# #####
# SCT_neuron_subsetted_scARC_Velmeshev_2023 <-
#   SCTransform(neuron_subsetted_scARC_Velmeshev_2023)
# var_neuron_subsetted_scARC_Velmeshev_2023 <-
#   FindVariableFeatures(neuron_subsetted_scARC_Velmeshev_2023,
#                        selection.method = "vst")

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
  x <- cell_3_types_2K_each[[i]]

  shared_gene_list <-
    merge(x = data.frame(Gene_Symbol = rownames(x)),
          y = ENSG_anno_gene_indexed,
          by = "Gene_Symbol")
  shared_gene_list <-
    shared_gene_list[!duplicated(shared_gene_list$Geneid), ]
  shared_gene_list <-
    shared_gene_list[!duplicated(shared_gene_list$Gene_Symbol), ]

  x <-
    x[rownames(x) %in% shared_gene_list$Gene_Symbol, ]
  x@assays$RNA@counts@Dimnames[[1]] <-
    shared_gene_list$Geneid
  x@assays$RNA@data@Dimnames[[1]] <-
    shared_gene_list$Geneid
  # x_backup <- x
  print("!")
  new_raw_x <-
    CreateSeuratObject(counts = as.matrix(x@assays$RNA@counts),
                       assay = "RNA",
                       meta.data = x@meta.data,
                       project = names(cell_3_types_2K_each)[i])
  print("!!")


  new_raw_x <- NormalizeData(new_raw_x,
                             normalization.method = "LogNormalize")
  new_raw_x <- ScaleData(new_raw_x)
  new_raw_x <-
    FindVariableFeatures(new_raw_x)
  # new_raw_x <- SCTransform(new_raw_x)
  Integrate_raw_df[[i]] <- new_raw_x

  transfer_anchor_list[[i]] <-
    FindTransferAnchors(reference = var_neuron_subsetted_scARC_Velmeshev_2023,
                        query = new_raw_x,
                        dims = 1:30,
                        normalization.method = "LogNormalize",
                        reference.reduction = "pca")
}

names(Integrate_raw_df) <-
  names(cell_3_types_2K_each)
names(transfer_anchor_list) <-
  names(cell_3_types_2K_each)





## preprocess Velmeshev #####
gene_list_Velmeshev_2023 <-
  merge(x = data.frame(Geneid = rownames(neuron_subsetted_scARC_Velmeshev_2023)),
        y = ENSG_anno_gene_indexed,
        by = "Geneid")
gene_list_Velmeshev_2023 <-
  gene_list_Velmeshev_2023[!duplicated(gene_list_Velmeshev_2023$Gene_Symbol), ]

neuron_subsetted_scARC_Velmeshev_2023 <-
  neuron_subsetted_scARC_Velmeshev_2023[rownames(neuron_subsetted_scARC_Velmeshev_2023) %in%
                                          gene_list_Velmeshev_2023$Geneid, ]
neuron_subsetted_scARC_Velmeshev_2023@assays$RNA@counts@Dimnames[[1]] <-
  gene_list_Velmeshev_2023$Gene_Symbol
neuron_subsetted_scARC_Velmeshev_2023@assays$RNA@data@Dimnames[[1]] <-
  gene_list_Velmeshev_2023$Gene_Symbol
rownames(neuron_subsetted_scARC_Velmeshev_2023)
neuron_subsetted_scARC_Velmeshev_2023_backup <-
  neuron_subsetted_scARC_Velmeshev_2023

# save.image("neuron_query_projection_Velmeshev_2023.RData")
# load("neuron_query_projection_Velmeshev_2023.RData")
# test make a Seurat object de novo
neuron_subsetted_raw <-
  CreateSeuratObject(counts = neuron_subsetted_scARC_Velmeshev_2023_backup@assays$RNA@counts,
                     data = neuron_subsetted_scARC_Velmeshev_2023_backup@assays$RNA@data,
                     project = "Velmeshev_2023_ref",
                     meta.data = neuron_subsetted_scARC_Velmeshev_2023_backup@meta.data)

neuron_subsetted_raw@assays$RNA@layers$counts@Dimnames[[1]] <-
  neuron_subsetted_raw@assays$RNA@features[[1]]
neuron_subsetted_raw@assays$RNA@layers$counts@Dimnames[[2]] <-
  neuron_subsetted_raw@assays$RNA@cells[[1]]

neuron_subsetted_raw@assays$RNA@layers$data@Dimnames[[1]] <-
  neuron_subsetted_raw@assays$RNA@features[[1]]
neuron_subsetted_raw@assays$RNA@layers$data@Dimnames[[2]] <-
  neuron_subsetted_raw@assays$RNA@cells[[1]]

# Run SCTransform
neuron_subsetted_scARC_Velmeshev_2023 <-
  PercentageFeatureSet(neuron_subsetted_raw,
                       pattern = "^MT-",
                       col.name = "percent.mt")
neuron_subsetted_scARC_Velmeshev_2023 <-
  SCTransform(neuron_subsetted_scARC_Velmeshev_2023,
              vars.to.regress = "percent.mt",
              verbose = T,
              seed.use = 42)

# neuron_subsetted_scARC_Velmeshev_2023 <-
#   NormalizeData(neuron_subsetted_raw,
#                 verbose = T)
# # neuron_subsetted_scARC_Velmeshev_2023@assays$RNA@scale.data <-
# #   neuron_subsetted_scARC_Velmeshev_2023_backup@assays$RNA@scale.data
# neuron_subsetted_scARC_Velmeshev_2023 <-
#   FindVariableFeatures(neuron_subsetted_scARC_Velmeshev_2023,
#                        verbose = T)
# neuron_subsetted_scARC_Velmeshev_2023 <-
#   ScaleData(neuron_subsetted_scARC_Velmeshev_2023,
#             verbose = T)
neuron_subsetted_scARC_Velmeshev_2023 <-
  RunPCA(neuron_subsetted_scARC_Velmeshev_2023,
         assay = "SCT",
         verbose = T)
neuron_subsetted_scARC_Velmeshev_2023 <-
  RunUMAP(neuron_subsetted_scARC_Velmeshev_2023,
          dims = 1:30,
          reduction = "pca",
          assay = "SCT",
          return.model = T,
          # reduction.name = "umap.original",
         verbose = T)
neuron_subsetted_scARC_Velmeshev_2023 <-
  FindNeighbors(neuron_subsetted_scARC_Velmeshev_2023,
                verbose = T)
neuron_subsetted_scARC_Velmeshev_2023 <-
  FindClusters(neuron_subsetted_scARC_Velmeshev_2023,
               verbose = T)

DimPlot(neuron_subsetted_scARC_Velmeshev_2023,
        group.by = "development_stage",
        label = T,
        repel = T)
DefaultAssay(neuron_subsetted_scARC_Velmeshev_2023)
Assays(neuron_subsetted_scARC_Velmeshev_2023)
# neuron_subsetted_scARC_Velmeshev_2023 <-
#   IntegrateLayers(object = neuron_subsetted_scARC_Velmeshev_2023,
#                   method = RPCAIntegration,
#                   normalization.method = "SCT",
#                   # assay = "SCT",
#                   # orig.reduction = "pca",
#                   # new.reduction = "harmony",
#                   verbose = T)

saveRDS(neuron_subsetted_scARC_Velmeshev_2023,
        file = "neuron_subsetted_scARC_Velmeshev_2023_assay5.RDs")

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

## !! Need to run SCTModel on the original object #####


# Integrate_raw_df <-
#   cell_3_types_2K_each
for (i in 1:length(cell_3_types_2K_each)) {
  x <- cell_3_types_2K_each[[i]]
  # x <- CreateSeuratObject(counts = x@assays$RNA@counts,
  #                         data = x@assays$RNA@data,
  #                         # project = "Velmeshev_2023_ref",
  #                         meta.data = x@meta.data)
  x <-
    x[rownames(x) %in%
        ENSG_anno_gene_indexed$Gene_Symbol, ]
  # rownames(x) <-
  #   ENSG_anno_gene_indexed$Geneid[ENSG_anno_gene_indexed$Gene_Symbol %in%
  #                                   rownames(x)]
  # x@assays$RNA@layers$counts@Dimnames[[1]] <-
  #   ENSG_anno_gene_indexed$Geneid[match(x = rownames(x),
  #                                       table = ENSG_anno_gene_indexed$Gene_Symbol)]
  #   # ENSG_anno_gene_indexed$Geneid[ENSG_anno_gene_indexed$Gene_Symbol %in%
  #   #                                 rownames(x)]
  # x@assays$RNA@layers$counts@Dimnames[[1]] <-
  #   ENSG_anno_gene_indexed$Geneid[match(x = rownames(x),
  #                                       table = ENSG_anno_gene_indexed$Gene_Symbol)]

  x@assays$RNA@counts@Dimnames[[1]] <-
    ENSG_anno_gene_indexed$Geneid[match(x = rownames(x),
                                        table = ENSG_anno_gene_indexed$Gene_Symbol)]
  # ENSG_anno_gene_indexed$Geneid[ENSG_anno_gene_indexed$Gene_Symbol %in%
  #                                 rownames(x)]
  # x@assays$RNA@layers$counts@Dimnames[[1]] <-
  #   ENSG_anno_gene_indexed$Geneid[match(x = rownames(x),
  #                                       table = ENSG_anno_gene_indexed$Gene_Symbol)]

      # ENSG_anno_gene_indexed$Geneid[ENSG_anno_gene_indexed$Gene_Symbol %in%
    #                                 rownames(x)]
  # x <- x[rownames(x) %in% gene_list_Velmeshev_2023$Gene_Symbol, ]
  # x <- ScaleData(x)
  x <- NormalizeData(x)
  Integrate_raw_df[[i]] <- x

  transfer_anchor_list[[i]] <-
    FindTransferAnchors(reference = neuron_subsetted_scARC_Velmeshev_2023,
                        query = x,
                        dims = 1:30,
                        # normalization.method = "SCT",
                        reference.reduction = "pca")
}

  length(ENSG_anno_gene_indexed$Geneid[ENSG_anno_gene_indexed$Gene_Symbol %in%
                                         rownames(x)])
length(!duplicated(ENSG_anno_gene_indexed$Geneid[ENSG_anno_gene_indexed$Gene_Symbol %in%
                                                   rownames(x)]))
length(!duplicated(ENSG_anno_gene_indexed$Gene_Symbol[ENSG_anno_gene_indexed$Gene_Symbol %in%
                                                   rownames(x)]))
sum(ENSG_anno_gene_indexed$Gene_Symbol %in%
      rownames(x))
length(ENSG_anno_gene_indexed$Gene_Symbol %in%
      rownames(x))
length(!duplicated(ENSG_anno_gene_indexed$Gene_Symbol))
length(!duplicated(ENSG_anno_gene_indexed$Geneid))
length(!duplicated(rownames(x)))

ENSG_anno_gene_indexed$Geneid[match(x = rownames(x),
       table = ENSG_anno_gene_indexed$Gene_Symbol)]
length(ENSG_anno_gene_indexed$Geneid[match(x = rownames(x),
                                           table = ENSG_anno_gene_indexed$Gene_Symbol)])

names(Integrate_raw_df) <-
  names(cell_3_types_2K_each)
names(transfer_anchor_list) <-
  names(cell_3_types_2K_each)
#
# Integrate_raw_df <-
#   lapply(X = Integrate_raw_df,
#          FUN = function(x) {
#            print(names(x))
#            x <- CreateSeuratObject(counts = x@assays$RNA@counts,
#                                     data = x@assays$RNA@data,
#                                     # project = "Velmeshev_2023_ref",
#                                     meta.data = x@meta.data)
#            x <- x[rownames(x) %in% gene_list_Velmeshev_2023$Gene_Symbol, ]
#            # x <- ScaleData(x)
#            x <- SCTransform(x)
#            # x <- FindVariableFeatures(object = x,
#            #                           selection.method = "vst",
#            #                           nfeatures = 3000)
#            x[["transfer_anchors"]] <-
#              FindTransferAnchors(reference = neuron_subsetted_scARC_Velmeshev_2023,
#                                  query = x,
#                                  dims = 1:30,
#                                  reference.reduction = "pca")
#          })

neuron_subsetted_scARC_Velmeshev_2023 <-
  RunUMAP(neuron_subsetted_scARC_Velmeshev_2023,
          dims = 1:30,
          reduction = "pca",
          return.model = T)
plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[1]],
           reference = neuron_subsetted_scARC_Velmeshev_2023,
           query = Integrate_raw_df[[1]],
           refdata = list(celltype = "development_stage"),
           reference.reduction = "pca",
           reduction.model = "umap")
plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[1]],
           reference = scARC_Velmeshev_2023,
           query = Integrate_raw_df[[1]],
           refdata = list(celltype = "development_stage"),
           reference.reduction = "pca",
           reduction.model = "umap")

plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[2]],
           reference = neuron_subsetted_scARC_Velmeshev_2023,
           query = Integrate_raw_df[[2]],
           refdata = list(celltype = "development_stage"),
           reference.reduction = "pca",
           reduction.model = "umap")

plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[3]],
           reference = neuron_subsetted_scARC_Velmeshev_2023,
           query = Integrate_raw_df[[3]],
           refdata = list(celltype = "development_stage"),
           reference.reduction = "pca",
           reduction.model = "umap")

DimPlot(plot_query,
        reduction = "ref.umap",
        group.by = "predicted.celltype")

DimPlot(neuron_subsetted_scARC_Velmeshev_2023,reduction = "pca",
        group.by = "seurat_clusters")
FeaturePlot(neuron_subsetted_scARC_Velmeshev_2023,
            features = "SLC17A7")
Assays(neuron_subsetted_scARC_Velmeshev_2023)
DefaultAssay(neuron_subsetted_scARC_Velmeshev_2023)
DefaultAssay(neuron_subsetted_scARC_Velmeshev_2023) <- "RNA"
FeaturePlot(neuron_subsetted_scARC_Velmeshev_2023,
            features = "TLE4")


## SCT on each of the list items ####
Integrate_raw_df <-
  vector(mode = "list",
         length = 3L)
transfer_anchor_list <-
  vector(mode = "list",
         length = 3L)

# Integrate_raw_df <-
#   cell_3_types_2K_each
for (i in 1:length(cell_3_types_2K_each)) {
  x <- cell_3_types_2K_each[[i]]
  x <- CreateSeuratObject(counts = x@assays$RNA@counts,
                          data = x@assays$RNA@data,
                          # project = "Velmeshev_2023_ref",
                          meta.data = x@meta.data)
  x <- x[rownames(x) %in% gene_list_Velmeshev_2023$Gene_Symbol, ]
  # x <- ScaleData(x)
  x <- NormalizeData(x)
  Integrate_raw_df[[i]] <- x

  transfer_anchor_list[[i]] <-
    FindTransferAnchors(reference = neuron_subsetted_scARC_Velmeshev_2023,
                        query = x,
                        dims = 1:30,
                        reference.reduction = "pca")
}

names(Integrate_raw_df) <-
  names(cell_3_types_2K_each)
names(transfer_anchor_list) <-
  names(cell_3_types_2K_each)


## !! Need to run SCTModel on the original object #####
# save.image("to_run_SCTModel_10Nov2023.RData")

neuron_subsetted_scARC_Velmeshev_2023 <-
  RunUMAP(neuron_subsetted_scARC_Velmeshev_2023,
          dims = 1:30,
          reduction = "pca",
          return.model = T)
DimPlot(neuron_subsetted_scARC_Velmeshev_2023)
unique(neuron_subsetted_scARC_Velmeshev_2023$development_stage)

# assemble developmental stages
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified <-
  "prenatal"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "[123]\\-month\\-old.*")] <-
  "First_trimester"
unique(neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified)
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "[456]\\-month\\-old.*")] <-
  "Second_trimester"
unique(neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified)
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'fourth LMP')] <-
  "Second_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'fifth LMP')] <-
  "Second_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'sixth LMP')] <-
  "Second_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'seventh LMP')] <-
  "Third_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'eighth LMP')] <-
  "Third_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'ninth LMP')] <-
  "Third_trimester"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = 'under-1-year-old')] <-
  "1-3_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "newborn.*")] <-
  "1-3_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^1-year-old.*")] <-
  "1-3_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^2-year-old.*")] <-
  "1-3_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^3-year-old.*")] <-
  "1-3_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^[4-9]-year-old.*")] <-
  "4-9_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^1[0-9]-year-old.*")] <-
  "10-19_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^2[0-9]-year-old.*")] <-
  "20-29_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^3[0-9]-year-old.*")] <-
  "more_than_30_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^4[0-9]-year-old.*")] <-
  "more_than_30_year_old"
neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified[str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                                       pattern = "^5[0-9]-year-old.*")] <-
  "more_than_30_year_old"
DimPlot(neuron_subsetted_scARC_Velmeshev_2023,
        group.by = "dev_stages_simplified")

#####
plot_query <-
  NormalizeData(cell_3_types_2K_each[[1]])
plot_query <-
  FindTransferAnchors(reference = neuron_subsetted_scARC_Velmeshev_2023,
                      query = plot_query,
                      reference.reduction = "pca",
                      dims = 1:30,
                      verbose = T)


plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[1]],
           reference = neuron_subsetted_scARC_Velmeshev_2023,
           query = Integrate_raw_df[[1]],
           refdata = list(celltype = "dev_stages_simplified"),
           reference.reduction = "pca",
           reduction.model = "umap")
plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[1]],
           reference = scARC_Velmeshev_2023,
           query = Integrate_raw_df[[1]],
           refdata = list(celltype = "dev_stages_simplified"),
           reference.reduction = "pca",
           reduction.model = "umap")

plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[2]],
           reference = neuron_subsetted_scARC_Velmeshev_2023,
           query = Integrate_raw_df[[2]],
           refdata = list(celltype = "development_stage"),
           reference.reduction = "pca",
           reduction.model = "umap")

plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[3]],
           reference = neuron_subsetted_scARC_Velmeshev_2023,
           query = Integrate_raw_df[[3]],
           refdata = list(celltype = "development_stage"),
           reference.reduction = "pca",
           reduction.model = "umap")

DimPlot(plot_query,
        reduction = "ref.umap",
        group.by = "predicted.celltype")

DimPlot(neuron_subsetted_scARC_Velmeshev_2023,reduction = "pca",
        group.by = "seurat_clusters")
FeaturePlot(neuron_subsetted_scARC_Velmeshev_2023,
            features = "SLC17A7")
Assays(neuron_subsetted_scARC_Velmeshev_2023)
DefaultAssay(neuron_subsetted_scARC_Velmeshev_2023)
DefaultAssay(neuron_subsetted_scARC_Velmeshev_2023) <- "RNA"
FeaturePlot(neuron_subsetted_scARC_Velmeshev_2023,
            features = "TLE4")
