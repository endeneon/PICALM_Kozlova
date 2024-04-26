# 18 Dec 2023 Siwei
# Relabel Velmeshev 2023 cells and convert to Assay5
# project 20000 cells, do not integrate

# init ####
{
  library(Seurat)
  # library(Signac)
  library(readr)
  library(future)
  library(parallel)
  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
  library(viridis)
  library(gridExtra)
  library(ggpubr)
}

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
options(Seurat.object.assay.version = "v3")
options(future.seed = T)
set.seed(42)


# saveRDS(neuron_subsetted_scARC_Velmeshev_2023,
#         file = "neuron_subsetted_scARC_Velmeshev_2023_assay5.RDs")

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

## convert neuron_subsetted to v3 object
v3_neuron_subsetted_scARC_Velmeshev_2023 <-
  neuron_subsetted_scARC_Velmeshev_2023

v3_neuron_subsetted_scARC_Velmeshev_2023[["RNA"]] <-
  as(object = neuron_subsetted_scARC_Velmeshev_2023[["RNA"]],
     Class = "Assay")

# process cells_3_types_2K
cell_3_types_2K_each <-
  readRDS("cell_3_types_from_76_lines_2K_each_4_Velmeshev.RDs")

# v3_neuron_subsetted_scARC_Velmeshev_2023_backup <-
#   v3_neuron_subsetted_scARC_Velmeshev_2023
v3_neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified <-
  factor(v3_neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified,
         levels = c("First_trimester",
                    "Second_trimester",
                    "Third_trimester",
                    "1-3_year_old",
                    "4-9_year_old",
                    "10-19_year_old",
                    "20-29_year_old",
                    "more_than_30_year_old"))


# [1] "10-19_year_old"        "20-29_year_old"        "4-9_year_old"
# [4] "more_than_30_year_old" "1-3_year_old"          "Third_trimester"
# [7] "First_trimester"       "Second_trimester"
## SCT on each of the list items ####
Integrate_raw_df <-
  vector(mode = "list",
         length = 3L)
transfer_anchor_list <-
  vector(mode = "list",
         length = 3L)

for (i in 1:length(cell_3_types_2K_each)) {
  x <- cell_3_types_2K_each[[i]]
  # x <- CreateSeuratObject(counts = x@assays$RNA@counts,
  #                         data = x@assays$RNA@data,
  #                         # project = "Velmeshev_2023_ref",
  #                         meta.data = x@meta.data)
  # x <-
  #   x[rownames(x) %in%
  #       ENSG_anno_gene_indexed$Gene_Symbol, ]
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

  # x@assays$RNA@counts@Dimnames[[1]] <-
  #   ENSG_anno_gene_indexed$Geneid[match(x = rownames(x),
  #                                       table = ENSG_anno_gene_indexed$Gene_Symbol)]
  # # ENSG_anno_gene_indexed$Geneid[ENSG_anno_gene_indexed$Gene_Symbol %in%
  # #                                 rownames(x)]
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
    FindTransferAnchors(reference = v3_neuron_subsetted_scARC_Velmeshev_2023,
                        query = x,
                        dims = 1:30,
                        normalization.method = "SCT",
                        reference.reduction = "pca")
}

names(Integrate_raw_df) <-
  names(cell_3_types_2K_each)
names(transfer_anchor_list) <-
  names(cell_3_types_2K_each)

# v3_neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified <-
#   factor(v3_neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified,
#          levels = c("First_trimester", #1
#                     "Second_trimester", #2
#                     "Third_trimester", #3
#                     "1-3_year_old", #4
#                     "4-9_year_old", #5
#                     "10-19_year_old", #6
#                     "20-29_year_old", #7
#                     "more_than_30_year_old")) #8

plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[1]],
           reference = v3_neuron_subsetted_scARC_Velmeshev_2023,
           query = Integrate_raw_df[[1]],
           refdata = list(celltype = "dev_stages_simplified"),
           reference.reduction = "pca",
           reduction.model = "umap")
plot_query$predicted.celltype <-
  factor(plot_query$predicted.celltype,
         levels = c("Second_trimester",
                    "1-3_year_old",
                    "4-9_year_old",
                    "10-19_year_old"))
plot_GABA <-
  DimPlot(plot_query,
          reduction = "ref.umap",
          group.by = "predicted.celltype",
          cols = plasma(n = length(unique(v3_neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified)))[c(2,4,5,6)]) +
  ggtitle(names(cell_3_types_2K_each)[1])


plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[2]],
           reference = v3_neuron_subsetted_scARC_Velmeshev_2023,
           query = Integrate_raw_df[[2]],
           refdata = list(celltype = "dev_stages_simplified"),
           reference.reduction = "pca",
           reduction.model = "umap")
plot_query$predicted.celltype <-
  factor(plot_query$predicted.celltype,
         levels = c("Second_trimester",
                    "Third_trimester",
                    "1-3_year_old",
                    "4-9_year_old",
                    "10-19_year_old"))
plot_nmglut <-
  DimPlot(plot_query,
          reduction = "ref.umap",
          group.by = "predicted.celltype",
          cols = plasma(n = length(unique(v3_neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified)))[c(2,3,4,5,6)]) +
  ggtitle(names(cell_3_types_2K_each)[2])

plot_query <-
  MapQuery(anchorset = transfer_anchor_list[[3]],
           reference = v3_neuron_subsetted_scARC_Velmeshev_2023,
           query = Integrate_raw_df[[3]],
           refdata = list(celltype = "dev_stages_simplified"),
           reference.reduction = "pca",
           reduction.model = "umap")
plot_query$predicted.celltype <-
  factor(plot_query$predicted.celltype,
         levels = c("Second_trimester",
                    "Third_trimester",
                    "1-3_year_old",
                    # "4-9_year_old",
                    "10-19_year_old"))
plot_npglut <-
  DimPlot(plot_query,
          reduction = "ref.umap",
          group.by = "predicted.celltype",
          cols = plasma(n = length(unique(v3_neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified)))[c(2,3,4,6)]) +
  ggtitle(names(cell_3_types_2K_each)[3])



DimPlot(plot_query,
        reduction = "ref.umap",
        group.by = "predicted.celltype")


# make reference plot
plot_ref <-
  DimPlot(v3_neuron_subsetted_scARC_Velmeshev_2023,
          reduction = "umap",
          group.by = "dev_stages_simplified",
          cols = plasma(n = length(unique(v3_neuron_subsetted_scARC_Velmeshev_2023$dev_stages_simplified)))) +
  ggtitle("Reference Developmental Stages")

# save.image("Ready_to_plot_18Dec2023.RData")

# plot_all <-

ggarrange(plot_ref,
            plot_GABA,
            plot_nmglut,
            plot_npglut,
          ncol = 2,
          nrow = 2)

grid.arrange(plot_ref,
            plot_GABA,
            plot_nmglut,
            plot_npglut,
          ncol = 2)


FeaturePlot(v3_neuron_subsetted_scARC_Velmeshev_2023,
            features = c("SLC17A6", "GAD2"))
FeaturePlot(v3_neuron_subsetted_scARC_Velmeshev_2023,
            features = c("PAX6", "TLE4"))
FeaturePlot(v3_neuron_subsetted_scARC_Velmeshev_2023,
            features = c("VIM", "CUX2",
                         "SLC17A6", "GAD2"))
