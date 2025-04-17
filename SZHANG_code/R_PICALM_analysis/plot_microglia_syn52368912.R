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

Immune_cells_meta <- # has projid and cell_type_high_resolution only
  Immune_cells@meta.data

brainsel_before_cluster_subset <-
  read_rds(file = "brainsel_before_cluster_subset.RDs")
brainsel_meta <- # index:subject, ROSMAP-xxxxx
  brainsel_before_cluster_subset@meta.data
rm(brainsel_before_cluster_subset)

# load ROSMAP lookup matrices ####
lookup_indivID_2_sunID <- # Rxxxxxxx - ROSMAP-xxxxx
  read.csv(file = "human_microglia_state_dynamics/metadata/ROSMAP_meta_syn52293430/MIT_ROSMAP_Multiomics_individual_metadata.csv")
lookup_indivID_2_sunID <-
  lookup_indivID_2_sunID[, c(1:4, 19)]
lookup_indivID_2_sunID <-
  lookup_indivID_2_sunID[!duplicated(lookup_indivID_2_sunID$subject), ]
lookup_indivID_2_sunID <-
  lookup_indivID_2_sunID[!duplicated(lookup_indivID_2_sunID$individualID), ]

# projid, from clinical
lookup_indivID_2_projID <- # projid - individualID
  read.csv(file = "human_microglia_state_dynamics/metadata/ROSMAP_meta_all/ROSMAP_clinical.csv")
lookup_indivID_2_projID <-
  lookup_indivID_2_projID[!(duplicated(lookup_indivID_2_projID$projid)), ]
lookup_indivID_2_projID <-
  lookup_indivID_2_projID[!(duplicated(lookup_indivID_2_projID$individualID)), ]

lookup_indivID_projID_sunID <-
  merge(lookup_indivID_2_sunID,
        lookup_indivID_2_projID,
        by = "individualID") # individualID: Rxxxxxxx

# sum(is.na(lookup_indivID_projID_sunID$))

Immune_cells_P2RY12 <-
  Immune_cells[, Immune_cells$cell_type_high_resolution == 'Mic P2RY12']
Immune_cells_P2RY12_meta <-
  Immune_cells_P2RY12@meta.data ## projid xxxxxxxx

Immune_cells_P2RY12_meta_diag <-
  merge(x = Immune_cells_P2RY12_meta[!duplicated(Immune_cells_P2RY12_meta$projid), ],
        y = lookup_indivID_projID_sunID,
        by.x = "projid",
        by.y = "projid",
        all.x = T)
# sum(is.na(Immune_cells_P2RY12_meta_diag$))

brainsel_meta_subject <-
  brainsel_meta[!duplicated(brainsel_meta$subject), ]

Immune_cells_P2RY12_meta_diag <-
  merge(x = Immune_cells_P2RY12_meta_diag,
        y = brainsel_meta_subject,
        by.x = "subject",
        by.y = "subject")
sum(is.na(Immune_cells_P2RY12_meta_diag$ADdiag3types))
Immune_cells_P2RY12_meta_diag$g_projid <-
  paste0("g",
         Immune_cells_P2RY12_meta_diag$projid)
length(unique(Immune_cells_P2RY12_meta_diag$g_projid))

samples_count_less_than_50 <-
  read.table("projid_less_than_50.txt")
samples_count_less_than_50 <-
  unlist(samples_count_less_than_50)

sum(!unique(Immune_cells_P2RY12_meta_diag$g_projid) %in%
      samples_count_less_than_50)

Immune_cells_P2RY12_meta_diag <-
  Immune_cells_P2RY12_meta_diag[!(Immune_cells_P2RY12_meta_diag$g_projid %in%
                                    samples_count_less_than_50), ]

# unsorted_counts <-
#   Immune_cells_P2RY12@assays$RNA@counts[, Immune_cells_P2RY12$projid %in%
#                                           unique(Immune_cells_P2RY12_meta_diag$projid)]
base_counts <-
  Immune_cells_P2RY12@assays$RNA@counts[, Immune_cells_P2RY12$projid %in%
                                          Immune_cells_P2RY12_meta_diag$projid]
base_meta <-
  Immune_cells_P2RY12@meta.data
base_meta <-
  base_meta[base_meta$projid %in% Immune_cells_P2RY12_meta_diag$projid, ]
base_meta <-
  merge(x = base_meta,
        y = Immune_cells_P2RY12_meta_diag,
        by = "projid",
        all.X = T)

Immune_cells_P2RY12_diag <-
  CreateSeuratObject(counts = base_counts,
                     project = "Immune_cells_P2RY12_diag",
                     meta.data = base_meta)

{
  Idents(Immune_cells_P2RY12_diag) <- "ADdiag3types"

  Immune_cells_P2RY12_diag <-
    Immune_cells_P2RY12_diag %>%
    NormalizeData(normalization.method = "LogNormalize") %>%
    FindVariableFeatures(selection.method = "vst",
                         nfeatures = 2000) %>%
    ScaleData()

  Idents(Immune_cells_P2RY12_diag) <- "batch"

}


{

  Immune_cells_P2RY12_diag_harmony <-
    Immune_cells_P2RY12_diag %>%
    RunPCA(features = VariableFeatures(Immune_cells_P2RY12_diag)) %>%
    harmony:::RunHarmony.Seurat(group.by.vars = "batch",
               plot_convergence = T) %>%
    # Embeddings('harmony') #%>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = 0.5) %>%
    identity()



  Idents(Immune_cells_P2RY12_diag_harmony) <- "ADdiag3types"

  Immune_cells_MAST_early <-
    FindMarkers(Immune_cells_P2RY12_diag_harmony,
                ident.1 = "earlyAD",
                ident.2 = "nonAD",
                logfc.threshold = 0,
                test.use = "MAST",
                verbose = T)

  Immune_cells_MAST_late <-
    FindMarkers(Immune_cells_P2RY12_diag_harmony,
                ident.1 = "lateAD",
                ident.2 = "nonAD",
                logfc.threshold = 0,
                test.use = "MAST",
                verbose = T)
}



hist(sort(table(Immune_cells_P2RY12_meta_diag$subject)), breaks = 20)

# Immune_cells_

Immune_cells_P2RY12_pseudobulk <-
  AggregateExpression(Immune_cells_P2RY12_diag,
                      assays = "RNA",
                      return.seurat = T,
                      group.by = "projid",
                      # scale.factor = 1e6,
                      normalization.method = "LogNormalize")
# Immune_cells_P2RY12_pseudobulk <-
#   Seurat:::PseudobulkExpression(Immune_cells_P2RY12_diag,
#                       assays = "RNA",
#                       return.seurat = T,
#                       group.by = "projid",
#                       method = "agg",
#                       normalization.method = "LogNormalize")



# df_2_writeout <-
#   data.frame(individualID = Immune_cells_P2RY12_pseudobulk$individualID,
#              cell_count = table(Immune_cells_P2RY12_pseudobulk$))


# try edgeR #####
Immune_cells_df_2_plot <-
  as.data.frame(as.matrix(Immune_cells_P2RY12_pseudobulk@assays$RNA@layers$counts))

# Immune_cells_df_2_plot <-
#   as.data.frame(as.matrix(Immune_cells_P2RY12_pseudobulk@assays$RNA@layers$scale.data))
# brainsel_df_2_plot <-
#   as.data.frame(as.matrix(brainsel_pseudobulk@assays$RNA@layers$scale.data))

rownames(Immune_cells_df_2_plot) <-
  rownames(Immune_cells_P2RY12_pseudobulk)
colnames(Immune_cells_df_2_plot) <-
  colnames(Immune_cells_P2RY12_pseudobulk)

#
# df_2_writeout <-
#   data.frame(PICALM_value = unlist(Immune_cells_df_2_plot[rownames(Immune_cells_df_2_plot) == "PICALM", ]),
#              projid = colnames(Immune_cells_df_2_plot),
#              cell_count = table(Immune_cells_P2RY12_meta_diag$projid))
#
# write.table(df_2_writeout,
#             file = "PICALM_P2RY12_normalised_CPM_projid_cellCount.tsv",
#             quote = F, sep = "\t",
#             row.names = F, col.names = T)
# head(df_2_writeout)
# head(table(Immune_cells_P2RY12_meta_diag$projid))

# Immune_cells_df_2_plot <-
#   Immune_cells_df_2_plot[, order(colnames(Immune_cells_df_2_plot))]
# head(colnames(Immune_cells_df_2_plot))

# Immune_cells_P2RY12_meta_diag_dedup <-
#   Immune_cells_P2RY12_meta_diag[!duplicated(Immune_cells_P2RY12_meta_diag_dedup$subject), ]
# Immune_cells_P2RY12_meta_diag_dedup$g_projid <-
#   paste0("g",
#          Immune_cells_P2RY12_meta_diag_dedup$projid)
# Immune_cells_P2RY12_meta_diag_dedup <-
#   Immune_cells_P2RY12_meta_diag_dedup[Immune_cells_P2RY12_meta_diag_dedup$g_projid %in% colnames(Immune_cells_df_2_plot), ]

# Immune_cells_P2RY12_meta_diag_dedup$g_projid <-
#   paste0("g",
#          Immune_cells_P2RY12_meta_diag_dedup$projid)



# Immune_cells_df_2_plot <-
nrow(Immune_cells_P2RY12_meta_diag)
Immune_cells_df_2_plot_meta <-
  Immune_cells_P2RY12_meta_diag[Immune_cells_P2RY12_meta_diag$projid %in% base_meta$projid, ]
  # base_meta[base_meta$]
Immune_cells_df_2_plot_meta <-
  Immune_cells_df_2_plot_meta[order(Immune_cells_df_2_plot_meta$projid), ]


Immune_cells_df_2_plot <-
  Immune_cells_df_2_plot[, order(colnames(Immune_cells_df_2_plot))]
#
#
# Immune_cells_P2RY12_meta_diag_dedup <-
#   Immune_cells_P2RY12_meta_diag_dedup[!(Immune_cells_P2RY12_meta_diag_dedup$g_projid %in%
#                                         samples_count_less_than_50), ]
# Immune_cells_P2RY12_meta_diag_dedup <-
#   Immune_cells_P2RY12_meta_diag_dedup[order(Immune_cells_P2RY12_meta_diag_dedup$g_projid), ]
#
# sum(colnames(Immune_cells_df_2_plot) %in% Immune_cells_P2RY12_meta_diag_dedup$g_projid)
# Immune_cells_df_2_plot <-
#   Immune_cells_df_2_plot[, (colnames(Immune_cells_df_2_plot) %in% Immune_cells_P2RY12_meta_diag_dedup$g_projid)]
# Immune_cells_df_2_plot <-
#   Immune_cells_df_2_plot[, order(colnames(Immune_cells_df_2_plot))]
# # Immune_cells_df_2_plot_agg <-
#   data.frame(samples = colnames(Immune_cells_df_2_plot),
#              aggregated_counts = colSums(Immune_cells_df_2_plot))
# write.table(Immune_cells_df_2_plot_agg,
#             file = "PICALM_P2RY12_aggCount.tsv",
#             quote = F, sep = "\t",
#             row.names = F, col.names = T)


pseudobulk_DGE <-
  DGEList(counts = as.matrix(Immune_cells_df_2_plot),
          genes = rownames(Immune_cells_df_2_plot),
          samples = colnames(Immune_cells_df_2_plot),
          group = Immune_cells_df_2_plot_meta$ADdiag3types,
          remove.zeros = T)

# length(Immune_cells_P2RY12_meta_diag_dedup$ADdiag3types[!(Immune_cells_P2RY12_meta_diag_dedup$g_projid %in% samples_count_less_than_50)])
#
# length(Immune_cells_P2RY12_meta_diag_dedup$ADdiag3types)
sum(is.na(as.matrix(pseudobulk_DGE$counts)))
ncol(pseudobulk_DGE)

sum(is.na(pseudobulk_DGE$samples$group))
GLM_model_df <-
  data.frame(AD3diag = factor(pseudobulk_DGE$samples$group))
GLM_model_df$AD3diag <-
  relevel(GLM_model_df$AD3diag,
          ref = "nonAD")
#
# pseudobulk_DGE <-
#   pseudobulk_DGE[!is.na(pseudobulk_DGE$counts), ]
# pseudobulk_DGE <-
#   pseudobulk_DGE[!is.na(rowSums(pseudobulk_DGE$counts)), ]
#
pseudobulk_DGE <-
  calcNormFactors(pseudobulk_DGE)
pseudobulk_DGE <-
  estimateDisp(pseudobulk_DGE)

nrow(stats::model.matrix(~ 0 + factor(AD3diag),
                         data = GLM_model_df))



# GLM_model_df <-
#   data.frame(AD3diag = factor(pseudobulk_DGE$samples$group))
# pseudobulk_DGE$samples$groupsum(is.na(as.matrix(Immune_cells_df_2_plot)))

QLF_test <-
  glmQLFit(pseudobulk_DGE,
           design = stats::model.matrix(~ 0 + factor(AD3diag),
                                        data = GLM_model_df))
QLF_test$coefficients
QLF_test_early <-
  glmQLFTest(QLF_test,
             contrast = c(-1, 1, 0))
View(QLF_test_early$table)
QLF_test_early$table[rownames(QLF_test_early$table) == "PICALM", ]

QLF_test_late <-
  glmQLFTest(QLF_test,
             contrast = c(-1, 0, 1))
View(QLF_test_late$table)
QLF_test_late$table[rownames(QLF_test_late$table) == "PICALM", ]


cpm_pseudobulk_DGE <-
  as.data.frame(edgeR::cpm(pseudobulk_DGE, log = T))
rownames(cpm_pseudobulk_DGE)
colnames(cpm_pseudobulk_DGE)

# df_writeout <-
#   data.frame(CPM = unlist(cpm_pseudobulk_DGE[rownames(cpm_pseudobulk_DGE) == "PICALM", ]),
#              projid = colnames(cpm_pseudobulk_DGE))
# write.table(df_writeout,
#             file = "PICALM_P2RY12_normalised_CPM_projid_cellCount_edgeR_more_than_50.tsv",
#             quote = F, sep = "\t",
#             row.names = F, col.names = T)
#
# write.table()

#
# cpm_pseudobulk_DGE <-
#   Immune_cells_P2RY12_pseudobulk@assays$RNA@layers$data
# rownames(cpm_pseudobulk_DGE) <-
#   rownames(Immune_cells_P2RY12_pseudobulk)
#
# colnames(cpm_pseudobulk_DGE) <-
#   colnames(Immune_cells_P2RY12_pseudobulk)

# cpm_pseudobulk_DGE <-
#   cpm_pseudobulk_DGE[, order(colnames(cpm_pseudobulk_DGE))]
# ncol(cpm_pseudobulk_DGE)

Sun_185_project_ids <-
  readRDS("~/backuped_space/Siwei_misc_R_projects/R_Alena_PICALM/Sun_185_project_ids.RDs")


df_2_plot <-
  data.frame(value = unlist(cpm_pseudobulk_DGE[rownames(cpm_pseudobulk_DGE) == "PICALM", ]),
             ADdiag3types = pseudobulk_DGE$samples$group)
df_2_plot <-
  df_2_plot[!is.na(df_2_plot$ADdiag3types), ]
# df_2_plot <-
#   df_2_plot[df_2_plot$projid %in% Sun_185_project_ids$bare_projid, ]
df_2_plot$ADdiag3types <-
  factor(df_2_plot$ADdiag3types,
         levels = c("nonAD",
                    "earlyAD",
                    "lateAD"))


ggerrorplot(df_2_plot,
            x = "ADdiag3types",
            y = "value",
            color = "ADdiag3types",
            # shape = "Clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = brewer.pal(n = 3,
                                          name = "Dark2"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = 0.95,
                     # label.x.npc = 0,
                     method = 'wilcox.test',
                     hide.ns = F,
                     ref.group = "nonAD",
                     paired = F) +
  geom_jitter(aes(colour = ADdiag3types),
              width = 0.1,
              size = 1) +
  # scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Normalized CPM") +
  scale_y_continuous(expand = c(0, 0.5)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("Sun et al. 425 individuals w/ genotype")

# glht_output <-
df_2_aov <-
  df_2_plot[!(is.nan(df_2_plot$value)), ]
df_2_aov <-
  df_2_plot[!(is.na(df_2_plot$value)), ]
df_2_aov <-
  df_2_plot[!(is.infinite(df_2_plot$value)), ]

summary(aov(value ~ ADdiag3types,
            data = df_2_aov))
summary(multcomp::glht(aov(value ~ ADdiag3types,
                           data = df_2_aov),
                       linfct = mcp(ADdiag3types = "Dunnett")))
table(ROSMAP.ImmuneCells.6regions.snRNAseq.meta$brainRegion)

Idents(ROSMAP.ImmuneCells.6regions.seurat.harmony) <-"RNA_snn_res.0.5"
UMAPPlot(ROSMAP.ImmuneCells.6regions.seurat.harmony)

table(ROSMAP.ImmuneCells.6regions.seurat.harmony$brainRegion)
FeaturePlot(ROSMAP.ImmuneCells.6regions.seurat.harmony,
            features = "P2RY12")
FeaturePlot(ROSMAP.ImmuneCells.6regions.seurat.harmony,
            features = "PICALM",
            alpha = 0.5)

UMAPPlot(Immune_cells)
table(Immune_cells$cell_type_high_resolution)
