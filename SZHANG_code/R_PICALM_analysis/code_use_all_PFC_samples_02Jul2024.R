# Siwei 02 Jul 2024
# Import scRNA-seq data of previous microglia to identify potential samples
# use all PFC samples, check PICALM-AD expression, disregard genotype

# init ####
{
  library(Seurat)
  library(Signac)

  library(edgeR)

  library(future)

  library(stringr)

  library(harmony)
  library(MAST)
  library(SingleCellExperiment)

  library(future)
  library(dplyr)

  library(readr)

  library(ggplot2)
  library(RColorBrewer)

  library(lsr)
}

plan("multisession", workers = 6)
options(mc.cores = 32)
set.seed(42)
options(future.globals.maxSize = 229496729600)

## read back Sun et all including all individuals####
processed_sun_et_al <-
  readRDS("Seurat_sun_PFC_cells_harmony.RDs")
unique(processed_sun_et_al$subject)

colnames(processed_sun_et_al@meta.data)
# [1] "orig.ident"             "nCount_RNA"             "nFeature_RNA"           "subject"
# [5] "brainRegion"            "batch"                  "barcode"                "percent.mt"
# [9] "age_death.x"            "msex.x"                 "pmi.x"                  "ADdiag3types"
# [13] "percent.rp"             "seurat_clusters"        "individualID"           "individualIdSource"
# [17] "species"                "sex"                    "projid"                 "Study"
# [21] "msex.y"                 "educ"                   "race"                   "spanish"
# [25] "apoe_genotype"          "age_at_visit_max"       "age_first_ad_dx"        "age_death.y"
# [29] "cts_mmse30_first_ad_dx" "cts_mmse30_lv"          "pmi.y"                  "braaksc"
# [33] "ceradsc"                "cogdx"                  "dcfdx_lv"               "nCount_SCT"
# [37] "nFeature_SCT"           "SCT_snn_res.0.5"

# split by subject id and re-merge by Harmony

# Split_seurat_sun_PFC <-
#   SplitObject(processed_sun_et_al,
#               split.by = "subject")

logNormalize_processed_sun_et_al <-
  processed_sun_et_al %>%
  NormalizeData(assay = "RNA",
                normalization.method = "LogNormalize",
                verbose = T) %>%
  FindVariableFeatures(selection.method = "vst",
                       binning.method = "equal_frequency",
                       verbose = T) %>%
  ScaleData(vars.to.regress = c("percent.mt",
                                "percent.rp",
                                "msex.x"),
            model.use = "poisson",
            use.umi = T,
            do.scale = T,
            do.center = T,
            verbose = T) %>%
  RunPCA(seed.use = 42,
         verbose = T) %>%
  FindNeighbors(reduction = "pca",
                verbose = T) %>%
  FindClusters(resolution = 0.5,
               algorithm = 1,
               random.seed = 42,
               verbose = T) %>%
  RunTSNE(reduction = "pca",
          seed.use = 42,
          tsne.method = "Rtsne",
          dims = 2:30) %>%
  RunUMAP(dims = 1:30,
          # umap.method = "umap-learn",
          seed.use = 42,
          verbose = T)
saveRDS(logNormalize_processed_sun_et_al,
        file = "log_normailsed_seurat_sun_PFC_all_samples.RDs")

# pseudo_exp_Sun <-
#   PseudobulkExpression(logNormalize_processed_sun_et_al,
#                        assays = "RNA",
#                        layer = "counts",
#                        group.by = "subject",
#                        verbose = T)

logNormalize_processed_sun_et_al <-
  log_normailsed_seurat_sun_PFC_all_samples

pseudo_exp_Sun <-
  AggregateExpression(logNormalize_processed_sun_et_al,
                      assays = "RNA",
                      group.by = "projid",
                      normalization.method = "CLR",
                      return.seurat = T,
                      verbose = T)

colnames(pseudo_exp_Sun)
rownames(pseudo_exp_Sun@assays$RNA@layers$data)
colnames(pseudo_exp_Sun@assays$RNA@layers$data)
hist(pseudo_exp_Sun@assays$RNA@layers$data)

df_2_plot <-
  data.frame(exp.value = pseudo_exp_Sun@assays$RNA@layers$data[rownames(pseudo_exp_Sun) == "PICALM", ],
             projid = colnames(pseudo_exp_Sun))

df_ADdiag3types <-
  data.frame(projid = str_c("g",
                            logNormalize_processed_sun_et_al$projid),
             ADdiag3types = logNormalize_processed_sun_et_al$ADdiag3types,
             brainRegion = logNormalize_processed_sun_et_al$brainRegion)
df_ADdiag3types <-
  df_ADdiag3types[!duplicated(df_ADdiag3types$projid), ]
df_ADdiag3types <-
  df_ADdiag3types[df_ADdiag3types$brainRegion == "PFC", ]

df_2_plot_final <-
  merge(x = df_2_plot,
        y = df_ADdiag3types,
        by.x = "projid",
        by.y = "projid")
df_2_plot_final$ADdiag3types <-
  factor(df_2_plot_final$ADdiag3types,
         levels = c("nonAD",
                    "earlyAD",
                    "lateAD"))


ggplot(df_2_plot_final,
       aes(x = ADdiag3types,
           y = exp.value,
           fill = ADdiag3types)) +
  geom_boxplot(colour = "black",
               outlier.size = 0,
               width = 0.5) +
  # geom_dotplot(binaxis = "y",
  #              stackdir = "center",
  #              dotsize = 0.1,
  #              binwidth = 0.25) +
  geom_jitter(size = 0.1,
              width = 0.1) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  # stat_summary(fun.data = mean_cl_normal) +
  stat_smooth(aes(group = 1),
              method = 'lm',
              linetype = 2,
              # formula = exp.value ~ ADdiag3types,
              se = F,
              colour = brewer.pal(n = 8,
                                  name = "Dark2")[4]) +
  # geom_violin(scale = "width") +
  theme_classic() +
  theme(axis.text = element_text(size = 12))

wilcox.test(x = df_2_plot_final$exp.value[df_2_plot_final$ADdiag3types == "nonAD"],
            y = df_2_plot_final$exp.value[df_2_plot_final$ADdiag3types == "earlyAD"],
            paired = F, exact = F)

## try edgeR and plot CPM

DEG_ADdiag3types <-
  DGEList(counts = pseudo_exp_Sun@assays$RNA@layers$counts,
          samples = pseudo_exp_Sun$projid,
          genes = rownames(pseudo_exp_Sun))
head(colnames(pseudo_exp_Sun))
head(pseudo_exp_Sun$projid)
head(df_2_plot_final)

DEG_ADdiag3types <-
  DEG_ADdiag3types[, order(DEG_ADdiag3types$samples$samples)]

df_2_plot_final <-
  df_2_plot_final[order(df_2_plot_final$projid), ]

DEG_ADdiag3types <-
  calcNormFactors(DEG_ADdiag3types)

cpm_ADdiag3types <-
  as.data.frame(edgeR::cpm(DEG_ADdiag3types))
colnames(cpm_ADdiag3types) <-
  df_2_plot_final$projid
head(colnames(cpm_ADdiag3types))
rownames(cpm_ADdiag3types)
head(df_2_plot_final)

df_2_plot_final$cpm.PICALM <-
  unlist(cpm_ADdiag3types[DEG_ADdiag3types$genes$genes == "PICALM", ])

ggplot(df_2_plot_final,
       aes(x = ADdiag3types,
           y = cpm.PICALM,
           fill = ADdiag3types)) +
  geom_boxplot(colour = "black",
               outlier.size = 0,
               width = 0.5) +
  # geom_dotplot(binaxis = "y",
  #              stackdir = "center",
  #              dotsize = 0.1,
  #              binwidth = 0.25) +
  geom_jitter(size = 0.1,
              width = 0.1) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  # stat_summary(fun.data = mean_cl_normal) +
  stat_smooth(aes(group = 1),
              method = 'lm',
              linetype = 2,
              # formula = exp.value ~ ADdiag3types,
              se = F,
              colour = brewer.pal(n = 8,
                                  name = "Dark2")[4]) +
  # geom_violin(scale = "width") +
  theme_classic() +
  theme(axis.text = element_text(size = 12)) +
  ggtitle("eta.sq = 0.0162; \np{EarlyVsNone}=0.0412; \np{LateVsNone}=0.0263")

wilcox.test(x = df_2_plot_final$cpm.PICALM[df_2_plot_final$ADdiag3types == "nonAD"],
            y = df_2_plot_final$cpm.PICALM[df_2_plot_final$ADdiag3types == "earlyAD"],
            paired = F, exact = F)
wilcox.test(x = df_2_plot_final$cpm.PICALM[df_2_plot_final$ADdiag3types == "nonAD"],
            y = df_2_plot_final$cpm.PICALM[df_2_plot_final$ADdiag3types == "lateAD"],
            paired = F, exact = F)

lm_output <-
  lm(formula = cpm.PICALM ~ ADdiag3types,
     data = df_2_plot_final)
coef(lm(formula = cpm.PICALM ~ ADdiag3types,
        data = df_2_plot_final))
summary(lm(formula = cpm.PICALM ~ ADdiag3types,
           data = df_2_plot_final)$r.squared)
etaSquared(lm_output)

anova_output <-
  aov(formula = cpm.PICALM ~ ADdiag3types,
     data = df_2_plot_final)
etaSquared(anova_output,
           anova = T)
