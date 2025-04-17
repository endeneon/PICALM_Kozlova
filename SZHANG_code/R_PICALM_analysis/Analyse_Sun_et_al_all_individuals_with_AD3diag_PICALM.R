# Extract Microglia only, of the 185 Sun et al samples.



{
  library(stringr)
  library(Seurat)

  library(parallel)
  library(future)


  library(glmGamPoi)

  library(edgeR)

  library(data.table)

  library(readr)

  plan("multisession", workers = 3)
  # options(mc.cores = 32)
  set.seed(42)
  options(future.globals.maxSize = 429496729600)

  setwd("~/backuped_space/Siwei_misc_R_projects/R_Alena_PICALM")
}


load("workspace_21Jun2024.RData")


Sun_all_cells_PFC_noNa <-
  Seurat_sun_et_al[, (Seurat_sun_et_al$brainRegion == "PFC") & !(is.na(Seurat_sun_et_al$ADdiag3types))]
UMAPPlot(Sun_all_cells_PFC_noNa)

length(unique(Sun_all_cells_PFC_noNa$projid))

Sun_all_pseudobulk <-
  AggregateExpression(Sun_all_cells_PFC_noNa,
                      assays = "RNA",
                      return.seurat = T,
                      group.by = "projid",
                      normalization.method = "CLR")
Sun_all_pseudobulk$projid ## use the order of this one
sort(unique(Sun_all_pseudobulk$projid))



lookup_projID_ADdiag3 <-
  Sun_all_cells_PFC_noNa@meta.data
lookup_projID_ADdiag3 <-
  lookup_projID_ADdiag3[!duplicated(lookup_projID_ADdiag3$projid), ]
lookup_projID_ADdiag3$g_projid <-
  str_c("g",
        lookup_projID_ADdiag3$projid,
        sep = "")
lookup_projID_ADdiag3$g_projid
# sort(unique(lookup_projID_ADdiag3$g_projid))

# lookup_projID_ADdiag3 <-
#   lookup_projID_ADdiag3[!duplicated(lookup_projID_ADdiag3$g_projid), ]
# lookup_projID_ADdiag3$g_projid

Sun_all_pseudobulk_meta <-
  Sun_all_pseudobulk@meta.data
# Sun_all_pseudobulk_meta$projid
# sort(unique(Sun_all_pseudobulk_meta$projid))

all(Sun_all_pseudobulk$projid == Sun_all_pseudobulk_meta$projid)

Sun_all_pseudobulk_meta_ADdiag3 <-
  Sun_all_pseudobulk_meta_ADdiag3[match(x = Sun_all_pseudobulk$projid,
                                        table = Sun_all_pseudobulk_meta_ADdiag3$projid), ]
all(Sun_all_pseudobulk_meta_ADdiag3$projid == Sun_all_pseudobulk$projid)
# sort(unique(Sun_all_pseudobulk_meta_ADdiag3$projid))

# Sun_all_pseudobulk_meta_ADdiag3 <-
#   Sun_all_pseudobulk_meta_ADdiag3[Sun_all_pseudobulk_meta$projid %in% Sun_all_pseudobulk_meta_ADdiag3$projid, ]
# colnames(Sun_all_pseudobulk)
#
# Sun_all_pseudobulk_ordered <-
#   Sun_all_pseudobulk[, order(Sun_all_pseudobulk$projid)]
# Sun_all_pseudobulk_meta_ADdiag3_ordered <-
#   Sun_all_pseudobulk_meta_ADdiag3[order(Sun_all_pseudobulk_meta_ADdiag3$projid), ]
# head(colnames(Sun_all_pseudobulk_ordered))
# head(Sun_all_pseudobulk_meta_ADdiag3_ordered$projid)
#
# all(colnames(Sun_all_pseudobulk_ordered),
#     Sun_all_pseudobulk_meta_ADdiag3_ordered$projid)

DGE_Sun_pseudobulk <-
  DGEList(counts = as.matrix(Sun_all_pseudobulk@assays$RNA@layers$counts),
          samples = colnames(Sun_all_pseudobulk),
          genes = rownames(Sun_all_pseudobulk),
          group = Sun_all_pseudobulk_meta_ADdiag3$ADdiag3types)

df_2_plot <-
  as.data.frame(cpm(DGE_Sun_pseudobulk,
                    log = F,
                    normalized.lib.sizes = T))
rownames(df_2_plot) <-
  rownames(Sun_all_pseudobulk)
colnames(df_2_plot) <-
  colnames(Sun_all_pseudobulk)

df_2_plot <-
  data.frame(value = unlist(df_2_plot[rownames(df_2_plot) == "PICALM", ]),
             Phenotype = Sun_all_pseudobulk_meta_ADdiag3$ADdiag3types)
df_2_plot$Phenotype <-
  factor(df_2_plot$Phenotype,
         levels = c("nonAD",
                    "earlyAD",
                    "lateAD"))
df_2_plot$value <-
  log2(df_2_plot$value)

ggerrorplot(df_2_plot,
            x = "Phenotype",
            y = "value",
            color = "Phenotype",
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
  geom_jitter(aes(colour = Phenotype),
              width = 0.1,
              size = 1) +
  # scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "CPM") +
  scale_y_continuous(expand = c(0, 0)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("Sun et al. all 425 individuals")



Immune_pseudobulk_logCPM_syn52368912 <-
  read.table("PICALM_Microglia_P2RY12_from_Immune_syn52368912.txt",
             sep = "\t",
             header = T, quote = "")

ggerrorplot(Immune_pseudobulk_logCPM_syn52368912,
            x = "ADdiag3types",
            y = "expression",
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
       y = "log2CPM") +
  scale_y_continuous(expand = c(0, 0)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1))
