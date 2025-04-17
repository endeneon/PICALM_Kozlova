# Extract Microglia only, of all 425 samples.



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


#
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
Immune_cells$projid

Immune_cells_meta <-
  Immune_cells@meta.data
# Immune_cells_meta$barcode <-
#   rownames(Immune_cells_meta)
# print(unique(str_split(Immune_cells_meta$barcode,
#                        pattern = '-',
#                        simplify = T)[, 3]))
print(length(unique(Immune_cells_meta$projid)))
Immune_lookup_table <-
  Immune_cells_meta[!duplicated(Immune_cells_meta$projid), ]

Sun_meta <-
  readRDS("sun_et_al/personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")
print(length(unique(Sun_meta$subject)))
# Sun_meta_425 <-
#   Sun_meta[!duplicated(Sun_meta$)]
sun_et_al_cells <-
  readRDS("human_microglia_state_dynamics/ROSMAP.ImmuneCells.6regions.snRNAseq.counts.rds")
colnames(sun_et_al_cells)


#####
processed_sun_et_al <-
  readRDS("Seurat_sun_PFC_cells_harmony.RDs")
UMAPPlot(processed_sun_et_al)
processed_sun_et_al$ADdiag3types
length(unique(processed_sun_et_al$subject))

# subset cluster 0:12 ####
Idents(processed_sun_et_al) <- "seurat_clusters"
unique(Idents(processed_sun_et_al))
# unique(processed_sun_et_al)


processed_sun_et_al_pseudobulk <-
  AggregateExpression(processed_sun_et_al,
                      assays = "RNA",
                      return.seurat = T,
                      group.by = "subject",
                      normalization.method = "CLR")
length(unique(processed_sun_et_al_pseudobulk$subject))

meta_processed_sun_et_al <-
  processed_sun_et_al@meta.data
meta_subject_ADdiag3 <-
  meta_processed_sun_et_al[!duplicated(meta_processed_sun_et_al$subject), ]

Sun_pseudobulk_matrix <-
  processed_sun_et_al_pseudobulk@assays$RNA@layers$counts
nrow(Sun_pseudobulk_matrix)
rownames(processed_sun_et_al_pseudobulk)
rownames(Sun_pseudobulk_matrix) <-
  rownames(processed_sun_et_al_pseudobulk)
colnames(Sun_pseudobulk_matrix) <-
  colnames(processed_sun_et_al_pseudobulk)

Sun_pseudobulk_matrix_reordered <-
  Sun_pseudobulk_matrix[, order(colnames(Sun_pseudobulk_matrix))]
meta_subject_ADdiag3_reordered <-
  meta_subject_ADdiag3[order(meta_subject_ADdiag3$subject), ]
all(colnames(Sun_pseudobulk_matrix_reordered) == meta_subject_ADdiag3_reordered$subject)


Sun_pseudobulk_DGE <-
  DGEList(counts = Sun_pseudobulk_matrix_reordered,
          sample = meta_subject_ADdiag3_reordered$subject,
          genes = rownames(Sun_pseudobulk_matrix_reordered),
          group = meta_subject_ADdiag3_reordered$ADdiag3types)

cpm_Sun_pseudobulk_DGE <-
  as.data.frame(edgeR::cpm(Sun_pseudobulk_DGE, log = T))
cpm_Sun_pseudobulk_PICALM <-
  cpm_Sun_pseudobulk_DGE[rownames(cpm_Sun_pseudobulk_DGE) == "PICALM", ]

df_2_plot <-
  data.frame(value = unlist(cpm_Sun_pseudobulk_PICALM),
             ADdiag3types = meta_subject_ADdiag3_reordered$ADdiag3types)


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
