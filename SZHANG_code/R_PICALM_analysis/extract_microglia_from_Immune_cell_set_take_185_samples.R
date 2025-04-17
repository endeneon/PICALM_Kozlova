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



# Immune_cells <-
#   readRDS(x)
Immune_cells <-
  readRDS(all_RDS_files[5])
Immune_cells_meta <-
  Immune_cells@meta.data
Immune_cells_meta$barcode <-
  rownames(Immune_cells_meta)
print(unique(str_split(Immune_cells_meta$barcode,
                       pattern = '-',
                       simplify = T)[, 3]))
print(length(unique(Immune_cells_meta$projid)))
# Immune_cells_meta$trimmed_barcode <-
#   str_c(str_split(string = Immune_cells_meta$barcode,
#                   pattern = '-',
#                   simplify = T)[, 1],
#         str_split(string = Immune_cells_meta$barcode,
#                   pattern = '-',
#                   simplify = T)[, 2],
#         sep = '-')

# Immune_cells_cluster_0 <-
#   Immune_cells_meta[str_detect(string = Immune_cells_meta$barcode,
#                                pattern = "\\-0$"), ]
# unique(Immune_cells_cluster_0$projid)



Immune_lookup_table <-
  Immune_cells_meta[!duplicated(Immune_cells_meta$projid), ]

Sun_meta <-
  readRDS("sun_et_al/personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")
print(length(unique(Sun_meta$subject)))

meta_indiv_index_from_immune <-
  readRDS("~/backuped_space/Siwei_misc_R_projects/R_Alena_PICALM/meta_indiv_index_from_immune.RDs")
Sun_185_project_ids <-
  readRDS("~/backuped_space/Siwei_misc_R_projects/R_Alena_PICALM/Sun_185_project_ids.RDs")

Immune_cells_subsetted <-
  Immune_cells[, Immune_cells$projid %in% Sun_185_project_ids$bare_projid]
Immune_cells_subsetted <-
  Immune_cells_subsetted[, Immune_cells_subsetted$cell_type_high_resolution == 'Mic P2RY12']
# combined_Immune_cells_meta_Sun_meta <-
#   merge(x = Immune_cells_meta,
#         y = Sun_meta,
#         by.x = "trimmed_barcode",
#         by.y = "barcode")
# print(length(unique(combined_Immune_cells_meta_Sun_meta$projid)))

# Sun_immune_lookup_table <-
#   combined_Immune_cells_meta_Sun_meta[!duplicated(combined_Immune_cells_meta_Sun_meta$projid), ]
# Sun_immune_lookup_table <-
#   Sun_immune_lookup_table[, c(2, 5, 15)]
#
# Immune_cells_meta_reassign <-
#   merge(x = Immune_cells_meta,
#         y = Sun_immune_lookup_table,
#         by.x = "projid",
#         by.y = "projid",
#         all.x = T)
#
# Immune_cells$ADdiag3types <-
#   Immune_cells_meta_reassign$ADdiag3types
Immune_cells_subsetted_meta <-
  Immune_cells_subsetted@meta.data
Immune_cells_subsetted_meta <-
  merge(x = Immune_cells_subsetted_meta,
        y = Sun_185_project_ids,
        by.x = "projid",
        by.y = "bare_projid",
        all.x = T)
Immune_cells_subsetted$ADdiag3types <-
  Immune_cells_subsetted_meta$ADdiag3types



print('=== 3 === ')

Immune_cells_subsetted <-
  UpdateSeuratObject(Immune_cells_subsetted)

Immune_cells_noNA <-
  Immune_cells_subsetted[, !is.na(Immune_cells_subsetted$ADdiag3types)]
# Immune_cells_noNA <-
#   Immune_cells_noNA[, Immune_cells_noNA$cell_type_high_resolution == 'Mic P2RY12']

print('=== assigned ADdiag3types === ')

saveRDS(Immune_cells_subsetted,
        file = "Immune_cells_subsetted_15Jan2025")



Immune_cells_noNA <-
  Immune_cells_noNA %>%
  NormalizeData(assay = "RNA",
                normalization.method = "LogNormalize",
                verbose = T) %>%
  FindVariableFeatures(selection.method = "vst",
                       binning.method = "equal_frequency",
                       verbose = T) %>%
  ScaleData(model.use = "poisson",
            use.umi = T,
            do.scale = T,
            do.center = T,
            verbose = T)
# Immune_cells_noNA <-
#   SCTransform(Immune_cells_noNA,
#               vst.flavor = "v2",
#               verbose = T)

Idents(Immune_cells_noNA) <- "projid"

Immune_pseudobulk <-
  AggregateExpression(Immune_cells_noNA,
                      assays = "RNA",
                      return.seurat = F,
                      group.by = "ident",
                      normalization.method = "LogNormalize")
Immune_pseudobulk <-
  Immune_pseudobulk$RNA

Immune_pseudobulk_DGE <-
  DGEList(counts = as.matrix(Immune_pseudobulk),
          sample = colnames(Immune_pseudobulk),
          genes = rownames(Immune_pseudobulk))
Immune_pseudobulk_logCPM <-
  as.data.frame(cpm(Immune_pseudobulk_DGE,
                    log = T))

Immune_pseudobulk_logCPM <-
  data.frame(expression = Immune_pseudobulk_logCPM[rownames(Immune_pseudobulk_logCPM) == "PICALM", ])
Immune_pseudobulk_logCPM$g_projid <-
  rownames(Immune_pseudobulk_logCPM)

Immune_lookup_table$g_projid <-
  str_c("g",
        Immune_lookup_table$projid,
        sep = "")
Sun_immune_lookup_table$g_projid <-
  str_c("g",
        Sun_immune_lookup_table$projid,
        sep = "")

Immune_pseudobulk_logCPM <-
  merge(x = Immune_pseudobulk_logCPM,
        y = Sun_immune_lookup_table,
        by = "g_projid")

print(Immune_pseudobulk_logCPM)

return(list(df = Immune_pseudobulk_logCPM
            # df2 = Immune_pseudobulk
            # nonADvslateAD = wilcox.test(x = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"],
            #             y = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "lateAD"],
            #             paired = F,
            #             alternative = "t"),
            #
            # nonADvsearlyAD = wilcox.test(x = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"],
            #             y = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "earlyAD"],
            #             paired = F,
            #             alternative = "t"),

            # mean_nonAD = mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"]),
            # mean_lateAD = mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "lateAD"]),
            # mean_earlyAD = mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "earlyAD"])))
))
#
