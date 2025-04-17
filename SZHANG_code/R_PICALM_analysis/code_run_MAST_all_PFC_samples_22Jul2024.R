# Siwei 02 Jul 2024
# Import scRNA-seq data of previous microglia to identify potential samples
# use all PFC samples, check PICALM-AD expression, disregard genotype
# run MAST for exact exp. values
# adapt from Lexi's scDE code

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
  library(data.table)

  library(readr)

  library(ggplot2)
  library(RColorBrewer)

  library(lsr)
}

plan("multisession", workers = 6)
options(mc.cores = 32)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# functions to use ####
## a simple freq() function for Seurat objects ####
filterGenes <-
  function(sobj, thres) {
  # calculate the percentage of cells expressing each gene in a Seurat object
  # input: sobj, a Seurat object
  ct_mat <- GetAssayData(object = sobj, assay = "SCT", slot = "data")
  print(nrow(ct_mat))
  passfilter <- (rowSums(ct_mat) / ncol(ct_mat)) > thres
  return(passfilter)
}

# load log normalised Seurat subject ####
logNormalize_processed_sun_et_al <-
  readRDS(file = "log_normailsed_seurat_sun_PFC_all_samples.RDs")

types <- "PFC"
lines <- sort(unique(logNormalize_processed_sun_et_al$individualID))

deg_counts <-
  array(dim = c(1, 10),
        dimnames = list("PFC",
                        c("total genes",
                          "passed filter_nonAD", "upregulated_nonAD", "downregulated_nonAD",
                          "passed filter_earlyAD", "upregulated_earlyAD", "downregulated_earlyAD",
                          "passed filter_lateAD", "upregulated_lateAD", "downregulated_lateAD")))

if (!dir.exists("./MAST_scDE")) {
  dir.create("./MAST_scDE")
}
if (!dir.exists("./MAST_scDE/volcano_plots")) {
  dir.create("./MAST_scDE/volcano_plots")
}


unique(logNormalize_processed_sun_et_al$ADdiag3types)

for (i in 1:length(types)) {
  type <- types[i]
  subobj <- logNormalize_processed_sun_et_al
  deg_counts[i, 1] <-
    nrow(subobj)


  # nonAD <-
  #   subset(subobj,
  #          ADdiag3types == "nonAD")
  earlyAD_vs_nonAD <-
    subset(subobj,
           ADdiag3types %in% c("earlyAD",
                               "nonAD"))
  lateAD_vs_nonAD <-
    subset(subobj,
           ADdiag3types %in% c("lateAD",
                               "nonAD"))

  timeobj_list <-
    list(earlyAD_vs_nonAD,
         lateAD_vs_nonAD)
  names(timeobj_list) <-
    c("earlyAD_vs_nonAD",
      "lateAD_vs_nonAD")

  # free up memory
  rm(earlyAD_vs_nonAD,
     lateAD_vs_nonAD)
  gc()

  print(type)

  for (j in 1:length(timeobj_list)) {
    if (i == 3 & j == 1) {
      print("npglut_0hr met")
    } else {
      # time <- times[j]
      print(names(timeobj_list)[j])
      time <- names(timeobj_list)[j]
      timeobj <- timeobj_list[[j]]
      feat_mat <- data.frame(geneid = timeobj@assays$SCT@data@Dimnames[[1]],
                             row.names = timeobj@assays$SCT@data@Dimnames[[1]])
      cell_mat <- timeobj@meta.data # has to be a data.frame, contains metadata for each cell
      stopifnot(ncol(timeobj@assays$SCT@data) == nrow(cell_mat),
                nrow(timeobj@assays$SCT@data) == nrow(feat_mat))

      scaRaw <-
        FromMatrix(exprsArray = as.array(timeobj@assays$SCT@data),
                   cData = cell_mat,
                   fData = feat_mat)
      save(scaRaw,
           file = paste0("./MAST_scDE/scaRAW_obj_", type, "_", time, ".RData"))
      expressed_genes <- freq(scaRaw) > 0.05 # proportion of nonzero values across all cells
      sca <- scaRaw[expressed_genes, ]
      rm(scaRaw)
      gc()

      deg_counts[i, 3 * j - 1] <- sum(expressed_genes)
      colData(sca)$ADdiag3types <- factor(colData(sca)$ADdiag3types)
      colData(sca)$ADdiag3types <- relevel(colData(sca)$ADdiag3types,
                                           ref = "nonAD")
      colData(sca)$sex <- factor(colData(sca)$sex)
      colData(sca)$batch <- factor(colData(sca)$batch)
      # colData(sca)$batch <- factor(colData(sca)$batch)

      zlmCond <-
        zlm(~ ADdiag3types +
              age_death.x +
              batch +
              sex,
            sca = sca,
            parallel = T, )
      save(zlmCond, file = paste0("./MAST_scDE/MAST_zlmCond_", type, "_", time, ".RData"))

      summaryCond <- summary(zlmCond, doLRT = 'affcontrol')
      rm(zlmCond)
      gc()

      save(summaryCond,
           file = paste0("./MAST_scDE/summaryCond_", type, "_", time, ".RData"))

      summaryDt <- summaryCond$datatable
      rm(summaryCond)
      gc()

      fcHurdle <-
        merge(summaryDt[contrast == 'affcontrol' &
                          component == 'H',
                        .(primerid, `Pr(>Chisq)`)], #hurdle P values
              summaryDt[contrast == 'affcontrol' &
                          component == 'logFC',
                        .(primerid, coef, ci.hi, ci.lo)],
              by = 'primerid')

      fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
      # fcHurdle <- merge(fcHurdle[fdr < 0.05],
      #                   as.data.table(mcols(sca)), by='primerid')
      setorder(fcHurdle, fdr)
      save(fcHurdle, file = paste0("./MAST_scDE/results_", type, "_", time, ".RData"))
      deg_counts[i, 3 * j] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef > 0)
      deg_counts[i, 3 * j + 1] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef < 0)

      # plot
      fcHurdle$significance <- "nonsignificant"
      fcHurdle$significance[fcHurdle$fdr < 0.05 & fcHurdle$coef > 0] <- "up"
      fcHurdle$significance[fcHurdle$fdr < 0.05 & fcHurdle$coef < 0] <- "down"
      unique(fcHurdle$significance)
      fcHurdle$significance <- factor(fcHurdle$significance, levels = c("up", "nonsignificant", "down"))
      fcHurdle$neg_log_pval <- (0 - log2(fcHurdle$`Pr(>Chisq)`))
      png(paste0("./MAST_scDE/volcano_plots/", type, "_", time,
                 "_casevsctrl_volcano_plot.png"), width = 700, height = 700)
      p <- ggplot(data = as.data.frame(fcHurdle),
                  aes(x = coef,
                      y = neg_log_pval,
                      color = significance)) +
        geom_point(size = 0.2) +
        scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
        theme_minimal() +
        xlim(c(-2, 2)) +
        ggtitle(paste0("Case vs Control DE - ", type, " ", time))
      print(p)
      dev.off()

      pdf(paste0("./MAST_scDE/volcano_plots/", type, "_", time,
                 "_casevsctrl_volcano_plot.pdf"))
      p <- ggplot(data = as.data.frame(fcHurdle),
                  aes(x = coef,
                      y = neg_log_pval,
                      color = significance)) +
        geom_point(size = 0.2) +
        scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
        theme_minimal() +
        xlim(c(-2, 2)) +
        ggtitle(paste0("Case vs Control DE - ", type, " ", time))
      print(p)
      dev.off()
    }
  }
}
