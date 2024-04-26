# Siwei 03 Jun 2022
# make PCA plot of AST, GA, MG, and NGN2-RGlut
# intervals used summit +/- 250 bp of 4 cell 
# type-specific peaks

# init
library(sva)
library(edgeR)

library(readr)

library(factoextra)
library(Rfast)

library(ggplot2)
library(ggrepel)
library(gplots)
library(RColorBrewer)

library(stringr)

# load the main data sheet
master_sheet_4_R <- 
  read_delim("data_4_PCA/output/master_sheet_4_R.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

batch_info <-
  factor(c(rep_len("2021", length.out = 12),
           rep_len("2021", length.out = 20),
           rep_len("2022", length.out = 14),
           rep_len("2021", length.out = 17),
           rep_len("2022", length.out = 11),
           rep_len("2021", length.out = 20),
           rep_len("2022", length.out = 21)))


cell_types <-
  factor(c(rep_len("AST", length.out = 12),
           rep_len("GA_2021", length.out = 20),
           rep_len("GA_2022", length.out = 14),
           rep_len("MG_2021", length.out = 17),
           rep_len("MG_2022", length.out = 11),
           rep_len("RGlut_2020", length.out = 20),
           rep_len("RGlut_2022", length.out = 21)))

# pheno <- 
#   data.frame(batch_info = factor(c(rep_len("AST", length.out = 12),
#                                    rep_len("GA_2021", length.out = 20),
#                                    rep_len("GA_2022", length.out = 14),
#                                    rep_len("MG_2021", length.out = 17),
#                                    rep_len("MG_2022", length.out = 11),
#                                    rep_len("RGlut_2020", length.out = 20),
#                                    rep_len("RGlut_2022", length.out = 21))))

# modcombat <- model.matrix(~1, data = pheno)
covar_mat <-
  as.matrix(data.frame(cells = factor(c(rep_len("AST", length.out = 12),
                                        rep_len("GA", length.out = 34),
                                        rep_len("MG", length.out = 28),
                                        rep_len("RGlut", length.out = 41)))))

# combat_adjusted_input <-
#   ComBat_seq(counts = as.matrix(master_sheet_4_R),
#              batch = batch_info,
#              # group = cell_types,
#              # covar_mod = covar_mat,
#              full_mod = T)

DGE_master_data <-
  DGEList(counts = as.matrix(master_sheet_4_R),
          samples = colnames(master_sheet_4_R), 
          group = cell_types,
          remove.zeros = T)
DGE_master_data <-
  calcNormFactors(DGE_master_data)
DGE_master_data <-
  estimateDisp(DGE_master_data)

DGE_master_data <-
  glmFit(DGE_master_data)

fitted_counts <- DGE_master_data$fitted.values
cpm_from_fitted <-
  apply(fitted_counts, 2, function(x)(x / sum(x)*1e6))
log_cpm_fitted <-
  log1p(cpm_from_fitted)

cpm_DGE <- cpm(DGE_master_data$counts + 1)

log_cpm_fitted <-
  log1p(cpm_DGE)

hist(log_cpm_fitted, breaks = 1000)
hist(rowSds(log_cpm_fitted), breaks = 1000)

log_cpm_fitted_plot <-
  log_cpm_fitted[rowSds(log_cpm_fitted) < 0.5, ]

## calculate PCA use prcomp
## ! plot the first two PCs use the x matrix,
## ! not the rotation matrix of the prcomp return
pca_results <- 
  prcomp(x = as.data.frame(t(cpm_DGE)),
         center = T,
         scale. = T,
         # retx = F,
         tol = 0)

fviz_pca_ind(pca_results,
             col.ind = "contrib",
             gradient.cols = c("darkblue", "grey50", "red"),
             repel = T)

pca_to_plot <- pca_results$x
pca_to_plot <- as.data.frame(pca_to_plot)
pca_to_plot_gg <- data.frame(x = pca_to_plot$PC1,
                             y = pca_to_plot$PC2,
                             label = colnames(cpm_DGE),
                             cell_types = cell_types,
                             stringsAsFactors = F)
pca_to_plot_gg$label <-
  str_replace_all(string = pca_to_plot_gg$label,
                  pattern = "Glut_rapid_neuron",
                  replacement = "NGN2_")

pca_to_plot_gg$label <-
  str_remove_all(string = pca_to_plot_gg$label,
                 pattern = "_new")

ggplot(pca_to_plot_gg, 
       aes(x = x, 
           y = y,
           label = label,
           color = cell_types)) +
  geom_point(size = 1) +
  xlab("PC1 = 36.3%") +
  ylab("PC2 = 17.0%") +
  scale_color_manual(values = brewer.pal(n = 8, name = "Dark2")) +
  theme_bw() +
  geom_text_repel(show.legend = F, 
                  force = 30,
                  min.segment.length = 0,
                  force_pull = 0)






cpm_DGE <- cpm(DGE_master_data)
