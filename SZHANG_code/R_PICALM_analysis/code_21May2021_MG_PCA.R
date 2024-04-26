# Siwei 21 May 2021
# Make PCA plot for 17 iMG samples


# init
library(readr)

library(factoextra)
library(Rfast)

library(ggplot2)
library(ggrepel)
library(gplots)
library(RColorBrewer)

library(stringr)

# load data
## load PCA raw data from featureCount output
### subsampled to 0.1 of the total reads
df_featureCount_raw <- 
  read_delim("summits_500bp_subsampled_0.1.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE, 
             skip = 1)

### NOT subsampled, scaled to the smaller sample 
### (send all bam files separately to macs2)
df_featureCount_raw <- 
  read_delim("summits_500bp_un_subsampled_summits_500bp_autosomes.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE, 
             skip = 1)

## section out reads for a separate df
df_read_counts <- df_featureCount_raw[, 7:ncol(df_featureCount_raw)]
rownames(df_read_counts) <- df_featureCount_raw$Geneid
sample_names <- str_remove_all(colnames(df_read_counts),
                               pattern = "_WASPed\\.bam")
colnames(df_read_counts) <- sample_names

## calculate PCA use prcomp
## ! plot the first two PCs use the x matrix,
## ! not the rotation matrix of the prcomp return
pca_results <- prcomp(x = as.data.frame(t(df_read_counts)),
                      center = T,
                      scale. = T,
                      # retx = F,
                      tol = 0)

plot(pca_results)
summary(pca_results)

fviz_pca_ind(pca_results,
             col.ind = "contrib",
             gradient.cols = c("darkblue", "grey50", "red"),
             repel = T)

# fviz_pca_contrib(pca_results,
#              col.ind = "cos2",
#              gradient.cols = c("red", "grey50", "darkblue"),
#              repel = T)

pca_to_plot <- pca_results$x
pca_to_plot <- as.data.frame(pca_to_plot)
pca_to_plot_gg <- data.frame(x = unlist(pca_to_plot[1 , ]),
                             y = unlist(pca_to_plot[2 , ]),
                             stringsAsFactors = F)

## make the plot
ggplot(data = pca_to_plot_gg, aes(x = x,
                                  y = y,
                                  label = colnames(df_read_counts))) +
  geom_point() +
  theme_classic() +
  geom_text_repel()
