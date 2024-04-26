# Siwei 20 Oct 2023
# convert logOR and P to Z-score

library(readr)
library(parallel)
library(future)
library(stringr)

df_input <-
  read_delim("UKB_MDD_2019_2_convert.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

clust_df_assembled <- makeCluster(type = "FORK", 20)
clusterExport(cl = clust_df_assembled,
              varlist = "df_input",
              envir = environment())
df_input$zscore <-
  parApply(cl = clust_df_assembled,
           X = df_input,
           MARGIN = 1,
           FUN = function(x)(as.numeric(x[4]) / abs(qnorm(p = as.numeric(x[5]) / 2))))
stopCluster(clust_df_assembled)
rm(clust_df_assembled)

se <- (4.627 - 1.248) / (2 * 1.96)
z <- 2.938 / se
Pval <- exp(z * (-0.717) - (0.416) * z ^ 2)
