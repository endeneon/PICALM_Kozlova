# Siwei 12 Jan 2024
# test Hanwen's FACS data with flowCell and ggCyto

# init ####
library(ggcyto)

# load data ####
raw_flowFrame <-
  read.FCS("293T_256_STOP1.fcs",
           transformation = F,
           alter.names = T)
