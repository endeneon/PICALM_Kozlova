library(readr)
library(dplyr)
library(stringr)


variants <-
  read_tsv(file = "sample_minqiao_SCZGWAS.zscore",
           col_names = F,
           col_types = 'ccd')
colnames(variants) = c('variant', 'ld', 'zscore')
var_info <- str_split(variants$variant, ":")
variants <- mutate(variants, chr = paste0("chr", sapply(var_info, function(x){x[1]})),
                   pos = sapply(var_info, function(x){x[2]})) %>%
  mutate(start = as.numeric(pos), stop=as.numeric(pos) + 1) %>%
  select(chr, start, stop, variant)
options(scipen=1000) # So that positions are always fully written out)
write.table(variants,
            file="output.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
