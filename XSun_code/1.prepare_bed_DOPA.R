library(rtracklayer)
library(GenomicRanges)

setwd("/project2/xinhe/xsun/neuron_simulation/2.torus")

df_snp <- readxl::read_xlsx("data/Bulk_DOPA_ASoC for TORUS.xlsx",sheet = 1)

gr_hg38 <- GRanges(seqnames = df_snp$CHROM, 
                   ranges = IRanges(names = df_snp$ID, start = as.numeric(df_snp$POS), end = as.numeric(df_snp$POS))
)
chain <- import.chain("~/xsun_xin/psych_analysis/1.torus/data/hg38ToHg19.over.chain")
result_hg19 <- as.data.frame(liftOver(gr_hg38,chain))
result_hg19 <- result_hg19[,-1]
colnames(result_hg19)[1:2] <- c("ID", "CHR")

index <- seq(2,nrow(result_hg19)*2, by =2)
chr <- unlist(strsplit(as.character(result_hg19$CHR), split = "chr"))[index]

bed <- cbind(chr,as.numeric(as.character(result_hg19$start)) - 1, as.numeric(as.character(result_hg19$end)))
write.table(bed,file = "./data/DOPA.bed",col.names = F, row.names = F,quote = F, sep = "\t")


