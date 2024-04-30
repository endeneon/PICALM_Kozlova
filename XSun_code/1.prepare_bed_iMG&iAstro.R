library(rtracklayer)
library(GenomicRanges)

setwd("/project2/xinhe/xsun/neuron_simulation/2.torus")

df_snp <- readxl::read_xlsx("data/SNPs_for_Torus.xlsx",sheet = 1)

gr_hg38 <- GRanges(seqnames = df_snp$CHROM, 
                   ranges = IRanges(names = df_snp$ID, start = as.numeric(df_snp$POS), end = as.numeric(df_snp$POS))
                   )
chain <- import.chain("~/xsun_xin/psych_analysis/1.torus/data/hg38ToHg19.over.chain")
result_hg19 <- as.data.frame(liftOver(gr_hg38,chain))
result_hg19 <- result_hg19[,-1]
colnames(result_hg19)[1:2] <- c("ID", "CHR")

#save(result_hg19, file = "./data/SNP_hg19.rdata")

index <- seq(2,nrow(result_hg19)*2, by =2)
chr <- unlist(strsplit(as.character(result_hg19$CHR), split = "chr"))[index]

bed <- cbind(chr,as.numeric(as.character(result_hg19$start)) - 1, as.numeric(as.character(result_hg19$end)))
write.table(bed,file = "./data/iAstro.bed",col.names = F, row.names = F,quote = F, sep = "\t")

df_snp <- readxl::read_xlsx("data/iMG_ASoC_QTL_SNPs.xlsx",sheet = 2)

gr_hg38 <- GRanges(seqnames = df_snp$CHROM, 
                   ranges = IRanges(names = df_snp$ID, start = as.numeric(df_snp$POS), end = as.numeric(df_snp$POS))
)
#chain <- import.chain("~/xsun_xin/psych_analysis/1.torus/data/hg38ToHg19.over.chain")
result_hg19 <- as.data.frame(liftOver(gr_hg38,chain))
result_hg19 <- result_hg19[,-1]
colnames(result_hg19)[1:2] <- c("ID", "CHR")

#save(result_hg19, file = "./data/SNP_hg19.rdata")

index <- seq(2,nrow(result_hg19)*2, by =2)
chr <- unlist(strsplit(as.character(result_hg19$CHR), split = "chr"))[index]

bed <- cbind(chr,as.numeric(as.character(result_hg19$start)) - 1, as.numeric(as.character(result_hg19$end)))
write.table(bed,file = "./data/iMG.bed",col.names = F, row.names = F,quote = F, sep = "\t")