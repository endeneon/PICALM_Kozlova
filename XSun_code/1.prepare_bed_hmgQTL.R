setwd("/project2/xinhe/xsun/neuron_simulation/2.torus")

df_snp <- readxl::read_xlsx("data/SNPs_for_Torus.xlsx",sheet = 3)
colnames(df_snp)[(ncol(df_snp)-2):ncol(df_snp)] <- c("chr","start","end")
df_snp <- df_snp[complete.cases(df_snp$chr),]


index <- seq(2,nrow(df_snp)*2, by =2)
chr <- unlist(strsplit(df_snp$chr, split = "chr"))[index]

bed <- cbind(chr,as.numeric(as.character(df_snp$start)), as.numeric(as.character(df_snp$end)))
write.table(bed,file = "./data/hMGcaQTL.bed",col.names = F, row.names = F,quote = F, sep = "\t")


df_snp <- readxl::read_xlsx("data/SNPs_for_Torus.xlsx",sheet = 4)
colnames(df_snp)[(ncol(df_snp)-2):ncol(df_snp)] <- c("chr","start","end")
df_snp <- df_snp[complete.cases(df_snp$chr),]

index <- seq(2,nrow(df_snp)*2, by =2)
chr <- unlist(strsplit(df_snp$chr, split = "chr"))[index]

bed <- cbind(chr,as.numeric(as.character(df_snp$start)), as.numeric(as.character(df_snp$end)))
write.table(bed,file = "./data/hMGeQTL.bed",col.names = F, row.names = F,quote = F, sep = "\t")
