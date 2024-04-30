library(mapgen)

setwd("/project2/xinhe/xsun/neuron_simulation/2.torus")

bigSNP <- bigsnpr::snp_attach(rdsfile = '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds')
data('Euro_LD_Chunks', package='mapgen')

folder_gwas_mapgen <- "/project2/xinhe/xsun/psych_analysis/1.torus/data/gwas/"
file_gwas <- list.files(folder_gwas_mapgen)

enrichment <- c()
for (i in 1:length(file_gwas)) {
  
  gwas_file <- paste0(folder_gwas_mapgen,file_gwas[i])
  if (grepl(gwas_file, pattern = "gz")) {
    trait <- unlist(strsplit(file_gwas[i],split= ".txt.gz"))
  }else {
    trait <- unlist(strsplit(file_gwas[i],split= ".txt"))
  }
  
  gwas <- process_gwas_sumstats(gwas_file,
                                chr='CHR', pos='BP',
                                beta='BETA', se='SE',
                                a0='REF', a1='ALT',
                                snp='SNP', pval='P',
                                bigSNP=bigSNP,
                                LD_Blocks=LD_Blocks)
  
  bed_annotations_dir <- "data/"
  annotation_bed_files <- "data//iAstro.bed"
  
  torus.files <- prepare_torus_input_files(gwas, annotation_bed_files, torus_input_dir = "./torus_input3")
  
  torus.result <- run_torus(torus.files$torus_annot_file, 
                            torus.files$torus_zscore_file,
                            option = "est-prior",
                            torus_path = "/home/sxt1229/software/torus/src/torus") # set the path to 'torus' executable.
  
  torus.enrich <- torus.result$enrich
  tmp1 <- torus.enrich[torus.enrich$term == "iAstro.bed.1", ]
  tmp1$snp <- "iAstro"
  tmp1[,1] <- trait
  
  enrichment <- rbind(enrichment,tmp1)
  save(enrichment,file = "enrichment_iAstro.rdata")
}

gwas_file <- paste0("/project2/xinhe/xsun/neuron_simulation/2.torus/gwas/Alzheimer_hg19_hm.txt.gz")
trait <- "Alzheimer_hg19_hm"

gwas <- process_gwas_sumstats(gwas_file,
                              chr='chr', pos='pos',
                              beta='beta', se='se',
                              a0='ref', a1='alt',
                              snp='rsid', pval='p',
                              bigSNP=bigSNP,
                              LD_Blocks=LD_Blocks)

bed_annotations_dir <- "data/"
annotation_bed_files <- "data//iAstro.bed"

torus.files <- prepare_torus_input_files(gwas, annotation_bed_files, torus_input_dir = "./torus_input3")

torus.result <- run_torus(torus.files$torus_annot_file, 
                          torus.files$torus_zscore_file,
                          option = "est-prior",
                          torus_path = "/home/sxt1229/software/torus/src/torus") # set the path to 'torus' executable.

torus.enrich <- torus.result$enrich
tmp1 <- torus.enrich[torus.enrich$term == "iAstro.bed.1", ]
tmp1$snp <- "iAstro"
tmp1[,1] <- trait

enrichment <- rbind(enrichment,tmp1)
save(enrichment,file = "enrichment_iAstro.rdata")