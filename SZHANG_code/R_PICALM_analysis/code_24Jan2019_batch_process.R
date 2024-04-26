# 24 Jan 2019 Siwei
# Write a loop to walk through all ready-made BAM files
# SplitCigarNReads skipped, should not cause major problem

# Libraries
library(zoo)
library(gplots)
library(stringr)

# init environment
######

k <- 1
source_file_list <-
  list.files(path = "eSNPKaryotyping/R",
             full.names = T)
for (k in 1:length(source_file_list)) {
  source(source_file_list[k])
}

bam_file_list <- list.files(path = "~/4TB_2/RNASeq_Jun2018/BAMs/sorted/MiNND_batch1", full.names = T, recursive = F, pattern = "*bam$")
Directory <- "/home/cpg/4TB_2/RNASeq_Jun2018/BAMs/sorted/working/"
# Directory <- ""
Picard_Path <- "/home/cpg/Database/Tools/"
Genome_Fa <- "/home/cpg/Database/Databases/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
GATK_Path <- "/home/cpg/Database/Tools/GATK37/"
Organism <- "Human"
plot_output_directory <- "/home/cpg/4TB_2/RNASeq_Jun2018/eSNP_Karotype/output/"
# library(zoo)
# library(gplots)
dev.off()

k <- 1
for (k in 1:length(bam_file_list)) {
  print(paste(k, bam_file_list[k], sep = " "))
  print(date())
  # set the file name and title for the PDF
  pdf_series_name <- gsub(pattern = "/home/cpg/4TB_2/RNASeq_Jun2018/BAMs/sorted/MiNND_batch1/",
                          replacement = "",
                          x = bam_file_list[k])
  pdf_series_name <- gsub(pattern = "\\.bam",
                          replacement = "",
                          x = pdf_series_name)

  CreateVCF(Directory = Directory, Genome_Fa = Genome_Fa,
            Picard_Path = Picard_Path, GATK_Path = GATK_Path, bam_file_name = bam_file_list[k])

  EditVCF(Directory = Directory, Organism = Organism)
  # read the table in
  VCF_table <- read.delim(file = paste(Directory, "variantTable.csv", sep = ""), header = T, sep = "\t",
                          quote = "", dec = ".")
  # or
  # VCF_table <-
  #   EditVCF(Directory = Directory, Organism = Organism)
  VCF_table$chr <- as.numeric(VCF_table$chr)
  VCF_table <- VCF_table[order(VCF_table$chr, VCF_table$position), ]
  VCF_table <- VCF_table[VCF_table$chr > 0, ]
  #### backup
  # VCF_table_backup <- VCF_table
  # VCF_table <- VCF_table_backup
  # VCF_table$chr <- paste("chr", VCF_table$chr, sep = "")


  ###### Return MajorMinorCalc results
  MajorMinorCal_results <- MajorMinorCalc(Table = VCF_table, minDP = 20, maxDP = 1000, minAF = 0.2)
  ##### Plot Allelic ratio along the genome for duplication detection
  pdf(file = paste(plot_output_directory, pdf_series_name, "_genome.pdf", sep = ""),
      paper = "USr")
  PlotGenome(orderedTable = MajorMinorCal_results, Window = 151, Ylim = 3, PValue = TRUE, Organism = Organism)
  dev.off()
  # save current working directory
  # current_working_directory <- getwd()
  #pdf(plot)
  # setwd(plot_output_directory) # send plot to plot output directory

  # setwd(current_working_directory)

  ##### intersect the observed SNP with common SNP from dbSNPs
  # Siwei: the Genome_Fa_dict is simply a waste of memory...
  # moreover, it did not consider the situation that if the .dict has additional header but not started with @SQ
  tbl_DeletionTable_output <- DeletionTable(Directory = Directory,Table = MajorMinorCal_results,
                                   dbSNP_Data_Directory = "/home/cpg/Database/Databases/hg38/common_snp150_by_chr/",
                                   dbSNP_File_Name = "Edited_Common_chr",
                                   Genome_Fa_dict = "/home/cpg/Database/Databases/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.dict",
                                   Organism = "Human")


  # 9. Plot each SNP, without any summarization
  # pdf(file = paste(plot_output_directory, pdf_series_name, "_zygosity_sinle.pdf", sep = ""),
  #     paper = "USr")
  # Plot_Zygosity_Sinle(Table = tbl_DeletionTable_output, Organism = "Human")
  # # Argument: 1. Table - The LOH table containing the output of the DeletionTable function
  # #           2. Organism - "Human" or "Mouse"
  # dev.off()
  # 
  # 


  # 10. Plot blocks of heterozygous and homozygous SNPs
  pdf(file = paste(plot_output_directory, pdf_series_name, "_zygosity_blocks.pdf", sep = ""),
      paper = "USr")
  Plot_Zygosity_Blocks(Table = tbl_DeletionTable_output,
                       Window = 1500000,
                       Max = 6,
                       Max2 = 20,
                       Organism = "Human")
  # Argument: 1. Table - The deletion table containing the output of the DeletionTable function
  #           2. window - the block size in bp, usually 1500000
  #           3. Max - How many Heterozygouse SNP need to be in a block to get the full color, usually 6
  #           4. Max2 - How many Homozygouse SNP need to be in a block to get the full color, usually 60
  #           5. Organism - "Human" or "Mouse"
  dev.off()
}

###################
# bam_file_list[1]
# paste("I=",Directory,bam_file_list[1])
# 
# Edit_dbSNP_Files(Directory = "~/Database/Databases/hg38/common_snp150_by_chr/",
#                  File_Name = "dbSNP_150_common_chr",
#                  Organism = "Human")

