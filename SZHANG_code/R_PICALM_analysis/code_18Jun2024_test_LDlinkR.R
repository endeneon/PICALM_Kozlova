# Siwei 18 Jun 2024
# test LDlinkR

# install.packages("LDlinkR")

# init ####
library(LDlinkR)
library(vcfR)

library(stringr)

rs10792832_proxy_hg38 <-
  LDproxy(snp = "rs10792832",
          pop = "CEU",
          token = "ed9e7f4a5d87",
          genome_build = "grch38", r2d = "r2")

# The master SNP genotype file is in hg38 ####!
master_SNP <-
  read.vcfR(file = "freshmicroglia_geno/FreshMicroglia_DeIDNames_noCK_genotyped_only_109ind.vcf",
            limit = 1e+8)
head(master_SNP@fix)

SNP_list <-
  master_SNP@fix
SNP_list <-
  as.data.frame(SNP_list)

rs10792832_proxy_SNP_id_hg38 <-
  rs10792832_proxy_hg38$RS_Number
rs10792832_proxy_SNP_id_hg38 <-
  unique(rs10792832_proxy_SNP_id_hg38)
rs10792832_proxy_SNP_id_hg38 <-
  rs10792832_proxy_SNP_id_hg38[order(rs10792832_proxy_SNP_id_hg38)]
rs10792832_proxy_SNP_id_hg38 <-
  rs10792832_proxy_SNP_id_hg38[-c(1)]

rs10792832_found <-
  SNP_list[SNP_list$ID %in% rs10792832_proxy_SNP_id_hg38, ]

rs10792832_gt_found <-
  master_SNP@gt[SNP_list$ID %in% rs10792832_proxy_SNP_id_hg38, ]
rs10792832_gt_found <-
  as.data.frame(rs10792832_gt_found)

extract_gt_column <-
  function(x) {
    gt <-
      str_split(x,
                pattern = "\\:",
                simplify = T)[, 1]
    return(gt)
  }

rs10792832_gt_only <-
  sapply(rs10792832_gt_found,
         FUN = extract_gt_column,
         simplify = T)
rs10792832_gt_only <-
  as.data.frame(rs10792832_gt_only)
colnames(rs10792832_gt_only) <-
  str_split(string = colnames(rs10792832_gt_only),
            pattern = "__",
            simplify = T)[, 1]

rs10792832_combined_hg38 <-
  cbind(rs10792832_found,
        rs10792832_gt_only)

MG17_SNP <-
  read.vcfR(file = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/Microglia/GATK_call/unprocessed_vcf/het_vcf/output/MG_merged_SNP_01Jun2021.vcf",
            limit = 1e+9)
MG17_SNP_list <-
  as.data.frame(MG17_SNP@fix)

MG17_SNP_rs10792832_gt_found <-
  MG17_SNP@gt[MG17_SNP_list$ID %in% c(rs10792832_combined_hg38$ID,
                                      "rs10792832"), ]
MG17_SNP_rs10792832_gt_found <-
  as.data.frame(MG17_SNP_rs10792832_gt_found)

MG17_SNP_rs10792832_fix_found <-
  as.data.frame(MG17_SNP@fix[MG17_SNP_list$ID %in% c(rs10792832_combined_hg38$ID,
                                                    "rs10792832"), ])

MG17_SNP_list_combined <-
  cbind(MG17_SNP_rs10792832_fix_found,
        MG17_SNP_rs10792832_gt_found)

MG17_SNP_list_combined <-
  MG17_SNP_list_combined[rowSums(!is.na(MG17_SNP_list_combined[, 10:26])) > 0, ]
# is.na(MG17_SNP_list_combined[, 10:26])

# rs10792832 (A/G) == rs3851179 (T/C), A=T, G=C
# rs3851179 R2=0.98

# rs10792832 (A/G) == rs561655 (G/A), A=G, G=A
# rs561655 R2=0.651

most_correlated_SNPs <-
  rs10792832_combined_hg38[rs10792832_combined_hg38$ID %in% c("rs3851179",
                                                              "rs561655"), ]

saveRDS(most_correlated_SNPs,
        file = "rs10792832_most_correlated_SNPs_hg38.RDs")
