# Siwei 14 Feb 2024
# load HiC Data RDS files and export links

# init ####
library(Chicago)

# load data ####
CD11_0hr_5kb <-
  readRDS("/data/FASTQ/HiC_RData/CD11-0_5000bp.Rds")
CD11_1hr_5kb <-
  readRDS("/data/FASTQ/HiC_RData/CD11-1_5000bp.Rds")
CD11_6hr_5kb <-
  readRDS("/data/FASTQ/HiC_RData/CD11-6_5000bp.Rds")


hist(CD11_0hr_5kb@x$score,
     breaks = 100,
     xlim = c(1, 10))
sum(CD11_0hr_5kb@x$score > 2)
sum(CD11_0hr_5kb@x$score > 3)
sum(CD11_0hr_5kb@x$score > 3 & !(is.na(CD11_0hr_5kb@x$distbin)) & !(is.na(CD11_0hr_5kb@x$refBinMean)))
# length(is.na(CD11_1hr_5kb@x$distbin))

sum(CD11_0hr_5kb@x$score > 4)
sum(CD11_0hr_5kb@x$score > 5)
nrow(CD11_0hr_5kb@x)

## export interactions at score = 3
CD11_list <-
  list(CD11_0hr_5kb,
       CD11_1hr_5kb,
       CD11_6hr_5kb)
names(CD11_list)
names(CD11_list) <-
  c("CD11_0hr_5kb",
    "CD11_1hr_5kb",
    "CD11_6hr_5kb")
head(CD11_0hr_5kb@x)#$distbin

bait_223 <-
  CD11_0hr_5kb@x[CD11_0hr_5kb@x$baitID == 223, ]
max(CD11_0hr_5kb@x$baitID)

for (i in 1:length(CD11_list)) {
  print(i)
  exportResults(cd = CD11_list[[i]],
                outfileprefix = paste0(names(CD11_list)[[i]],
                                       "_"),
                cutoff = 3,
                scoreCol = "score",
                format = "washU_text",
                order = "position")
}

