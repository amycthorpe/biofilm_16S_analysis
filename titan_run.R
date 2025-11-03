# biofilm analysis for manuscripts - 16S
# titan on R bash


# load libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(TITAN2)


set.seed(123)


path <- "path/to/data_directory"


# import data
ASV_tab <- data.frame(fread(file.path(path, "data/ASV_counts_16S_processed.csv"), sep = ",", header = TRUE, data.table = FALSE), check.names=F, row.names=1)
tax_tab <- read.table(file.path(path, "data/ASV_tax_16S_processed.csv"), header=T, row.names=1, check.names=F, sep=",")
meta_all <- read.table(file.path(path, "data/seq_metadata_all_amplicons.csv"), header=T, check.names=F, sep=",")

meta <- meta_all[meta_all$`16S_ID` %in% rownames(ASV_tab), ] # subset metadata to seq ID (change depending on amplicon)
rownames(meta) <- NULL
meta <- column_to_rownames(meta, var="16S_ID") # make seq ID rownames

meta <- meta %>% # add season columm
  mutate(
    season = case_when(
      month(Sample_date) %in% c(03, 04, 05, 06, 07, 08) ~ "Spring",
      month(Sample_date) %in% c(09, 10, 11, 12, 01, 02) ~ "Autumn"
    )
  )

meta$Sample_date <- as.Date(meta$Sample_date, format = "%d/%m/%Y")  # change format to date


# convert to relative abundance
ASV_rel <- ASV_tab %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(across(everything(), ~ . / sum(.))) %>%
  t(.) # samples as rows


# environmental metadata
meta_merged <- read.csv(file.path(path, "data/meta_merged.csv"), header = T, row.names = 1, check.names = F) %>% 
  select(Alkalinity_combined, Nitrate_N, DO, Water_temp)

meta_merged <- as.data.frame(na.omit(meta_merged)) # remove rows with NA

variables <- "Nitrate_N" # subset metadata to numeric variables of interest
meta_sub <- meta_merged[variables]

ASV_sub <- as.data.frame(ASV_rel[rownames(ASV_rel) %in% rownames(meta_sub), ]) # subset to remaining samples
meta_sub <- as.data.frame(meta_sub[intersect(rownames(meta_sub), rownames(ASV_rel)), , drop=F]) # subset to remaining samples
ASV_sub <- ASV_sub[match(rownames(meta_sub), rownames(ASV_sub)), ] # match order of rows


# test a subset
# meta_sub2 <- meta_sub[sample(rownames(meta_sub), 20), , drop = FALSE]
# ASV_sub2 <- ASV_sub[rownames(ASV_sub) %in% rownames(meta_sub2), ]
# ASV_sub2 <- ASV_sub2[match(rownames(meta_sub2), rownames(ASV_sub2)), ] # match order of rows
# ASV_sub3 <- ASV_sub2[, colSums(ASV_sub2 > 0) >= 3] # minimum occurrence frequency for titan is 3
# 
# start <- date()
# test <- titan(meta_sub2, ASV_sub3, minSplt=5, numPerm=500, boot=T, nBoot=500, imax=F, ivTot=F, pur.cut=0.95, rel.cut=0.95, ncpus=52)
# end <- date()
#
# start
# end
# 
# saveRDS(test, file = file.path(path, "output/titan_test.rds"))


# filter
# ASV_sub4 <- ASV_sub[, colSums(ASV_sub > 0) >= 78] # increase minimum occurrence frequency to present in >5% sites to reduce dataset size and run time - 24067 to 2651 ASVs
ASV_sub4 <- ASV_sub[, colSums(ASV_sub > 0) >= 156] # increase minimum occurrence frequency to present in >10% sites to reduce dataset size and run time - 1287 ASVs
# ASV_sub4 <- ASV_sub[, colSums(ASV_sub > 0) >= 312] # increase minimum occurrence frequency to present in >20% sites to reduce dataset size and run time - 526 ASVs
# ASV_sub4 <- ASV_sub[, colSums(ASV_sub > 0) >= 467] # increase minimum occurrence frequency to present in >30% sites to reduce dataset size and run time - 278 ASVs
# ASV_sub4 <- ASV_sub4[, colSums(ASV_sub4 > 0) >= 3] # minimum occurrence frequency for titan is 3


# run
start <- date()
titan_alk <- titan(meta_sub, ASV_sub4, minSplt=5, numPerm=500, boot=T, nBoot=500, imax=F, ivTot=F, pur.cut=0.95, rel.cut=0.95, ncpus=52)
end <- date()

saveRDS(titan_alk, file = file.path(path, "output/titan/titan_nitrate_10perc.rds"))

start
end
# done
