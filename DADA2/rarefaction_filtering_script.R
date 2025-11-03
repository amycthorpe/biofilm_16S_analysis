# filtering and rarefaction of repeated DADA2 16S EA biofilms for manuscript


# load libraries, set wd
library(plyr)
library(dplyr)
library(vegan)
library(seqinr)

setwd("/hdd0/amytho/EA_biofilm_nextseqruns/EA_16S_all_DADA2_repeated_for_manuscript/")

set.seed(123)


# import data (change file paths)
tax_tab<-read.table("output/taxonomy_16S_merged.tsv",header=T,row.names=1,check.names=F,sep="\t")
counts_tab<-read.table("output/counts_16S_merged.tsv",header=T,row.names=1,sep="\t",check.names=FALSE)
meta<-read.table("seq_metadata_all_amplicons.csv",header=T,sep=",",fill=T,check.names=F)

meta <- meta[!is.na(meta$`16S_ID`) & meta$`16S_ID` != "", ] # remove non-16S samples

# correct taxonomy headings (18S only)
#colnames(tax_tab) <- c("Domain", "Supergroup", "Division", "Phylum", "Class", "Order", "Family", "Genus", "Species")


# subset dataset
counts_tab <- counts_tab[, colnames(counts_tab) %in% meta$`16S_ID`]
counts_tab <- counts_tab[rowSums(counts_tab[])>0, ]
tax_tab <- tax_tab[rownames(tax_tab) %in% rownames(counts_tab), ]


# ASVs per sample (good to check number of reads for any blanks or samples which did not work well)
ASVs_per_sample<-as.data.frame(colSums(counts_tab))
ASVs_per_sample_meta<-merge(ASVs_per_sample, meta, by.x="row.names", by.y="16S_ID")
write.csv(ASVs_per_sample,"output/filtering_rarefaction/ASVs_per_sample.csv")


# filtering
# filter based on taxonomy
tax_tab[is.na(tax_tab)] <- "unassigned" # NA must be converted to unassigned first otherwise o_NAs are removed alongside Chloroplast
tax_tab_filt<-subset(tax_tab,Kingdom=="Bacteria") # eukaryota = 2847 ASVs, archaea = 2600 ASVs, unassigned = 5987 ASVs
tax_tab_filt<-subset(tax_tab_filt,Order!="Chloroplast") # chloroplast = 6212 ASVs
tax_tab_filt<-subset(tax_tab_filt,Phylum!="unassigned") # unassigned = 70314

# subset counts based on filtered tax (so both tax and counts tables have the same samples in after filtering)
counts_tab_filt<-counts_tab[rownames(counts_tab)%in%rownames(tax_tab_filt),]

# filter based on abundance
counts_tab_filt2<-counts_tab_filt[rowSums(counts_tab_filt)>5,] # remove ASVs with <X reads
counts_tab_filt3<-counts_tab_filt2[rowSums(counts_tab_filt2>0)>=5,] # remove ASVs that only appear in <X samples


# rarefaction
# check row sums (samples need to be rows so transpose first)
counts_tab_filt3<-t(counts_tab_filt3)
rowsums<-as.data.frame(rowSums(counts_tab_filt3))
plot(sort(rowsums[,1]))

min(rowsums[,1])
mean(rowsums[,1])
hist(rowSums(counts_tab_filt3))

# rarefaction curve
plot<-rarecurve(counts_tab_filt3, step = 1000, col = "blue", cex = 0.7, subsample = 1000, label = FALSE, xlab = "Number of Sequences",
                ylab = "Observed Species Richness",
                main = "Rarefaction Curve")
plot

# rarefy
counts_tab_rare<-rrarefy(counts_tab_filt3, 10000)
counts_tab_rare<-as.data.frame(counts_tab_rare)

# remove samples below rarefied depth (if any)
counts_tab_rare<-counts_tab_rare[rowSums(counts_tab_rare[])==10000,]

# check rarefied rowsums - min and mean should now be the same
rowsums<-as.data.frame(rowSums(counts_tab_rare))
plot(sort(rowsums[,1]))

min(rowsums[,1]) 
mean(rowsums[,1]) 

# remove 0 abundance ASVs from counts after rarefying (if samples below rarefaction depth were removed, these are ASVs only in those samples)
counts_tab_rare <- counts_tab_rare[rownames(counts_tab_rare) != "EA-AH12-Blank17-6", ]
counts_tab_rare<-counts_tab_rare[,colSums(counts_tab_rare[])>0]

# subset tax to new rarefied ASV list (only if you removed samples)
counts_tab_rare<-as.data.frame(t(counts_tab_rare))
tax_tab_rare<-tax_tab[rownames(tax_tab)%in%rownames(counts_tab_rare),]

# subset meta to rarefied sample list (only if you removed samples)
counts_tab_rare<-as.data.frame(t(counts_tab_rare))
meta_rare<-meta[meta$`16S_ID` %in% rownames(counts_tab_rare), ]


# save filtered and rarefied tables (change file paths)
write.csv(counts_tab_rare,"output/filtering_rarefaction/ASV_counts_16S_processed.csv")
write.csv(tax_tab_rare,"output/filtering_rarefaction/ASV_tax_16S_processed.csv")
write.csv(meta_rare,"output/filtering_rarefaction/meta_16S_processed.csv")
