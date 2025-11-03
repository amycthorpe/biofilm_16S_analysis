# DADA2 for big data
# EA biofilms 16S - merging runs
# repeated for consistency for manuscript 24.03.25


library(dada2); packageVersion("dada2")

set.seed(100)

output <- "/hdd0/amytho/EA_biofilm_nextseqruns/EA_16S_all_DADA2_repeated_for_manuscript"


# load seq tabs
seqtab.optim <- read.table("/hdd0/amytho/EA_biofilm_nextseqruns/EA_16S_all_DADA2_repeated_for_manuscript/run_optimisation/output/seqtab_16S.csv",header=T,row.names=1,sep=",",check.names=FALSE)
seqtab.2021_22 <- read.table("/hdd0/amytho/EA_biofilm_nextseqruns/EA_16S_all_DADA2_repeated_for_manuscript/run_2021_22/output/seqtab_16S.csv",header=T,row.names=1,sep=",",check.names=FALSE)
seqtab.2023 <- read.table("/hdd0/amytho/EA_biofilm_nextseqruns/EA_16S_all_DADA2_repeated_for_manuscript/run_2023/output/seqtab_16S.csv",header=T,row.names=1,sep=",",check.names=FALSE)

# transpose so samples are rows
seqtab.optim <- t(seqtab.optim)
seqtab.2021_22 <- t(seqtab.2021_22)
seqtab.2023 <- t(seqtab.2023)

# subset to only relevant samples - zymo extraction kit
seqtab.optim <- seqtab.optim[grep("^Z", rownames(seqtab.optim)), ]


# merge runs
seqtab.all <- mergeSequenceTables(seqtab.optim, seqtab.2021_22, seqtab.2023)


# remove chimeras after merge
seqtab.all.nochim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=TRUE, minFoldParentOverAbundance=4)

dim(seqtab.all.nochim) # check table dimensions
sum(seqtab.all.nochim)/sum(seqtab.all) # check proportion non-chimera


# transpose and save merged seqtab
seqtab.all.nochim_t <- t(seqtab.all.nochim)
write.csv(seqtab.all.nochim_t, file.path(output, "seqtab16S_merged.csv"))


# assign taxonomy
taxa <- assignTaxonomy(seqtab.all.nochim, "/hdd0/amytho/training datasets/silva_nr99_v138.2_toGenus_trainset.fa.gz", minBoot = 80) # downloaded from https://zenodo.org/records/14169026

taxa.print <- taxa # removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) # check taxonomy

asv_seqs <- colnames(seqtab.all.nochim)
asv_headers <- vector(dim(seqtab.all.nochim)[2], mode="character")

for (i in 1:dim(seqtab.all.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# get fasta
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(output, "ASVs_16S_merged.fa"))

# get counts tab
asv_tab <- t(seqtab.all.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file.path(output, "counts_16S_merged.tsv"), sep="\t", quote=F, col.names=NA)

# get tax tab
rownames(taxa) <- gsub(pattern=">", replacement="", x=asv_headers) 
write.table(taxa, file.path(output, "taxonomy_16S_merged.tsv"), sep = "\t", quote=F, col.names=NA)

# species
rownames(taxa) <- getSequences(seqtab.all.nochim)
taxa_sp <- addSpecies(taxa, "/hdd0/amytho/training datasets/silva_v138.2_assignSpecies.fa.gz") # add species

taxa.print <- taxa_sp # removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) # check taxonomy

rownames(taxa_sp) <- gsub(pattern=">", replacement="", x=asv_headers) 
write.table(taxa_sp, file.path(output, "taxonomy_species_16S_merged.tsv"), sep = "\t", quote=F, col.names=NA)


