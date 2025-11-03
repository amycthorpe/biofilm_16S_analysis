# DADA2 for big data
# EA biofilms 16S - 2023 samples
# repeated for consistency for manuscript 24.03.25


library(dada2); packageVersion("dada2")


# File parsing
# cd /hdd0/amytho/EA_optimisation_samples/16S/Data/cutadapt/
# mkdir cutadapt_F
# mkdir cutadapt_R
# cp *R1_001.fastq.gz cutadapt_F -vrf
# cp *R2_001.fastq.gz cutadapt_R -vrf

pathF <- "/hdd0/amytho/EA_biofilm_nextseqruns/EA_biofilm_16S_ITS_2023/16S/Data/cutadapt/cutadapt_F" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/hdd0/amytho/EA_biofilm_nextseqruns/EA_biofilm_16S_ITS_2023/16S/Data/cutadapt/cutadapt_R" # CHANGE ME ...
list.files(pathF)
list.files(pathR)

filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...

fastqFs <- sort(list.files(pathF, pattern="_R1_001.fastq")) # sort ensures forward/reverse reads are in same order
fastqRs <- sort(list.files(pathR, pattern="_R2_001.fastq"))

sample.names <- sapply(strsplit(fastqFs, "_"), `[`, 1)

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")


# Filtering: optimise parameters by sequencing run
fastqFs_p <- sort(list.files(pathF, pattern="_R1_001.fastq", full.names=T)) # sort ensures forward/reverse reads are in same order
fastqRs_p <- sort(list.files(pathR, pattern="_R2_001.fastq", full.names=T))

plotQualityProfile(fastqFs_p[1:10])  # add full.names=TRUE to fastqFs and Rs lines for this to work
plotQualityProfile(fastqRs_p[1:10]) 

out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                     rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                     truncLen=c(230,220), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

write.csv(out,file="/hdd0/amytho/EA_biofilm_nextseqruns/EA_16S_all_DADA2_repeated_for_manuscript/run_2023/output/filtered_out.csv") # check reads out


# File parsing
filtpathF <- "/hdd0/amytho/EA_biofilm_nextseqruns/EA_biofilm_16S_ITS_2023/16S/Data/cutadapt/cutadapt_F/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/hdd0/amytho/EA_biofilm_nextseqruns/EA_biofilm_16S_ITS_2023/16S/Data/cutadapt/cutadapt_R/filtered" # CHANGE ME ...

filtFs <- list.files(filtpathF, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(100)


# Infer sequence variants - error rates are being learned from a subset of the data. Learning error rates is computationally intensive, as it requires multiple iterations of the core algorithm. As a rule of thumb, a million 100nt reads (or 100M total bases) is more than adequate to learn the error rates.
# Learn the error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]], verbose = TRUE)
  ddF <- dada(derepF, err=errF, multithread=TRUE, verbose = TRUE)
  derepR <- derepFastq(filtRs[[sam]], verbose = TRUE)
  ddR <- dada(derepR, err=errR, multithread=TRUE, verbose = TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose = TRUE)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)


# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # check table dimensions
table(nchar(getSequences(seqtab))) # inspect distribution of sequence lengths 

#remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE, minFoldParentOverAbundance=4) #higher minFold value reduces rate of false positive chimeras
dim(seqtab.nochim) # check table dimensions
sum(seqtab.nochim)/sum(seqtab) # check proportion non-chimera (0.72 with trimleft suggesting some primers remain, 0.84 with cutadapt so improved but still relatively high for 16S so included minFold to reduce FP)

seqtab16S=t(seqtab.nochim)
write.csv(seqtab16S, file = "/hdd0/amytho/EA_biofilm_nextseqruns/EA_16S_all_DADA2_repeated_for_manuscript/run_2023/output/seqtab_16S.csv")

#track number of reads that made it through each step in the pipeline. This is a great place to do a last sanity check. Outside of filtering (depending on how stringent you want to be) there should no step in which a majority of reads are lost. If a majority of reads failed to merge, you may need to revisit the  truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon. If a majority of reads failed to pass the chimera check, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
rownames(track) <- sample.names
head(track)

write.csv(track,"/hdd0/amytho/EA_biofilm_nextseqruns/EA_16S_all_DADA2_repeated_for_manuscript/run_2023/output/track_reads_16S_cutadapt.csv")