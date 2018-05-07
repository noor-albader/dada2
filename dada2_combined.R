#! /usr/bin/env Rscript

##Filter and trim
##truncLen set to 240 for both pairs for now, can adjust as needed
##maxEE set to 2 for F and 5 for R, can adjust as needed

##Outputs are filtered fastq file and table of read lengths before and after filter

library("dada2")

pathF <- "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/adapter_trimmed_data/R1_trimmed"
pathR <- "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/adapter_trimmed_data/R2_trimmed"
fnFs <- sort(list.files(pathF, pattern="_R1_001.fastq.gz.trimmed", full.names = TRUE))
fnRs <- sort(list.files(pathR, pattern="_R2_001.fastq.gz.trimmed", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# File parsing
pathF <- "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/adapter_trimmed_data/R1_trimmed" 
pathR <- "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/adapter_trimmed_data/R2_trimmed" 
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") 
fastqFs <- sort(list.files(pathF, pattern="R1_001.fastq.gz.trimmed"))
fastqRs <- sort(list.files(pathR, pattern="R2_001.fastq.gz.trimmed"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering
reads <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(240,240), maxEE=c(2,5), truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
write.table(x = reads, "filtered_table.tsv", append = FALSE, sep = "\t")


##Learn error rates, sample inference (dereplicate), and merge reads
##produces sequence table

##Output is ASV sequence table in .rds format

# File parsing
filtpathF <- "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/adapter_trimmed_data/R1_trimmed/filtered" 
filtpathR <- "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/adapter_trimmed_data/R2_trimmed/filtered" 
filtFs <- list.files(filtpathF, pattern="R1_001.fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="R2_001.fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) 
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100) #not sure what this does...

# Learn forward error rates
errF <- learnErrors(filtFs, nread=1e6, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nread=1e6, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/seqtab.rds")
foo <- readRDS("/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/seqtab.rds")
#write.csv(foo,"") #if we want to export this matrix
#head(foo) #view seqtab



##Removes chimeras and assigns taxonomy
##Uses Silva database
## maintained and formatted by dada2
## https://benjjneb.github.io/dada2/training.html
##Output is .csv and .rds files of assigned taxonomy
## plus .csv file of read counts throughout pipeline

##At this point there are commands to merge multiple runs together, if we decide to do this
## But we shouldn't need to do this because we will only have one seqtab table from one "run"

st.all <- readRDS("/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/seqtab.rds")

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
#track <- cbind(reads, sapply(ddF, getN), sapply(ddR, getN), sapply(merger, getN), rowSums(seqtab))
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
track <- cbind(reads, rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "merged", "nochim")
rownames(track) <- sample.names
#head(track)
write.csv(track,file="read_track.csv")

# Assign taxonomy to genus level, using training set
tax <- assignTaxonomy(seqtab.nochim, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# Add species level assignment to taxonomy table
tax <- addSpecies(tax, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/silva_species_assignment_v132.fa.gz")

# Write to disk
saveRDS(seqtab.nochim, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/tax_final.rds") 
write.csv(seqtab,file="seqtab_final.csv")
write.csv(tax,file="tax_final.csv")
