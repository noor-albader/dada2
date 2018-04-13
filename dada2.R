#This is a script for running dada2 on 16s sequence data
#Resultant output will be an ASV (amplicon sequence variant) table, plus desired figures from phyloseq
#modified from tutorial found at "https://benjjneb.github.io/dada2/tutorial.html"

#Need to change commands to write output into files

#! /usr/bin/env Rscript

#import packages
library(dada2)
library(phyloseq)
library(ggplot2)

#set variable for file names. This should be included in the R command line in Makefile "Rscript path/dada2.txt $(SEQID)"
seqid <- commandArgs(trailingOnly=TRUE)

#Define path (might not need this???)
path <- "" #insert directory containing sequence fasta files

#Prepare for filter and trimming
#Create a list of forward and reverse fastq files in matched order
fnFs <- sort(list.files(path, pattern="seqid_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="seqid_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Plot the quality profile of the reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#Filter and trim reads
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(out)

#Learn error rates and plot them
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Dereplicate filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Infer sequence variants
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#Merge the forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

#Make sequence table (comparable to OTU table)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Final check that tracks the number of reads at each step of the pipeline
#If large decrease in sequences at one step, might need to adjust settings
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy
#need to adjust to the database we decide to use
taxa <- assignTaxonomy(seqtab.nochim, "Training/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "Training/silva_species_assignment_v128.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Evaluate accuracy
#I don't think we will do this part because we will not have a reference community, right?

#Create plots using phyloseq
#construct data frame
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

#construct phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps

#construct alpha-diversity plots
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When") + theme_bw()

#construct ordinate (NMDS plot)
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
plot_ordination(ps, ord.nmds.bray, color="When", title="Bray NMDS")

#construct bar plot of community abundances
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")