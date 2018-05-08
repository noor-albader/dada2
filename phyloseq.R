source('http://bioconductor.org/biocLite.R')
biocLite()
#biocLite('phyloseq', type = "sourse")
biocLite('phyloseq')
packageVersion('phyloseq')
library(phyloseq)
sessionInfo()
taxa <- readRDS("/Users/nooral-bader/Desktop/tax_final.rds")
seqtab <- readRDS("/Users/nooral-bader/Desktop/seqtab_final.rds")
samples.out <- rownames(seqtab)
samples.out
subject <- sapply(strsplit(samples.out, "-"), "[", 10)
subject
time <- sapply(strsplit(samples.out, "-"), `[`, 9)
time
df <- data.frame(Subject=subject, Time=time, stringsAsFactors = FALSE)
df$treatement <- c("rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal")
rownames(df) <-samples.out
head(df)
ps <-phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), sample_data(df), tax_table(taxa))
ps
library("ggplot2")
library("ape")
plot_richness(ps, x="Time", measures=c("Shannon", "Simpson"), color="Time")

