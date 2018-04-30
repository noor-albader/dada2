##Filter and trim
##truncLen set to 240 for both pairs for now, can adjust as needed
##maxEE set to 2 for F and 5 for R, can adjust as needed

#! /usr/bin/env Rscript

library(dada2)

# File parsing
pathF <- "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/adapter_trimmed_data" 
pathR <- "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/adapter_trimmed_data" 
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") 
fastqFs <- sort(list.files(pathF, pattern="fastq.gz.trimmed"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz.trimmed"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(240,240), maxEE=c(2,5), truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
