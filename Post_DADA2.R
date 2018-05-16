##READ ME:
#start at the begining 
#run everything in order unless otherwise specified
#optional steps will be noted as optional/can change
#the files produced should match the files in the google drive folder
#some of the tutorials I pulled theses out of (the list is long so  am placing the two used most extensively)
http://master.bioconductor.org/packages/release/workflows/html/rnaseqGene.html#time
http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#log-fold-change-shrinkage-for-visualization-and-ranking

#please add more to the list if you use the tutorial/code/explanation 
##########################################################################
#prep data frame for phyloseq
source('http://bioconductor.org/biocLite.R')
biocLite()
biocLite('phyloseq', type = "sourse")
biocLite('phyloseq')
packageVersion('phyloseq')
library(phyloseq)
taxa <- readRDS("/Users/nooral-bader/Desktop/tax_final.rds")
head(taxa)
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
str(df)
##########################################################################
#transform to phyloseq object
ps <-phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), sample_data(df), tax_table(taxa))
ps
##########################################################################
# DESeq2 conversion
source('http://bioconductor.org/biocLite.R')
biocLite()
biocLite('DESeq2')
library("DESeq2")
packageVersion("DESeq2")
# phyloseq object to deseq2 object
# from here on out we will use DESeq2 object to run normalization
diagdds = phyloseq_to_deseq2(ps, ~ Time + treatement)
diagdds
#function for
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local") # hmm what does local means, why are we estimating size factors, are we lowkey normalizing without knowing it
# add levels to treatment so DESeq2 know what to compare to
# add more levels if 
diagdds$treatement %<>% relevel("basal")
diagdds$treatement
#### now you have a DESeq2 object called diadds ####
##########################################################################

## Optional/can change or play around with below:
# There are 5 ways of normalization done
# 1. trimmed means and replace outliers (I think that this is a normalization type, found in some DESeq2 manuals)
# 2. log2 +1 (no csv print out,  just meanSD plot)
# 3. regulerized log2
# 4. variance stabilizing 
# 5. least likelihood ratio (time course specific normalization) # Sarah if you could focus on understanding this part, that would be grand!
# 6. proportions (this is only done in the very end when producing Bray-Curtis distances)


##########################################################################

####printing out csv files pre-notmalization 
#optional/can change or play around with

##One way to get results of Deseq2 in table format 
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab) #does not look like example in tutorial because OTU so long in row names 
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
head(posigtab)
write.csv(posigtab, file="DESeq2-geomeans-estimatesizefactor-results.csv")# only ~400 OTUs present

##different way of making a table #much better
resdata <- merge(as.data.frame(res), as.data.frame(counts(diagdds,normalized=T)), by='row.names',sort=F)
write.csv(resdata, file="DESeq2-geomeans-estimatesizefactor-results-way2.csv") # ~4,000 OTUs present 

# produce DataFrame of results of statistical tests
# could way to record experimental design
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = "DESeq2-test-conditions.csv")
##########################################################################
# 1. replacing outlier value with estimated value as predicted by distribution using "trimmed mean" approach. 
# recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates
# A way to normalize using trimmed means
ddsClean <- replaceOutliersWithTrimmedMean(diagdds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(diagdds)$padj < 0.1,
             cleaned = results(ddsClean)$padj < 0.1)
addmargins(tab)
#print trimmed means
write.csv(as.data.frame(tab),file = 'DESeq2-replaceoutliers.csv')
# 
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
write.csv(as.data.frame(resClean),file = 'DESeq2-replaceoutliers-results.csv')
##########################################################################
# DESeq2's normalization methods: (3. and 4.)
rld <- rlogTransformation(diagdds, blind=T) #long time for >50 samples
vsd <- varianceStabilizingTransformation(diagdds, blind=T)#much faster than rlog transformation 
# save normalized results in csv files
write.table(as.data.frame(assay(rld)), file = 'DESeq2-rlog-transformed-counts.txt', sep = '\t')
write.table(as.data.frame(assay(vsd)), file = 'DESeq2-vst-transformed-counts.txt', sep = '\t')
##########################################################################
# DESeq2's time series normalizaton method: (5.)
# must start from phyloseq object and add interaction term (Time:treatement) to design option
psdeseq2 = phyloseq_to_deseq2(ps, ~ Time + treatement + Time:treatement)
geoMeans = apply(counts(psdeseq2), 1, gm_mean)
psdeseq2 = estimateSizeFactors(psdeseq2, geoMeans = geoMeans)
ddsTC <- DESeq(psdeseq2, test="LRT", reduced = ~ Time + treatement)
resTC <- results(ddsTC)
head(resTC[order(resTC$padj),], 4)
#
fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("Time","treatement"), returnData = TRUE)
#
ggplot(fiss,
       aes(x =Time, y = count, color = treatement, group = treatement)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10()
#loess= the smoothing method, using less than 1,00 observation, based on the size of the largest grp depending on all grps (?)
resultsNames(ddsTC)
betas <- coef(ddsTC)
colnames(betas)
#
topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)
##########################################################################

# Differnt plotting (dispersion, Mean vs SD, MA plot and PcoA plots) scripts (examples)
# the libraries that need to be called should be per block, so there will also be repeats. 
# optional/can change or play around with:


##########################################################################
##bar plot showing the log2-fold-change, showing Genus and Phylum
#optional/can change or play around with
library("ggplot2")
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
head(sigtabgen)
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggplot(sigtabgen, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
##########################################################################
# plot ordinante plots and pull out alpha diversity
library("ggplot2")
library("ape") #for phylogeny tree (?)
library("plyr") #for ordinate plots (?)
###ordinate plots pre transformation
# visualize alpha-diversity 
plot_richness(ps, x="Time", measures=c("Shannon", "Simpson"), color="Time")
# transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Time", title="Bray NMDS")
### Ordination did not pick-out clear seperation between timepoints
##########################################################################
#MA plot
#estimated size factor with geo means plots:
pdf(file="/Users/nooral-bader/Desktop/DESeq2_MAplot_.pdf")
plotMA(diagdds,ylim=c(-2,2),main="DESeq2")
dev.off()
##########################################################################
# MA plot
pdf(file="DESeq2_variance_stabilizing.pdf")
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()
##The x axis is the square root of variance over the mean for all samples, so this will naturally included variance due to the treatment. 
##The goal here is to flatten the curve so that there is consistent variance across the read counts, and that is what we got.
##########################################################################
# Principal components plot
# will show additional clustering of samples
# showing basic PCA function in R from DESeq2 package
# this lacks sample IDs and only broad sense of sample clustering
# its not nice - but it does the job
# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
pdf(file="/Users/nooral-bader/Desktop/DESeq2_PCA_initial_analysis.pdf")
print(plotPCA(diagdds, intgroup = c("Time")))
dev.off()
##########################################################################
#MeanSD plots
# this is the three transformations from DESeq 2 tutorial
#biocLite()
#biocLite('vsn')
library("vsn")
#biocLite('hexbin')
library("hexbin")
meanSdPlot(cts, ranks = FALSE)
meanSdPlot(diagdds.filtered, ranks=FALSE)
##logrithmicly transformed 
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)
## this gives log2(n + 1)
ntd <-normTransform(diagdds.filtered)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))
##########################################################################
# plot to show effect of transformation:
# NOTE:  first simply using the log2 function (after adding 1, to avoid taking the log of zero), and then using the VST and rlog-transformed values. 
# For the log2 approach, we need to first estimate size factors to account for sequencing depth, and then specify normalized=TRUE. 
# Sequencing depth correction is done automatically for the vst and rlog.
#dds <- estimateSizeFactors(dds) #did already; done automatically is vsd and rlog
biocLite('dplyr')
library("dplyr")
library("ggplot2")
df <- bind_rows(
  as_data_frame(log2(counts(diagdds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 
##########################################################################
#Scatterplot of transformed counts from two samples. Shown are scatterplots using the log2 transform of normalized counts (left), using the VST (middle), and using the rlog (right). 
#While the rlog is on roughly the same scale as the log2 counts, the VST has a upward shift for the smaller values. 
#It is the differences between samples (deviation from y=x in these scatterplots) which will contribute to the distance calculations and the PCA plot.
#We can see how genes with low counts (bottom left-hand corner) seem to be excessively variable on the ordinary logarithmic scale, while the VST and rlog compress differences for the low count genes for which the data provide little information about differential expression.
#
#R function dist to calculate the Euclidean distance between samples. 
#To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
#We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.
sampleDists <- dist(t(assay(rld)))
sampleDists
biocLite('pheatmap')
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatement, rld$Time, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
#Note that we have changed the row names of the distance matrix to contain treatement type and time instead of sample/lane ID, so that we have all this information in view when looking at the heatmap.
##The x axis is the square root of variance over the mean for all samples, so this will naturally included variance due to the treatment. 
##The goal here is to flatten the curve so that there is consistent variance across the read counts, and that is what we got.
##########################################################################
# Principal components plot
# will show additional clustering of samples
# showing basic PCA function in R from DESeq2 package
# this lacks sample IDs and only broad sense of sample clustering
# its not nice - but it does the job
pdf(file="/Users/nooral-bader/Desktop/DESeq2_PCA_rld.pdf")
print(plotPCA(rld, intgroup = c("Time")))
dev.off()
#
pdf(file="/Users/nooral-bader/Desktop/DESeq2_PCA_vsd.pdf")
print(plotPCA(vsd, intgroup = c("Time")))
dev.off()
#
pdf(file="/Users/nooral-bader/Desktop/DESeq2_PCA_rld_2.pdf")
print(plotPCA(rld, intgroup = c("Time","treatement")))
dev.off()
#
pca_rld_2<- plotPCA(rld, intgroup = c("Time","treatement"), returnData=TRUE)
percentVar <- round(100 * attr(pca_rld_2, "percentVar"))
ggplot(pca_rld_2, aes(PC1, PC2, color=Time, shape=treatement)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#
pdf(file="/Users/nooral-bader/Desktop/DESeq2_PCA_vsd_2.pdf")
print(plotPCA(vsd, intgroup = c("Time","treatement")))
dev.off()
#
pca_vsd_2<- plotPCA(vsd, intgroup = c("Time","treatement"), returnData=TRUE)
percentVar <- round(100 * attr(pca_vsd_2, "percentVar"))
ggplot(pca_vsd_2, aes(PC1, PC2, color=Time, shape=treatement)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
