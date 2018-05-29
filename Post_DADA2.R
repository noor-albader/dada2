library(phyloseq)
library(DESeq2)
library(vegan)
##########################################################################
#prep data frame for phyloseq
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
df$treatment <- c("rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","rice","rice","basal","basal","rice","rice","basal","basal","basal","basal","basal","basal","rice","basal","basal","rice","rice","rice","rice","basal","basal")
df$sex <- c("Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female","Female","Female","Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female","Female","Female","Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female","Female","Female","Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female","Female","Female","Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female","Female","Female","Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Male",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female",	"Female","Female","Female")
rownames(df) <-samples.out
head(df)
str(df)
##########################################################################
#transform to phyloseq object and subset
ps <-phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), sample_data(df), tax_table(taxa))
ps
rank_names(ps)
sample_variables(ps)
otu_table(ps)[1:5, 1:5]
taxa_names(ps)
ps_time1 <- subset_samples(ps, time =="Timepoint1")
ps_time2 <- subset_samples(ps, time == "Timepoint2")
ps_time3 <- subset_samples(ps, time== "Timepoint3")
ps_time4 <- subset_samples(ps, time== "Timepoint4")
ps_time5 <- subset_samples(ps, time== "Timepoint5")
ps_time6 <- subset_samples(ps, time== "Timepoint6")
##########################################################################
#transform subsetted phyloseq object's t DESeq2 objects 
dds_all = phyloseq_to_deseq2(ps, ~treatment)
dds = phyloseq_to_deseq2(ps_time1, ~treatment)
dds_2 =phyloseq_to_deseq2(ps_time2, ~treatment)
dds_3 =phyloseq_to_deseq2(ps_time3, ~treatment)
dds_4 =phyloseq_to_deseq2(ps_time4, ~treatment)
dds_5 =phyloseq_to_deseq2(ps_time5, ~treatment)
dds_6 =phyloseq_to_deseq2(ps_time6, ~treatment)
##########################################################################
#DESeq2 - part 1: estimate size factors
ps_dds_all <- estimateSizeFactors(dds_all, type ="poscounts")
ps_dds_t1 <- estimateSizeFactors(dds, type = "poscounts")
ps_dds_t2 <- estimateSizeFactors(dds_2, type = "poscounts")
ps_dds_t3 <- estimateSizeFactors(dds_3, type = "poscounts")
ps_dds_t4 <- estimateSizeFactors(dds_4, type = "poscounts")
ps_dds_t5 <- estimateSizeFactors(dds_5, type = "poscounts")
ps_dds_t6 <- estimateSizeFactors(dds_6, type = "poscounts")
##########################################################################
#DESeq - part 2: get normalized results
dds_norm <- DESeq(ps_dds_all,fitType="parametric", betaPrior = FALSE) 
#####
dds_2_norm <- DESeq(ps_dds_t2,fitType="parametric", betaPrior = FALSE) 
dds_2_norm$treatement <- relevel(dds_2_norm$treatement, ref = "basal")
res = results(dds_2_norm)
res$padj[is.na(res$padj)] = 1
sig <- res[res$padj <.05,]

sigtab = cbind(as(sig, "data.frame"), as(tax_table(ps_time2)[rownames(sig), ], "matrix"))
write.csv(sigtab, file="/Users/nooral-bader/Desktop/DESeq2-results-t2.csv") 

tab = cbind(as(res, "data.frame"), as(tax_table(ps_time2)[rownames(res), ], "matrix"))
write.csv(tab, file="/Users/nooral-bader/Desktop/DESeq2-unfiltered-results-t2.csv") 
#####
dds_3_norm <- DESeq(ps_dds_t3,fitType="parametric", betaPrior = FALSE) 
dds_3_norm$treatement <- relevel(dds_3_norm$treatement, ref = "basal")
res = results(dds_3_norm)
res$padj[is.na(res$padj)] = 1
sig <- res[res$padj <.05,]

sigtab = cbind(as(sig, "data.frame"), as(tax_table(ps_time3)[rownames(sig), ], "matrix"))
write.csv(sigtab, file="/Users/nooral-bader/Desktop/DESeq2-results-t3.csv") 

tab = cbind(as(res, "data.frame"), as(tax_table(ps_time3)[rownames(res), ], "matrix"))
write.csv(tab, file="/Users/nooral-bader/Desktop/DESeq2-unfiltered-results-t3.csv") 
#####
dds_4_norm <- DESeq(ps_dds_t4,fitType="parametric", betaPrior = FALSE) 
dds_4_norm$treatement <- relevel(dds_4_norm$treatement, ref = "basal")
res = results(dds_4_norm)
res$padj[is.na(res$padj)] = 1
sig <- res[res$padj <.05,]

sigtab = cbind(as(sig, "data.frame"), as(tax_table(ps_time4)[rownames(sig), ], "matrix"))
write.csv(sigtab, file="/Users/nooral-bader/Desktop/DESeq2-results-t4.csv") 

tab = cbind(as(res, "data.frame"), as(tax_table(ps_time4)[rownames(res), ], "matrix"))
write.csv(tab, file="/Users/nooral-bader/Desktop/DESeq2-unfiltered-results-t4.csv") 
#####
dds_5_norm <- DESeq(ps_dds_t5,fitType="parametric", betaPrior = FALSE) 
dds_5_norm$treatement <- relevel(dds_5_norm$treatement, ref = "basal")
res = results(dds_5_norm)
res$padj[is.na(res$padj)] = 1
sig <- res[res$padj <.05,]

sigtab = cbind(as(sig, "data.frame"), as(tax_table(ps_time5)[rownames(sig), ], "matrix"))
write.csv(sigtab, file="/Users/nooral-bader/Desktop/DESeq2-results-t5.csv") 

tab = cbind(as(res, "data.frame"), as(tax_table(ps_time5)[rownames(res), ], "matrix"))
write.csv(tab, file="/Users/nooral-bader/Desktop/DESeq2-unfiltered-results-t5.csv") 
#####
dds_6_norm <- DESeq(ps_dds_t6,fitType="parametric", betaPrior = FALSE) 
dds_6_norm$treatement <- relevel(dds_6_norm$treatement, ref = "basal")
res = results(dds_6_norm)
res$padj[is.na(res$padj)] = 1
sig <- res[res$padj <.05,]

tab = cbind(as(res, "data.frame"), as(tax_table(ps_time6)[rownames(res), ], "matrix"))
write.csv(tab, file="/Users/nooral-bader/Desktop/DESeq2-unfiltered-results-t6.csv") 

sigtab = cbind(as(sig, "data.frame"), as(tax_table(ps_time6)[rownames(sig), ], "matrix"))
write.csv(sigtab, file="/Users/nooral-bader/Desktop/DESeq2-results-t6.csv") 
##########################################################################
#DESeq - part 2 alternative: use estimated size factors to (1) estimate sispersions and (2) normalize with vst
ps_dds_all <- estimateDispersions(ps_dds_all)
abund_all <- getVarianceStabilizedData(ps_dds_all)
abund_all <- abund_all + abs(min(abund_all)) #don't allow deseq to return negative counts
#####
ps_dds_t1 <- estimateDispersions(ps_dds_t1)
abund_t1 <- getVarianceStabilizedData(ps_dds_dis)
abund_t1 <- abund_t1 + abs(min(abund_t1)) #don't allow deseq to return negative counts
#####
ps_dds_t2 <- estimateDispersions(ps_dds_t2)
abund_t2 <- getVarianceStabilizedData(ps_dds_tmp)
abund_t2 <- abund + abs(min(abund)) #don't allow deseq to return negative counts
#####
ps_dds_t3 <- estimateDispersions(ps_dds_t3)
abund_t3 <- getVarianceStabilizedData(ps_dds_t3)
abund_t3 <- abund_t3 + abs(min(abund_t3)) #don't allow deseq to return negative counts
#####
ps_dds_t4  <- estimateDispersions(ps_dds_t4)
abund_t4 <- getVarianceStabilizedData(ps_dds_t4)
abund_t4 <- abund_t4 + abs(min(abund_t4)) #don't allow deseq to return negative counts
#####
ps_dds_t5  <- estimateDispersions(ps_dds_t5)
abund_t5 <- getVarianceStabilizedData(ps_dds_t5)
abund_t5 <- abund_t5 + abs(min(abund_t5)) #don't allow deseq to return negative counts
#####
ps_dds_t6  <- estimateDispersions(ps_dds_t6)
abund_t6 <- getVarianceStabilizedData(ps_dds_t6)
abund_t6 <- abund_t6 + abs(min(abund_t6)) #don't allow deseq to return negative counts

##########################################################################
#Back to to phyloseq with normalized abundance 
ps_deSeq_all <- phyloseq(otu_table(abund_all, taxa_are_rows = T), sample_data(ps), tax_table(ps))
ps_deSeq_1 <- phyloseq(otu_table(abund_t1, taxa_are_rows = T), sample_data(ps_time1), tax_table(ps_time1))
ps_deSeq_2 <- phyloseq(otu_table(abund, taxa_are_rows = T), sample_data(ps_time2), tax_table(ps_time2))
ps_deSeq_3 <- phyloseq(otu_table(abund_t3, taxa_are_rows = T), sample_data(ps_time3), tax_table(ps_time3))
ps_deSeq_4 <- phyloseq(otu_table(abund_t4, taxa_are_rows = T), sample_data(ps_time4), tax_table(ps_time4))
ps_deSeq_5 <- phyloseq(otu_table(abund_t5, taxa_are_rows = T), sample_data(ps_time5), tax_table(ps_time5))
ps_deSeq_6 <- phyloseq(otu_table(abund_t6, taxa_are_rows = T), sample_data(ps_time6), tax_table(ps_time6))

##########################################################################
## PCoA plots: ###
#non constrained
ps_pcoa_non<-ordinate(
  physeq =ps_deSeq_all,
  method='MDS',
  distance="bray"
)

#constrained
ps_pcoa_constraint <-ordinate(
  physeq =ps_deSeq_all ,
  method='CAP',
  distance="bray",
  formula = ~ treatment + sex
)
### Plot the PCoA plots: ###

plot_ordination(
  physeq =ps_deSeq_all,
  ordination = ps_pcoa_constraint,
  color = "treatment",
  shape = "sex",
  axes= c(1,2),
  title = "Constraint PCoA"
)

plot_ordination(
  physeq =ps_deSeq_6,
  ordination = ps_pcoa_non,
  color = "treatment",
  shape = "sex",
  axes= c(1,2),
  title = "Nonconstraint PCoA"
)
##########################################################################
