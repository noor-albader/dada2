#
# Makefile
# Date: 12 Apr 2018
# This makefile is for running the dada2 R script to process 16s sequences

#Set Variables
#Modify SEQID once we know how our sequence files will be named
DIRPATH = `pwd`
#PROCS =
SGE OPTIONS = -q 'stevens'
#SEQID=`ls -1 *R1_001.fastq* | sed -r 's/(.*)BSL(.*)/BSL\2/' | sed -r 's/(.*)-([0-9]+)_(.*)/\1-\2/g'`
#ALLDIR=`ls -1 | grep -v Makefile | grep Sample`
###instead of running all of the pipline with one make command, I propose this following
all: 
	@echo " This directory is used to process 16s sequences"
	@echo " The current directory is: $(DIRPATH)
	@echo ""
	@echo " make file and path : runs R script for reading in fastq files and string manipulation for forwas and reverse files"
	@echo " make quality profiles : runs R script to plot raw qaulity for F and R file"
	@echo " make filter	: run R script for filtering"
	@echo " make error rate  : run R script to run paramertic error model for error rates"
	@echo " make error profile : run R script to plot estimated error rates"
	@echo " make dereplication : run R script to combine identical sequencing and calculate abundance"
	# not sure what this step does conceptually 
	@echo " make sample inference : run R script to the core sequence-variant inference" 
	# might want to inspect object creating to look at denoising results fro F and R 
	@echo " make merged paired reads : run R script to merge overlapping reads
	@echo " make sequence table : run R script to create seq table of seq variants"
	@echo " make remove chimeras : run R script to remove chimeras"
	@echo " make track reads : run R script to visualize read filtering thus far"
	@echo " make taxonomy	: run R script to assign taxonomy"
	@echo " make accuracy	 : run R script if comparing to known communities is desired"
	###EXTRA VISUALIZATION
	@echo " make alpha diversity : run R script to plot alpha diversity"
	@echo " make ordinate : run R script to transform and standarize results"
	@echo " make ordinate plot : run R script to plot transformed results"
	@echo " make family plot " run R script to plot top 20 taxonomic group distribution across sample"
dada2
	@echo "Use the R script dada2.R to process 16s sequences using the dada2 pipeline"
	@echo "Output files will be ASV table and phyloseq graphs"
	@echo for dir in $(ALLDIR); do \
	cd $${dir}; \
	SGE_Batch -c "Rscript /path/dada2.R $(SEQID)" -r dada2-$(SEQID)-StdOut $(SGE_OPTIONS);\
	echo ""; \
	cd ..; \
	done ;

