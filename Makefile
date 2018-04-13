7#
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

all: 
	@echo " This directory is used to process 16s sequences"
	@echo " The current directory is: $(DIRPATH)
	@echo ""
	@echo " make dada2	: runs R script of dada2 pipeline on sequence files"

dada2
	@echo "Use the R script dada2.R to process 16s sequences using the dada2 pipeline"
	@echo "Output files will be ASV table and phyloseq graphs"
	@echo for dir in $(ALLDIR); do \
	cd $${dir}; \
	SGE_Batch -c "Rscript /path/dada2.R $(SEQID)" -r dada2-$(SEQID)-StdOut $(SGE_OPTIONS);\
	echo ""; \
	cd ..; \
	done ;

