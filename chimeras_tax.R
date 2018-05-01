##Removes chimeras and assigns taxonomy
##Uses Silva database
## maintained and formatted by dada2
## https://benjjneb.github.io/dada2/training.html

#! /usr/bin/env Rscript

library(dada2)

##At this point there are commands to merge multiple runs together, if we decide to do this

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy to genus level, using training set
tax <- assignTaxonomy(seqtab, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# Add species level assignment to taxonomy table
taxa <- addSpecies(taxa, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/silva_species_assignment_v132.fa.gz")

# Write to disk
saveRDS(seqtab, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/tax_final.rds") 
