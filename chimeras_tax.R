##Removes chimeras and assigns taxonomy
##will need to get Silva db
##See https://www.arb-silva.de/projects/ssu-ref-nr/ , under downloads, want archive fastq file
## believe found here https://www.arb-silva.de/no_cache/download/archive/release_132/Exports/
## titled SILVA_132_SSURef_Nr99_tax_silva.fasta.gz

#! /usr/bin/env Rscript

library(dada2)

##At this point there are commands to merge multiple runs together, if we decide to do this

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/SILVA_132_SSURef_Nr99_tax_silva.fasta.g", multithread=TRUE)
# Write to disk
saveRDS(seqtab, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "/ACTF/Course/mb599_bds_s18/data/share/group3_Rothenberg/tax_final.rds") 