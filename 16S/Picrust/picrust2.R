#Picrust analysis 

#loadingss
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("biomformat")

library(phyloseq)
library(Biostrings)
library(biomformat)


phylobiom <- physeq_16S_nit
phylobiom@otu_table <- phylobiom@otu_table[,colSums(phylobiom@otu_table) != 0]

#Write a fasta file
writeXStringSet(refseq(phylobiom), filepath = "F:/16S_HOLOGREEN/core_microbiome/physeq_16S_nit.fasta", format = "fasta")



#Write a biom file 
phylobiom <- t(as.matrix(otu_table(phylobiom)))
write_biom(make_biom(phylobiom), "F:/16S_HOLOGREEN/core_microbiome/physeq_16S_nit.biom")
list.files(pattern = "*.biom")
