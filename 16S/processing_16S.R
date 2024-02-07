#######TRAITEMENT OF 16S rRNA SEQUENCE FROM HOLOGREEN PROJECT######


#1. Loadings----
library(dada2)
library(phyloseq)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ShortRead)
library(microDecon)
library(DECIPHER)
library(phangorn)
library(ymd)

#setup the working directory of the project 
#setwd("working_directory") #not necessary when you work with project

#setup the directory of the location of the raw sequences
path <- "raw_data/sequences_16S/raw_seq/"
pathn1<- "raw_data/sequences_16S"
pathobj <- "preprocessing/preprocessing_final/objets/"


#2. Preprocessing----

##2.1 Samples loadings----

#Check if every samples are listed 
list.files(path)

#Load the raw sequences in R1 (forward) and R2 (reverse) vector 
fnR1s <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnR2s <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

#Setup the samples names into the "samples.names" variable
sample.names <- sapply(strsplit(list.files(path, pattern="_R1"),"---IDT_i5_"), "[",2)
sample.names <- sapply(strsplit(sample.names,"_16S"), "[",1)
sample.names <- sapply(strsplit(sample.names,"[.]"), "[",2)
saveRDS(sample.names, file="preprocessing/preprocessing_final/objets/sample.names")
#In this case the samplenames 

#Subset of the samples by samples names
negatif <- c("NEGATIF_F","NEGATIF_A","BlankpcrCES_R1", "NEGATIF_PCR_1", "NEGATIF_PCR_2")
filter_EDM <- grep("EDM_F", sample.names, perl = TRUE, value = TRUE)
filter_ENR <- grep("ENR_F", sample.names, perl = TRUE, value = TRUE)
algae_ENR <- grep("EDM_A_A|EDM_B_A|EDM_C_A", sample.names, value = TRUE)
algae_EDM <- grep("ENR_A_A|ENR_B_A|ENR_C_A", sample.names, value = TRUE)

##2.2 Check the quality of raw sequences----
###2.2.1 Quality check of the R1 (forward)----
#The following command will create a pdf file with the quality summary of R1 sequences of each  samples 
p=plotQualityProfile(fnR1s[1:length(sample.names)])
pdf("preprocessing/preprocessing_final/objets/Qualityprofile_R1_16S.pdf", height=30, width=30)
print(p)
dev.off()

###2.2.2 Quality check of the R2 (reverse)----
#The following command will create a pdf file with the quality summary of R1 sequences of each  samples 
p2=plotQualityProfile(fnR2s[1:length(sample.names)]) 
pdf("preprocessing/preprocessing_final/objets/Qualityprofile_R2_16S.pdf", height=30, width=30)
print(p2)
dev.off()


##2.3 Primer removing------
#Using cutadapt software based on Python 

#Amorce sequences:
#Multiple forward primers  
FWD4 <- "CAACCTACGGGNGGCWGCAG"
FWD3 <-  "ACCCTACGGGNGGCWGCAG"
FWD2 <-   "TCCTACGGGNGGCWGCAG"
FWD1 <-    "CCTACGGGNGGCWGCAG"
FWD <- c(FWD1, FWD2,FWD3,FWD4)

#Reverse primer 
REV <- "CMGGGTATCTAATCCKGTT"

#Check if its the right primers sequences and if its the right orientation on the sequence:

#Function   
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

#Make an vector of every orientation of the primers
FWD.orients <- c(allOrients(FWD1), allOrients(FWD2), allOrients(FWD3),allOrients(FWD4))

names(FWD.orients) <- c("Forward FP1", "Complement FP1", "Reverse FP1", "RevComp FP1","Forward FP2", "Complement FP2", "Reverse FP2", "RevComp FP2", 
                        "Forward FP3", "Complement FP3", "Reverse FP3", "RevComp FP3", "Forward FP4", "Complement FP4", "Reverse FP4", "RevComp FP4")

REV.orients <- allOrients(REV)


#Ambiguous bases (Ns) interfer with the mapping of short sequences. Ns bases should be thrimmed as follow:

#Creation of an sub-folder "filtN", and place the filtered file 
path.filtN <- file.path(pathn1, "filtN_seq")
if(!dir.exists(path.filtN)) dir.create(path.filtN)
fnR1s.filtN <- file.path(path.filtN, basename(fnR1s)) # Put N-filterd files in filtN/ subdirectory
fnR2s.filtN <- file.path(path.filtN, basename(fnR2s))

#Filtration of ambiguous bases (Ns) 
filterAndTrim(fnR1s, fnR1s.filtN, fnR2s, fnR2s.filtN, maxN = 0, multithread = TRUE)  

#Function to  count the presence of the primer in every orientation
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#See the number of the primer hit in each orientation for 1 sample 
numsam <- 1 #sample number 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnR1s.filtN[[numsam]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnR2s.filtN[[numsam]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnR1s.filtN[[numsam]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnR2s.filtN[[numsam]]))
#It expect that almost the globality of the primer is in the right sens and some in reverse complement. 


#Cutadapt tool to remove primers

#Install it on your machine using python on a terminal prompt with this command line:
#python3 -m pip install --user --upgrade cutadapt

#inform the path to the cutadapt tool 
cutadapt <- "XXXXXX" # CHANGE ME to the cutadapt path on your machine
#Run Shell command to check if cutadapt is well installed and if its the good version 
system2(cutadapt, args = "--version") # Run shell commands from R

#Creation of an sub-folder "cutadapt_seq", and place the primers thrimmed file   
path.cut <- file.path(pathn1, "cutadapt_seq")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnR1s.cut <- file.path(path.cut, basename(fnR1s))
fnR2s.cut <- file.path(path.cut, basename(fnR2s))

#Inform the reverse complement of each primers
FWD.RC <- dada2:::rc(FWD) #forward 
REV.RC <- dada2:::rc(REV) #reverse 

#Trim the forward and the reverse-complement of the reverse primer off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim the reverse and the reverse-complement of the forward primer off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnR1s)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnR1s.cut[i], "-p", fnR2s.cut[i], # output files
                             fnR1s.filtN[i], fnR2s.filtN[i], # input files
                             "--json", "filename.cutadapt.json")) # log sheet with output
}

#Check the output filename.cutadapt.json

#Check if cutadapt has remove the primers in the sequences in 1 sample
numsam <- 1 #sample number 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnR1s.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnR2s.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnR1s.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnR2s.cut[[1]]))


##2.4 General filtration----

###2.4.1 Filtration----   
#Filtration of the sequences based on their quality

#Attribute sample names to the future filtered fastq.gz file
#Creation of an sub-folder "filtered", and place the filtered file
path.filtered <- file.path(pathn1, "filtered_seq")
if(!dir.exists(path.filtered)) dir.create(path.filtered)
filtR1s <- file.path(path.filtered,  paste0(sample.names, "_R1"))  # Forward
filtR2s <- file.path(path.filtered,  paste0(sample.names, "_R2"))   # Reverse

#Save every object in case 
saveRDS(file=paste0(pathobj,"filtR1s"), filtR1s)
saveRDS(file=paste0(pathobj,"filtR2s"), filtR2s)



##Use of the filterAndTrim() function
out <- filterAndTrim(fwd=fnR1s.cut,filt= filtR1s, rev=fnR2s.cut, filt.rev=filtR2s, maxN = 0, maxEE = c(5,5), matchIDs=TRUE,
                     truncQ = 2, truncLen = c(250,215),minLen = 150, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

#Save object
saveRDS(file=paste0(pathobj,"out"), out)


#Quality check after filtration

###2.4.2 Quality check of the R1 (forward)----
#The following command will create a pdf file with the quality summary of R1 sequences of each  samples 
p=plotQualityProfile(filtR1s[1:length(sample.names)])
pdf(paste0(pathobj,"Qualityprofile_filtered_R1_16S.pdf"), height=30, width=30)
print(p)
dev.off()

###2.4.3 Quality check of the R2 (reverse)----
#The following command will create a pdf file with the quality summary of R1 sequences of each  samples 
p2=plotQualityProfile(filtR1s[1:length(sample.names)]) 
pdf(paste0(pathobj,"Qualityprofile_filtered_R2_16S.pdf"), height=30, width=30)
print(p2)
dev.off()

###2.4.4 Dereplication----
#Dereplication of the same sequences to reduce the calculation time 
derepR1s <- derepFastq(filtR1s, verbose=TRUE)#forward
names(derepR1s) <- sample.names
derepR2s <- derepFastq(filtR2s, verbose=TRUE)#reverse
names(derepR2s) <- sample.names

#Save object
saveRDS(file=paste0(pathobj,"derepR1s"), derepR1s)
saveRDS(file=paste0(pathobj,"derepR2s"), derepR2s)

##2.5 Correction of sequencing errors----

###2.5.1 Learn the errors rate model----  
#Learn the errors rate 
errF <- learnErrors(filtR1s)#Forward
errR <- learnErrors(filtR2s)#Reverse 

#Save object
saveRDS(file=paste0(pathobj,"errF"), errF)
saveRDS(file=paste0(pathobj,"errR"), errR)


#Graphic representation
plotErrors(errF, nominalQ=TRUE)#Forward  
plotErrors(errR, nominalQ=TRUE)#Reverse 

#Save the graphic representation in pdf
#Forward 
pdf(paste0(pathobj,"errF.pdf"))
plotErrors(errF, nominalQ=TRUE)
dev.off()
#Reverse 
pdf(paste0(pathobj,"errR.pdf"))
plotErrors(errR, nominalQ=TRUE)
dev.off()

###2.5.2 Sequencing error correction---- 

#Forward
dadaR1s <- dada(derepR1s, err=errF, multithread=TRUE)

#Reverse
dadaR2s <- dada(derepR2s, err=errR, multithread=TRUE)

#Save object
saveRDS(file=paste0(pathobj,"dadaR1s"), dadaR1s)
saveRDS(file=paste0(pathobj,"dadaR2s"), dadaR2s)


##2.6 Merge the forward and reverse---- 
mergers <- mergePairs(dadaR1s, filtR1s, dadaR2s, filtR2s, verbose=TRUE)

#Save object
saveRDS(file=paste0(pathobj,"mergers"), mergers)

#Sequence table construction   
seqtable <- makeSequenceTable(mergers)

#Save object
saveRDS(file=paste0(pathobj,"seqtable"), seqtable)

#Check the distribution of the sequence length
pdf(paste0(pathobj,"seqlengthdistribution.pdf"))
plot(table(nchar(getSequences(seqtable))))
plot(table(dim(seqtable)))
dev.off()

##2.7 Bimera remove----
seqtab.nochim <- removeBimeraDenovo(seqtable, method="consensus", multithread=TRUE, verbose=TRUE)

#Save object
saveRDS(file=paste0(pathobj,"seqtab.nochim"), seqtab.nochim)

#Check the distribution of the sequence length with no chim    
pdf(paste0(pathobj, "seqlengthdistribution_nonchim.pdf"))
plot(table(nchar(getSequences(seqtab.nochim))))
plot(table(dim(seqtab.nochim)))
dev.off()

##2.8 Summary of preprocessing----

#Construction of the table  
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaR1s, getN), sapply(dadaR2s, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

#Save the summary table
write.table(track,paste0(pathobj,"summary_preprocessing_16S.txt"),sep=";",quote=FALSE)

#3 Taxonomic assignment-----

#Assign the taxonomy from kingdom to the genus level
taxa1 <- assignTaxonomy(seqtab.nochim, "raw_data/silva_nr99_v138.1_train_set.fa.gz", 
                        taxLevels= c("Kingdom", "Phylum","Class", "Order", "Family","Genus", "Species"), multithread=TRUE, minBoot=50)

#Assign the taxonomy to the species level
taxa1 <- addSpecies(taxa1, "raw_data/silva_species_assignment_v138.fa.gz",allowMultiple=FALSE) 

#Save object
saveRDS(file=paste0(pathobj,"taxa1"), taxa1)

#4 Removal of contaminations----
#Removal of contaminants with the microDecon package 

#Transform the seqtab.nochim to the right format for the microDecon package
#Table for microDecon should have 
# ASV in row and the first column name is "ASV_ID"   
# Blank sample is in the first column
# Other samples are in the following columns  
# Last column contains taxonomy: "K_Bacteria; P_Phylum; C_Class; O_Order; F_Family; G_Genus; S_Species"
#Here is an example:
# example <- cbind.data.frame(c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6"),
#                                                          c(0,200,1000,50,0,25),
#                                                          c(0,220,800,30,0,10),
#                                                          c(0,180,1300,70,0,30),
#                                                          c(60,660,1440,70,2400,30),
#                                                          c(64,520,1000,48,1900,20),
#                                                          c(40,480,700,35,2100,15),
#                                                          c("K_Bacteria; P_Phylum; C_Class; O_Order; F_Family; G_Genus; S_Species"))
# colnames(example) <- c("OTU_ID","Blank1","Blank2","Blank3","Pop1_Sample1","Pop1_Sample2","Pop2_Sample3","Taxa")

#Creation of the desired taxonomy column
taxa2 <- taxa1
for (i in 2:6){
  taxa2[,i] <- paste(taxa2[,i-1] , taxa2[,i], sep = ";")
}

#Save object
saveRDS(file=paste0(pathobj,"taxa2"), taxa2)

#Creation of the data frame 
df_decon <- t(seqtab.nochim)%>%as.data.frame()
mbind<- function(...){
  Reduce(function(x,y){cbind(x,y[match(row.names(x),row.names(y)),])}, list(...) )
} 
df_decon <- mbind(x=df_decon, y=taxa2)
df_decon <- df_decon[, c(1:93, 99)] 
#ASV ID column
df_decon <- cbind(rownames(df_decon), df_decon)
#Name the column
colnames(df_decon) <- c("OTU_ID", sample.names, "taxonomy")

#Sort the columns to have the blanks first and then samples and then taxonomy 
df_decon_1 <- df_decon[,c("OTU_ID",negatif, algae_EDM, algae_ENR,filter_EDM, filter_ENR, "taxonomy")] 

#Aplication of the decon package to remove contamiants detected in blank samples 
decontaminated <- decon(df_decon_1, numb.blanks = 5, numb.ind = c(27,27,17,17), taxa=T, run=2)

#Check which OTU are removed 
decontaminated$OTUs.removed
#Check which reads are removed 
decontaminated$reads.removed


#Recovery of the ASV table with the removed contaminants 
decontaminated_1<-as.data.frame(select(decontaminated$decon.table, -c("Mean.blank", "taxonomy")))

#Modification of the format of the generated decon ASV table to look like the seqtab.nochim
seqtab.microdecon <- decontaminated_1

#Save names (ASVs sequences) into a vector 
save.rows <- seqtab.microdecon[,"OTU_ID"] 

#Replace the introduced NAs by 0 values
seqtab.microdecon[is.na(seqtab.microdecon)] <- 0 
seqtab.microdecon<- apply(seqtab.microdecon, 2, as.integer)

#Rename the ASVs with their sequences 
rownames(seqtab.microdecon) <- save.rows
#Get rid of the first column  
seqtab.microdecon <- t(seqtab.microdecon[,c(-1)] )


#Get off the removed ASV from the initial taxonomy table  
OTU_rmv_decon <- as.vector(decontaminated$OTUs.removed$OTU_ID)
taxa3 <- taxa1[ -match(OTU_rmv_decon, table= rownames(taxa1)), ] 


#Check if the dimension correspond with tuhe number of remove ASVs
dim(df_decon)
dim(decontaminated$decon.table)
dim(seqtab.microdecon)

#Save object
saveRDS(file=paste0(pathobj,"df_decon"), df_decon)
saveRDS(file=paste0(pathobj,"df_decon_1"), df_decon_1)
saveRDS(file=paste0(pathobj,"decontaminated"), decontaminated)
saveRDS(file=paste0(pathobj,"seqtab.microdecon"), seqtab.microdecon)
saveRDS(file=paste0(pathobj,"taxa3"), taxa3)


#5 Phylogenetic tree----

#Get the ASV sequences
seqs <-getSequences(seqtab.microdecon)

#Alignment by DECIPHER package
aln <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

#See the alignment
BrowseSeqs(aln, highlight=0)

#Conversion among Sequence Formats 
phang.align <- phyDat(as(aln, "matrix"), type="DNA")

#Change label (leaves) of the tree
dm <- dist.ml(phang.align)

#Distance Matrix
treeNJ <- NJ(dm)

#Root the phylogetenic tree
treeNJ <- phangorn::midpoint(tree = treeNJ)

#Change the name of the label to built the phyloseq object
taxa_names(treeNJ) <- taxa_names(otu_table(seqtab.microdecon, taxa_are_rows = FALSE))
#Save object
saveRDS(file=paste0(pathobj,"treeNJ"), treeNJ)


#6 Phyloseq object----
#Build the phyloseq object  

#Integrate the environnemental parameters 
MAPFILE <-import_qiime_sample_data("Sam_data.txt")


#Construction of the phyloseq object
Final_16S <- phyloseq(otu_table(seqtab.microdecon, taxa_are_rows = FALSE), sample_data(MAPFILE) ,tax_table(as.matrix(taxa3)), treeNJ)


#ADD ASV sequences to your phyloseq object
#add the class sequence (refseq) to the object final1
dna <- DNAStringSet(taxa_names(Final_16S))
names(dna) <- taxa_names(Final_16S)
Final1_16S <- merge_phyloseq(Final_16S, dna)

#Remove the Mitochondria and Chloroplast sequences
Final1_16S <- subset_taxa( Final1_16S , Family != "Mitochondria")
Final1_16S <- subset_taxa( Final1_16S , Order != "Chloroplast")

#See the proportions of each
sum(otu_table(subset_taxa( Final_16S , Family == "Mitochondria")))/sum(otu_table(Final_16S))*100
sum(otu_table(subset_taxa( Final_16S , Order == "Chloroplast")))/sum(otu_table(Final_16S))*100


#Remove singletons
#keep only sequences present at least 2 times and seen in 2 samples 
Final1_16S <- phyloseq::filter_taxa(Final1_16S, function(x) sum(x > 2) > round(((2/88)*length(x))), TRUE)

#Modification of the date in sam_data
Final1_16S@sam_data$date <- ymd(Final1_16S@sam_data$date)

#Change the ASV names from a sequence to an ID, ex : ASV1
taxa_names(Final1_16S) <- paste0("ASV", seq(ntaxa(Final1_16S)))


#Check the final object 
Final1_16S@otu_table
Final1_16S@tax_table
Final1_16S@sam_data
Final1_16S@phy_tree
Final1_16S@refseq


#Save object
saveRDS(file=paste0(pathobj,"Final_16S"), Final_16S)
saveRDS(file=paste0(pathobj,"Final1_16S"), Final1_16S)




