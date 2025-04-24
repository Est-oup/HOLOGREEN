#######TRAITEMENT OF 18S rRNA SEQUENCE FROM HOLOGREEN PROJECT######

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

#setup the directory of the location of the raw sequences
path <- "raw_data/sequences_18S/raw_seq/"
pathn1<- "raw_data/sequences_18S"
pathobj <- "preprocessing/preprocessing_final/objets/"



# #ou trouver les seq brutes
# fnR1s <- sort(list.files(paste0(path,"/seq_brute"), pattern="_R1.fastq.gz", full.names = TRUE))   # Forward
# fnR2s <- sort(list.files(paste0(path,"/seq_brute"), pattern="_R2.fastq.gz", full.names = TRUE))   # Reverse
# 
# 
# 
# # # Ici on stocke les noms des échantillons dans la variable "sample.names"
# sample.names <- sapply(strsplit(list.files(path, pattern="_R1"),"---IDT_i5_"), "[",2)
# sample.names <- sapply(strsplit(sample.names,"TAReuk"), "[",1)
# sample.names <- sapply(strsplit(sample.names,"[.]"), "[",2)
# 
# sample.namesR1 <- paste0(sample.names, "_R1")
# sample.namesR2 <- paste0(sample.names, "_R2")
# #on donne le nom des ech 
# names(fnR1s) <- sample.namesR1
# names(fnR2s) <- sample.namesR2
# 
# #Faire des subset algues filtres négatifs
# dfsampnam <- as.data.frame(sample.names)
# algae <- subset( dfsampnam[grep(pattern = "_A", x= sample.names),], dfsampnam[grep(pattern = "_A", x= sample.names),] != "NEGATIF_A" )
# algae_ENR <- na.omit(algae[grep(pattern= "ENR", x=algae)])
# algae_EDM <- na.omit(algae[grep(pattern= "EDM", x=algae)]) 
# saveRDS(algae, file="algae")
# 
# filter<- subset( dfsampnam[grep(pattern = "_F", x= sample.names),], dfsampnam[grep(pattern = "_F", x= sample.names),] != "NEGATIF_F")
# filter_ENR <- na.omit(filter[grep(pattern= "ENR", x=filter)])
# filter_EDM <- na.omit(filter[grep(pattern= "EDM", x=filter)]) 
# saveRDS(filter, file="filter")
# 
# 
# negatif <- c(dfsampnam[grep(pattern = "NEG", x= sample.names),], dfsampnam[grep(pattern = "Blank", x= sample.names),])
# negatif_A<- subset(negatif, negatif!="NEGATIF_F")
# negatif_F<- subset(negatif, negatif!="NEGATIF_F")
# saveRDS(negatif, file="negatif")
# 
# 
# names_type_enrich_sample <- list(algae_EDM, algae_ENR, filter_EDM, filter_ENR)
# 
# #Subset des dates d'échantillonage des algues 
# date_algae <- c("210309", "210330" ,"210402" ,"210414" ,"210427" ,"210511" ,"210525" ,"210622" ,"210706")
# 
# ##II.2.Quality check---- 
# 
# #Vérification de la qualité et de la taille desséquences avant analyse 
# #
# #Profil qualité des séquences brutes
# p=plotQualityProfile(fnR1s[1:93])  # Forward
# pdf("Qualityprofile_fnR1s_18S.pdf", height=30, width=30)
# print(p)
# dev.off()
# p2=plotQualityProfile(fnR2s[1:93])  # Reverse
# pdf("Qualityprofile_fnR2s_18S.pdf", height=30, width=30)
# print(p2)
# dev.off()
# 
# 
# ##II.3.Cutadapt-retrait amorces---------
# #Prise en compte de la composition en nucléotides - cutadapt
# #+ d'info à https://benjjneb.github.io/dada2/ITS_workflow.html
# #####Identifier les amorces
# FWD <- "CCAGCASCYGCGGTAATTCC"
# REV <- "ACTTTCGTTCTTGATYRA"

# #Pour être sûr d'avoir les bonnes amorces et l'orientation correcte des amorces sur les séquences,
# #nous allons vérifier la présence et l'orientation des amorces dans les données.
# 
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = Biostrings::complement(dna), 
#                Reverse = reverse(dna), RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
#
# FWD.orients <-allOrients(FWD)
# REV.orients <-allOrients(REV)
# 
# #La présence de bases ambigues (Ns) dans les séquences rendent difficile le mapping des courtes séquences d'amorces
# #Maintenant, nous allons pré-filtrer les séquences pour enlever ces Ns mais sans faire d'autre forme de filtration.
# fnR1s.filtN <- file.path(path, "filtN", basename(fnR1s)) # Put N-filterd files in filtN/ subdirectory
# fnR2s.filtN <- file.path(path, "filtN", basename(fnR2s))
# filterAndTrim(fnR1s, fnR1s.filtN, fnR2s, fnR2s.filtN, maxN = 0, multithread = TRUE)
# 
# #compter le nombre de fois les amorces apparaissent dans les séquences forward et reverse, 
# 
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnR1s.filtN[[1]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnR2s.filtN[[1]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnR1s.filtN[[1]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnR2s.filtN[[1]]))
# 
# #Ces amorces peuvent être enlevées à l'aide d'un outil spécialisé pour ça: cutadapt 
# #Instructions pour le télécharger, l'installer et l'utiliser ici: 
# #http://cutadapt.readthedocs.io/en/stable/index.html
# #On a terminal
# # python3 -m pip install --user --upgrade cutadapt
# #Pour mettre à jour python: python3 -m pip install --upgrade pip --user

# # cutadapt <- "path to cutadapt"
# system2(cutadapt, args = "--version") # Run shell commands from R
# path.cut <- file.path(path, "cutadapt")
# if(!dir.exists(path.cut)) dir.create(path.cut)
# fnR1s.cut <- file.path(path.cut, basename(fnR1s))
# fnR2s.cut <- file.path(path.cut, basename(fnR2s))
# 
# FWD.RC <- dada2:::rc(FWD)
# REV.RC <- dada2:::rc(REV)
# # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
# R1.flags <- paste("-g", FWD, "-a", REV.RC)
# # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
# R2.flags <- paste("-G", REV, "-A", FWD.RC)
# 
# names(FWD)
# names(FWD.RC)
# 
# 
# # Run Cutadapt 
# 
# for(i in seq_along(fnR1s)) {
#   system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
#                              "-o", fnR1s.cut[i], "-p", fnR2s.cut[i], # output files
#                              fnR1s.filtN[i], fnR2s.filtN[i], # input files
#                              "> cutadapt.log", # log sheet with output
#                              "-j", "22")) #multithreading 
# }
# 
# #Check sample by sample
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnR1s.cut[[5]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnR2s.cut[[5]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnR1s.cut[[5]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnR2s.cut[[5]]))
# 

# #Séquences avec les primers retirés par cutadapt 
# 
# #ou trouver les sequences filtrés N
# fnR1s.filtN <- file.path(path, "filtN", basename(fnR1s))
# fnR2s.filtN <- file.path(path, "filtN", basename(fnR2s))
# 
# #Profil qualité des séquences filtrés N
# # p=plotQualityProfile(fnR1s.filtN[1:93])  # Forward
# # pdf("qualityprofile/Qualityprofile_fnR1s.filtN_16S.pdf", height=30, width=30)
# # print(p)
# # dev.off()
# # p2=plotQualityProfile(fnR2s.filtN[1:93])  # Reverse
# # pdf("qualityprofile/Qualityprofile_fnR2s.filtN_16S.pdf", height=30, width=30)
# # print(p2)
# # dev.off()
# 
# ####ou trouver les seq coupées par le package cutadapt 
# path.cut <- file.path(path, "cutadapt")
# if(!dir.exists(path.cut)) dir.create(path.cut)
# fnR1s.cut <- file.path(path.cut, basename(fnR1s))
# fnR2s.cut <- file.path(path.cut, basename(fnR2s))
# 
# #subset des ech qui ont trop peu de séquences
# fnR1s.cut.sub <- fnR1s.cut[!fnR1s.cut %in% c("~/HOLOGREEN/compress/cutadapt/MI.M05812_0313.001.FLD_ill_189_i7---IDT_i5_2.210706_ENR_B_A_16S_R1.fastq.gz",
#                                    "~/HOLOGREEN/compress/cutadapt/MI.M05812_0313.001.FLD_ill_188_i7---IDT_i5_2.210706_ENR_A_A_16S_R1.fastq.gz")] 
# 
# fnR2s.cut.sub <- fnR2s.cut[!fnR2s.cut %in% c("~/HOLOGREEN/compress/cutadapt/MI.M05812_0313.001.FLD_ill_188_i7---IDT_i5_2.210706_ENR_A_A_16S_R2.fastq.gz",
#                                    "~/HOLOGREEN/compress/cutadapt/MI.M05812_0313.001.FLD_ill_189_i7---IDT_i5_2.210706_ENR_B_A_16S_R2.fastq.gz" )] 
# 
# 
# 
# # #Profil qualité des séquences cut amorces
# # p=plotQualityProfile(fnR1s.cut[1:93])  # Forward
# # pdf("qualityprofile/Qualityprofile_fnR1s.cut_16S.pdf", height=30, width=30)
# # print(p)
# # dev.off()
# # p2=plotQualityProfile(fnR2s.filtN[1:93])  # Reverse
# # pdf("qualityprofile/Qualityprofile_fnR2s.cut_16S.pdf", height=30, width=30)
# # print(p2)
# # dev.off()
# 
# # Attribuer les noms de fichier aux fichiers fastq.gz filtrés
# # Placer les fichiers filtrés dans le sous-dossier "filtered" (création de ce sous-dossier)
# filtR1s <- file.path(path, "filtered", sample.namesR1)  # Forward
# filtR2s <- file.path(path, "filtered", sample.namesR2)   # Reverse
# 
# # On redonne aux fichiers forward R1/reverse R2 leurs noms d'échantillons
# sample.names <- sapply(strsplit(list.files(paste0(path,"/seq_brute"), pattern="_R1"),"---IDT_i5_"), "[",2)
# sample.names <- sapply(strsplit(sample.names,"TAReuk"), "[",1)
# sample.names <- sapply(strsplit(sample.names,"[.]"), "[",2)
# names(filtR1s) <- sample.names
# names(filtR2s) <- sample.names
# 
# 
# # Attribuer les noms de fichier aux fichiers fastq.gz filtrés sur les fichiers sub
# # Placer les fichiers filtrés dans le sous-dossier "filtered" (création de ce sous-dossier)
# filtR1s.sub <- file.path(path, "filtered", sample.namesR1[-c(89,90)] )  # Forward
# filtR2s.sub <- file.path(path, "filtered", sample.namesR2[-c(89,90)])   # Reverse
# 
# # On redonne aux fichiers forward R1/reverse R2 leurs noms d'échantillons
# names(filtR1s.sub) <- sample.names[-c(89,90)]
# names(filtR2s.sub) <- sample.names[-c(89,90)]
# 
# 
# sample.namesR1 <- readRDS(file = "sample.namesR1")
# sample.namesR2 <- readRDS(file = "sample.namesR2")
# 
# algae <-readRDS(file="algae")
# filter <- readRDS(file="filter")
# negatif <- readRDS(file="negatif")


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

seqtab.nochim <- readRDS(file = "G:/18S_HOLOGEEN/preprocess_18S/data_preprocess/test88/var/seqtab.nochim88")

##II.5.Taxonomic Assignment-----
taxa1 <- assignTaxonomy(seqtab.nochim, "pr2_version_5.0.0_SSU_dada2.fasta.gz", 
                        taxLevels= c("Kingdom", "Phylum","Class", "Order", "Family","Genus", "Species"), multithread=TRUE, minBoot=50)

#pour mettre l'entièreté des nom des taxons sur chaque case 
taxa2 <- taxa1
for (i in 2:6){
  taxa2[,i] <- paste(taxa2[,i-1] , taxa2[,i], sep = ";")
}


##II.6.Retrait contaminations -------
#Dans notre cas aucune séquences des blanc n'est passé après le filtre
#on retire les blanc qui ont pas de séquences 
rowSums(seqtab.nochim[c("NEGATIF_A","NEGATIF_PCR_1"),])

seqtab.nochim <- seqtab.nochim[rownames(seqtab.nochim) != c("NEGATIF_A","NEGATIF_PCR_1"), ] 


#Contextual data
#tableau de données envrionnementales 
MAPFILE <-import_qiime_sample_data("Sam_data.txt")


#contruction de l'objet phyloseq
Final_18S <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(MAPFILE), tax_table(as.matrix(taxa1))) 


#ADD ASV sequences to your phyloseq object
#add the class sequence (refseq) to the object final1

dna <- Biostrings::DNAStringSet(taxa_names(Final_18S))
names(dna) <- taxa_names(Final_18S)


Final1_18S <- merge_phyloseq(Final_18S, dna)

#Modifictaion de la date  
Final1_18S@sam_data$date <- ymd(Final1_18S@sam_data$date)

#Probleme --> the name of ASV is the sequence --> change to ASV1 ...
taxa_names(Final1_18S) <- paste0("ASV", seq(ntaxa(Final1_18S)))



##II.9.Retrait Ulve----
colnames(Final1_18S@tax_table) <- c("Kingdom", "Phylum",  "Class" ,  "Order",   "Family",  "Genus",   "Species" ,"NA.", "NA..")
Final1_18S_noulva <- subset_taxa( Final1_18S , NA. != "Ulva")

#Quantifier le % de read d'Ulves dans les échantillons
round(mean(sample_sums(subset_taxa(subset_samples(Final1_18S, type %in% "biofilm"), NA. == "Ulva"))/ sample_sums(subset_samples(Final1_18S, type %in% "biofilm"))) *100,1)
round(mean(sample_sums(subset_taxa(subset_samples(Final1_18S, type %in% "filtre"), NA. == "Ulva"))/ sample_sums(subset_samples(Final1_18S, type %in% "filtre"))) *100,1)

#nombre de reands non ulva 
mean(sample_sums(subset_samples(Final1_18S_noulva, type %in% "biofilm")))
mean(sample_sums(subset_samples(Final1_18S_noulva, type %in% "filtre"))) 

##II.9.Retrait singletons----

#1. Only 0.2-3, >3, centrifugation
#retrait  des samples qui sont présent dans au mn deux échantillons aumoin une fois 
physeq_18S = phyloseq::filter_taxa(Final1_18S_noulva, function(x) sum(x > 1) > round(((2/88)*length(x))), TRUE)

physeq_18S <- prune_samples(!sample_sums(physeq_18S) %in% c("210427_EDM_C_A","210427_EDM_F","210511_ENR_C_A","210525_ENR_C_A","210601_ENR_F","210706_ENR_A_A"), physeq_18S)
physeq_18S <- prune_samples(sample_sums(physeq_18S) > 0, physeq_18S)

#transformer en pourcentage  
physeq_18S_p <- transform_sample_counts(physeq_18S, function(x) x/sum(x)*100)


saveRDS(physeq_18S, "physeq_18S.rds")