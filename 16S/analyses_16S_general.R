#Script analyses diversité microbienne 16S 

#Loadings----
library(phyloseq)
library(dplyr)
library(ymd)

Final1_16S <- readRDS(file = "F:/16S_HOLOGREEN/preprocessing/preprocessing_final/objets/Final1_16S")
physeq_16S <- Final1_16S



#Traiment of the phyloseq object----


##Modification of the mapfile----

sample_data(physeq_16S) <- sample_data(import_qiime_sample_data("Sam_data.txt"))
#Modification of the date in sam_data
physeq_16S@sam_data$date <- ymd(physeq_16S@sam_data$date)


#retrait des archée
physeq_16S <- subset_taxa(physeq_16S, Kingdom !="Archaea")

#transformer en pourcentage  
physeq_16S_p <- transform_sample_counts(physeq_16S, function(x) x/sum(x)*100)


#Normalisation---- 
total_rm_sing= median(sample_sums(physeq_16S))
standf_rm_sing = function(x, t=total_rm_sing) round(t * (x / sum(x)))
physeq2_16S = transform_sample_counts(physeq_16S, standf_rm_sing)

#Le nombre de lectures utilisées pour la normalisation est de `r total`.
total_rm_sing


# https://vaulot.github.io/tutorials/Phyloseq_tutorial.html





##II.10.Subsampling général----


##Subsampling ech-- 

# Pour faire un subset des dates non souhaités / par enrichissements  
# physeq_16S_nit <- subset_samples(physeq_16S, date != "210601" & date != "210616" & date !="210622" & date !="210629" & date !="210706" & date !="210713" )
physeq_16S_nit <- subset_samples(physeq_16S, date < "2021-06-01" )

# physeq2_16S_nit <- subset_samples(physeq2_16S, date != "210601" & date != "210616" & date !="210622" & date !="210629" & date !="210706" & date !="210713" )
physeq2_16S_nit <- subset_samples(physeq2_16S, date < "2021-06-01" )


#transformer en % 
physeq2_16S_p <- transform_sample_counts(physeq2_16S, function(x) x/sum(x)*100) 
physeq_16S_nit_p  <- transform_sample_counts(physeq_16S_nit, function(x) x/sum(x)*100)


physeq2_16S_nit_p <- transform_sample_counts(physeq2_16S_nit, function(x) x/sum(x)*100)
physeq2_16S_nit_biofilm_p <- transform_sample_counts(physeq2_16S_nit_biofilm, function(x) x/sum(x)*100)
physeq2_16S_nit_filter_p <- transform_sample_counts(physeq2_16S_nit_filter, function(x) x/sum(x)*100)


