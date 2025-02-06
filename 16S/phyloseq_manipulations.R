#Becareful obsolete




# #########Phyloseq object manipulations
# 
# Final1_16S <- readRDS(file="preprocessing/preprocessing_final/objets/Final1_16S")
# 
# 
# 
# ##Retrait singletons----
# 
# #1. Only 0.2-3, >3, centrifugation
# #retrait  des samples qui sont présent dans aumoin deux échantillons aumoin une fois 
# physeq_16S = phyloseq::filter_taxa(Final1_16S, function(x) sum(x > 1) > round(((2/88)*length(x))), TRUE)
# 
# 
# #transformer en pourcentage  
# physeq_16S_p <- transform_sample_counts(physeq_16S, function(x) x/sum(x)*100)
# 
# ##II.10.Normalisation---- 
# total_rm_sing= median(sample_sums(physeq_16S))
# standf_rm_sing = function(x, t=total_rm_sing) round(t * (x / sum(x)))
# physeq2_16S = transform_sample_counts(physeq_16S, standf_rm_sing)
# 
# #Le nombre de lectures utilisées pour la normalisation est de `r total`.
# total_rm_sing
# 
# 
# # https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
# 
# 
# 
# 
# 
# ##II.10.Subsampling général----
# 
# 
# ##Subsampling ech-- 
# 
# # Pour faire un subset des dates non souhaités / par enrichissements  
# # physeq_16S_nit <- subset_samples(physeq_16S, date != "210601" & date != "210616" & date !="210622" & date !="210629" & date !="210706" & date !="210713" )
# physeq_16S_nit <- subset_samples(physeq_16S, date < "2021-06-01" )
# save(physeq_16S_nit_p,file= "physeq_16S_nit_p.RData")
# physeq_16S_nitur <- subset_samples(physeq_16S, date !="2021-07-13" )
# 
# # physeq2_16S_nit <- subset_samples(physeq2_16S, date != "210601" & date != "210616" & date !="210622" & date !="210629" & date !="210706" & date !="210713" )
# physeq2_16S_nit <- subset_samples(physeq2_16S, date < "2021-06-01" )
# 
# physeq2_16S_nitur <- subset_samples(physeq2_16S, date !="2021-07-13" )
# 
# # Subset des Biofilm/filter
# physeq2_16S_nit_biofilm <- subset_samples(physeq2_16S_nit, type=="biofilm")
# physeq2_16S_nit_filter <- subset_samples(physeq2_16S_nit, type=="filtre")
# 
# #transformer en % 
# physeq2_16S_p <- transform_sample_counts(physeq2_16S, function(x) x/sum(x)*100) 
# physeq_16S_nit_p <- transform_sample_counts(physeq_16S_nit, function(x) x/sum(x)*100)
# physeq_16S_nitur_p <- transform_sample_counts(physeq_16S_nitur, function(x) x/sum(x)*100)
# physeq2_16S_nit_p <- transform_sample_counts(physeq2_16S_nit, function(x) x/sum(x)*100)
# physeq2_16S_nitur_p <- transform_sample_counts(physeq2_16S_nitur, function(x) x/sum(x)*100)
# physeq2_16S_nit_biofilm_p <- transform_sample_counts(physeq2_16S_nit_biofilm, function(x) x/sum(x)*100)
# physeq2_16S_nit_filter_p <- transform_sample_counts(physeq2_16S_nit_filter, function(x) x/sum(x)*100)
# 
# 
# #merge d'objet phyloseq selon une var environnementales 
# merged_physeq2_16S_p <- merge_samples(physeq2_16S_p, "type_enrich")
# 
# #Quand le nombre d'ech n'est pas égal entre catégorie , ifaut rajouter une étape 
# merged_physeq2_16S_nit_p <- transform_sample_counts(merge_samples(physeq2_16S_nit_p, "type_enrich"), function(x) x/sum(x)*100)
# merged_physeq2_16S_nit_p@sam_data$type_enrich <- c("biofilm_SW", "biofilm_ENR", "filtre_SW", "filtre_ENR")
# merged_physeq2_16S_nit_p@sam_data$type <- c("biofilm", "biofilm", "filtre", "filtre")
# 
# merged_physeq2_16S_nitur_p <- transform_sample_counts(merge_samples(physeq2_16S_nitur_p, "type_enrich"), function(x) x/sum(x)*100)
# 
# #Merge de l'objet phyloseq  
# pyseq_merged_type_enrich <- merge_samples(Final2_16S_Ts_P,"type_enrich")
# otu_table(pyseq_merged_type_enrich) <- otu_table(as.data.frame(rbind(pyseq_merged_type_enrich@otu_table["biofilm_EDM",]/nsamples(subset_samples(Final2_16S_Ts_P, type_enrich =="biofilm_EDM")),
#                                                                      pyseq_merged_type_enrich@otu_table["biofilm_ENR",]/nsamples(subset_samples(Final2_16S_Ts_P, type_enrich =="biofilm_ENR")),
#                                                                      pyseq_merged_type_enrich@otu_table["filtre_EDM",]/nsamples(subset_samples(Final2_16S_Ts_P, type_enrich =="filtre_EDM")),
#                                                                      pyseq_merged_type_enrich@otu_table["filtre_ENR",]/nsamples(subset_samples(Final2_16S_Ts_P, type_enrich =="filtre_ENR")))), 
#                                                  taxa_are_rows = F)
# 
# pyseq_merged_type_enrich_nit <- subset_samples(pyseq_merged_type_enrich, date != "210616" & date !="210622" & date !="210629" & date !="210706" & date !="210713" )
# pyseq_merged_type_enrich_nitur <- subset_samples(pyseq_merged_type_enrich, date !="210713" )
# 
# 
