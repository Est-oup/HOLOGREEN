#venn Diagram

#loadings 

library("VennDiagram")
library("venn")
library(ggVennDiagram)
library(eulerr)
library(microbiome)

library(microbiomeutilities)



#valeurs pour diagramme de venn
valeurs.venn <- function(physeq){
  physeq@otu_table <- physeq@otu_table[,colSums(physeq@otu_table) != 0]
  
  biofilm_SW <- subset_samples(physeq_16S_nit_p, type_enrich=="biofilm_SW")
  biofilm_SW@otu_table <- biofilm_SW@otu_table[,colSums(biofilm_SW@otu_table) != 0]
  biofilm_ENR <- subset_samples(physeq_16S_nit_p, type_enrich =="biofilm_ENR")
  biofilm_ENR@otu_table <- biofilm_ENR@otu_table[,colSums(biofilm_ENR@otu_table) != 0]
  biofilm <- subset_samples(physeq_16S_nit_p, type=="biofilm")
  biofilm@otu_table <- biofilm@otu_table[,colSums(biofilm@otu_table) != 0]
  
  water_SW <- subset_samples(physeq_16S_nit_p, type_enrich=="water_SW")
  water_SW@otu_table <- water_SW@otu_table[,colSums(water_SW@otu_table) != 0]
  water_ENR <- subset_samples(physeq_16S_nit_p, type_enrich =="water_ENR")
  water_ENR@otu_table <- water_ENR@otu_table[,colSums(water_ENR@otu_table) != 0]
  water <- subset_samples(physeq_16S_nit_p, type=="water")
  water@otu_table <- water@otu_table[,colSums(water@otu_table) != 0]
  
  
  SW <- subset_samples(physeq_16S_nit_p, enrich=="SW")
  SW@otu_table <- SW@otu_table[,colSums(SW@otu_table) != 0]
  ENR <- subset_samples(physeq_16S_nit_p, enrich=="ENR")
  ENR@otu_table <- ENR@otu_table[,colSums(ENR@otu_table) != 0]
  
  ASV_biofilm_SW <- taxa_names(biofilm_SW)
  ASV_biofilm_ENR <- taxa_names(biofilm_ENR)
  ASV_water_SW <- taxa_names(water_SW)
  ASV_water_ENR <- taxa_names(water_ENR)
  
  nb_ASV_tot_cond_SW <- length(taxa_names(SW))
  nb_ASV_tot_cond_ENR <- length(taxa_names(ENR))
  nb_ASV_tot_type_biof <- length(taxa_names(biofilm))
  nb_ASV_tot_type_water <- length(taxa_names(water))
  

  #calculs 
  
  #biofilm 
  unique_biof_SW <- setdiff(taxa_names(biofilm_SW),taxa_names(biofilm_ENR))
  unique_biof_ENR <- setdiff(taxa_names(biofilm_ENR),taxa_names(biofilm_SW))
  sharedbiof <- intersect(taxa_names(biofilm_SW),taxa_names(biofilm_ENR))
  

  #biofilm abondance 
  abund_unique_biof_SW <- round(mean(rowSums(phyloseq::prune_taxa(unique_biof_SW, phyloseq::subset_samples(physeq, type =="biofilm"))@otu_table)),1)
  abund_unique_biof_ENR <- round(mean(rowSums(phyloseq::prune_taxa(unique_biof_ENR, phyloseq::subset_samples(physeq, type =="biofilm"))@otu_table)),1)
  abund_shared_biof <- round(mean(rowSums(phyloseq::prune_taxa(sharedbiof, phyloseq::subset_samples(physeq, type =="biofilm"))@otu_table)),1)
  
  
  #water 
  unique_water_SW <- setdiff(taxa_names(water_SW),taxa_names(water_ENR))
  unique_water_ENR <- setdiff(taxa_names(water_ENR),taxa_names(water_SW))
  sharedwater <- intersect(taxa_names(water_SW),taxa_names(water_ENR))
  
  #water abondance 
  abund_unique_water_SW <- round(mean(rowSums(phyloseq::prune_taxa(unique_water_SW, phyloseq::subset_samples(physeq, type =="water"))@otu_table)),1)
  abund_unique_water_ENR <- round(mean(rowSums(phyloseq::prune_taxa(unique_water_ENR, phyloseq::subset_samples(physeq, type =="water"))@otu_table)),1)
  abund_shared_water <- round(mean(rowSums(phyloseq::prune_taxa(sharedwater, phyloseq::subset_samples(physeq, type =="water"))@otu_table)),1)
  
  
  
  
  #compartment 
  #SW  
  unique_biof_SW_comp <- setdiff(taxa_names(biofilm_SW),taxa_names(water_SW))
  unique_water_SW_comp <- setdiff(taxa_names(water_SW),taxa_names(biofilm_SW))
  sharedSW <- intersect(taxa_names(biofilm_SW),taxa_names(water_SW))
  
  #sw abondance 
  abund_unique_biof_SW_comp <- round(mean(rowSums(phyloseq::prune_taxa(unique_biof_SW_comp, phyloseq::subset_samples(physeq, enrich =="SW"))@otu_table)),1)
  abund_unique_water_SW_comp <- round(mean(rowSums(phyloseq::prune_taxa(unique_water_SW_comp, phyloseq::subset_samples(physeq, enrich =="SW"))@otu_table)),1)
  abund_sharedSW <- round(mean(rowSums(phyloseq::prune_taxa(sharedSW, phyloseq::subset_samples(physeq, enrich =="SW"))@otu_table)),1)
  
  #ENR
  unique_biof_ENR_comp <- setdiff(taxa_names(biofilm_ENR),taxa_names(water_ENR))
  unique_water_ENR_comp <- setdiff(taxa_names(water_ENR),taxa_names(biofilm_ENR))
  sharedENR <- intersect(taxa_names(biofilm_ENR),taxa_names(water_ENR))
  
  #enr abondance 
  abund_unique_biof_ENR_comp <- round(mean(rowSums(phyloseq::prune_taxa(unique_biof_ENR_comp, phyloseq::subset_samples(physeq, enrich =="ENR"))@otu_table)),1)
  abund_unique_water_ENR_comp <- round(mean(rowSums(phyloseq::prune_taxa(unique_water_ENR_comp, phyloseq::subset_samples(physeq, enrich =="ENR"))@otu_table)),1)
  abund_sharedENR <- round(mean(rowSums(phyloseq::prune_taxa(sharedENR, phyloseq::subset_samples(physeq, enrich =="ENR"))@otu_table)),1)
  
  
  
  
  # 
  print("figure article compartment")
  ##% d'asv biof
  print( paste0("% ASV unique biofilm_SW dans type biof: ", round(length(unique_biof_SW)/nb_ASV_tot_type_biof*100)),2)
  print( paste0("% ASV unique biofilm_ENR dans type biof: ", round(length(unique_biof_ENR)/nb_ASV_tot_type_biof*100)),2)
  print( paste0("% ASV shared dans type biof: ", round(length(sharedbiof)/nb_ASV_tot_type_biof*100)),2)
  ##% d'abondance biof
  print( paste0("% d'abondance ASV unique biofilm_SW dans type biof: ", abund_unique_biof_SW))
  print( paste0("% d'abondance ASV unique biofilm_ENR dans type biof: ", abund_unique_biof_ENR))
  print( paste0("% d'abondance ASV shared dans type biof: ", abund_shared_biof))
  
  #% asv water
  print( paste0("% ASV unique waterilm_SW dans type water: ", round(length(unique_water_SW)/nb_ASV_tot_type_water*100)),2)
  print( paste0("% ASV unique waterilm_ENR dans type water: ", round(length(unique_water_ENR)/nb_ASV_tot_type_water*100)),2)
  print( paste0("% ASV shared dans type water: ", round(length(sharedwater)/nb_ASV_tot_type_water*100)),2)
  ##% d'abondance water
  print( paste0("% d'abondance ASV unique waterilm_SW dans type water: ", abund_unique_water_SW))
  print( paste0("% d'abondance ASV unique waterilm_ENR dans type water: ", abund_unique_water_ENR))
  print( paste0("% d'abondance ASV shared dans type water: ", abund_shared_water))
  
         
  
  print("figure supp compartment")
  ##% d'asv SW
  print( paste0("% ASV unique biofilm_SW dans cond SW: ", round(length(unique_biof_SW_comp)/nb_ASV_tot_cond_SW*100)),2)
  print( paste0("% ASV unique water_SW dans cond SW: ", round(length(unique_water_SW_comp)/nb_ASV_tot_cond_SW*100)),2)
  print( paste0("% ASV shared dans cond SW: ", round(length(sharedSW)/nb_ASV_tot_cond_SW*100)),2)
  ##% d'abondance SW
  print( paste0("% d'abondance ASV unique biofilm_SW dans cond SW: ", abund_unique_biof_SW_comp))
  print( paste0("% d'abondance ASV unique water_SW dans cond SW: ", abund_unique_water_SW_comp))
  print( paste0("% d'abondance ASV unique water_SW dans cond SW: ", abund_sharedSW))
  
  ##% d'asv ENR
  print( paste0("% ASV unique biofilm_ENR dans cond ENR: ", round(length(unique_biof_ENR_comp)/nb_ASV_tot_cond_ENR*100)),2)
  print( paste0("% ASV unique water_ENR dans cond ENR: ", round(length(unique_water_ENR_comp)/nb_ASV_tot_cond_ENR*100)),2)
  print( paste0("% ASV shared dans cond ENR: ", round(length(sharedENR)/nb_ASV_tot_cond_ENR*100)),2)
  ##% d'abondance ENR
  print( paste0("% d'abondance ASV unique biofilm_ENR dans cond ENR: ", abund_unique_biof_ENR_comp))
  print( paste0("% d'abondance ASV unique water_ENR dans cond ENR: ", abund_unique_water_ENR_comp))
  print( paste0("% d'abondance ASV unique water_ENR dans cond ENR: ", abund_sharedENR))
  
  
}


valeurs.venn(physeq_16S_nit_p)


















##Diagramme toutes conditions confondues
pseq= subset_samples(physeq_16S_nit, enrich ==" SW")
table(meta(pseq)$type, useNA = "always")
pseq.rel <- microbiome::transform(pseq, "compositional")
list_core <- c() # an empty object to store information
for (n in pseq@sam_data$type){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, type == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g
                         detection = 0.00, # 0.001 in atleast 90% samples
                         prevalence = 0.0)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
venn(list_core,zcolor = c("#999999", "#E69F00", "#56B4E9", "red") , opacity = 0.3 ,plotsize = 15, box = FALSE, ilcs = 1.5, sncs = 0)
venn(list_core) 


ggVennDiagram(list_core, label_alpha = 0, category.names =c("Biofilm eau de mer", "Biofilm engrais","Filtre eau de mer", "Filtre engrais") )+ ggplot2::scale_fill_gradient(low="skyblue", high="pink")+
  scale_x_continuous(expand=expansion(mult=.2))

ggVennDiagram(list_core, label_alpha = 0, category.names =c("Biofilm eau de mer", "Filtre eau de mer") )+ ggplot2::scale_fill_gradient(low="skyblue", high="pink")+
  scale_x_continuous(expand=expansion(mult=.2))



##Diagramme toutes conditions confondues
pseq= subset_samples(physeq_16S_nit_p, type == "water")

list_core <- c() # an empty object to store information
for (n in pseq@sam_data$type_enrich){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq, type_enrich == n) # Choose sample from DiseaseState by n
  
  core_m <- microbiome::core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g
                         detection = 0.00, # 0.001 in atleast 90% samples
                         prevalence = 0.0)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
venn(list_core,zcolor = c("#999999", "#E69F00", "#56B4E9", "red") , opacity = 0.3 ,plotsize = 15, box = FALSE, ilcs = 1.5, sncs = 0)
venn(list_core) 


ggVennDiagram(list_core, label_alpha = 0, category.names =c("Biofilm eau de mer", "Biofilm engrais","Filtre eau de mer", "Filtre engrais") )+ ggplot2::scale_fill_gradient(low="skyblue", high="pink")+
  scale_x_continuous(expand=expansion(mult=.2))
