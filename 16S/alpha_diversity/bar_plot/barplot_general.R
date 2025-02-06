

library(dplyr)
library(phyloseq)
library(ggplot2)
library(tidyverse)



#Functions----
bar.plot.data.glom <- function(physeq, percent, env_var){
  
  #agglomeration de l'objet phyloseq 
  
  #merge d'objet phyloseq selon une var environnementales 
  merged_physeq <- merge_samples(physeq, env_var)
  
  #Quand le nombre d'ech n'est pas égal entre catégorie , ifaut rajouter une étape 
  merged_physeq <- transform_sample_counts(merge_samples(physeq, env_var), function(x) x/sum(x)*100)
  merged_physeq@sam_data$type_enrich <- c("biofilm_ENR","biofilm_SW","filter_ENR","filter_SW")
  merged_physeq@sam_data$type <- c("biofilm", "biofilm", "filter", "filter")
  
  
  
  
  
  #Version améliorée avec les 4 variables 
  glom <- tax_glom(merged_physeq, taxrank = 'Genus', NArm = FALSE, bad_empty = c(""))
  
  data_glom <- psmelt(glom)
  
  for (i in 27:32) {  # Boucle de la colonne Genus à la colonne Kingdom
    data_glom[, i] <- ifelse(is.na(data_glom[, i]), 
                             paste("unassigned", data_glom[, (i-1)], sep = " "), 
                             data_glom[, i])
  }
  
  
  data_glom$Genus <- as.character(data_glom$Genus)
  data_glom$Class <- as.character(data_glom$Class)
  data_glom <-data.frame(data_glom %>% dplyr::group_by(OTU,Sample,Abundance,type_enrich,Kingdom,Phylum,Class,Order,Family,Genus)%>% summarize(Abundance = sum(Abundance)))
  
  
  for (c in unique(data_glom$Class)) {
    tmp <- data_glom[data_glom$Class == c, ]
    other_genus <- tmp[, "Genus"]
    if (max(tmp$Abundance) < 4) {
      data_glom[rownames(data_glom) %in% rownames(unique(tmp)),][data_glom[rownames(data_glom) %in% rownames(unique(tmp)),]$Abundance < 4,]$Genus <- "Others <4%"
    }
    if (max(tmp$Abundance) >= 4){
      data_glom[rownames(data_glom) %in% rownames(unique(tmp)),][data_glom[rownames(data_glom) %in% rownames(unique(tmp)),]$Abundance < 4,]$Genus <- paste0("Others <4% ",c)
    }
  }
  
  
  
  
  # Créer une nouvelle colonne "title" en utilisant les conditions demandées
  data_glom$title <- ifelse(data_glom$Genus != "Others <4%", data_glom$Class, "Others")
  
  
  # Réordonner le tableau en fonction de l'ordre obtenu
  data_glom_ordered <- data_glom
  
  #calcul des moyennes d'abondances
  data_glom_ordered <- data_glom_ordered %>% 
    group_by(Class, Genus) %>% 
    summarize(Abundance = sum(Abundance)) %>% 
    ungroup()
  
  
  
  #réordonner suivant la classe
  
  #on extrait Other < 2% on fait la moyenne abund par colass on ordonne et on récupère les noms 
  
  class_order <- (subset(data_glom_ordered, !grepl("^Other", Genus)) %>%
                    group_by(Class) %>%
                    summarize(Abundance = sum(Abundance)) %>%
                    ungroup() %>% data.frame()%>%
                    arrange(desc(Abundance)))$Class
  
  data_glom_ordered <- data_glom_ordered %>%
    arrange(factor(Class, levels = class_order))
  
  #réordonner suivant le genre
  
  #pour que dans chaque classe les genres soient ordonnés par abondance avec other a la fin 
  order_data_glom <- NA
  for (i in class_order){
    data <- subset(data_glom_ordered, Class==i)#extraire les données pour chaque classe
    data1 <- subset(data, grepl("^Other", Genus))#extraire Other 
    data2 <- subset(data, !grepl("^Other", Genus))#extraire tout le reste
    #ON réarange par ordre d'abondance 
    data2 <- (data2 %>%
                group_by(Genus) %>%
                summarize(Abundance = sum(Abundance)) %>%
                ungroup() %>% data.frame()%>%
                arrange(desc(Abundance)))
    #on rajoute a la fin le other 
    order_data_glom <- c(order_data_glom, data2$Genus, data1[,-1]$Genus)#La on fait le vecteur de nom de genus pour avoir l'ordre souhaité 
  }
  #elever le NA de début
  order_data_glom <- order_data_glom[-1]
  #ajout du Other 2%
  order_data_glom <- c(order_data_glom, "Others <4%")
  #reordonner le tableau final 
  data_glom$Genus <- factor(data_glom$Genus, levels = order_data_glom)
  
  
  return(data_glom)
}



final.bar.plot.data.glom <- function(data_glom){
  #déduplication des lignes des ASv qui ont le meme genre pour le format svg 
  data_glom <- data_glom %>%
    group_by(Sample, type_enrich, Genus) %>%
    summarize(Sum_Abundance = sum(Abundance))
  return(data_glom)
}


























#Calculs----


data_glom_general1 <- bar.plot.data.glom(physeq_16S_nit_p, 2, "type_enrich")

data_glom_general <- final.bar.plot.data.glom(data_glom_general1)

oder_glom_general <- order.vec.barplot(data_glom_general1)



#Graphiques-----
#Sélection des couleurs 

#palette couleur Gammaproteobacteria
pal1 <- colorRampPalette(c("#E7D4E8", "#C2A5CF", "#9970B5", "#762A83"))

#nb de gamma 
length(unique(subset(data_glom, Class== "Gammaproteobacteria")$Genus)) # = 4 


#palette couleur Alphaproteobacteria 

# pal3 <- colorRampPalette(c("#663300", "#994C00", "#CC6600", "#FF8000", "#FF9933", "#FFB266", "#FFCC99", "#FFE5CC", "#FFE0E0", "#FFCCCC", "#FFBFBF", "#FFAFAF", "#FF9F9F", "#FF8F8F", "#FF7F7F"))
pal2 <- colorRampPalette(c("#F1E5D5", "#E5C29F", "#D9A46A", "#CD8740", "#C16A17", "#A95611", "#8C450B", "#703508", "#572B06", "#3E2103", "#251800"))

#nb de Alpha 
length(unique(subset(data_glom, Class== "Alphaproteobacteria")$Genus)) # = 15 


#palette couleur Bacteroidia
pal3 <- colorRampPalette(c("#EFF3FF", "#C3DAF6", "#96C1ED", "#6AA8E4", "#3D8FD9", "#1176D0", "#0F63BA", "#0E509E", "#0D3D82", "#0C2A66"))

#nb de bactero 
length(unique(subset(data_glom, Class== "Bacteroidia")$Genus)) # = 10 


#palette couleur Acidimicrobiia
pal4 <- colorRampPalette(c("#CCEBC5", "#8FD68E", "#5EBF57"))

#nb de Acidi 
length(unique(subset(data_glom, Class== "Acidimicrobiia")$Genus)) # = 3 

#palette couleur Deinococci
pal5 <- colorRampPalette(c("#FEE5D9", "#FCAE91"))

#nb de Deinococci 
length(unique(subset(data_glom, Class== "Deinococci")$Genus)) # = 2 

#palette Others
pal8 <- c("#E0E0E0")


#Palette générale :
pal_fan <- c(pal3(7),pal1(4),pal2(11), pal4(2), pal5(2) ,pal6(2),pal8)



#Plot 

ggplot(data = data_glom_general1, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  scale_fill_manual(values = pal_fan, breaks = oder_glom_general) +
  guides(fill = guide_legend(ncol = 2))





#Ecriture article ----

#ptites features en chiffre de l'article 

#moyenne d'une classe sur le jeu de données global 
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR"), Class == "Bacteroidia")@otu_table)),)
mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type == "biofilm"), Class == "Alphaproteobacteria")@otu_table))
mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "filter_SW"), Class == "Alphaproteobacteria")@otu_table))
mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "filter_ENR"), Class == "Alphaproteobacteria")@otu_table))

round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type == "water"), Class == "Bacteroidia")@otu_table)),1)


round(mean(rowSums(subset_taxa(physeq_16S_nit_p, Class %in% c("Alphaproteobacteria", "Bacteroidia", "Gammaproteobacteria"))@otu_table)),1)
round(mean(rowSums(subset_taxa(physeq_16S_nit_p, Class == "Alphaproteobacteria")@otu_table)),1)
round(mean(rowSums(subset_taxa(physeq_16S_nit_p, Class == "Bacteroidia")@otu_table)),1)
round(mean(rowSums(subset_taxa(physeq_16S_nit_p, Class == "Gammaproteobacteria")@otu_table)),1)

#Regarder au type 
mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type == "biofilm"), Class == "Alphaproteobacteria")@otu_table))
mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type == "biofilm"), Class == "Gammaproteobacteria")@otu_table))

round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Class == "Alphaproteobacteria")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Class == "Gammaproteobacteria")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Class == "Bacteroidia")@otu_table)),1)


mean(rowSums(subset_taxa(physeq_16S_nit_p, Class %in% c("Alphaproteobacteria"))@otu_table))
mean(rowSums(subset_taxa(physeq_16S_nit_p, Class %in% c("Bacteroidia"))@otu_table))
mean(rowSums(subset_taxa(physeq_16S_nit_p, Class %in% c("Gammaproteobacteria"))@otu_table))

#moyenne du genre 
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_SW"), Genus == "Granulosicoccus")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_ENR"), Genus == "Granulosicoccus")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_SW"), Genus == "Glaciecola")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_ENR"), Genus == "Glaciecola")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_SW"), Genus == "Agaribacterium")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_ENR"), Genus == "Agaribacterium")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_SW"), Genus == "Jannaschia")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_ENR"), Genus == "Jannaschia")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_SW"), Genus == "Croceitalea")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_ENR"), Genus == "Croceitalea")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_SW"), Genus == "Maribacter")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_ENR"), Genus == "Maribacter")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_SW"), Genus == "Winogradskyella")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_ENR"), Genus == "Winogradskyella")@otu_table)),1)

round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Genus == "Aurantivirga")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR"), Genus == "Aurantivirga")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Genus == "Polaribacter")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR"), Genus == "Polaribacter")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Genus == "Marivita")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR"), Genus == "Marivita")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Genus == "Nereida")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR"), Genus == "Nereida")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Genus == "Yoonia-Loktanella")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR"), Genus == "Yoonia-Loktanella")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Genus == "NS3a marine group")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR"), Genus == "NS3a marine group")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Genus == "Ponticoccus")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR"), Genus == "Ponticoccus")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_SW"), Genus == "Polaribacter")@otu_table)),1)
round(mean(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR"), Genus == "Polaribacter")@otu_table)),1)

#>2% analysis 
#Number of genera > 2% in type
unique(subset(data_glom, type_enrich == "biofilm_ENR" )$Genus)
unique(subset(data_glom, type_enrich == "biofilm_SW" )$Genus)
unique(subset(data_glom, type_enrich == "filter_ENR" )$Genus)
unique(subset(data_glom, type_enrich %in% c("biofilm_ENR", "biofilm_SW") )$Genus)


#Mean of the abundance of genere in the >X% analysis
#by the support 

#biofilm
data.frame(subset(data_glom, type_enrich %in% c("biofilm_ENR", "biofilm_SW"))) %>% 
                    dplyr::group_by(OTU,Kingdom,Phylum,Class,Order,Family,Genus)%>% 
                    summarize(Abundance = mean(Abundance)) %>%
                    dplyr::group_by(Genus)%>% 
                    summarize(Abundance = sum(Abundance)) 
     
#filter
data.frame(subset(data_glom_general, type_enrich == "filter_SW")) %>% 
  dplyr::group_by(Genus)%>% 
  summarize(Abundance = sum(Sum_Abundance))

#Biofilm SW
data.frame(subset(data_glom_general, type_enrich == "biofilm_SW")) %>% 
  dplyr::group_by(Genus)%>% 
  summarize(Abundance = sum(Sum_Abundance))

#Biofilm ENR 
data.frame(subset(data_glom, type_enrich == "biofilm_ENR")) %>% 
  dplyr::group_by(OTU,Kingdom,Phylum,Class,Order,Family,Genus)%>% 
  summarize(Abundance = mean(Abundance)) %>%
  dplyr::group_by(Genus)%>% 
  summarize(Abundance = sum(Abundance)) 

#Filter SW
data.frame(subset(data_glom_general1, type_enrich == "filter_SW")) %>% 
  dplyr::group_by(OTU,Kingdom,Phylum,Class,Order,Family,Genus)%>% 
  summarize(Abundance = mean(Abundance)) %>%
  dplyr::group_by(Genus)%>% 
  summarize(Abundance = sum(Abundance))   

#Filter ENR
data.frame(subset(data_glom_general1, type_enrich == "filter_ENR")) %>% 
  dplyr::group_by(Class)%>% 
  summarize(Abundance = sum(Abundance))%>%
  data.frame()

#connaitre les genre < 2% 
subset_taxa(glom, Genus == "Croceitalea")@otu_table

  
#savoir quelles sont l'ordre des class 
(data_glom[, c("OTU","Abundance", "type_enrich","Class")] %>% group_by(OTU,Class) %>% summarize(Abundance = mean(Abundance))) %>% group_by(Class) %>% summarize(Abundance = sum(Abundance))%>%
  arrange(desc(Abundance)) %>% data.frame()




subset(data_glom, type_enrich=="biofilm_ENR") %>% group_by(OTU,Class)
#Nombre de classes différentes dans le jeu de données: 
unique(physeq_16S_nit@tax_table[,"Class"])








#Supplementary table ----

supp.data.glom <- function(physeq, percent, env_var){
  
  #agglomeration de l'objet phyloseq 
  
  #merge d'objet phyloseq selon une var environnementales 
  merged_physeq <- merge_samples(physeq, env_var)
  
  #Quand le nombre d'ech n'est pas égal entre catégorie , ifaut rajouter une étape 
  merged_physeq <- transform_sample_counts(merge_samples(physeq, env_var), function(x) x/sum(x)*100)
  merged_physeq@sam_data$type_enrich <- c("biofilm_ENR","biofilm_SW","filter_ENR","filter_SW")
  merged_physeq@sam_data$type <- c("biofilm", "biofilm", "filter", "filter")
  
  
  #Version améliorée avec les 4 variables 
  glom <- tax_glom(merged_physeq, taxrank = 'Genus', NArm = FALSE, bad_empty = c(""))
  
  
  #Modification de la taxonomie
  glom_tax <- matrix(data = NA, ncol = 1, nrow =nrow(glom@tax_table), dimnames = list(NULL,"taxonomy"))
  rownames(glom_tax) <- rownames(glom@tax_table)
  tmp_tax <- data.frame(tax_table(glom))    
  
  
  for (i in 2:6) {  # Boucle de la colonne Genus à la colonne Kingdom
    tmp_tax[, i] <- ifelse(is.na(tmp_tax[, i]), 
                           paste("Unassigned", tmp_tax[, (i-1)], sep = " "), 
                           tmp_tax[, i])
  }
  
  
  
  #Then replace the names in the new dataset
  for (i in 1:nrow(glom@tax_table)){
    glom_tax[i,1] <- paste0( "r__",tmp_tax[i,"Kingdom"],
                             ";p__",tmp_tax[i,"Phylum"],
                             ";c__",tmp_tax[i,"Class"],
                             ";o__",tmp_tax[i,"Order"],
                             ";f__",tmp_tax[i,"Family"],
                             ";g__",tmp_tax[i,"Genus"]
    )
  }
  
  
  
  
  #Then replace the names of ASV on the otu table 
  final_tab <- t(data.frame(glom@otu_table))
  
  
  for (row_name in rownames(final_tab)) {
    new_row_name <- glom_tax[rownames(glom_tax) == row_name, "taxonomy"]
    rownames(final_tab)[rownames(final_tab) == row_name] <- new_row_name
  }
  
  #finally reinteger the phyloseq object 
  glom@otu_table <- otu_table(t(final_tab), taxa_are_rows = F)
  
  
  
  data_glom <- psmelt(glom)
  
  for (i in 27:32) {  # Boucle de la colonne Genus à la colonne Kingdom
    data_glom[, i] <- ifelse(is.na(data_glom[, i]), 
                             paste("unassigned", data_glom[, (i-1)], sep = " "), 
                             data_glom[, i])
  }
  
  
  data_glom$Genus <- as.character(data_glom$Genus)
  data_glom$Class <- as.character(data_glom$Class)
  data_glom <-data.frame(data_glom %>% dplyr::group_by(OTU,Sample,Abundance,type_enrich,Kingdom,Phylum,Class,Order,Family,Genus)%>% summarize(Abundance = sum(Abundance)))
  
  
  for (c in unique(data_glom$Class)) {
    tmp <- data_glom[data_glom$Class == c, ]
    other_genus <- tmp[, "Genus"]
    if (max(tmp$Abundance) < 2) {
      data_glom[rownames(data_glom) %in% rownames(unique(tmp)),][data_glom[rownames(data_glom) %in% rownames(unique(tmp)),]$Abundance < 2,]$Genus <- "Other bacterial classes <2%"
    }
    if (max(tmp$Abundance) >= 2){
      data_glom[rownames(data_glom) %in% rownames(unique(tmp)),][data_glom[rownames(data_glom) %in% rownames(unique(tmp)),]$Abundance < 2,]$Genus <- paste0("Other <2% ",c)
    }
  }
  
  
  # Réordonner le tableau en fonction de l'ordre obtenu
  data_glom_ordered <- data_glom
  
  #calcul des moyennes d'abondances
  data_glom_ordered <- data_glom_ordered %>% 
    group_by(Class, Genus) %>% 
    summarize(Abundance = sum(Abundance)) %>% 
    ungroup()
  
  #réordonner suivant la classe
  
  #on extrait Other < 2% on fait la moyenne abund par colass on ordonne et on récupère les noms 
  
  class_order <- (subset(data_glom_ordered, !grepl("^Other", Genus)) %>%
                    group_by(Class) %>%
                    summarize(Abundance = sum(Abundance)) %>%
                    ungroup() %>% data.frame()%>%
                    arrange(desc(Abundance)))$Class
  
  data_glom_ordered <- data_glom_ordered %>%
    arrange(factor(Class, levels = class_order))
  
  #réordonner suivant le genre
  
  #pour que dans chaque classe les genres soient ordonnés par abondance avec other a la fin 
  order_data_glom <- NA
  for (i in class_order){
    data <- subset(data_glom_ordered, Class==i)#extraire les données pour chaque classe
    data1 <- subset(data, grepl("^Other", Genus))#extraire Other 
    data2 <- subset(data, !grepl("^Other", Genus))#extraire tout le reste
    #ON réarange par ordre d'abondance 
    data2 <- (data2 %>%
                group_by(Genus) %>%
                summarize(Abundance = sum(Abundance)) %>%
                ungroup() %>% data.frame()%>%
                arrange(desc(Abundance)))
    #on rajoute a la fin le other 
    order_data_glom <- c(order_data_glom, data2$Genus, data1[,-1]$Genus)#La on fait le vecteur de nom de genus pour avoir l'ordre souhaité 
  }
  #elever le NA de début
  order_data_glom <- order_data_glom[-1]
  #ajout du Other 2%
  order_data_glom <- c(order_data_glom, "Other bacterial classes <2%")
  
  
  # Réorganiser le tableau en fonction du vecteur order_data_glom
  data_glom <- data_glom %>% 
    arrange(factor(Genus, levels = order_data_glom))
  
  # Dernières modifications du tableau en pivotant
  data_glom <- data_glom[,c("OTU","Abundance","type_enrich","Genus","Class")]
  data_glom <- data_glom %>%
    pivot_wider(names_from = type_enrich, values_from = Abundance)
  
  colnames(data_glom) <- c("Taxa","Graphical denomination","Class", "Biofilm ENR", "Biofilm SW", "Water ENR", "Water SW")
  
  return(data_glom)
}

supp_data_glom_general <- supp.data.glom(physeq_16S_nit_p, 2, "type_enrich")



write.table(supp_data_glom_general, "alpha_diversity/bar_plot/objets/supp_data_glom_general.txt", sep = "\t", row.names= FALSE, quote =FALSE)






  
  