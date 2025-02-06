#barplot_detail

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)



bar.plot.data <- function(physeq, percent, env_var, env_var_val){
  
  #subset of the pyloseq object
  # Créer un masque logique pour filtrer les échantillons
  mask <- rownames(physeq@sam_data[physeq@sam_data$type_enrich %in% env_var_val,])
  # Sous-ensemble de l'objet phyloseq en utilisant le masque
  physeq1 <- prune_samples(mask, physeq)
  
  
  
  
  glom <- tax_glom(physeq1 , taxrank = 'Genus', NArm = FALSE, bad_empty = c(""))
  
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
  
  
  #Ajout de la colonne date pour le plot 
  data_glom <- data_glom %>% mutate(Date = sub("(.{6}).*?_([A-Z])_A", "\\1_\\2", Sample))
  data_glom <- data_glom %>% mutate(Date = ifelse(grepl("_F$", Date), substr(Date, 1, 6), Date))
  

}





order.vec.barplot <- function(data_glom){
  # Réordonner le tableau en fonction de l'ordre obtenu
  data_glom_ordered <- data_glom

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
    order_data_glom <- c(order_data_glom, as.character(data2$Genus), as.character(data1$Genus[1]))#La on fait le vecteur de nom de genus pour avoir l'ordre souhaité 
  }
  #elever le NA de début
  order_data_glom <- order_data_glom[-1]
  #ajout du Other 2%
  order_data_glom <- c(order_data_glom, "Others <4%")
  
return(order_data_glom)
}



final.bar.plot.data <- function(data_glom){
  #déduplication des lignes des ASv qui ont le meme genre pour le format svg 
  data_glom <- data_glom %>%
    group_by(Sample, type_enrich, Genus, Date) %>%
    summarize(Sum_Abundance = sum(Abundance))
  return(data_glom)
}




sd.calculation <- function(physeq, tax_level, env_var, env_var_val){
  #subset of the pyloseq object
  # Créer un masque logique pour filtrer les échantillons
  mask <- rownames(physeq@sam_data[physeq@sam_data$type_enrich %in% env_var_val,])
  # Sous-ensemble de l'objet phyloseq en utilisant le masque
  physeq1 <- prune_samples(mask, physeq)
  
  
  
  
  glom <- tax_glom(physeq1 , taxrank = "Genus", NArm = FALSE, bad_empty = c(""))
  
  data_glom <- psmelt(glom)
  
  for (i in 27:32) {  # Boucle de la colonne Genus à la colonne Kingdom
    data_glom[, i] <- ifelse(is.na(data_glom[, i]), 
                             paste("unassigned", data_glom[, (i-1)], sep = " "), 
                             data_glom[, i])
  }
  
  data_glom$Genus <- as.character(data_glom$Genus)
  data_glom$Class <- as.character(data_glom$Class)
  data_glom <-data.frame(data_glom %>% dplyr::group_by(OTU,Sample,Abundance,type_enrich,Kingdom,Phylum,Class,Order,Family,Genus)%>% summarize(Abundance = sum(Abundance)))
  
  
  #sd calculation 
  if(tax_level == "Class"){data_glom <- data_glom %>%
    group_by(Class) %>%  # Regroupez les donnéess
    summarize(SD_Abundance = sd(Abundance)) 
  }
  if(tax_level == "Genus"){data_glom <-data_glom %>%
      group_by(Genus) %>%  # Regroupez les données
      summarize(SD_Abundance = sd(Abundance)) 
  }
  return(data_glom)
}




#Biofilm plot----

data_glom_biof1 <- bar.plot.data(physeq_16S_nit_p, env_var = "type_enrich", env_var_val= c("biofilm_SW", "biofilm_ENR"))
data_glom_biof <- final.bar.plot.data(data_glom_biof1)
oder_glom_biof <- order.vec.barplot(data_glom_biof1)
  
#Sélection des couleurs

#palette couleur Gammaproteobacteria
pal1biof <- colorRampPalette(c("#E7D4E8", "#C2A5CF", "#9970B5", "#762A83"))

#nb de gamma
length(unique(subset(data_glom_biof1, Class== "Gammaproteobacteria")$Genus)) # = 4

#palette couleur Bacteroidia
pal2biof <- colorRampPalette(c("#EFF3FF", "#C3DAF6", "#96C1ED", "#6AA8E4", "#3D8FD9", "#1176D0", "#0F63BA", "#0E509E", "#0D3D82", "#0C2A66"))

#nb de bactero
length(unique(subset(data_glom_biof1, Class== "Bacteroidia")$Genus)) # = 10

#palette couleur Alphaproteobacteria
pal3biof <-colorRampPalette(c("#F1E5D5",  "#D9A46A",  "#C16A17", "#8C450B",  "#572B06","#FFDFDF","#FF7F7F","#B22222","#960A00","#251800"))


#nb de Alha
length(unique(subset(data_glom_biof1, Class== "Alphaproteobacteria")$Genus)) # = 15

#palette couleur Acidimicrobiia
pal4biof <- colorRampPalette(c("#CCEBC5", "#8FD68E", "#5EBF57"))

#nb de Acidi
length(unique(subset(data_glom_biof1, Class== "Acidimicrobiia")$Genus)) # = 3

#palette couleur Deinococci
pal5biof <- colorRampPalette(c("#FEE5D9", "#FCAE91"))

#nb de Deinococci
length(unique(subset(data_glom_biof1, Class== "Deinococci")$Genus)) # = 2

#palette couleur Actinobacteria
pal6biof <- colorRampPalette(c("#FFF4D2", "#FFE89B", "#FFDE63"))

#nb de Actinobacteria
length(unique(subset(data_glom_biof1, Class== "Actinobacteria")$Genus)) # = 3

#palette couleur Rhodothermia
pal7biof <- colorRampPalette( c("#7FFFD4", "#40E0D0"))

#nb de Rhodothermia
length(unique(subset(data_glom_biof1, Class== "Rhodothermia")$Genus)) # = 2

#palette Others
pal8biof <- c("#E0E0E0")

#Palette générale :
pal_fan_biof <- c(pal1biof(length(unique(subset(data_glom_biof1, Class== "Gammaproteobacteria")$Genus))),
             pal3biof(length(unique(subset(data_glom_biof1, Class== "Alphaproteobacteria")$Genus))),
             pal2biof(length(unique(subset(data_glom_biof1, Class== "Bacteroidia")$Genus))),
             pal4biof(length(unique(subset(data_glom_biof1, Class== "Acidimicrobiia")$Genus))),
             pal5biof(length(unique(subset(data_glom_biof1, Class== "Deinococci")$Genus))),
             pal8biof)


#Plot
ggplot(data = data_glom_biof, aes(x = Date, y = Sum_Abundance, fill = Genus)) +
  facet_wrap(~type_enrich, scales = "free_x",
             labeller = labeller(type_enrich = c("biofilm_SW" = "biofilm SW", "biofilm_ENR" = "biofilm ENR"))) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  scale_fill_manual(values = pal_fan_biof, breaks = oder_glom_biof) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(y = "Relative sequence abundance") +
  theme(legend.text = element_text(face = "italic", family = "Calibri"), text = element_text(family = "Calibri"))






# Water plot ----

data_glom_water1 <- bar.plot.data(physeq_16S_nit_p, env_var = "type_enrich", env_var_val= c("water_SW", "water_ENR"))
data_glom_water <- final.bar.plot.data(data_glom_water1)
oder_glom_water <- order.vec.barplot(data_glom_water1)

#Sélection des couleurs

#palette couleur Gammaproteobacteria
pal1wat <- colorRampPalette(c("#E7D4E8", "#762A83"))

#nb de gamma
length(unique(subset(data_glom_water1, Class== "Gammaproteobacteria")$Genus)) # = 4

#palette couleur Bacteroidia
pal2wat <- colorRampPalette(c("#EFF3FF", "#C3DAF6", "#96C1ED", "#1176D0",  "#0E509E",  "#0C2A66"))

#nb de bactero
length(unique(subset(data_glom_water1, Class== "Bacteroidia")$Genus)) # = 10

#palette couleur Alphaproteobacteria
pal3wat <-colorRampPalette(c("#F1E5D5",  "#D9A46A",  "#C16A17", "#8C450B",  "#572B06","#FFDFDF","#FF7F7F","#B22222","#960A00","#251800"))

#nb de Alha
length(unique(subset(data_glom_water1, Class== "Alphaproteobacteria")$Genus)) # = 15

#palette couleur Acidimicrobiia
pal4wat <- colorRampPalette(c("#CCEBC5", "#5EBF57"))

#nb de Acidi
length(unique(subset(data_glom_water1, Class== "Acidimicrobiia")$Genus)) # = 3

#palette couleur Deinococci
pal5wat <- colorRampPalette(c("#FEE5D9"))

#nb de Deinococci
length(unique(subset(data_glom_water1, Class== "Deinococci")$Genus)) # = 2

#palette couleur Actinobacteria
pal6wat <- colorRampPalette(c("#FFF4D2", "#FFE89B", "#FFDE63"))

#nb de Actinobacteria
length(unique(subset(data_glom_water1, Class== "Actinobacteria")$Genus)) # = 3

#palette couleur Rhodothermia
pal7wat <- colorRampPalette( c("#7FFFD4", "#40E0D0"))

#nb de Rhodothermia
length(unique(subset(data_glom_water1, Class== "Rhodothermia")$Genus)) # = 2

#palette Others
pal8wat <- c("#E0E0E0")


#Palette générale :
pal_fan_wat <- c(pal2wat(length(unique(subset(data_glom_water1, Class== "Bacteroidia")$Genus))),
             pal3wat(length(unique(subset(data_glom_water1, Class== "Alphaproteobacteria")$Genus))),
             pal1wat(length(unique(subset(data_glom_water1, Class== "Gammaproteobacteria")$Genus))),
             pal6wat(length(unique(subset(data_glom_water1, Class== "Actinobacteria")$Genus))),
             pal4wat(length(unique(subset(data_glom_water1, Class== "Acidimicrobiia")$Genus))),
             pal7wat(length(unique(subset(data_glom_water1, Class== "Rhodothermia")$Genus))),
             pal8wat)


#Plot

ggplot(data = data_glom_water, aes(x = Date, y = Sum_Abundance, fill = Genus)) +
  facet_wrap(~type_enrich, scales = "free_x",
             labeller = labeller(type_enrich = c("water_SW" = "water SW", "water_ENR" = "water ENR"))) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  scale_fill_manual(values = pal_fan_wat, breaks = oder_glom_water) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(y = "Relative sequence abundance") +
  theme(legend.text = element_text(face = "italic", family = "Calibri"), text = element_text(family = "Calibri"))






#Ecriture article---- 




#Etude des fluctuastion sur données plot----
subset(subset(data_glom_water, type_enrich == "water_SW")%>%group_by(Genus) %>%
         ungroup() %>% data.frame(), Genus =="Aurantivirga")

#highest class fluctuations 
max(subset(subset(data_glom_water, type_enrich == "water_SW")%>%group_by(Class,Sample)%>%
         summarize(Abundance = sum(Sum_Abundance)) %>%
         ungroup() %>% data.frame(), Class =="Bacteroidia")[,"Abundance"])- 
min(subset(subset(data_glom_water, type_enrich == "water_SW")%>%group_by(Class,Sample)%>%
         summarize(Abundance = sum(Sum_Abundance)) %>%
         ungroup() %>% data.frame(), Class =="Bacteroidia")[,"Abundance"])

#fluctuantion of GEnus and genera for water 
subset(subset(data_glom_water, type_enrich == "water_ENR")%>%group_by(Genus,Sample)%>%
         summarize(Abundance = sum(Sum_Abundance)) %>%
         ungroup() %>% data.frame(), Genus =="Nereida")


subset(subset(data_glom_water1, type_enrich == "water_SW")%>%group_by(Class,Sample)%>%
         summarize(Abundance = sum(Abundance)) %>%
         ungroup() %>% data.frame(), Class =="Gammaproteobacteria")



##Fluctuantion of Genus biofilm ----
#fonctionne en complémentarité avec les fonction dessous pour uniquement les unassigned etc etc sinon préferer l'autre méthode
data_glom_biof1_abundview <- data_glom_biof1
data_glom_biof1_abundview$Date <- substr(data_glom_biof1_abundview$Date, 1, 6)

#Pour rhodo
  data_glom_biof1_abundview_rodo <- rbind(data_glom_biof1_abundview,c("ASV1730", "210309_ENR_C_A", 0.00, "biofilm_ENR", "Bacteria", "Proteobacteria", "Alphaproteobacteria", "Rhodobacterales", "Rhodobacteraceae", "unassigned Rhodobacteraceae","Alphaproteobacteria", 210309))
  data_glom_biof1_abundview_rodo$Abundance <- as.numeric(data_glom_biof1_abundview_rodo$Abundance)
  data_glom_biof1_abundview_rodo$Date <- substr(data_glom_biof1_abundview_rodo$Date, 1, 6)
  
  
  test <- data_glom_biof1_abundview_rodo %>% subset(type_enrich == "biofilm_SW" & Genus =="unassigned Rhodobacteraceae") %>%
    group_by(Genus,Sample,Date) %>% summarize(Abundance = sum(Abundance)) %>%  ungroup() %>% data.frame() %>%
    group_by(Genus, Date) %>% summarize(Abundance = mean(Abundance), SD = sd(Abundance))
  test$Abundance <- round(test$Abundance,1)
  test
  test <- data_glom_biof1_abundview_rodo %>% subset(type_enrich == "biofilm_ENR" & Genus =="unassigned Rhodobacteraceae") %>%
    group_by(Genus,Sample,Date) %>% summarize(Abundance = sum(Abundance)) %>%  ungroup() %>% data.frame() %>%
    group_by(Genus, Date) %>% summarize(Abundance = sd(Abundance))
  test$Abundance <- round(test$Abundance,1)
  test
  
  subset(data_glom_biof1_abundview_rodo, type_enrich == "biofilm_ENR" & Genus =="unassigned Rhodobacteraceae") %>%
    group_by(Sample, Genus,Date) %>% summarise(Abundance=sum(Abundance))%>%
    group_by(Genus) %>% summarise(Abundance=mean(Abundance))

  subset(test, type_enrich == "biofilm_ENR" & Genus =="unassigned Rhodobacteraceae") 
  
  
#etudider les >X% a une date 
subset(data_glom_water, Sample == "210507_ENR_F")%>%group_by(Genus,Sample)%>%
         summarize(Abundance = sum(Sum_Abundance)) %>%
         ungroup() %>% data.frame()

#etudier les classes dominantes
subset(data_glom_water, type_enrich == "water_ENR")%>%group_by(Class, Sample) %>%
  summarize(Abundance = sum(Sum_Abundance)) %>%
  ungroup() %>% data.frame()%>%
  arrange(desc(Abundance))



#le genre le plus abondant 
unique((subset(data_glom_water, type_enrich == "water_ENR") %>%
          group_by(Sample) %>%
          slice_max(Sum_Abundance) %>%
          ungroup())$Genus)

#version sans les others 
unique((subset(data_glom_biof %>%
                 filter(!grepl("^Others", as.character(Genus)))
               , type_enrich == "biofilm_ENR") %>%
          group_by(Sample) %>%
          slice_max(Sum_Abundance) %>%
          ungroup())$Genus)



#Etude globale article---- 



##Etude genus----
study.genus.RSA <- function(genus, samp_type){
  genus_name <- as.character(substitute(genus))
  
  # Abondance RSA d'un genre type d'ech en %
  abondance <- eval(bquote(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == .(samp_type)), Genus == .(genus_name))@otu_table)))
  
  # Récupérer les noms des échantillons
  sample_names <- eval(bquote(rownames(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == .(samp_type)), Genus == .(genus_name))@otu_table)))
  
  # Si le type d'échantillon est un biofilm, calculer la moyenne et l'écart type des triplicats
  if(samp_type == "biofilm_SW" | samp_type == "biofilm_ENR"){
    mean_abondance <- sapply(seq(1, length(abondance), by=3), function(i){
      mean(abondance[i:(i+2)])
    })
    sd_abondance <- sapply(seq(1, length(abondance), by=3), function(i){
      sd(abondance[i:(i+2)])
    })
    sample_names <- sample_names[seq(1, length(sample_names), by=3)]
  } else {
    mean_abondance <- abondance
    sd_abondance <- sd(mean_abondance)  # Pas de triplicats, donc pas de SD calculable
  }
  
  # min, max, moyenne globale et SD globale
  min_abondance <- min(mean_abondance, na.rm = TRUE)
  max_abondance <- max(mean_abondance, na.rm = TRUE)
  global_mean_abondance <- mean(mean_abondance, na.rm = TRUE)
  global_sd_abondance <- mean(sd_abondance, na.rm = TRUE)
  
  # Afficher les résultats
  cat("Abondance par date (Moyenne ± SD):", "\n")
  for (i in 1:length(mean_abondance)) {
    # cat(sample_names[i], ":", round(mean_abondance[i], 1), "±", round(sd_abondance[i], 1), "\n")
    cat(sample_names[i], ":", mean_abondance[i], "±", round(sd_abondance[i], 1), "\n")
  }
  cat("Min:", round(min_abondance, 1), "\n")
  cat("Max:", round(max_abondance, 1), "\n")
  cat("Moyenne globale:", round(global_mean_abondance, 1), "\n")
  cat("SD globale:", round(global_sd_abondance, 1), "\n")
}

study.genus.RSA("Nereida","water_SW")
study.genus.RSA("Nereida","water_ENR")
study.genus.RSA("Jannaschia","water_SW")
study.genus.RSA("Jannaschia","water_ENR")
study.genus.RSA("Yoonia-Loktanella","water_SW")
study.genus.RSA("Yoonia-Loktanella","water_ENR")
study.genus.RSA("Aurantivirga","water_SW")
study.genus.RSA("Aurantivirga","water_ENR")
study.genus.RSA("Polaribacter","water_SW")
study.genus.RSA("Polaribacter","water_ENR")
study.genus.RSA("Ponticoccus","water_SW")
study.genus.RSA("Ponticoccus","water_ENR")
study.genus.RSA("NS3a marine group","water_SW")
study.genus.RSA("NS3a marine group","water_ENR")
study.genus.RSA("Marivita","water_SW")
study.genus.RSA("Marivita","water_ENR")
study.genus.RSA("Nereida","water_ENR")
study.genus.RSA("Nereida","water_SW")
study.genus.RSA("Glaciecola","water_SW")
study.genus.RSA("Glaciecola","water_ENR")

study.genus.RSA("Glaciecola","biofilm_ENR")
study.genus.RSA("Glaciecola","biofilm_SW")
study.genus.RSA("Granulosicoccus","biofilm_SW")
study.genus.RSA("Granulosicoccus","biofilm_ENR")
study.genus.RSA("Jannaschia","biofilm_SW")
study.genus.RSA("Jannaschia","biofilm_ENR")
study.genus.RSA("Agaribacterium","biofilm_SW")
study.genus.RSA("Agaribacterium","biofilm_ENR")
study.genus.RSA("Croceitalea","biofilm_SW")
study.genus.RSA("Croceitalea","biofilm_ENR")
study.genus.RSA("Maribacter","biofilm_SW")
study.genus.RSA("Maribacter","biofilm_ENR")
study.genus.RSA("Winogradskyella","biofilm_SW")
study.genus.RSA("Winogradskyella","biofilm_ENR")
study.genus.RSA("Croceitalea","biofilm_SW")
study.genus.RSA("Croceitalea","biofilm_ENR")
study.genus.RSA("Litorimonas","biofilm_SW")
study.genus.RSA("Litorimonas","biofilm_ENR")
study.genus.RSA("Hyphomonas","biofilm_SW")
study.genus.RSA("Hyphomonas","biofilm_ENR")
study.genus.RSA("Roseobacter","biofilm_SW")
study.genus.RSA("Roseobacter","biofilm_ENR")

test <- tax_glom(physeq_16S_nit_p , taxrank = 'Genus', NArm = FALSE, bad_empty = c(""))
subset(psmelt(test), Family == "Cryomorphaceae" & is.na(Genus))%>% 
  group_by(type_enrich, date)%>%
  summarise(Abundance =mean(Abundance))%>%
  data.frame()



  

##Etude des genres----

study.genera.RSA <- function(physeq, samp_type, n_top_genera = 15){
  # Subset the physeq object based on the sample type
  physeq_sub <- eval(bquote(phyloseq::subset_samples(physeq, type_enrich == .(samp_type))))
  
  # Convert physeq object to dataframe
  df <- psmelt(physeq_sub)
  
  # Create a new column that contains either the Genus or the ASV ID if Genus is NA
  df$Genus_or_ASV <- ifelse(is.na(df$Genus), as.character(df$OTU), as.character(df$Genus))
  print(df)
  # Your existing code for biofilm or other sample types, but use Genus_or_ASV for grouping
  if(samp_type == "biofilm_SW" | samp_type == "biofilm_ENR"){
    summarized_df <- df %>%
      group_by(Sample, Genus_or_ASV) %>%
      summarise(Abundance = sum(Abundance))
    summarized_df$Sample <- paste0(substr(summarized_df$Sample, 1, 10),substr(summarized_df$Sample, 13, 14))
    summarized_df <- summarized_df %>%
      group_by(Sample, Genus_or_ASV) %>%
      summarise(Abundance = mean(Abundance))
  } else {
    summarized_df <- df %>%
      group_by(Sample, Genus_or_ASV) %>%
      summarise(Abundance = sum(Abundance))
  }
  
  # Rest of your code remains largely the same, but use Genus_or_ASV
  final_df <- summarized_df %>%
    spread(key = Genus_or_ASV, value = Abundance, fill = 0)
  
  # Order columns by colSums in descending order
  ordered_cols <- order(-colSums(final_df[, -1]))
  top_genera_cols <- c(1, ordered_cols[1:n_top_genera] + 1)
  final_df <- final_df[, top_genera_cols]
  
  print(final_df)
  # Calcul des statistiques descriptives pour chaque genre
  stats_list <- lapply(final_df[-1], function(x) c(mean = mean(x, na.rm = TRUE),
                                                   sd = sd(x, na.rm = TRUE),
                                                   max = max(x, na.rm = TRUE),
                                                   min = min(x, na.rm = TRUE)))
  
  stats_df <- as.data.frame(do.call(cbind, stats_list))
  print(stats_df)
  # Ajout d'une colonne "Sample" pour les statistiques
  stats_df$Sample <- c("mean", "sd", "max", "min")
  
  # Ajout des statistiques descriptives en dessous du tableau final
  final_df_with_stats <- data.frame(rbind(final_df, stats_df))
  
  rownames(final_df_with_stats) <- final_df_with_stats[,"Sample"]
  final_df_with_stats <- round(final_df_with_stats[,-1],1)
  # Affichage du tableau final avec les statistiques
  print(final_df_with_stats)
  
}




# Example usage:
study.genera.RSA(physeq_16S_nit_p, "water_SW")
study.genera.RSA(physeq_16S_nit_p, "water_SW")
study.genera.RSA(physeq_16S_nit_p, "biofilm_SW")
study.genera.RSA(physeq_16S_nit_p, "biofilm_ENR")

#voir rappidement la dominannce
for (i in 1:(nrow(test)-4)){
  print(colnames(test)[which.max(test[i,])])
}


##Etude family----

study.family.RSA <- function(family, samp_type){
  family_name <- as.character(substitute(family))
  
  # Abondance RSA d'un genre type d'ech en %
  abondance <- eval(bquote(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == .(samp_type)), family == .(family_name))@otu_table)))
  
  # Récupérer les noms des échantillons
  sample_names <- eval(bquote(rownames(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == .(samp_type)), family == .(family_name))@otu_table)))
  
  # min
  min_abondance <- min(abondance)
  
  # max
  max_abondance <- max(abondance)
  
  # SD d'un genre dans un type d'ech
  sd_abondance <- sd(abondance)
  
  # Moyenne d'un genre dans un type d'ech 
  mean_abondance <- mean(abondance)
  
  # Afficher les résultats
  cat("Abondance:", "\n")
  for (i in 1:length(abondance)) {
    cat(sample_names[i], ":", round(abondance[i],1), "\n")
  }
  cat("Min:", round(min_abondance,1), "\n")
  cat("Max:", round(max_abondance,1), "\n")
  cat("SD:", round(sd_abondance,1), "\n")
  cat("Moyenne:", round(mean_abondance,1), "\n")
}


##Etude by class----

study.class.RSA <- function(class, samp_type){
  class_name <- as.character(substitute(class))
  
  # Abondance RSA d'un genre type d'ech en %
  abondance <- eval(bquote(rowSums(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == .(samp_type)), Class == .(class_name))@otu_table)))
  
  # Récupérer les noms des échantillons
  sample_names <- eval(bquote(rownames(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == .(samp_type)), Class == .(class_name))@otu_table)))
  # Si le type d'échantillon est un biofilm, calculer la moyenne et l'écart type des triplicats
  if(samp_type == "biofilm_SW" | samp_type == "biofilm_ENR"){
    mean_abondance <- sapply(seq(1, length(abondance), by=3), function(i){
      mean(abondance[i:(i+2)])
    })
    sd_abondance <- sapply(seq(1, length(abondance), by=3), function(i){
      sd(abondance[i:(i+2)])
    })
    sample_names <- sample_names[seq(1, length(sample_names), by=3)]
  } else {
    mean_abondance <- abondance
    sd_abondance <- sd(mean_abondance, na.rm =T)  # Pas de triplicats, donc pas de SD calculable
  }
  
  # min, max, moyenne globale et SD globale
  min_abondance <- min(mean_abondance, na.rm = TRUE)
  max_abondance <- max(mean_abondance, na.rm = TRUE)
  global_mean_abondance <- mean(mean_abondance, na.rm = TRUE)
  global_sd_abondance <- mean(sd_abondance, na.rm = TRUE)
  
  # Afficher les résultats
  cat("Abondance par date (Moyenne ± SD):", "\n")
  for (i in 1:length(mean_abondance)) {
    cat(sample_names[i], ":", round(mean_abondance[i], 1), "±", round(sd_abondance[i], 1), "\n")
  }
  cat("Min:", round(min_abondance, 1), "\n")
  cat("Max:", round(max_abondance, 1), "\n")
  cat("Moyenne globale:", round(global_mean_abondance, 1), "\n")
  cat("SD globale:", round(global_sd_abondance, 1), "\n")
}


study.class.RSA("Gammaproteobacteria","biofilm_SW")
study.class.RSA("Gammaproteobacteria","biofilm_ENR")
study.class.RSA("Alphaproteobacteria","biofilm_SW")
study.class.RSA("Alphaproteobacteria","biofilm_ENR")
study.class.RSA("Bacteroidia","biofilm_SW")
study.class.RSA("Bacteroidia","biofilm_ENR")


study.class.RSA("Bacteroidia","water_SW")
study.class.RSA("Bacteroidia","water_ENR")
study.class.RSA("Alphaproteobacteria","water_SW")
study.class.RSA("Alphaproteobacteria","water_ENR")
study.class.RSA("Gammaproteobacteria","water_SW")
study.class.RSA("Gammaproteobacteria","water_ENR")



##Etude genera by class----

study.class.bygenera.RSA <- function(class, samp_type){
  class_name <- as.character(substitute(class))
  
  # Filter samples by specified sample type and taxa by the specified class
  filtered_physeq <- eval(bquote(subset_taxa(subset_samples(physeq_16S_nit_p, type_enrich == .(samp_type)), Class == .(class_name))))
  
  # Extract the OTU table
  glom_physeq <- tax_glom(filtered_physeq, taxrank = 'Genus', NArm = FALSE, bad_empty = c(""))
  
  otu_data <- as.data.frame(glom_physeq@otu_table)
  
  # Sum the abundance for each genus
  genera_abundance <- colSums(otu_data)
  
  # Sort the genera by abundance and take top 10
  top_genera <- names(sort(genera_abundance, decreasing = TRUE)[1:15])
  
  # Extract the taxonomy table
  tax_table_data <- as.data.frame(tax_table(glom_physeq))
  
  # Obtain genus names for the top_genera ASVs
  top_genera_names <- tax_table_data[top_genera, "Genus"]
  
  # Filter OTU data for these top genera
  top_genera_data <- otu_data[, top_genera, drop = FALSE]
  
  # Set the column names of top_genera_data to genus names
  colnames(top_genera_data) <- top_genera_names
  
  # If the sample type is a biofilm, average the triplicates
  if(samp_type == "biofilm_SW" | samp_type == "biofilm_ENR"){
    top_genera_data <- sapply(seq(1, ncol(top_genera_data)), function(j){
      sapply(seq(1, nrow(top_genera_data), by=3), function(i){
        mean(top_genera_data[i:(i+2), j])
      })
    })
    rownames(top_genera_data) <- unique(rownames(otu_data))[seq(1, nrow(otu_data), by=3)]
    colnames(top_genera_data) <- top_genera_names  # Set column names to genus names
  }
  
  # Compute statistics and add to the dataframe
  means <- colMeans(top_genera_data, na.rm = TRUE)
  sds <- apply(top_genera_data, 1, sd, na.rm = TRUE)
  mins <- apply(top_genera_data, 1, min, na.rm = TRUE)
  maxs <- apply(top_genera_data, 1, max, na.rm = TRUE)
  
  top_genera_data <- rbind(top_genera_data, means, sds, mins, maxs)
  rownames(top_genera_data)[(nrow(top_genera_data)-3):(nrow(top_genera_data))] <- c("Mean", "SD", "Min", "Max")
  
  # Round all numbers in the dataframe to one decimal places
  top_genera_data <- round(top_genera_data, 1)
  
  # Display the results
  print(top_genera_data)
}

study.class.bygenera.RSA("Gammaproteobacteria","biofilm_SW")
study.class.bygenera.RSA("Gammaproteobacteria","biofilm_ENR")
study.class.bygenera.RSA("Alphaproteobacteria","biofilm_SW")
study.class.bygenera.RSA("Alphaproteobacteria","biofilm_ENR")
study.class.bygenera.RSA("Bacteroidia","biofilm_SW")
study.class.bygenera.RSA("Bacteroidia","biofilm_ENR")

study.class.bygenera.RSA("Gammaproteobacteria","water_SW")
study.class.bygenera.RSA("Gammaproteobacteria","water_ENR")
study.class.bygenera.RSA("Alphaproteobacteria","water_SW")
study.class.bygenera.RSA("Alphaproteobacteria","water_ENR")
study.class.bygenera.RSA("Bacteroidia","water_SW")
study.class.bygenera.RSA("Bacteroidia","water_ENR")





##Variation----
#Etudier les variations entres classes et genres 
sd_water_SW_genus <-sd.calculation(physeq_16S_nit_p, tax_level ="Genus", env_var = "type_enrich", env_var_val= c("water_SW"))
sd_water_SW_class <- sd.calculation(physeq_16S_nit_p, tax_level = "Class", env_var = "type_enrich", env_var_val= c("water_SW"))

sd_water_ENR_genus <-sd.calculation(physeq_16S_nit_p, tax_level ="Genus", env_var = "type_enrich", env_var_val= c("water_ENR"))
sd_water_ENR_class <- sd.calculation(physeq_16S_nit_p, tax_level ="Class", env_var = "type_enrich", env_var_val= c("water_ENR"))

sd_biofilm_SW_genus <-sd.calculation(physeq_16S_nit_p, tax_level ="Genus", env_var = "type_enrich", env_var_val= c("biofilm_SW"))
sd_biofilm_SW_class <- sd.calculation(physeq_16S_nit_p, tax_level ="Class", env_var = "type_enrich", env_var_val= c("biofilm_SW"))

sd_biofilm_ENR_genus <-sd.calculation(physeq_16S_nit_p, tax_level ="Genus", env_var = "type_enrich", env_var_val= c("biofilm_ENR"))
sd_biofilm_ENR_class <- sd.calculation(physeq_16S_nit_p, tax_level ="Class", env_var = "type_enrich", env_var_val= c("biofilm_ENR"))


data.frame(sd_water_SW_class)%>% arrange(desc(SD_Abundance))
data.frame(sd_water_SW_genus)%>% arrange(desc(SD_Abundance))

data.frame(sd_water_ENR_class)%>% arrange(desc(SD_Abundance))
data.frame(sd_water_ENR_genus)%>% arrange(desc(SD_Abundance))

data.frame(sd_biofilm_SW_class)%>% arrange(desc(SD_Abundance))
data.frame(sd_biofilm_SW_genus)%>% arrange(desc(SD_Abundance))

data.frame(sd_biofilm_ENR_class)%>% arrange(desc(SD_Abundance))
data.frame(sd_biofilm_ENR_genus)%>% arrange(desc(SD_Abundance))


#quel écosysteme est le plys dynamique ? 
data.frame(sd_water_SW_class)$SD_Abundance %>% sum()
data.frame(sd_water_SW_class)$SD_Abundance %>% mean()

data.frame(sd_water_SW_genus)$SD_Abundance %>% sum()
data.frame(sd_water_SW_genus)$SD_Abundance %>% mean()

data.frame(sd_water_ENR_class)%>% arrange(desc(SD_Abundance))
data.frame(sd_water_ENR_genus)%>% arrange(desc(SD_Abundance))














#Supplemntary table ----


supp.data.glom.d <- function(physeq, percent, env_var, env_var_val){
  
  
  #subset of the pyloseq object
  # Créer un masque logique pour filtrer les échantillons
  mask <- rownames(physeq@sam_data[physeq@sam_data$type_enrich %in% env_var_val,])
  # Sous-ensemble de l'objet phyloseq en utilisant le masque
  physeq1 <- prune_samples(mask, physeq)
  
  
  
  
  glom <- tax_glom(physeq1 , taxrank = 'Genus', NArm = FALSE, bad_empty = c(""))
  
  
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
                             ",p__",tmp_tax[i,"Phylum"],
                             ",c__",tmp_tax[i,"Class"],
                             ",o__",tmp_tax[i,"Order"],
                             ",f__",tmp_tax[i,"Family"],
                             ",g__",tmp_tax[i,"Genus"]
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
    if (max(tmp$Abundance) < 4) {
      data_glom[rownames(data_glom) %in% rownames(unique(tmp)),][data_glom[rownames(data_glom) %in% rownames(unique(tmp)),]$Abundance < 4,]$Genus <- "Other bacterial classes <4%"
    }
    if (max(tmp$Abundance) >= 4){
      data_glom[rownames(data_glom) %in% rownames(unique(tmp)),][data_glom[rownames(data_glom) %in% rownames(unique(tmp)),]$Abundance < 4,]$Genus <- paste0("Other <4% ",c)
    }
  }
  
  
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
  order_data_glom <- c(order_data_glom, "Other bacterial classes <4%")
  
  print("ok")
  # Réorganiser le tableau en fonction du vecteur order_data_glom
  data_glom <- data_glom %>% 
    arrange(factor(Genus, levels = order_data_glom))
  # Dernières modifications du tableau en pivotant
  data_glom <- data_glom[,c("OTU","Abundance","Sample","Genus","Class")]
  data_glom <- data_glom %>%
    pivot_wider(names_from = Sample, values_from = Abundance)
  colnames(data_glom[,1:3]) <- c("Taxa","Graphical denomination","Class" )
  print("ok")
  return(data_glom)
}



supp_data_glom_d_b <- supp.data.glom.d(physeq_16S_nit_p, env_var = "type_enrich", env_var_val= c("biofilm_SW", "biofilm_ENR"))

write.table(supp_data_glom_d_b, "dynamics/objets/supp_data_glom_d_b.txt", sep = ";", row.names = FALSE, quote = FALSE)
write.csv2(data.frame(supp_data_glom_d_b), "dynamics/objets/supp_data_glom_d_b.csv")




supp_data_glom_d_w <- supp.data.glom.d(physeq_16S_nit_p, env_var = "type_enrich", env_var_val= c("water_SW", "water_ENR"))

write.table(data.frame(supp_data_glom_d_w), "dynamics/objets/supp_data_glom_d_w.txt", sep = ";", row.names = FALSE, quote = FALSE)
write.csv2(data.frame(supp_data_glom_d_w), "dynamics/objets/supp_data_glom_d_w.csv")
