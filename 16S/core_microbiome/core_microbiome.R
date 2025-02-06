library(cowplot)


extract.core <- function(physeq, type_val, sampmax){
  #récupérer le nombre de lignes
  nsamp <- nrow(otu_table(eval(bquote(subset_samples(physeq, type == .(type_val))))))
  #extraire le nom des ASV selon la condition
  core_asv <- names(colSums(otu_table(eval(bquote(subset_samples(physeq, type == .(type_val)))))>0)>=nsamp-sampmax)[
    colSums(otu_table(eval(bquote(subset_samples(physeq, type == .(type_val)))))>0)>=nsamp-sampmax]
  # Extraire les ASV souhaités de l'objet phyloseq
  physeq_core <- prune_taxa(core_asv, eval(bquote(subset_samples(physeq, type == .(type_val)))))
  print(physeq_core)
}



nb.asv.core.ech <- function(physeq, type_val, sampmax){
  # Récupérer le nombre de lignes
  nsamp <- nrow(otu_table(eval(bquote(subset_samples(physeq, type == .(type_val))))))
  # Initialiser un tableau pour stocker les résultats
  result_table <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(result_table) <- c("Echantillons", "Nombre_ASV_core")
  
  for(i in 1: nsamp){
    core_asv <- 0
    # Extraire le nom des ASV selon la condition
    core_asv <- names(colSums(otu_table(eval(bquote(subset_samples(physeq, type == .(type_val))))) > 0) >= i)[
      colSums(otu_table(eval(bquote(subset_samples(physeq, type == .(type_val))))) > 0) >= i]
    # Ajouter le résultat au tableau
    result_table <- rbind(result_table, data.frame(Nb_Echantillons_min = i, Nombre_ASV_core = length(core_asv)))
  }
  
  
  plot <- ggplot(data = result_table, aes(x= Nb_Echantillons_min, y=Nombre_ASV_core))+
    geom_bar(stat = "identity")+
    theme_minimal()
  
  # Afficher le tableau des résultats
  print(plot)
}


 



#Core micro biofilm----
nb.asv.core.ech(physeq_16S_nit_p, "biofilm")

core_biofilm <- extract.core(physeq_16S_nit_p, "biofilm",0)
core_biofilm_p <- transform_sample_counts(core_biofilm, function(x) x/sum(x)*100)


#Core micro water ----
nb.asv.core.ech(physeq_16S_nit_p, "water")
list.core(physeq_16S_nit, "water")

core_water <- extract.core(physeq_16S_nit_p, "water", 0)
core_water_p <- transform_sample_counts(core_water, function(x) x/sum(x)*100)




#Taxonmy of core----


##Biofilm----

###Plot----
# Préparation des données
core_biofilm_melt <- psmelt(core_biofilm_p)

# Convertir la colonne Abundance en numérique
core_biofilm_melt$Abundance <- as.numeric(core_biofilm_melt$Abundance)

# Remplacer les valeurs NA par "unassigned" et le niveau taxonomique supérieur
for (i in 27:32) {
  core_biofilm_melt[, i] <- ifelse(is.na(core_biofilm_melt[, i]), 
                                   paste("unassigned", core_biofilm_melt[, (i-1)], sep = " "), 
                                   core_biofilm_melt[, i])
}

# Calculer l'abondance relative
core_biofilm_melt <- core_biofilm_melt %>% 
  group_by(Sample, Class, Genus) %>% 
  summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance / sum(Abundance, na.rm = TRUE)) %>% 
  ungroup()

# Calculez l'abondance totale par classe
class_abundance <- core_biofilm_melt %>%
  group_by(Class) %>%
  summarise(TotalClassAbundance = sum(Abundance, na.rm = TRUE))

# Joignez avec le dataframe principal
core_biofilm_melt <- left_join(core_biofilm_melt, class_abundance, by = "Class")

# Réordonnez par abondance de classe puis par abondance relative de genre
core_biofilm_melt <- core_biofilm_melt %>%
  arrange(-desc(TotalClassAbundance), -desc(RelativeAbundance))

# Créer une variable ordonnée pour les genres
# core_biofilm_melt$Genus <- factor(core_biofilm_melt$Genus, levels = unique(core_biofilm_melt$Genus))
core_biofilm_melt$Genus <- factor(core_biofilm_melt$Genus, levels = factor(c("unassigned Microtrichaceae","unassigned Saprospiraceae","unassigned NS9 marine group",
                                                                             "Maribacter","Croceitalea","unassigned Hyphomonadaceae","Roseobacter","unassigned Rhizobiaceae",
                                                                             "Litorimonas","unassigned Rhodobacteraceae","Jannaschia","Methylotenera","Granulosicoccus")))

# Définir les palettes de couleurs
pal1biof <- colorRampPalette(c("#762A83","#E7D4E8")) #, "#C2A5CF", "#9970B5", ))
pal2biof <- colorRampPalette(c("#C3DAF6","#EFF3FF", "#96C1ED"))# , "#6AA8E4", "#3D8FD9", "#1176D0", "#0F63BA", "#0E509E", "#0D3D82", "#0C2A66"))
pal3biof <- colorRampPalette(c("#8C450B","#FFDFDF","#FF7F7F", "#C16A17", "#D9A46A","#F1E5D5"))#,   "#572B06","#B22222","#960A00","#251800"))
pal4biof <- colorRampPalette(c("#CCEBC5" ))#, "#8FD68E", "#5EBF57"))
pal5biof <- colorRampPalette(c("#FEE5D9", "#FCAE91"))
pal8biof <- c("#E0E0E0")

# Créer une liste de toutes les classes uniques
unique_classes <- unique(core_biofilm_melt$Class)

# Initialiser un vecteur de couleurs vide
final_color_vector <- c()

# Remplir le vecteur de couleurs
for(class in unique_classes) {
  current_palette <- switch(class,
                            "Gammaproteobacteria" = pal1biof,
                            "Alphaproteobacteria" = pal3biof,
                            "Bacteroidia" = pal2biof,
                            "Acidimicrobiia" = pal4biof,
                            "Deinococci" = pal5biof,
                            pal8biof) # default palette for "Others" and any other classes
  
  num_genus <- length(unique(subset(core_biofilm_melt, Class == class)$Genus))
  final_color_vector <- c(final_color_vector, current_palette(num_genus))
}

# Associer les couleurs aux genres uniques
color_mapping <- setNames(final_color_vector, unique(core_biofilm_melt$Genus))

# Déduplication des lignes des ASV qui ont le même Genus pour le format SVG
core_biofilm_melt <- core_biofilm_melt %>%
  group_by(Genus) %>%
  summarize(Mean_Abundance = mean(Abundance))




# Créer le graphique à secteurs
ggplot() +
  geom_bar(data= core_biofilm_melt, aes(x = "", y = Mean_Abundance, fill = Genus),width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Abondance relative des différents genres bactériens") +
  theme_minimal() +
  theme(legend.position = "bottom")



###Writing----


core_biof_tab <- psmelt(core_biofilm_p)
# Remplacer les valeurs NA par "unassigned" et le niveau taxonomique supérieur
for (i in 27:32) {
  core_biof_tab[, i] <- ifelse(is.na(core_biof_tab[, i]), 
                               paste("unassigned", core_biof_tab[, (i-1)], sep = " "), 
                               core_biof_tab[, i])
}
core_biof_tab[,c("Sample", "Class", "Abundance")]%>%
  group_by(Sample, Class)%>%
  summarise(sum_Abundance = sum(Abundance))%>%
  group_by(Class)%>%
  summarise(mean_Abundance = mean(sum_Abundance))

core_biof_tab[,c("Sample", "Genus", "Abundance")]%>%
  group_by(Sample, Genus)%>%
  summarise(sum_Abundance = sum(Abundance))%>%
  group_by(Genus)%>%
  summarise(mean_Abundance = mean(sum_Abundance))

#importance of core in total dataset
round(mean(data.frame(rowSums(extract.core(physeq_16S_nit_p, "biofilm", 0)@otu_table))[,1]),1)



##Water----

###Plot----
# Préparation des données
core_water_melt <- psmelt(core_water_p)

# Convertir la colonne Abundance en numérique
core_water_melt$Abundance <- as.numeric(core_water_melt$Abundance)

# Remplacer les valeurs NA par "unassigned" et le niveau taxonomique supérieur
for (i in 27:32) {
  core_water_melt[, i] <- ifelse(is.na(core_water_melt[, i]), 
                                 paste("unassigned", core_water_melt[, (i-1)], sep = " "), 
                                 core_water_melt[, i])
}

# Calculer l'abondance relative
core_water_melt <- core_water_melt %>% 
  group_by(Sample, Class, Genus) %>% 
  summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance / sum(Abundance, na.rm = TRUE)) %>% 
  ungroup()

# Calculez l'abondance totale par classe
class_abundance_water <- core_water_melt %>%
  group_by(Class) %>%
  summarise(TotalClassAbundance = sum(Abundance, na.rm = TRUE))

# Joignez avec le dataframe principal
core_water_melt <- left_join(core_water_melt, class_abundance_water, by = "Class")

# Réordonnez par abondance de classe puis par abondance relative de genre
core_water_melt <- core_water_melt %>%
  arrange(desc(TotalClassAbundance), desc(RelativeAbundance))

# Créer une variable ordonnée pour les genres
core_water_melt$Genus <- factor(core_water_melt$Genus, levels = unique(core_water_melt$Genus))



# Définir les palettes de couleurs
pal1water <- colorRampPalette(c( "#762A83"))#,"#E7D4E8", "#C2A5CF","#9970B5", ))
pal2water <- colorRampPalette(c("#B0C4DE", "#0E509E" ))#"#EFF3FF", "#C3DAF6", "#96C1ED", "#6AA8E4", "#3D8FD9", "#1176D0", "#0F63BA", "#0E509E", "#0D3D82", ))
pal3water <- colorRampPalette(c("#F1E5D5",  "#D9A46A",  "#C16A17"))#, "#8C450B",  "#572B06","#FFDFDF","#FF7F7F","#B22222","#960A00","#251800"))
pal4water <- colorRampPalette(c("#CCEBC5", "#8FD68E", "#5EBF57"))
pal5water <- colorRampPalette(c("#FEE5D9", "#FCAE91"))
pal8water <- c("#E0E0E0")
c("#C3DAF6","#EFF3FF", "#96C1ED"    , "#6AA8E4", "#3D8FD9", "#1176D0", "#0F63BA", "#0E509E", "#0D3D82", "#0C2A66", "#C3DAF6", "#B0C4DE")
# Créer une liste de toutes les classes uniques
unique_classes_water <- unique(core_water_melt$Class)

# Initialiser un vecteur de couleurs vide
final_color_vector_water <- c()

# Remplir le vecteur de couleurs
for(class in unique_classes_water) {
  current_palette <- switch(class,
                            "Gammaproteobacteria" = pal1water,
                            "Alphaproteobacteria" = pal3water,
                            "Bacteroidia" = pal2water,
                            "Acidimicrobiia" = pal4water,
                            "Deinococci" = pal5water,
                            pal8water) # default palette for "Others" and any other classes
  
  num_genus <- length(unique(subset(core_water_melt, Class == class)$Genus))
  final_color_vector_water <- c(final_color_vector_water, current_palette(num_genus))
}

# Associer les couleurs aux genres uniques
color_mapping_water <- setNames(final_color_vector_water, unique(core_water_melt$Genus))

# Déduplication des lignes des ASV qui ont le même Genus pour le format SVG
core_water_melt <- core_water_melt %>%
  group_by(Genus) %>%
  summarize(Mean_Abundance = mean(Abundance))

# Créer le graphique à secteurs
ggplot() +
  geom_bar(data= core_water_melt, aes(x = "", y = Mean_Abundance, fill = Genus),width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = color_mapping_water) +
  labs(title = "Abondance relative des différents genres bactériens dans l'eau") +
  theme_minimal() +
  theme(legend.position = "bottom")


###Writing----

core_water_tab <- psmelt(core_water_p)
# Remplacer les valeurs NA par "unassigned" et le niveau taxonomique supérieur
for (i in 27:32) {
  core_water_tab[, i] <- ifelse(is.na(core_water_tab[, i]), 
                               paste("unassigned", core_water_tab[, (i-1)], sep = " "), 
                               core_water_tab[, i])
}
core_water_tab[,c("Sample", "Class", "Abundance")]%>%
  group_by(Sample, Class)%>%
  summarise(sum_Abundance = sum(Abundance))%>%
  group_by(Class)%>%
  summarise(mean_Abundance = mean(sum_Abundance))

core_water_tab[,c("Sample", "Genus", "Abundance")]%>%
  group_by(Sample, Genus)%>%
  summarise(sum_Abundance = sum(Abundance))%>%
  group_by(Genus)%>%
  summarise(mean_Abundance = mean(sum_Abundance))


#importance of core in total dataset
mean(data.frame(rowSums(extract.core(physeq_16S_nit_p, "water", 0)@otu_table))[,1])


#writing for ASV abundance in whole dataset 

mean(rowSums(subset_samples(subset_taxa(prune_taxa(taxa_names(core_water), physeq_16S_nit_p), Class =="Bacteroidia"), type =="water")@otu_table))
mean(rowSums(subset_samples(subset_taxa(prune_taxa(taxa_names(core_water), physeq_16S_nit_p), Genus =="Nereida"), type == "water")@otu_table))





#supplementary table 

psmelt(core_biofilm_p) %>% group_by(OTU,Kingdom,Phylum,Class,Order,Family,Genus) %>%summarise(Abundance = mean(Abundance))

#tab taxonomy biofilm core 
core_tab_taxonomybiof <- psmelt(core_biofilm_p) %>% 
  mutate(Concatenated_Taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = ";"))%>%
  select(OTU,Sample_ID_E, Concatenated_Taxonomy, Genus, Abundance)

core_tab_taxonomybiof <- core_tab_taxonomybiof %>% 
  group_by(OTU,Concatenated_Taxonomy,Genus) %>%
  summarise(Abundance = mean(Abundance))

write.table(core_tab_taxonomybiof, file= "/16S_HOLOGREEN/core_microbiome/objets/core_tab_taxonomybiof.txt", sep = ",")

#tab taxonomy water core 
core_tab_taxonomywater <- psmelt(core_water_p) %>% 
  mutate(Concatenated_Taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = ";"))%>%
  select(OTU,Sample_ID_E, Concatenated_Taxonomy, Genus, Abundance)

core_tab_taxonomywater <- core_tab_taxonomywater %>% 
  group_by(OTU,Concatenated_Taxonomy,Genus) %>%
  summarise(Abundance = mean(Abundance))

write.table(core_tab_taxonomywater, file= "/16S_HOLOGREEN/core_microbiome/objets/core_tab_taxonomywater.txt", sep = ",")


#Values article----
round(sum(colMeans(core_biofilm@otu_table)),1)
round(sum(colMeans(core_water@otu_table)),1)

round(sum(colMeans(subset_taxa(core_water, Class =="Bacteroidia")@otu_table)),2)

round(sum(colMeans(subset_taxa(core_water, Genus =="Polaribacter")@otu_table)),1)

round(sum(colMeans(prune_taxa(c("ASV1730","ASV1854","ASV1892"), core_biofilm)@otu_table)),1)


#supp table  
write.table(
data.frame(rbind(biofilm = core_biof_tab%>% 
        group_by(OTU,Phylum,Class,Order,Family,Genus)%>%
        summarise(mean_Abund = mean(Abundance)), 
      water = core_water_tab%>% 
        group_by(OTU,Phylum,Class,Order,Family,Genus)%>%
        summarise(mean_Abund = mean(Abundance)))),
file = "F:/16S_HOLOGREEN/core_microbiome/objets/core_table.txt"
)







# #Differential ASV whole dataset 
# 
# ##Biofilm---
# 
# ###Aldex Maeva
# 
# 
# #Aldex2 package, voir tuto https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.pdf
# #clr_aldex2
# aldex_clr_biof<- aldex.clr( t(otu_table(subset_samples(physeq_16S_nit, type == "biofilm"))),
#                    sample_data(subset_samples(physeq_16S_nit, type == "biofilm"))$enrich,
#                    mc.samples=1000, denom="all", verbose=F)
# aldex_kw_biof <- aldex.kw(aldex_clr_biof)    #étape qui met bcp de temps (ici environ 2-3h je crois, si bcp d'ASV je laisse tourner sur la nuit). retourne un tableau avec ASV (rownames) et pvalue associée (4 colonnes, différentes méthodes).
# write.table(aldex_kw_biof,"aldex_kw_biof",sep = "\t")    #je sauve le output au cas où il y ait un bug, pour ne pas avoir à refaire tourner
# aldex_kw_biof_pval <- aldex_kw_biof[aldex_kw_biof$glm.eBH < 0.05,]  #garder seulement les ASV avec pvalue corrigée BH < 0.05
# 
# 
# diff_biofilm <- prune_taxa(rownames(aldex_kw_biof_pval),subset_samples(physeq_16S_nit_p, type=="biofilm"))
#  
# 
# ###Others tests
# # ancombc_biof <- ancombc.m(subset_samples(physeq_16S_nit, type =="biofilm"))
# # lefse_biof <- lefse.m(subset_samples(physeq_16S_nit, type =="biofilm"))
# # 
# # summary(ancombc_biof$feature %in% row.names(test5))
# # summary(lefse_biof$feature %in% row.names(test5))
# 
# 
# 
# 
# ##Water
# 
# 
# ###Aldex Maeva
# 
# 
# #Aldex2 package, voir tuto https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.pdf
# #clr_aldex2
# aldex_clr_water<- aldex.clr( t(otu_table(subset_samples(physeq_16S_nit, type == "water"))),
#                             sample_data(subset_samples(physeq_16S_nit, type == "water"))$enrich,
#                             mc.samples=1000, denom="all", verbose=F)
# aldex_kw_water <- aldex.kw(aldex_clr_water)    #étape qui met bcp de temps (ici environ 2-3h je crois, si bcp d'ASV je laisse tourner sur la nuit). retourne un tableau avec ASV (rownames) et pvalue associée (4 colonnes, différentes méthodes).
# write.table(aldex_kw_water,"aldex_kw_water",sep = "\t")    #je sauve le output au cas où il y ait un bug, pour ne pas avoir à refaire tourner
# aldex_kw_water_pval <- aldex_kw_water[aldex_kw_water$glm.eBH < 0.05,]  #garder seulement les ASV avec pvalue corrigée BH < 0.05
# 
# diff_water <- prune_taxa(rownames(aldex_kw_water_pval),subset_samples(physeq_16S_nit_p, type=="water"))
# 
# #comparaison
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# test <- names(colSums(subset_samples(physeq_16S_nit_p,type_enrich=="biofilm_SW")@otu_table)[colSums(subset_samples(physeq_16S_nit_p,type_enrich=="biofilm_SW")@otu_table)>0]) 
# test1 <- names(colSums(subset_samples(physeq_16S_nit_p,type_enrich=="biofilm_ENR")@otu_table)[colSums(subset_samples(physeq_16S_nit_p,type_enrich=="biofilm_ENR")@otu_table)>0])
# 
# test %in% test1
# mean(rowSums(prune_taxa(intersect(test, test1), subset_samples(physeq_16S_nit_p, type =="biofilm"))@otu_table))
# 

# #Differential abundances
# 
# ##ancomBC
# 
# ancombc.m <- function(physeq){
#   mm_ancombc <- run_ancombc_patched(
#     physeq,
#     group = "enrich",
#     taxa_rank = "none",
#     pvalue_cutoff = 0.001,
#     p_adjust = "fdr"
#   )
#   
#   mm_ancombc_table <- data.frame(mm_ancombc@marker_table)
#   return(mm_ancombc_table)
#   
# }
# 
# 
# # 
# # an_ef <- microbiomeMarker::plot_ef_bar(mm_ancombc)
# # y_labs <- ggplot_build(an_ef)$layout$panel_params[[1]]$y$get_labels()
# # an_abd <- microbiomeMarker::plot_abundance(mm_ancombc, group = "enrich") +
# #   scale_y_discrete(limits = y_labs)
# # gridExtra::grid.arrange(an_ef, an_abd, nrow = 1)
# 
# 
# 
# 
# ##LEFSE
# 
# lefse.m <- function(physeq){
#   mm_lefse <- microbiomeMarker::run_lefse(physeq, norm = "CPM",
#                                           wilcoxon_cutoff = 0.01,
#                                           group = "enrich",
#                                           taxa_rank = "none",
#                                           kw_cutoff = 0.01,
#                                           multigrp_strat = TRUE,
#                                           lda_cutoff = 4)
#   
#   mm_lefse_table <- data.frame(mm_lefse@marker_table)
#   return(mm_lefse_table)
# }
# 
# 
# 
# p_LDAsc <- microbiomeMarker::plot_ef_bar(mm_lefse)
# y_labs <- ggplot_build(p_LDAsc)$layout$panel_params[[1]]$y$get_labels()
# p_abd <- microbiomeMarker::plot_abundance(mm_lefse, group = "enrich") +
#   scale_y_discrete(limits = y_labs)
# gridExtra::grid.arrange(p_LDAsc, p_abd, nrow = 1)
# 
# 
# 
# ##ALDEX
# 
# aldex.m <- function(physeq){
#   mm_aldex <- microbiomeMarker::run_aldex(physeq, group = "enrich",
#                                           norm = "CPM",
#                                           taxa_rank = "none",
#                                           p_adjust = "fdr")
#   
#   mm_aldex_table <- data.frame(mm_aldex@marker_table)
#   return(mm_aldex_table)
# }
# 
# 
# 
# ##Biofilm
# aldex_biof <- aldex.m(physeq_16S_nit)
# ancombc_biof <- ancombc.m(physeq_16S_nit)
# lefse_biof <- lefse.m(physeq_16S_nit)
# 
# ancombc_biof$feature %in% lefse_biof$feature 
# ancombc_biof$feature %in% aldex_biof$feature 
# aldex_biof$feature %in% lefse_biof$feature 
# 
# 
# # Renommer les colonnes pour indiquer la méthode
# #aldex_biof <- rename(aldex_biof, ef_aldex = ef_aldex, pvalue_aldex = pvalue, padj_aldex = padj)
# # Renommer les colonnes pour indiquer la méthode
# ancombc_biof <- rename(ancombc_biof, ef_ancombc = ef_W, pvalue_ancombc = pvalue, padj_ancombc = padj)
# lefse_biof <- rename(lefse_biof, ef_lefse = ef_lda, pvalue_lefse = pvalue, padj_lefse = padj)
# 
# #fusionner la data_frame
# diff_abund_biof_table <- ancombc_biof %>%
#   inner_join(lefse_biof, by = c("feature", "enrich_group"))
# 
# # Fusionner les données d'abondance et les données taxonomiques
# core_biof_marker <- psmelt(prune_taxa(diff_abund_biof_table$feature, physeq_16S_nit_p))
# 
# # Calculer l'abondance moyenne pour chaque ASV dans chaque groupe
# core_biof_marker %>%
#   group_by(OTU) %>%
#   summarise(mean_abundance = mean(Abundance)) %>%
#   ungroup()%>%
#   arrange(desc(mean_abundance))
# 
# # Créer un vecteur d'ordre pour les ASV
# core_biof_marker_order <- c("ASV1446","ASV683", "ASV1805", "ASV680","ASV1958")
# core_biof_marker_order <- c("ASV1446","ASV1062","ASV683", "ASV1805", "ASV680")
# # Convertir la colonne OTU en un facteur ordonné
# core_biof_marker$OTU <- factor(core_biof_marker$OTU, levels = rev(core_biof_marker_order))
# 
# 
# # Créer le graphique
# ggplot(data = core_biof_marker, aes(y = OTU, x = Abundance, fill = enrich)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("SW" = "#33CC33", "ENR" = "#99FF66")) +
#   labs(title = "",y = "",x = "") +
#   theme_classic()
# 
# 
# ##Water
# aldex_water <- aldex.m(core_water)
# ancombc_water <- ancombc.m(core_water)
# lefse_water <- lefse.m(core_water)
# 
# ancombc_water$feature %in% lefse_water$feature 
# ancombc_water$feature %in% aldex_water$feature 
# aldex_water$feature %in% lefse_water$feature 
# 
# # Renommer les colonnes pour indiquer la méthode
# #aldex_water <- rename(aldex_water, ef_aldex = ef_aldex, pvalue_aldex = pvalue, padj_aldex = padj)
# # Renommer les colonnes pour indiquer la méthode
# ancombc_water <- rename(ancombc_water, ef_ancombc = ef_W, pvalue_ancombc = pvalue, padj_ancombc = padj)
# #lefse_water <- rename(lefse_water, ef_lefse = ef_lda, pvalue_lefse = pvalue, padj_lefse = padj)
# 
# 
# #fusionner la data_frame
# diff_abund_water_table <- ancombc_water # %>%
# # inner_join(lefse_biof, by = c("feature", "enrich_group"))
# 
# # Fusionner les données d'abondance et les données taxonomiques
# core_water_marker <- psmelt(prune_taxa(diff_abund_water_table$feature, physeq_16S_nit_p))
# 
# # Calculer l'abondance moyenne pour chaque ASV dans chaque groupe
# core_water_marker %>%
#   group_by(OTU) %>%
#   summarise(mean_abundance = mean(Abundance)) %>%
#   ungroup()
# 
# # Créer un vecteur d'ordre pour les ASV
# core_water_marker_order <- c("ASV853","ASV1742", "ASV849", "ASV683")
# 
# # Convertir la colonne OTU en un facteur ordonné
# core_water_marker$OTU <- factor(core_water_marker$OTU, levels = rev(core_water_marker_order))
# 
# # Créer le graphique
# ggplot(data = core_water_marker, aes(y = OTU, x = Abundance, fill = enrich)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("SW" = "#0099CC", "ENR" = "#33CCFF")) +
#   labs(title = "",y = "",x = "") +
#   theme_classic()