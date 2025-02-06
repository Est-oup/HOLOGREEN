#differential abundance test



diff.abundance.plot <- function(physeq_core, sample_type){
  #run de la fonction
  mm_aldex <- microbiomeMarker::run_aldex(physeq_core, group = "enrich",
                                          norm = "CPM",
                                          taxa_rank = "none",
                                          p_adjust = "fdr")
  
  mm_aldex_table <- data.frame(mm_aldex@marker_table)
  
  # Ordonner le dataframe selon la colonne ef_aldex
  mm_aldex_table <- merge(mm_aldex_table, cbind("feature" = noquote(row.names(physeq_core@tax_table)),as.data.frame(physeq_core@tax_table)), by= "feature")
  mm_aldex_table <- arrange(mm_aldex_table, ef_aldex)
  
  
  #etiquette des données 
  names_diff <- paste(mm_aldex_table$Genus, mm_aldex_table$feature, sep = " ")
  
  # Choix des couleurs en fonction du type d'échantillon
  color_values <- if(sample_type == "biofilm") {
    c("SW" = "#33CC33", "ENR" = "#99FF66")
  } else if(sample_type == "water") {
    c("SW" = "#0099CC", "ENR" = "#33CCFF")
  } else {
    stop("Invalid sample_type")
  }
  
  #Plot
  aldex2_plot <- ggplot(mm_aldex_table , aes(x = ef_aldex, y = feature, fill = enrich_group)) +
    geom_bar(stat = "identity") +
    xlab("ef_aldex") + ylab("Feature") +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_discrete(limits=arrange(mm_aldex_table, ef_aldex)$feature, labels= names_diff)+
    scale_fill_manual(values = color_values)+
    theme_classic()
  
  print(aldex2_plot)
}




diff.abundance.tab <- function(physeq_core, sample_type){
  #run de la fonction
  mm_aldex <- microbiomeMarker::run_aldex(physeq_core, group = "enrich",
                                          norm = "CPM",
                                          taxa_rank = "none",
                                          p_adjust = "fdr")
  
  mm_aldex_table <- data.frame(mm_aldex@marker_table)
  
  # Ordonner le dataframe selon la colonne ef_aldex
  mm_aldex_table <- merge(mm_aldex_table, cbind("feature" = noquote(row.names(physeq_core@tax_table)),as.data.frame(physeq_core@tax_table)), by= "feature")
  mm_aldex_table <- arrange(mm_aldex_table, ef_aldex)
  
  print(mm_aldex_table)
}




list.core <- function(physeq, type_val){
  physeq_sub <- eval(bquote(subset_samples(physeq, type == .(type_val))))
  
  phy_SW <- subset_samples(physeq_sub, enrich == "SW")
  phy_ENR <- subset_samples(physeq_sub, enrich == "ENR")
  
  SW <- names(colSums(otu_table(phy_SW))>0)[colSums(otu_table(phy_SW))>0]
  ENR <- names(colSums(otu_table(phy_ENR))>0)[colSums(otu_table(phy_ENR))>0]
  SW_ENR <- intersect(SW,ENR)
  
  # Créer une liste avec des noms
  list_core <- list("SW" = SW, "ENR" = ENR)
  
  # Définir les couleurs en fonction de type_val
  if (type_val == "biofilm") {
    colors <- c("#33CC33","#99FF66")
  } else if (type_val == "water") {
    colors <- c("#0099CC", "#33CCFF")
  } else {
    colors <- NULL  # ou une valeur par défaut
  }
  
  plot <- ggvenn(list_core, fill_color = colors, 
                 text_size = 6, 
                 stroke_alpha = 0,  # Enlever les lignes de contour
                 stroke_size = 0) +  # Taille du trait à 0
    theme(legend.position = "none")  # Enlever la légende
  return(plot)
}



heat.dynamics <- function(physeq_obj, aldex_tab) {
  # Obtenir les données de physeq_obj
  data <- psmelt(physeq_obj) %>% group_by(enrich, OTU) %>% summarise(mean_Abundance = mean(Abundance))
  
  # Fusionner avec aldex_tab
  colnames(aldex_tab)[colnames(aldex_tab) == "feature"] <- "OTU"
  data <- merge(data, aldex_tab, by = "OTU")
  
  # Obtenir l'ordre des OTU de aldex_tab
  ordered_otu <- aldex_tab$OTU
  
  # Réorganiser les données de 'data' pour qu'elles correspondent à l'ordre des OTU dans 'aldex_tab'
  data <- data %>%
    mutate(OTU = factor(OTU, levels = ordered_otu)) %>%
    arrange(desc(OTU))
  
  # Créer la heatmap
  p_dynamic <- ggplot(data = data, aes(x = enrich, y = OTU, fill = mean_Abundance)) +
    geom_tile(aes(fill = mean_Abundance), colour = "white") +
    scale_fill_gradient(low = "grey", high = "#FF0000") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y=element_blank()) +
    theme_minimal()
  
  return(p_dynamic)
}



bubble.dynamics <- function(physeq_obj, aldex_tab, type_val) {
  # Obtenir les données de physeq_obj
  data <- psmelt(physeq_obj) %>% group_by(enrich, OTU) %>% summarise(mean_Abundance = mean(Abundance))
  
  # Fusionner avec aldex_tab
  colnames(aldex_tab)[colnames(aldex_tab) == "feature"] <- "OTU"
  data <- merge(data, aldex_tab, by = "OTU")
  
  # Obtenir l'ordre des OTU de aldex_tab
  ordered_otu <- aldex_tab$OTU
  
  # Réorganiser les données de 'data' pour qu'elles correspondent à l'ordre des OTU dans 'aldex_tab'
  data <- data %>%
    mutate(OTU = factor(OTU, levels = ordered_otu)) %>%
    arrange(desc(OTU))
  
  # Définir les couleurs en fonction de type_val
  if (type_val == "biofilm") {
    colors <- c("SW" = "#33CC33", "ENR" = "#99FF66")
  } else if (type_val == "water") {
    colors <- c("SW" = "#0099CC", "ENR" = "#33CCFF")
  } else {
    colors <- NULL  # ou une valeur par défaut
  }
  
  # Créer le bubble chart
  p_dynamic <- ggplot(data = data, aes(x = enrich, y = OTU)) +
    geom_point(aes(size = mean_Abundance, color = enrich)) +
    geom_text(aes(label = paste0(round(mean_Abundance, 1), "")), vjust = -0.5) +  # Ajout des pourcentages
    scale_size(range = c(0, 15)) +  # Augmenter la plage de tailles des bulles
    scale_color_manual(values = colors) +
    scale_x_discrete(limits = c("SW", "ENR")) +  # Ordre des niveaux sur l'axe des x
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y=element_blank()) +
    theme_void() +
    guides(size = guide_legend(override.aes = list(stroke = 0)))  # Empiler les cercles de la légende
  
  return(p_dynamic)
}
