#Script to evaluate the dynamics of the core microbiome 


#Library

library(dplyr)
library(phyloseq)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(Bio)




heat.map.phy.glom <- function(physeq, percent, env_var, env_var_val){
  
  #subset of the pyloseq object
  # Créer un masque logique pour filtrer les échantillons
  mask <- rownames(physeq@sam_data[physeq@sam_data$type_enrich %in% env_var_val,])
  # Sous-ensemble de l'objet phyloseq en utilisant le masque
  physeq1 <- prune_samples(mask, physeq)
  
  
  #Calculation to have the vision of genus abundance 
  glom <- tax_glom(physeq1, taxrank = 'Genus', NArm = FALSE, bad_empty = c(""))
  data_glom <- psmelt(glom)
  
  # Rename missing Genus values with "unknown"
  for (i in 27:32) {  # Boucle de la colonne Genus à la colonne Kingdom
    data_glom[, i] <- ifelse(is.na(data_glom[, i]), 
                             paste("unknown", data_glom[, (i-1)], sep = " "), 
                             data_glom[, i])
  }
  
  # Clean data
  data_glom$Genus <- as.character(data_glom$Genus)
  data_glom$Class <- as.character(data_glom$Class)
  data_glom <- data_glom[, c("OTU", "Sample", "Abundance", "type_enrich", "Class", "Genus")]
  data_glom <- data.frame(data_glom %>% dplyr::group_by(OTU, Sample, type_enrich,Class, Genus) %>% summarize(Abundance = sum(Abundance)))
  
  # Handle "Others" category based on abundance
  for (c in unique(data_glom$Class)) {
    tmp <- data_glom[data_glom$Class == c, ]
    if (max(tmp$Abundance) < percent) {
      tmp[tmp$Abundance < percent, "Genus"] <- paste0("Others <", percent, "%")
    }
    if (max(tmp$Abundance) >= percent) {
      tmp[tmp$Abundance < percent, "Genus"] <- paste0("Others <", percent, "% ", c)
    }
    # Fusionner les résultats dans le tableau d'origine
    data_glom[data_glom$Class == c, "Genus"] <- tmp$Genus
  }
  
  
  
  #Reorder the dataset 
  # Reorder the table based on the obtained order
  data_glom_ordered <- data_glom %>% 
    group_by(Class, Genus) %>% 
    summarize(Abundance = sum(Abundance)) %>% 
    ungroup()
  
  
  # Order by Class
  #extract first "Others", then apply the mean of the abundance by class, order and then take the names
  class_order <- (subset(data_glom_ordered, !grepl("^Other", Genus)) %>%
                    group_by(Class) %>%
                    summarize(Abundance = sum(Abundance)) %>%
                    ungroup() %>% data.frame() %>%
                    arrange(desc(Abundance)))$Class
  
  data_glom_ordered <- data_glom_ordered %>%
    arrange(factor(Class, levels = class_order))
  
  # Order by Genus within each Class
  #to have the Genera ordered in each classes and the other at the end systematically
  order_data_glom <- as.character()
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
  #ajout du Other X%
  order_data_glom <- c(order_data_glom, paste0("Others <", percent, "%"))
  
  
  #reordonner le tableau final
  data_glom <- arrange(data_glom, match(Genus, order_data_glom))
  data_glom$Genus <- factor(data_glom$Genus, levels = order_data_glom)
  
  
  # Prepare data for heatmaps
  data_glom_heat <- data.frame(data_glom[, c("Sample", "Genus", "Abundance")])
  data_glom_heat <- as.data.frame(dcast(data_glom_heat, Sample ~ Genus, value.var = "Abundance", fun.aggregate = mean, fill = 0))
  rownames(data_glom_heat) <- data_glom_heat[,"Sample"]
  data_glom_heat <- data_glom_heat[, -which(names(data_glom_heat) == "Sample")]
  data_glom_heat <- t(data_glom_heat)
  
  
  
  # Order the samples
  heat1 <-  rownames(physeq1@sam_data[(physeq1@sam_data[,env_var] == env_var_val[1]),env_var])
  heat2 <-  rownames(physeq1@sam_data[(physeq1@sam_data[,env_var] == env_var_val[2]),env_var])
  
  # Find the maximum value
  max_value <- max(data_glom_heat)

  #create a color pal
  my_palette <- c("#A6CEE3", "#FFFFCC", "#FF0000")
  
  #create the heatmap 
  
  
  
  #first heatmap on the SW condition
  heatmap1 <- ggplot(melt(data_glom_heat[,heat1]), aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = my_palette, na.value = "white",
                         limits = c(0, max_value), breaks = c(0,2.5,5,10,15,20,25,30,45,60,max_value),
                         trans = "log1p", guide = "legend") +
    labs(x = "", y = "", fill = "Abundance") +
    geom_hline(yintercept = seq(0.5, nrow(data_glom_heat)+0.5, by = 1)-1, color = "black", linetype = "dotted", linewidth = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #modification layout 
  if (identical(env_var_val,  c("biofilm_SW", "biofilm_ENR"))){
    #ajouter les lignes verticales 
    heatmap1 <- heatmap1 + geom_vline(xintercept =seq(3, ncol(data_glom_heat)/2, by = 3)[-length(seq(3, ncol(data_glom_heat)/2, by = 3))]+0.5, color = "white", linetype = "dashed", linewidth = 1.5)
    
    #définir les étiquetttes
    dates_names <- as.character()
    dates_names <- str_sub(colnames(data_glom_heat[, heat1]), start = 1, end = 6)[duplicated(str_sub(colnames(data_glom_heat[, heat1]), start = 1, end = 6)) == 0]
    dates_names <- paste(substr(dates_names, 3, 4), substr(dates_names, 5, 6), sep="-")
    
    dates_names2 <- rep(" ", length(dates_names) * 3)
    dates_names2[seq(1, length(dates_names2), by = 3)] <- dates_names
    dates_names2 <- c( " ", dates_names2, " ")
      
      #réarranger avec la date : 
      heatmap1 <- heatmap1 + scale_x_discrete(labels = dates_names2)
      
  }else{
    #ajouter les lignes verticales 
    heatmap1 <- heatmap1 + geom_vline(xintercept =seq(1, ncol(data_glom_heat)/2, by = 1)[-length(seq(1, ncol(data_glom_heat)/2, by = 1))]+0.5, color = "white", linetype = "dashed", linewidth = 1.5)
    
    dates_names <- str_sub(rownames(otu_table(physeq1)), start = 1, end=6)
    dates_names <- paste(substr(dates_names, 3, 4), substr(dates_names, 5, 6), sep="-")
    
    #réarranger : 
    heatmap1 <- heatmap1 + scale_x_discrete(labels = dates_names)
  }
  
  
  
  #Second heatmap on ENR condition
  heatmap2 <- ggplot(melt(data_glom_heat[,heat2]), aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = my_palette, na.value = "white",
                         limits = c(0, max_value), breaks = c(0,2.5,5,10,15,20,25,30,45,60,max_value),
                         trans = "log1p", guide = "legend") +
    labs(x = "",y = "", fill = "Abundance") +
    geom_hline(yintercept = seq(0.5, nrow(data_glom_heat)+0.5, by = 1)-1, color = "black", linetype = "dotted", linewidth = 0.7)+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  #modification layout 
  if (identical(env_var_val,  c("biofilm_SW", "biofilm_ENR"))){
    #ajouter les lignes verticales 
    heatmap2 <- heatmap2 + geom_vline(xintercept =seq(3, ncol(data_glom_heat)/2, by = 3)[-length(seq(3, ncol(data_glom_heat)/2, by = 3))]+0.5, color = "white", linetype = "dashed", linewidth = 1.5)
    heatmap2 <- heatmap2 + geom_vline(xintercept =c(6,12,15)+0.5, color = "red", linetype = "dashed", linewidth = 1.5)
    
    #définir les étiquetttes
    dates_names <- as.character()
    dates_names <- str_sub(colnames(data_glom_heat[, heat2]), start = 1, end = 6)[duplicated(str_sub(colnames(data_glom_heat[, heat2]), start = 1, end = 6)) == 0]
    dates_names <- paste(substr(dates_names, 3, 4), substr(dates_names, 5, 6), sep="-")
    dates_names2 <- rep(" ", length(dates_names) * 3)
    dates_names2[seq(1, length(dates_names2), by = 3)] <- dates_names
    dates_names2 <- c( " ", dates_names2, " ")
    
    #réarranger avec la date : 
    heatmap2 <- heatmap2 + scale_x_discrete(labels = dates_names2)
    
  }else{
    #ajouter les lignes verticales
    heatmap2 <- heatmap2 + geom_vline(xintercept =seq(1, ncol(data_glom_heat)/2, by = 1)[-length(seq(1, ncol(data_glom_heat)/2, by = 1))]+0.5, color = "white", linetype = "dashed", linewidth = 1.5)
    heatmap2 <- heatmap2 + geom_vline(xintercept =c(2,5,7)+0.5, color = "red", linetype = "dashed", linewidth = 1.5)
    
    
    
    dates_names <- str_sub(rownames(otu_table(physeq1)), start = 1, end=6)
    dates_names <- paste(substr(dates_names, 3, 4), substr(dates_names, 5, 6), sep="-")
    
    #réarranger : 
    heatmap2 <- heatmap2 + scale_x_discrete(labels = dates_names)
  }
  
  
  
  
  
  
  
  # Organiser les heatmaps dans une figure 2x2
  heatmap_arrange <- plot_grid(
    heatmap1 +  guides(fill = FALSE),
    heatmap2 + theme(axis.text.y = element_blank()),
    nrow = 1, ncol = 2,
    labels = c("", "", "", ""),
    align = "hv",
    rel_heights = c(1, 1),
    rel_widths = c(1, 1)
  )
  
  # Ajouter le titre à la figure
  heatmap_arrange <- ggdraw() +
    draw_plot(heatmap_arrange) +
    draw_label("Heatmaps", size = 14, fontface = "bold", x = 0.5, y = 1.02)
  
  # Afficher les heatmaps arrangés
  print(heatmap_arrange)
  
}

#Plot pour les algues :  
heatmap_a <- heat.map.phy.glom(physeq_16S_nit_p, 4, "type_enrich", c("biofilm_SW", "biofilm_ENR"))

#Plot pour les filtres :
heatmap_f <- heat.map.phy.glom(physeq_16S_nit_p, 4, "type_enrich", c("water_SW", "water_ENR"))


