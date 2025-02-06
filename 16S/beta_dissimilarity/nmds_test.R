####NMDS----


##Loadings----
library(phyloseq)
library(zCompositions)
library(ggplot2)
library(vegan)




#transformation des données----
# Charger les données phyloseq

physeq <- physeq_16S_nit

#Enelever les colonnes avec des 0 
physeq@otu_table <- physeq@otu_table[,colSums(physeq@otu_table) != 0]

# we first replace the zeros using
# the Count Zero Multiplicative approach
tmp <- zCompositions::cmultRepl(physeq@otu_table,
                                method = "CZM",
                                label = 0)

# generate the centered log-ratio transformed. ASVs are in rows!!!!!
physeq_clr_asv <- apply(tmp, 1, function(x) log(x) - mean(log(x)))

#create a new phyloseq object with CLR tranformed counts
physeq_clr <- physeq
otu_table(physeq_clr) <- otu_table(t(physeq_clr_asv),
                                   taxa_are_rows = FALSE)










#Analyse NMDS----
physeq_clr_pre_nmds <- physeq_clr
#physeq_clr_pre_nmds@sam_data <- sample_data(meta_RDA)

#distance matrix calculation
physeq_clr_dist <- phyloseq::distance(physeq_clr_pre_nmds, method = "euclidean")


#NMDS plot on Aitchison distance
physeq_clr_nmds <- vegan::metaMDS(physeq_clr_dist, k=2, trymax=100) #Aitchison distance


vegan::stressplot(physeq_clr_nmds)

nmds_coord <- data.frame(physeq_clr_nmds$points)

#Data frame for hull
hull <- data.frame("Axis.1" = nmds_coord[,1],
                   "Axis.2" = nmds_coord[,2],
                   "sample" = as.data.frame(sample_data(physeq_clr_pre_nmds@sam_data)))

hull$sample.date <- substring(hull$sample.date, 6)


biofilm_SW_NMDS <- hull[hull$sample.type_enrich  == "biofilm_SW", ][chull(hull[hull$sample.type_enrich =="biofilm_SW", c("Axis.1", "Axis.2")]), ]  # hull values for biofilm_SW
biofilm_ENR_NMDS <- hull[hull$sample.type_enrich  == "biofilm_ENR", ][chull(hull[hull$sample.type_enrich =="biofilm_ENR", c("Axis.1", "Axis.2")]), ]  # hull values for biofilm_SW
water_SW_NMDS <- hull[hull$sample.type_enrich  == "water_SW", ][chull(hull[hull$sample.type_enrich =="water_SW", c("Axis.1", "Axis.2")]), ]  # hull values for water_SW
water_ENR_NMDS <- hull[hull$sample.type_enrich  == "water_ENR", ][chull(hull[hull$sample.type_enrich =="water_ENR", c("Axis.1", "Axis.2")]), ]  # hull values for water_SW


#Vector of color for hulls
hull_col <- c("#33CC33", "#99FF66", "#0099CC", "#33CCFF")
names(hull_col) <- c("biofilm_SW" , "biofilm_ENR", "water_SW", "water_ENR")



hull_data <- hull %>%
  dplyr::group_by(sample.type_enrich) %>%
  dplyr::slice(chull(Axis.1,Axis.2)) %>%
  dplyr::mutate(color = hull_col[sample.type_enrich])



###Graphic representaion ----
NMDS_plot <- ggplot(hull,aes(x = Axis.1, y = Axis.2)) +
  geom_hline(yintercept = 0, colour = "lightgrey", linetype = 2) + 
  geom_vline(xintercept = 0, colour = "lightgrey", linetype = 2) +
  scale_fill_manual(values = c("biofilm_SW" = "#33CC33", "biofilm_ENR" = "#99FF66", "water_SW" = "#0099CC", "water_ENR" = "#33CCFF")) +
  geom_point(data = hull,size = 4,aes(shape = sample.type_enrich, fill = sample.type_enrich)) +
  scale_shape_manual(values = c("biofilm_SW" = 21, "biofilm_ENR" = 22, "water_SW" = 23, "water_ENR" = 24))+
  scale_color_manual(values =  c("#33CC33", "#99FF66", "#0099CC", "#33CCFF")) +
  geom_polygon(data = hull_data,
               aes(group = sample.type_enrich,
                   fill = sample.type_enrich),
               alpha = 0.3) + # add the convex hulls)+
  geom_text(data = hull_data, x = -50, y = 50, label = paste("Stress =", round(physeq_clr_nmds$stress, 2)),colour = "Black",size = 4)  +
  ggrepel::geom_text_repel(data = hull,
                           aes(x=Axis.1, y=Axis.2, label = as.character(sample.date)), size= 2.5)+
  #stat_ellipse(aes(group = sample.type_enrich, fill = sample.type_enrich), geom = "polygon", alpha = 0.2) +  # Utiliser stat_ellipse avec aes(group = ...) et d'autres paramètres ++ suppose une distribution normale alors c'est mieux enveloppe convexe 
  xlab(paste("MDS1")) +
  ylab(paste("MDS2")) +
  guides(shape = guide_legend(title = NULL, color = "black"),
         fill = guide_legend(title = NULL))+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_equal() +
  theme(axis.title.x = element_text(size=14), # remove x-axis labels
        axis.title.y = element_text(size=14), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.position = "bottom")



#Dispersion

#de chaque groupe 
betadisper <- vegan::betadisper(physeq_clr_dist, hull$sample.type_enrich)

boxplot(betadisper)

permutest(betadisper)

# Calcul de la dispersion
betadisper <- vegan::betadisper(physeq_clr_dist, hull$sample.type_enrich)

# Visualisation de la dispersion
boxplot(betadisper$distances ~ hull$sample.type_enrich, main="Dispersion des échantillons", ylab="Distance au centroïde")

# Test de la dispersion
anova(betadisper)



#calcul dispersion  des triplicats 
#colonne de id de trilplicats 
df_tripli_disp <- subset(hull, sample.type_n== "biofilm")
                         
df_tripli_disp$tripli_ID <- substr(df_tripli_disp$sample.Sample_ID, start=1, stop=10)

#calcul de distance 
betadisper <- vegan::betadisper(physeq_clr_dist, df_tripli_disp$tripli_ID)


boxplot(betadisper)



colnames(physeq_clr_dist)



#PERMANOVA----

#PERMANOVA générale 
metadata <- data.frame(sample_data(physeq_clr))
results_permanova <- vegan::adonis2(physeq_clr_dist ~ type_enrich,
                                    data = metadata,
                                    perm = 1000)
results_permanova


#PERMANOVA deux a deux 

# Extraction des métadonnées
metadata <- data.frame(sample_data(physeq_clr))

# Obtenir les niveaux uniques de type_enrich
groups <- unique(metadata$type_enrich)

# Initialiser un data.frame pour stocker les résultats
results_permanova_df <- data.frame(
  Group1 = character(0),
  Group2 = character(0),
  R2 = numeric(0),
  p_value = numeric(0),
  stringsAsFactors = FALSE
)

# Boucle sur chaque combinaison de deux groupes et effectuer PERMANOVA
combinations <- combn(groups, 2, simplify = FALSE)
for (combo in combinations) {
  # Filtrer les données pour la combinaison actuelle de groupes
  subset_physeq <- subset_samples(physeq_clr, type_enrich %in% combo)
  subset_physeq_clr_dist <- distance(subset_physeq, method = "euclidean")
  
  # Exécuter PERMANOVA pour la combinaison actuelle de groupes
  result <- vegan::adonis2(subset_physeq_clr_dist ~ type_enrich, data = data.frame(sample_data(subset_physeq)), perm = 1000)
  print(combo)
  print(result)
  # Extraire les résultats pertinents
  R2 <- result$R2[1] # Prendre la première valeur de R2
  p_value <- result$`Pr(>F)`[1]
  
  # Ajouter les résultats au data.frame
  results_permanova_df <- rbind(results_permanova_df, data.frame(
    Group1 = combo[1],
    Group2 = combo[2],
    R2 = R2,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
  
  
}
# Appliquer la correction des p-valeurs
# Vous pouvez choisir la méthode qui vous convient le mieux
results_permanova_df$adjusted_p_value <- p.adjust(results_permanova_df$p_value, method = "BH")

# Afficher le data.frame des résultats
results_permanova_df




subset(physeq_16_nit_p, type_enrich =="biofilm_SW")


#SD triplicates----

#calcul de la moyenne des distance des ech face au centroide 

hull %>%
  group_by(sample.type_enrich) %>%
  summarize(mean_Axis1 = mean(Axis.1), mean_Axis2 = mean(Axis.2))


hull %>%
  group_by(sample.type_enrich) %>%
  mutate(distance_to_centroid = sqrt((Axis.1 - mean(Axis.1))^2 + (Axis.2 - mean(Axis.2))^2)) %>%
  summarise(mean_distance_to_centroid = mean(distance_to_centroid))



#pour les triplicats biofilms 

subset(hull, sample.type_enrich =="biofilm_SW") %>%
  group_by(sample.date) %>%
  summarize(mean_Axis1 = mean(Axis.1), mean_Axis2 = mean(Axis.2))



subset(hull, sample.type_enrich =="biofilm_SW") %>%
  group_by(sample.date) %>%
  summarize(sd_Axis1 = sd(Axis.1), sd_Axis2 = sd(Axis.2))



# Stocker le résultat du summarize dans un objet
sd_values <- subset(hull, sample.type_enrich =="biofilm_SW") %>%
  group_by(sample.date) %>%
  summarize(sd_Axis1 = mean(Axis.1), sd_Axis2 = mean(Axis.2))

# Calculer la moyenne des écarts-types
mean_sd_values <- rowMeans(select(sd_values, sd_Axis1, sd_Axis2), na.rm = TRUE)


# Fonction pour calculer la moyenne et l'écart-type
calculate.mean.sd <- function(hull_data, type_enrich_value) {
  filtered_data <- subset(hull_data, sample.type_enrich == type_enrich_value)
  
  # Calculer la moyenne et l'écart-type pour chaque date
  stats_values <- filtered_data %>%
    group_by(sample.date) %>%
    summarize(
      mean_Axis1 = mean(Axis.1, na.rm = TRUE),
      mean_Axis2 = mean(Axis.2, na.rm = TRUE),
      sd_Axis1 = sd(Axis.1, na.rm = TRUE),
      sd_Axis2 = sd(Axis.2, na.rm = TRUE)
    )
  
  # Ajouter une colonne pour le type_enrich pour référencer plus tard
  stats_values$type_enrich <- type_enrich_value
  
  return(stats_values)
}

# Utilisation de la fonction pour biofilm_SW et biofilm_ENR
stats_biofilm_SW <- calculate.mean.sd(hull, "biofilm_SW")
stats_biofilm_ENR <- calculate.mean.sd(hull, "biofilm_ENR")

# Fusionner les deux data.frames en un seul
final_stats_values <- full_join(stats_biofilm_SW, stats_biofilm_ENR, by = "sample.date", suffix = c("_SW", "_ENR"))

# Afficher le résultat final
print(final_stats_values)

#PCA----
#######PCA 

pca_results <- dudi.pca(as.matrix(otu_table(physeq_clr)), scannf = FALSE, nf = 2)
plot(pca_results)
