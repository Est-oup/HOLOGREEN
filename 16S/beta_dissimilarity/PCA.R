library(phyloseq)
library(ggplot2)
library(microViz)


physeq_clr %>% ord_calc() %>%
  ord_plot(color = "ibd", shape = "DiseaseState", size = 2) +
  scale_colour_brewer(palette = "Dark2")


#vegan::
  
  
# Effectuer la PCA
pca_result <- prcomp(otu_table(physeq_clr), scale. = TRUE)

# Accéder aux valeurs propres (eigenvalues)
eigenvalues <- pca_result$sdev^2

# Calculer le pourcentage de variance expliquée par chaque composante principale
variance_explained <- eigenvalues / sum(eigenvalues) * 100

# Créer un dataframe pour les valeurs propres et les pourcentages de variance expliquée
eigen_df <- data.frame(Component = 1:length(variance_explained), VarianceExplained = variance_explained)

# Visualiser les pourcentages de variance expliquée sous forme de graphique à barres
ggplot(eigen_df, aes(x = Component, y = VarianceExplained)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  labs(x = "Composante Principale", y = "Pourcentage de Variance Expliquée",
       title = "Pourcentage de Variance Expliquée par Composante Principale") +
  theme_minimal()

# Extraire les coordonnées des échantillons pour les deux premières composantes principales
pca_data <- as.data.frame(pca_result$x[, 1:2])

# Associer correctement la variable "type_enrich" aux échantillons dans pca_data
pca_data$type_enrich <- physeq_clr@sam_data$type_enrich[match(rownames(pca_data), rownames(physeq_clr@sam_data))]

# Ajouter les noms d'échantillons à la DataFrame
pca_data$Sample <- rownames(pca_data)

# Visualiser les résultats
ggplot(pca_data, aes(x = PC1, y = PC2, color = type_enrich)) +
  geom_point(size = 3) +
  labs(x = "Composante Principale 1", y = "Composante Principale 2",
       title = "PCA avec Coloration par type_enrich") +
  theme_minimal() +
  scale_color_manual(values = c("biofilm_SW" = "#33CC33", "biofilm_ENR" = "#99FF66", "filter_SW" = "#0099CC", "filter_ENR" = "#33CCFF"))



# a partir de résulttas de RDA 

# Extraire les coordonnées des échantillons à partir des résultats RDA
pca_data <- as.data.frame(scores.rda(rda_result$CCA)[, 1:2])

# Ajouter les noms d'échantillons à la DataFrame
pca_data$Sample <- rownames(pca_data)

# Visualiser les résultats
ggplot(pca_data, aes(x = Dim1, y = Dim2)) +
  geom_point(size = 3) +
  labs(x = "Composante 1", y = "Composante 2",
       title = "PCA basée sur les résultats de RDA") +
  theme_minimal()
rda_result$
  
  
  
  #PCOA
  

