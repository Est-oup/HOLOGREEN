########Alpha diversity


#loadings
library(ggplot2)
library(microbiome)
library(gridExtra)
library(car)
library(coin)
library(knitr)
library(kableExtra)




source()

##Calculs indices-----
physeq = physeq_16S_nit

####Statistical mean comparison between alpha indices

#Test richesse spécifique 
alpha_indices_nit <- data.frame(cbind(microbiome::alpha(physeq, index = "all"), physeq@sam_data[,c(1,2,4,5,6)] ))


####Partie plot  

data_alpha_indices <-alpha_indices_nit

# Définir l'ordre des niveaux de la variable "type_enrich"
data_alpha_indices$type_enrich <- factor(data_alpha_indices$type_enrich, levels = c("biofilm_SW", "biofilm_ENR", "water_SW", "water_ENR"))


graphA <- ggplot(data=data_alpha_indices, aes(x=type_enrich, y= observed, fill=type_enrich)) + 
  geom_boxplot(width=0.7, show.legend = F, outlier.size = 0.25, color = c("#33CC33", "#99FF66", "#0099CC", "#33CCFF"), alpha= 0.4) +
  geom_point(show.legend = F, size= 1)+
  # geom_text(aes(x=1, y=530, label = "a"), size= 5, hjust= -0.6,vjust=2, color= "#33CC33", family= "Times")+
  # geom_text(aes(x=2, y= 360, label = "a"), size= 5, hjust= -0.6,vjust=1, color= "#99FF66", family= "Times")+
  # geom_text(aes(x=3, y= 345, label = "a"), size= 5, hjust= -0.6,vjust=2, color= "#0099CC", family= "Times")+
  # geom_text(aes(x=4, y= 430, label = "a"), size= 5, hjust= -0.6,vjust=2, color= "#33CCFF", family= "Times")+
  labs(title="Specific richness", x="", y="")+
  scale_fill_manual(values=c("biofilm_SW" = "#33CC33", "biofilm_ENR" = "#99FF66", "water_SW" = "#0099CC", "water_ENR" = "#33CCFF"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, hjust = 1))
  
graphB <- ggplot(data=data_alpha_indices, aes(x=type_enrich, y= diversity_shannon, fill=type_enrich)) + 
  geom_boxplot(width=0.7, show.legend = F, outlier.size = 0.25, color = c("#33CC33",  "#99FF66",  "#0099CC", "#33CCFF"), alpha= 0.4) +
  geom_point(show.legend = F, size= 1)+
  # geom_text(aes(x=1, y= 3.98 , label = "a"), size= 5, hjust= -0.6,vjust=0, color= "#33CC33", family= "Times")+
  # geom_text(aes(x=2, y= 3.98 , label = "b"), size= 5, hjust= -0.6,vjust=0, color= "#99FF66", family= "Times")+
  # geom_text(aes(x=3, y= 2.9 , label = "c"), size= 5, hjust= -0.6,vjust=0, color= "#0099CC", family= "Times")+
  # geom_text(aes(x=4, y= 4.05 , label = "ab"), size= 5, hjust= 0,vjust=0, color= "#33CCFF", family= "Times")+
  labs(title="Shannon diversity", x="", y="")+
  scale_fill_manual(values=c("biofilm_SW" = "#33CC33", "biofilm_ENR" = "#99FF66", "water_SW" = "#0099CC", "water_ENR" = "#33CCFF"))+
    theme_minimal()+
  theme(axis.text.x = element_text(angle=90, hjust = 1))



  # Définir l'ordre des niveaux de la variable "type_enrich"
  data_alpha_indices$type_enrich <- factor(data_alpha_indices$type_enrich, levels = c("biofilm_SW", "biofilm_ENR", "water_SW", "water_ENR"))

  
  
  #plot 
  grid.arrange(graphA,
               graphB,
               nrow = 1)
   
  

####Test statistiques 2 à 2 

#test de la normalité 

data_alpha_indices %>%dplyr::select(observed, diversity_shannon)%>% indices_normality(nrow=1, ncol=2)

#Si c'est normal
#Shannon c'est des données normales : 
# Check homogeneity of variance between groups
# (avoid bias in ANOVA result & keep the power of the test)
# H0= equality of variances in the different populations
stats::bartlett.test(diversity_shannon ~ type_enrich, data_alpha_indices)

aov_observed <- stats::aov(diversity_shannon ~ type_enrich, data_alpha_indices)
summary(aov_observed)

signif_pairgroups <- stats::TukeyHSD(aov_observed, method = "bh")
signif_pairgroups
  

#nos données ne sont pas normales alors :  

#Utilisation de test non paramétriques pour + de >3 groupes  
res_test_kruscal_obs <- kruskal.test(observed ~ type_enrich, data=data_alpha_indices)
res_test_kruscal_sha <- kruskal.test(diversity_shannon ~ type_enrich, data=data_alpha_indices)
res_test_kruscal_simp <- kruskal.test(diversity_gini_simpson ~ type_enrich, data=data_alpha_indices)
res_test_kruscal_piel <- kruskal.test(evenness_pielou ~ type_enrich, data=data_alpha_indices)

res_test_kruscal_obs 
res_test_kruscal_sha 
res_test_kruscal_simp
res_test_kruscal_piel


#pour la parwaise comparaison : POST HOC DUNN TEST  
res_test_post_hoc_obs <- FSA::dunnTest(observed~ type_enrich, data= data_alpha_indices, method="bh")
res_test_post_hoc_sha  <- FSA::dunnTest(diversity_shannon~ type_enrich, data= data_alpha_indices, method="bh")
res_test_post_hoc_simp  <- FSA::dunnTest(diversity_gini_simpson~ type_enrich, data= data_alpha_indices, method="bh")
res_test_post_hoc_piel <- FSA::dunnTest(evenness_pielou~ type_enrich, data= data_alpha_indices, method="bh")

res_test_post_hoc_obs
res_test_post_hoc_sha
res_test_post_hoc_simp
res_test_post_hoc_piel








#test of sd comparison 
# Supposons que "data_alpha_indices" est votre data frame contenant les colonnes "observed", "diversity_shannon" et "type_enrich".

# Charger le package stats si ce n'est pas déjà fait
library(stats)

# Fonction pour effectuer le test de Fligner-Killeen entre deux variables données
perform_fligner_test <- function(variable_name, group1, group2) {
  formula <- as.formula(paste(variable_name, "~ type_enrich"))
  test_result <- fligner.test(formula, data = subset(data_alpha_indices, type_enrich %in% c(group1, group2)))
  return(test_result)
}

# Fonction pour afficher les résultats des tests de Fligner-Killeen
print_fligner_results <- function(results_df) {
  cat("Comparison\tStatistic\tdf1\tdf2\tp_value\tp_adjusted\n")
  for (i in 1:nrow(results_df)) {
    cat(paste(results_df$Comparison[i], "\t",
              results_df$Statistic[i], "\t",
              results_df$df1[i], "\t",
              results_df$df2[i], "\t",
              results_df$p_value[i], "\t",
              results_df$p_adjusted[i], "\n"))
  }
}

# Données : data_alpha_indices$observed et data_alpha_indices$diversity_shannon sont les variables à comparer
# data_alpha_indices$type_enrich est la variable catégorielle des groupes

# Obtenir les niveaux uniques de "type_enrich"
unique_groups <- unique(data_alpha_indices$type_enrich)

# Créer une liste pour stocker les résultats du test de Fligner-Killeen pour toutes les comparaisons par paires
fligner_test_results <- list()

# Effectuer le test de Fligner-Killeen pour toutes les comparaisons par paires et pour chaque variable
for (variable_name in c("observed", "diversity_shannon")) {
  for (i in 1:(length(unique_groups)-1)) {
    for (j in (i + 1):length(unique_groups)) {
      group1 <- unique_groups[i]
      group2 <- unique_groups[j]
      test_result <- perform_fligner_test(variable_name, group1, group2)
      comparison_name <- paste(variable_name, "-", group1, "-", group2)
      fligner_test_results[[comparison_name]] <- test_result
    }
  }
}

# Créer le DataFrame avec les résultats des tests de Fligner-Killeen
fligner_results_df <- data.frame(Comparison = character(),
                                 Variable = character(),
                                 Statistic = numeric(),
                                 df1 = numeric(),
                                 df2 = numeric(),
                                 p_value = numeric(),
                                 stringsAsFactors = FALSE)

# Créer des vecteurs pour stocker les résultats bruts des p-valeurs
raw_p_values <- c()

# Remplir le DataFrame avec les résultats des tests de Fligner-Killeen et les informations de comparaison
for (comparison in names(fligner_test_results)) {
  test_result <- fligner_test_results[[comparison]]
  comparison_parts <- strsplit(comparison, "-")
  variable_name <- comparison_parts[[1]][1]
  group1 <- comparison_parts[[1]][2]
  group2 <- comparison_parts[[1]][3]
  fligner_results_df <- rbind(fligner_results_df,
                              data.frame(Comparison = comparison,
                                         Variable = variable_name,
                                         Statistic = round(test_result$statistic, 1),
                                         df1 = test_result$parameter[1],
                                         df2 = test_result$parameter[2],
                                         p_value = round(test_result$p.value, 5)))
  raw_p_values <- c(raw_p_values, test_result$p.value)
}

# Calculer les p-valeurs ajustées avec la méthode de Benjamini-Hochberg (BH)
fligner_results_df$p_adjusted <- round(p.adjust(raw_p_values, method = "BH"),5)

# Afficher les résultats des tests de Fligner-Killeen
print_fligner_results(fligner_results_df)







#Stat table final----

stat_final <- fligner_results_df[, c(1,7)]
rownames(stat_final) <- fligner_results_df[,1]
stat_final <- stat_final[-1]


stat_final <- cbind(stat_final, c(round(res_test_post_hoc_obs$res$P.adj, 5), round(res_test_post_hoc_sha$res$P.adj, 5)))

colnames(stat_final) <- c("fligner test", "Post Hoc test")

pdf("stat_final.pdf")



##Article writing , information about ASV----

#number of ASv 
ncol(physeq_16S_nit_p@otu_table)

#minimum of ASV in samples
min(alpha_indices_nit$observed)
max(alpha_indices_nit$observed)


#mean of observed value in type_enrich
alpha_indices_nit %>% group_by(type_enrich) %>% summarize(observed2 = mean(observed))
#median of observed value in type_enrich
alpha_indices_nit %>% group_by(type_enrich) %>% summarize(observed2 = median(observed))





#alpha indice along time 

alpha_dynamics <-alpha_indices_nit 
alpha_dynamics$date <- paste0(substr(alpha_dynamics$Sample_ID, 3,4), "-",substr(alpha_dynamics$Sample_ID, 5,6)) 




# Calcul de la moyenne et de l'écart type par groupe
alpha_dynamics <- alpha_dynamics %>%
  group_by(type_enrich, date) %>%
  summarize(
    mean_observed = mean(observed),
    sd_observed = sd(observed),
    mean_shannon = mean(diversity_shannon),
    sd_shannon = sd(diversity_shannon)
    
  )

#test de comparaison 
#Utilisation de test non paramétriques pour + de >3 groupes  
res_test_kruscal_dynamics_obs <- kruskal.test(sd_observed  ~ type_enrich, data=alpha_dynamics)

#pour la parwaise comparaison : POST HOC DUNN TEST  
res_test_post_dynamics_obs <- FSA::dunnTest(sd_observed~ type_enrich, data= alpha_dynamics, method="bh")



ggplot(alpha_dynamics) +
  geom_boxplot(aes(x = type_enrich, y = mean_shannon, group = type_enrich, color = type_enrich)) +
  geom_point(aes(x = type_enrich, y = mean_shannon, color = type_enrich)) +
  labs(x = "Date", y = "Moyenne Shannon Value") +
  scale_color_manual(values = c("biofilm_SW" = "#33CC33", "biofilm_ENR" = "#99FF66", "water_SW" = "#0099CC", "water_ENR" = "#33CCFF")) +
  theme_minimal()



#brouillon 

for(x in unique(physeq_16S_nit_p@sam_data[,"type_enrich"])){
  for (i in 1:nrow(physeq_16S_nit_p@otu_table)){
      print(sd(subset_samples(physeq_16S_nit_p, type_enrich==x)@otu_table[,i]))
      }
}

sd.ASV <- function(physeq){
  library(phyloseq)
  
  # Initialisation
  types_enrich <- unique(sample_data(physeq)$type_enrich)
  results <- matrix(0, nrow = n_otus, ncol = length(types_enrich))
  rownames(results) <- rownames(otu_table(physeq))
  colnames(results) <- types_enrich
  
  # Calcul des écarts-types
  for(x in types_enrich){
    physeq_sub <- subset_samples(physeq, type_enrich == x)
    for (i in 1: nrow(otu_table(physeq_sub))){
      results[i, x] <- sd(physeq_sub@otu_table[,i])
    }
  }
  
  # Affichage des résultats
  print(results)
  
}

sd.ASV(physeq_16S_nit_p)

sd.ASV <- function(physeq){
  
types_enrich <- unique(sample_data(physeq)$type_enrich)
results <- matrix(0, nrow = ncol(otu_table(test_physeq_sub)), ncol = length(types_enrich))

for(x in types_enrich){
  test_physeq_sub <- subset_samples(physeq, type_enrich == x)
  for (i in 1: ncol(otu_table(test_physeq_sub))){
   #test <- sd(test_physeq_sub@otu_table[,i])
   results[i, x] <- sd(test_physeq_sub@otu_table[,i])
  }
}

# Affichage des résultats
print(results)

}




sd.ASV <- function(physeq){
  
  # Initialisation
  types_enrich <- unique(sample_data(physeq)$type_enrich)
  n_otus <- ncol(otu_table(physeq))
  results <- matrix(0, nrow = n_otus, ncol = length(types_enrich))
  rownames(results) <- taxa_names(physeq) # Utilisez taxa_names pour obtenir les noms des OTUs
  colnames(results) <- types_enrich
  
  # Calcul des écarts-types
  for(x in types_enrich){
    physeq_sub <- subset_samples(physeq, type_enrich == x)
    for (i in 1: ncol(otu_table(physeq_sub))){
      otu_name <- taxa_names(physeq_sub)[i]
      results[which(rownames(results) == otu_name), which(colnames(results) == x)] <- sd(physeq_sub@otu_table[,i])
    }
  }
  


# Exemple d'utilisation :
sd.ASV(physeq_16S_nit_p)


sd.ASV <- function(physeq){
  for(x in types_enrich){
    physeq_sub <- subset_samples(physeq, type_enrich == x)
    
    subset_samples(physeq_sub, type_enrich ==  x)@otu_table %>%
      data.frame %>%
      summarise_all(sd)
    
  }  
  # Renvoyer la matrice des résultats
  return(results)
}


library(dplyr)
library(phyloseq)

sd.ASV <- function(physeq){
  
  # Initialisation
  types_enrich <- unique(sample_data(physeq)$type_enrich)
  results_list <- list()
  
  for(x in types_enrich){
    physeq_sub <- subset_samples(physeq, type_enrich == x)
    
    # Calculer l'écart-type pour chaque ASV
    temp_result <- subset_samples(physeq_sub, type_enrich ==  x)@otu_table %>%
      data.frame() %>%
      summarise_all(sd)
    
    # Ajouter le résultat à la liste
    results_list[[x]] <- temp_result
  }  
  
  # Renvoyer la matrice des résultats
  return(results_list)
}

# Exemple d'utilisation :
sd.ASV_results <- sd.ASV(physeq_16S_nit_p)


#Utilisation de test non paramétriques pour + de >3 groupes  
kruskal.test(sd.ASV_results ~ type_enrich, data=data_alpha_indices)



#pour la parwaise comparaison : POST HOC DUNN TEST  
FSA::dunnTest(sd.ASV_results~ type_enrich, data= data_alpha_indices, method="bh")






sd.ASV <- function(physeq){
  
  # Initialisation
  types_enrich <- unique(sample_data(physeq)$type_enrich)
  results_df <- data.frame()
  
  for(x in types_enrich){
    physeq_sub <- subset_samples(physeq, type_enrich == x)
    
    # Calculer l'écart-type pour chaque ASV
    temp_result <- otu_table(physeq_sub) %>%
      data.frame() %>%
      summarise_all(sd) %>%
      mutate(type_enrich = x) # Ajout d'une colonne avec le nom du groupe
    
    # Ajouter le résultat au dataframe
    results_df <- rbind(results_df, temp_result)
  }  
  # Transpose the dataframe
  results_df <- t(results_df)
  
  # Set the column names based on the 'type_enrich' column
  colnames(results_df) <- results_df["type_enrich",]
  
  # Remove the 'type_enrich' row
  results_df <- data.frame(results_df[-which(rownames(results_df) == "type_enrich"),])
  
  # Transformation du tableau
  long_format <- sd.ASV_results %>%
    rownames_to_column("ASV") %>%
    gather(key = "type_enrich", value = "value", -ASV)
  
  
  # Renvoyer le dataframe des résultats
  return(results_df)
}

# Exemple d'utilisation :
sd.ASV_results <- sd.ASV(physeq_16S_nit_p)

# Utilisation de test non paramétriques pour + de >3 groupes  
kruskal_test_result <- kruskal.test(sd.ASV_results ~ type_enrich, data=sd.ASV_results)
print(kruskal_test_result)

# Pour la comparaison par paires : POST HOC DUNN TEST  
dunn_test_result <- FSA::dunnTest(sd.ASV_results[,1] ~ type_enrich, data= sd.ASV_results, method="bh")
print(dunn_test_result)

colnames(sd.ASV_results)
res_test_kruscal_obs <- kruskal.test(observed ~ type_enrich, data=data_alpha_indices)









sd.ASV <- function(physeq){
  
  # Initialisation
  types_enrich <- unique(sample_data(physeq)$type_enrich)
  results_list <- list()
  
  for(x in types_enrich){
    physeq_sub <- subset_samples(physeq, type_enrich == x)
    
    # Calculer l'écart-type pour chaque ASV
    temp_result <- otu_table(physeq_sub) %>%
      data.frame() %>%
      summarise_all(sd)
    
    # Ajouter le résultat à la liste
    results_list[[x]] <- temp_result
  }
print(results_list)
}

sd.ASV_results <- sd.ASV(physeq_16S_nit_p)



#Supplementary table----

write.table(alpha_indices_nit, file="/16S_HOLOGREEN/alpha_diversity/alpha_indices/objets/alpha_indices_nit.txt", sep=";")
