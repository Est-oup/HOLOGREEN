##STUDY OF THE ASV variance 


#Loadings----




#
sd_water_SW_asv <- data.frame(apply((subset_samples(physeq_16S_nit_p, type_enrich == "water_SW")@otu_table), 2, sd))
colnames(sd_water_SW_asv) <- "water_SW"
sd_water_ENR_asv <- apply((subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR")@otu_table), 2, sd)
sd_biofilm_SW_asv <- apply((subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_SW")@otu_table), 2, sd)
sd_biofilm_ENR_asv <- apply((subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_ENR")@otu_table), 2, sd)


ggplot(data = sd_water_SW_asv)+ geom_boxplot(y="water_SW")




# Convertir les vecteurs d'écarts-types en dataframes
sd_water_SW_asv <- data.frame(value = apply((subset_samples(physeq_16S_nit_p, type_enrich == "water_SW")@otu_table), 2, sd), type = "water_SW")
sd_water_ENR_asv <- data.frame(value = apply((subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR")@otu_table), 2, sd), type = "water_ENR")
sd_biofilm_SW_asv <- data.frame(value = apply((subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_SW")@otu_table), 2, sd), type = "biofilm_SW")
sd_biofilm_ENR_asv <- data.frame(value = apply((subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_ENR")@otu_table), 2, sd), type = "biofilm_ENR")

# Combiner tous les dataframes en un seul dataframe en format long
all_sd_asv <- rbind(sd_water_SW_asv, sd_water_ENR_asv, sd_biofilm_SW_asv, sd_biofilm_ENR_asv)



#Calculer les mean des ASV de manièere globales 
mean_water_SW_asv <- data.frame(value = apply((subset_samples(physeq_16S_nit_p, type_enrich == "water_SW")@otu_table), 2, mean), type = "water_SW")
mean_water_ENR_asv <- data.frame(value = apply((subset_samples(physeq_16S_nit_p, type_enrich == "water_ENR")@otu_table), 2, mean), type = "water_ENR")
mean_biofilm_SW_asv <- data.frame(value = apply((subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_SW")@otu_table), 2, mean), type = "biofilm_SW")
mean_biofilm_ENR_asv <- data.frame(value = apply((subset_samples(physeq_16S_nit_p, type_enrich == "biofilm_ENR")@otu_table), 2, mean), type = "biofilm_ENR")

# Combiner tous les dataframes en un seul dataframe en format long
all_mean_asv <- rbind(mean_water_SW_asv, mean_water_ENR_asv, mean_biofilm_SW_asv, mean_biofilm_ENR_asv)



# Créer un boxplot avec ggplot2
ggplot(data = all_sd_asv, aes(x = type, y = value)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 0.0155)) +  # Ajustez ces valeurs en fonction de votre besoin
  ggtitle("Distribution of Standard Deviations for Different Sample Types") +
  xlab("Sample Type") +
  ylab("Standard Deviation")




sd_total <- data.frame(c(mean(sd_water_SW_asv$value),mean(sd_water_ENR_asv$value),mean(sd_biofilm_SW_asv$value),mean(sd_biofilm_ENR_asv$value)))
colnames(sd_total) <- c("sd_mean")
sd_total$type_enrich<- c("water_SW","water_ENR","Biofilm_SW","Biofilm_ENR")




#clacul du test 


#nos données ne sont pas normales alors :  

#Utilisation de test non paramétriques pour + de >3 groupes  
kruskal.test(sd_value ~ type, all_sd_asv)

#pour la parwaise comparaison : POST HOC DUNN TEST  
FSA::dunnTest(sd_value ~ type, all_sd_asv, method="bh")




#supplementary table


# Ajouter une colonne 'rowname' aux deux dataframes
all_sd_asv$rowname <- rownames(all_sd_asv)
all_mean_asv$rowname <- rownames(all_mean_asv)

# Renommer les colonnes pour qu'elles soient plus descriptives
colnames(all_sd_asv)[colnames(all_sd_asv) == "value"] <- "sd_value"
colnames(all_mean_asv)[colnames(all_mean_asv) == "value"] <- "mean_value"

# Fusionner les dataframes en utilisant les colonnes 'rowname' et 'type'
merged_df <- merge(all_mean_asv, all_sd_asv, by = c("rowname", "type"))

# Vous pouvez également supprimer la colonne 'rowname' si vous n'en avez plus besoin
#merged_df$rowname <- NULL

write.table(merged_df,"/16S_HOLOGREEN/alpha_diversity/alpha_indices/objets/all_mean_sd_asv.txt", sep=";")

