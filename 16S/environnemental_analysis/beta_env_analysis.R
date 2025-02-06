#env_analysis

# Assurez-vous d'installer et de charger les packages requis
#install.packages("vegan")
library(vegan)
#install.packages("phyloseq")
library(phyloseq)
library(ggplot2)





#plot(vegan::rda(otu_table(physeq_16S_nit) ~ date+pH+silicates+nitrites+phosphates+nitrates+ammonium+T.bassins, data = data.frame(sample_data(physeq_16S_nit))[c(7,12:16,21:22)]))
#colnames(sample_data(physeq_16S_nit))[c(6,7,12:16,21:22)]


#normalisation des données 

data_env_norm <- microbiome::transform(physeq_16S_nit, 'clr')





#RDA redundan analysis----   


library(phyloseq)
library(zCompositions)

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

#Environmental values 
#choisir les variables environnementales voulu
meta_RDA <- data.frame(sample_data(physeq_clr))[, c(12:16,21:22)]
colnames(meta_RDA)[colnames(meta_RDA) == "T.bassins"] <- "T°C"

#normalizaton of the meta table (z-score)
#meta_RDA[,-c(1,2)] <- scale(meta_RDA[,-c(1,2)])
#meta_RDA[,-c(1,2)] <- vegan::decostand(meta_RDA[,-c(1,2)], method = "standardize")
meta_RDA <- vegan::decostand(meta_RDA, method = "standardize")

#place the date as factor
#meta_RDA[,2] <- as.factor(meta_RDA[,2])
#

# RDA of the Aitchinson distance
# constrained by all the environmental variables
# contained in metadata
#
# Observe the shortcut formula
rda_result <- vegan::rda(t(physeq_clr_asv) ~ ., meta_RDA)


# Forward selection of explanatory variables using vegan's ordiR2step()
rda_step_fwd <- vegan::rda(t(physeq_clr_asv) ~ 1,
                           data = meta_RDA)
step_forward <- vegan::ordiR2step(rda_step_fwd,
                                  scope = formula(rda_result),
                                  direction = "forward",
                                  pstep = 1000)


#keep only significant variables

meta_RDA <- meta_RDA[, c("phosphates","pH", "nitrates","T°C","nitrites","silicates")]


# RDA of the Aitchinson distance
# constrained by all the environmental variables
# contained in metadata
rda_result <- vegan::rda(t(physeq_clr_asv) ~ ., meta_RDA)




#voir les résultats 
print(rda_result)
head(summary(rda_result))  # Scaling 2 (default)



# Unadjusted R^2 retrieved from the rda object
R2 <- vegan::RsquareAdj(rda_result)$r.squared
R2

# Adjusted R^2 retrieved from the rda object
R2adj <- vegan::RsquareAdj(rda_result)$adj.r.squared
R2adj


## Global test of the RDA result
anova_rda <- anova(rda_result, step = 1000)

# Tests of all canonical axes
anova_rda_all <- anova(rda_result, by = "axis", step = 1000)


#You can select variables that is significant 
# Variance inflation factors (VIF)
vif_res <- vegan::vif.cca(rda_result)



# Parsimonious RDA
spe_rda_pars <- vegan::rda(t(physeq_clr_asv) ~ nitrites, data = meta_RDA[,-1])
anova_spe_rda_pars <- anova(spe_rda_pars, step = 1000)
anova_spe_rda_pars_all <- anova(spe_rda_pars, step = 1000, by = "axis")

R2adj_pars <- vegan::RsquareAdj(spe_rda_pars)$adj.r.squared

# Compare variance inflation factors
vegan::vif.cca(rda_result)





#plot
#extraction des données 
res_plot <- summary(rda_result)

#fabrication des points pour les sites 
#choix des sites et du nombre de rda
res_st <- res_plot$sites[,1:2]

#rajout des samples data 
res_st <- merge(res_st, data.frame(sample_data(physeq_clr))[,c("type_enrich", "date")], by = "row.names")
colnames(res_st) <- c("Sample", "RDA1", "RDA2", "type_enrich", "date")

#fabrication  des flèches
res_arw <- as.data.frame(res_plot$biplot[, 1:2]) * 6

#eigenvalues : explained variance 
eigen_values <- format(100 * res_plot$cont[[1]][2,], digits = 4)

# Extraction des chiffres de "03" à "09" de la colonne "date"
res_st$date <- substr(res_st$date, start = 6, stop = 10)  # Extration des deux premiers caractères


#plot
RDA_plot <- ggplot() +
  #ajouts des points des sites
  geom_point(data = res_st,
             size = 4,
             aes(x = RDA1, y = RDA2,
                 shape = type_enrich, fill = type_enrich)) +
  scale_shape_manual(values = c("biofilm_SW" = 21, "biofilm_ENR" = 22, "water_SW" = 23, "water_ENR" = 24))+
  ggrepel::geom_text_repel(data = res_st,
                           aes(x=RDA1, y=RDA2, label = as.character(date)), size= 2.5)+
  #ajout des flèches des var environnementales 
  geom_segment(data = res_arw,
               arrow = arrow(angle = 22.5,
                             length = unit(0.35, "cm"),
                             type = "closed"),
               linetype = 1, size = 0.6, colour = "red",
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2)) +
  ggrepel::geom_text_repel(data = res_arw, aes(RDA1, RDA2, label=row.names(res_arw)), colour="red",fontface = "bold",
                           nudge_x = c(0,-1,1,1,1,1,1,0.5),
                           nudge_y = c(0.5,-0.5,0,-0.4,0,0,-0.5,0))+
  labs(x = paste("RDA 1 (", round(as.numeric(eigen_values[1])), "%)", sep = ""),
       y = paste("RDA 2 (", round(as.numeric(eigen_values[2])), "%)", sep = "")) +
  geom_hline(yintercept = 0, linetype = 3, size = 1) + 
  geom_vline(xintercept = 0, linetype = 3, size = 1) +
  guides(shape = guide_legend(title = NULL, color = "black"),
         fill = guide_legend(title = NULL)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  # Coloration manuelle
  scale_fill_manual(values = c("biofilm_SW" = "#33CC33", "biofilm_ENR" = "#99FF66", "water_SW" = "#0099CC", "water_ENR" = "#33CCFF"))

print(RDA_plot)


#on peut rajouter des ellipses s'il le faut 
#stat_ellipse(data = res_st, aes(x = RDA1, y = RDA2, group = type_enrich, fill = type_enrich), geom = "polygon", alpha = 0.2) +
  






###for further analysis to the paper 16S 
# Sélection des variables de saisonnalité et d'enrichissement
variables_seasonality <- c("T.bassins", "date") # Modifier avec les noms réels de vos variables
variables_enrichment <- c("pH", "nitrates", "nitrites","phosphates","ammonium") # Modifier avec les noms réels de vos variables

# Exemple de test de comparaison de moyennes (test t de Student)
seasonality_data <- meta_RDA[, variables_seasonality]
enrichment_data <- meta_RDA[, variables_enrichment]

# Test de comparaison de moyennes entre les groupes
t_test_results <- lapply(1:length(variables_seasonality), function(i) {
  t.test(seasonality_data[, i], enrichment_data[, i])
})

# Affichage des résultats des tests
print("Résultats des tests de comparaison de moyennes :")
for (i in 1:length(variables_seasonality)) {
  variable_name <- variables_seasonality[i]
  t_test_result <- t_test_results[[i]]
  
  print(paste("Variable :", variable_name))
  print(t_test_result)
}








####PCA 

#prepare the ASV table to add taxonomy
tax_CLR <-  as.data.frame(tax_table(physeq_clr)) #get taxnomic tablehttp://127.0.0.1:17243/graphics/plot_zoom_png?width=1745&height=925
#concatene ASV with Family & Genus names
ASVname <- paste(rownames(tax_CLR), tax_CLR$Family, tax_CLR$Genus,sep="_")
#apply 
rownames(physeq_clr_asv) <- ASVname
p <- PCAtools::pca(physeq_clr_asv,
                   metadata = data.frame(sample_data(physeq_clr)))
PCAtools::screeplot(p, axisLabSize = 18, titleLabSize = 22)

#variance explained by each PC
p$variance


#Horn’s parallel analysis (Horn 1965) (Buja and Eyuboglu 1992)
horn <- PCAtools::parallelPCA(physeq_clr_asv)
horn$n

#elbow method
elbow <- PCAtools::findElbowPoint(p$variance)
elbow


# #Plotting the PCA
# PCAtools::biplot(
#   p,
#   lab = p$metadata$SampName,
#   colby = "type_enrich",
#   pointSize = 5,
#   hline = 0, vline = 0,
#   legendPosition = "right"
# )


#modification de la colonne date pourcorrespondre 
p$metadata$date <- substring(p$metadata$date , 6)

# Plotting the PCA based on rotated scores and metadata
PCA_plot <- ggplot() +
  geom_point(data = p$rotated,
             aes(x = PC1, y = PC2,
                 shape = p$metadata$type_enrich, fill = p$metadata$type_enrich),
             size = 4) +
  scale_shape_manual(values = c("biofilm_SW" = 21, "biofilm_ENR" = 22, "water_SW" = 23, "water_ENR" = 24))+
  ggrepel::geom_text_repel(data = p$rotated,
                           aes(x = PC1, y = PC2, label = as.character(p$metadata$date)),
                           size = 2.5) +
  stat_ellipse(data = p$rotated, aes(x = PC1, y = PC2, group = p$metadata$type_enrich, fill = p$metadata$type_enrich), geom = "polygon", alpha = 0.2) +
  labs(x = paste("PC1 (", round(as.numeric(p$variance[1])), "%)", sep = ""),
       y = paste("PC2 (", round(as.numeric(p$variance[2])), "%)", sep = "")) +
  guides(shape = guide_legend(title = NULL, color = "black"),
         fill = guide_legend(title = NULL)) +
  geom_hline(yintercept = 0, linetype = 3, size = 1) + 
  geom_vline(xintercept = 0, linetype = 3, size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = c("biofilm_SW" = "#33CC33", "biofilm_ENR" = "#99FF66", "water_SW" = "#0099CC", "water_ENR" = "#33CCFF"))

print(PCA_plot)



#######calculation of the dispertion on PCA ----



# Supposons que vous avez déjà effectué votre PCA et obtenu les groupes de type enrich
# Supposons que "p" contient vos résultats de PCA

# Créez un vecteur de groupes à partir de la variable type_enrich
groupes <- p$metadata$type_enrich

# Créez un vecteur de données pour chaque groupe
donnees_groupes <- list()
for (g in unique(groupes)) {
  donnees_groupes[[as.character(g)]] <- p$rotated[groupes == g, c("PC1", "PC2")]  # Remplacez PC1 et PC2 par vos dimensions de PCA
}

# Calculez la variance intra-cluster pour chaque groupe
variance_intra_cluster <- sapply(donnees_groupes, function(data) {
  mean(apply(data, 2, var))  # Calcule la variance de chaque dimension (PC1, PC2, etc.) pour chaque groupe
})
