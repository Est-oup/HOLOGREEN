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

# 
# 
# #Stasticical analysis 
# # Créer un tableau de données vide
# combined_data_pca <- data.frame()
# 
# # Pour chaque groupe de données, ajouter la colonne "type_enrich" avec le nom du groupe
# for (group_name in names(donnees_groupes)) {
#   group_data_pca <- donnees_groupes[[group_name]]
#   group_data_pca$type_enrich <- group_name
#   combined_data_pca <- rbind(combined_data_pca, group_data_pca)
# }
# 
# # Afficher le tableau de données combinées
# print(combined_data_pca)
# apply(combined_data_pca, 2, var)
# 
# 
# 
# kruskal.test(observed ~ type_enrich, data=combined_data_pca)
# FSA::dunnTest(observed~ type_enrich, data= data_alpha_indices, method="bh")



PCA_plot + RDA_plot
# 
# #####For extrem  values 
# 
# # Spécifier les sites extrêmes à notifier ou enlever
# sites_extremes <- c("210407_EDM_F", "210309_EDM_A_A", "210309_EDM_B_A", "210309_EDM_C_A")
# 
# # Filtrer les données pour exclure les sites extrêmes (facultatif)
# res_st_filtered <- res_st[!res_st$Sample %in% sites_extremes, ]
# 
# # Extraction des chiffres de "03" à "09" de la colonne "date"
# res_st_filtered$date <- substr(res_st_filtered$date, start = 6, stop = 10)  # Extration des deux premiers caractères
# 
# # Inversion des chiffres
# res_st_filtered$date <- sapply(strsplit(res_st_filtered$date, "-"), function(x) paste0(rev(x), collapse = "-"))
# 
# 
# #plot
# p <- ggplot() +
#   #ajouts des points des sites
#   geom_point(data = res_st_filtered,
#              size = 4,
#              aes(x = RDA1, y = RDA2,
#                  shape = type_enrich, fill = type_enrich)) +
#   scale_shape_manual(values = c(21:25)) +
#   ggrepel::geom_text_repel(data = res_st_filtered,
#                            aes(x=RDA1, y=RDA2, label = as.character(date)), size= 2.5)+
#   #ajout des flèches des var environnementales 
#   geom_segment(data = res_arw,
#                arrow = arrow(angle = 22.5,
#                              length = unit(0.35, "cm"),
#                              type = "closed"),
#                linetype = 1, size = 0.6, colour = "red",
#                aes(x = 0, y = 0, xend = RDA1, yend = RDA2)) +
#   ggrepel::geom_text_repel(data = res_arw, aes(RDA1, RDA2, label=row.names(res_arw)), colour="red",fontface = "bold",
#                            nudge_x = c(0,-1,1,1,1,1,1,0.5),
#                            nudge_y = c(0.5,-0.5,0,-0.4,0,0,-0.5,0))+
#   labs(x = paste("RDA 1 (", round(as.numeric(eigen_values[1])), "%)", sep = ""),
#        y = paste("RDA 2 (", round(as.numeric(eigen_values[2])), "%)", sep = "")) +
#   geom_hline(yintercept = 0, linetype = 3, size = 1) + 
#   geom_vline(xintercept = 0, linetype = 3, size = 1) +
#   guides(shape = guide_legend(title = NULL, color = "black"),
#          fill = guide_legend(title = NULL)) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   # Coloration manuelle
#   scale_fill_manual(values = c("biofilm_SW" = "#33CC33", "biofilm_ENR" = "#99FF66", "filter_SW" = "#0099CC", "filter_ENR" = "#33CCFF"))
# 
# print(p)






######################################For more fun-------


# 
# 
# 
# 
# # we first replace the zeros using
# # the Count Zero Multiplicative approach
# tmp <- zCompositions::cmultRepl(physeq@otu_table,
#                                 method = "CZM",
#                                 label = 0)
# 
# # generate the centered log-ratio transformed. ASVs are in rows!!!!!
# physeq_clr_asv <- apply(tmp, 1, function(x) log(x) - mean(log(x)))
# 
# 
# 
# meta_RDA <- data.frame(sample_data(physeq))[, c(7,12:16,21:22)]
# 
# 
# #create a new phyloseq object with CLR tranformed counts
# physeq_clr <- transform(physeq,'clr')
# otu_table(physeq_clr) <- otu_table(t(physeq_clr_asv),
#                                    taxa_are_rows = FALSE)
# data.frame(physeq_clr@otu_table@.Data[1:5,1:10])
# 
# 
# 
# 
# ###tennnnnnn 
# library(PCAtools)
# 
# #prepare the ASV table to add taxonomy
# tax_CLR <-  as.data.frame(tax_table(physeq_clr)) #get taxnomic table
# #concatene ASV with Family & Genus names
# ASVname <- paste(rownames(tax_CLR), tax_CLR$Family, tax_CLR$Genus,sep="_")
# #apply 
# rownames(physeq_clr_asv) <- ASVname
# p <- PCAtools::pca(physeq_clr_asv,
#                    metadata =meta_RDA)
# PCAtools::screeplot(p, axisLabSize = 18, titleLabSize = 22)
# 
# 
# 
# 
# PCAtools::biplot(
#   p, 
#   # loadings parameters
#   showLoadings = TRUE,
#   lengthLoadingsArrowsFactor = 1.5,
#   sizeLoadingsNames = 3,
#   colLoadingsNames = 'red4',
#   ntopLoadings = 3,
#   # other parameters
#   #lab = p$metadata$Sample_ID,
#   #colby = "type_enrich",
#   hline = 0, vline = 0,
#   legendPosition = "right"
# )
# ##############################
# 
# 
# 
# 
# # Préparer les données de métadonnées pour l'analyse RDA
# #meta_RDA <- data.frame(sample_data(physeq_clr))
# 
# # Effectuer l'analyse RDA avec les variables environnementales
# rda_result <- vegan::rda(t(physeq_clr_asv) ~ .,
#                          meta_RDA)
# 
# # Afficher le screeplot
# PCAtools::screeplot(rda_result, axisLabSize = 18, titleLabSize = 22)
# 
# # Afficher le biplot avec les variables environnementales
# PCAtools::biplot(
#   rda_result,
#   showLoadings = TRUE,
#   lengthLoadingsArrowsFactor = 1.5,
#   sizeLoadingsNames = 3,
#   colLoadingsNames = 'red4',
#   ntopLoadings = 3,
#   colby = "type_enrich",
#   hline = 0, vline = 0,
#   legendPosition = "right"
# )
# 
# 
# 
# 




#########################
# 
# #Horn’s parallel analysis (Horn 1965) (Buja and Eyuboglu 1992)
# horn <- PCAtools::parallelPCA(physeq_clr_asv)
# horn$n
# 
# 
# PCAtools::eigencorplot(
#   p,
#   components = PCAtools::getComponents(p, 1:horn$n),
#   metavars = c("silicates","nitrites","phosphates","nitrates","ammonium","Phaeo","t_ext","T.bassins","pH"),
#   col = c('white', 'cornsilk1', 'gold',
#                  'forestgreen', 'darkgreen'),
#                  cexCorval = 1.2,
#   fontCorval = 2,
#   posLab = "all",
#   rotLabX = 45,
#   scale = TRUE,
#   main = bquote(PC ~ Spearman ~ r^2 ~ environmental ~ correlates),
#   plotRsquared = TRUE,
#   corFUN = "spearman",
#   corUSE = "pairwise.complete.obs",
#   corMultipleTestCorrection = 'BH',
#   signifSymbols = c("****", "***", "**", "*", ""),
#   signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
# )
# 
# 
# 
# 
# 
# meta_RDA[, c(7,12:16,21:22)]
# 
# 
# # RDA of the Aitchinson distance
# # constrained by all the environmental variables
# # contained in metadata
# #
# # Observe the shortcut formula
# spe_rda <- vegan::rda(t(physeq_clr_asv) ~ .,
#                       meta_RDA)
# head(summary(spe_rda))  # Scaling 2, (de:fault)
# 
# 
# # Unadjusted R^2 retrieved from the rda object
# R2 <- vegan::RsquareAdj(spe_rda)$r.squared
# R2
# # Adjusted R^2 retrieved from the rda object
# R2adj <- vegan::RsquareAdj(spe_rda)$adj.r.squared
# R2adj
# 
# # Global test of the RDA result
# anova(spe_rda, step = 1000)
# 
# # Tests of all canonical axes
# anova(spe_rda, by = "axis", step = 1000)
# 
# # Variance inflation factors (VIF)
# vegan::vif.cca(spe_rda)
# 
# # Forward selection of explanatory variables using vegan's ordiR2step()
# step_forward <- vegan::ordiR2step(vegan::rda(t(physeq_clr_asv) ~ 1,
#                                              data = meta_RDA),
#                                   scope = formula(spe_rda),
#                                   direction = "forward",
#                                   pstep = 1000)
# 
# 
# # Parsimonious RDA
# spe_rda_pars <- vegan::rda(t(physeq_clr_asv) ~ date, data = meta_RDA)
# anova(spe_rda_pars, step = 1000)
# 
# anova(spe_rda_pars, step = 1000, by = "axis")
# 
# R2adj_pars <- vegan::RsquareAdj(spe_rda_pars)$adj.r.squared
# 
# # Compare variance inflation factors
# vegan::vif.cca(spe_rda)
# 
# vegan::vif.cca(spe_rda_pars)
# 
# # Preparation of the data for the plot
# #
# # View analysis results
# ii <- summary(spe_rda_pars)
# 
# # Depending on the drawing result
# # the drawing data can be enlarged or
# # reduced to a certain extent, as follows
# sp <- as.data.frame(ii$species[, 1:2]) * 2
# sp_top <- sp[order(abs(sp$RDA1), decreasing = TRUE), ][1:6, ]
# 
# st <- as.data.frame(ii$sites[, 1:2])
# st <- merge(st,
#             meta_RDA,
#             by = "row.names")
# 
# yz <- t(as.data.frame(ii$biplot[, 1:2]))
# row.names(yz) <- "Salinity"
# yz <- as.data.frame(yz)
# 
# eigen_values <- format(100 *ii$cont[[1]][2,], digits=4)
# 
# #plot
# ggplot() +
#   geom_point(data = st, size = 4,
#              aes(x = RDA1, y = PC1,
#                  shape = "type_enrich", fill = "type_enrich")) +
#   scale_shape_manual(values = c(21:25)) +
#   geom_segment(data = sp_top,
#                arrow = arrow(angle = 22.5,
#                              length = unit(0.35, "cm"),
#                              type = "closed"),
#                linetype = 1, size = 0.6, colour = "red",
#                aes(x = 0, y = 0, xend = RDA1, yend = PC1)) +
#   ggrepel::geom_text_repel(data = sp_top,
#                            aes(x = RDA1, y = PC1, label = row.names(sp_top))) +
#   geom_segment(data = yz,
#                arrow = arrow(angle = 22.5,
#                              length = unit(0.35,"cm"),
#                              type = "closed"),
#                linetype = 1, size = 0.6, colour = "blue",
#                aes(x = 0, y = 0, xend = RDA1, yend = PC1)) +
#   ggrepel::geom_text_repel(data = yz, aes(RDA1, PC1, label=row.names(yz)))+
#   labs(x = paste("RDA 1 (", eigen_values[1], "%)", sep = ""),
#        y = paste("PC 1 (", eigen_values[2], "%)", sep = ""))+
#   geom_hline(yintercept = 0,linetype = 3,size = 1) + 
#   geom_vline(xintercept = 0,linetype = 3,size = 1)+
#   guides(shape = guide_legend(title = NULL,
#                               color = "black"),
#          fill = guide_legend(title = NULL))+
#   theme_bw() +
#   theme(panel.grid = element_blank())









# ####Tutoriel Github
# https://chiliubio.github.io/microeco_tutorial/explainable-class.html
# 
# # Chapter 7 Explainable class
# # The trans_env and trans_func classes are placed into the section ‘Explainable class’, as environmental factors and microbial functions can be generally applied to explain microbial community structure and assembly.
# # 
# # 7.1 trans_env class
# # There may be some NA (missing value) in the user’s env data. If so, please add complete_na = TRUE for the interpolation when creating a trans_env object.
# # 
# # 7.1.1 Example
# # Creating trans_env object has at least two ways. The following is using additional environmental data which is not in the microtable object.
# # 
# # # add_data is used to add the environmental data
# t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
# ## Env data is stored in object$data_env ...
# Maybe a more general way is to directly use the data from sample_table of your microtable object. To show this operation, we first merge additional table into sample_table to generate a new microtable object.
# 
# new_test <- clone(dataset)
# new_test$sample_table <- data.frame(new_test$sample_table, env_data_16S[rownames(new_test$sample_table), ])
# # now new_test$sample_table has the whole data
# new_test
# ## microtable-class object:
# ## sample_table have 90 rows and 15 columns
# ## otu_table have 12766 rows and 90 columns
# ## tax_table have 12766 rows and 7 columns
# ## phylo_tree have 12766 tips
# ## Taxa abundance: calculated for Kingdom,Phylum,Class,Order,Family,Genus,Species 
# ## Alpha diversity: calculated for Observed,Chao1,se.chao1,ACE,se.ACE,Shannon,Simpson,InvSimpson,Fisher,Coverage 
# ## Beta diversity: calculated for bray,jaccard
# Now let’s use env_cols to select the required columns from sample_table in the microtable object.
# 
# t1 <- trans_env$new(dataset = new_test, env_cols = 8:15)
# ## Env data is stored in object$data_env ...
# Generally, it is beneficial to analyze environmental variables in order to better use more methods. The cal_diff function is used to test the significance of variables across groups like we have shown in trans_alpha and trans_diff class parts.
# 
# # use Wilcoxon Rank Sum Test as an example
# t1$cal_diff(group = "Group", method = "wilcox")
# head(t1$res_diff)
# ## The result is stored in object$res_diff ...
# Comparison	Measure	Group	P.adj	Significance
# IW - CW	Temperature	CW	3.283e-10	***
#   IW - TW	Temperature	TW	0.0001958	***
#   CW - TW	Temperature	CW	3.283e-10	***
#   IW - CW	Precipitation	CW	5.239e-09	***
#   IW - TW	Precipitation	IW	0.07758	ns
# CW - TW	Precipitation	CW	4.002e-07	***
#   IW - CW	TOC	IW	0.4687	ns
# Let’s perform the ANOVA and show the letters in the box plot. We use list to store all the plots for each factor and plot them together.
# 
# t1$cal_diff(method = "anova", group = "Group")
# # place all the plots into a list
# tmp <- list()
# for(i in colnames(t1$data_env)){
#   tmp[[i]] <- t1$plot_diff(measure = i, add_sig_text_size = 5, xtext_size = 12) + theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm"))
# }
# plot(gridExtra::arrangeGrob(grobs = tmp, ncol = 3))
# 
# 
# From v0.12.0, trans_env class supports the differential test of groups within each group by using the by_group parameter in cal_diff function.
# 
# t1$cal_diff(group = "Type", by_group = "Group", method = "anova")
# t1$plot_diff(measure = "pH", add_sig_text_size = 5)
# 
# 
# Then we show the autocorrelations among variables.
# 
# # require GGally package to be installed
# t1$cal_autocor()
# 
# 
# For different groups, please use group parameter to show the distributions of variables and the autocorrelations across groups.
# 
# t1$cal_autocor(group = "Group")
# 
# 
# Then let’s show the RDA analysis (db-RDA and RDA).
# 
# # use bray-curtis distance for dbRDA
# t1$cal_ordination(method = "dbRDA", use_measure = "bray")
# # t1$res_rda is the result list stored in the object
# t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
# # t1$res_rda_trans is the transformed result for plotting
# t1$plot_ordination(plot_color = "Group")
# 
# 
# From v0.14.0, the function cal_ordination_anova is implemented to check the significance of the ordination model instead of the encapsulation in cal_ordination. Furthermore, the function cal_ordination_envfit can be used to get the contribution of each variables to the model.
# 
# t1$cal_ordination_anova()
# t1$cal_ordination_envfit()
# Then, let’s try to do RDA at the Genus level.
# 
# # use Genus
# t1$cal_ordination(method = "RDA", taxa_level = "Genus")
# # As the main results of RDA are related with the projection and angles between different arrows,
# # we adjust the length of the arrow to show them clearly using several parameters.
# t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
# # t1$res_rda_trans is the transformed result for plot
# t1$plot_ordination(plot_color = "Group")
# 
# 
# For more visualization styles, run the following examples.
# 
# t1$plot_ordination(plot_color = "Group", plot_shape = "Group")
# t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
# t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "centroid"))
# t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull"))
# t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse", "centroid"))
# t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull", "centroid"), add_sample_label = "SampleID")
# t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = "centroid", centroid_segment_alpha = 0.9, centroid_segment_size = 1, centroid_segment_linetype = 1)
# t1$plot_ordination(plot_color = "Type", plot_type = c("point", "centroid"), centroid_segment_linetype = 1)
# Mantel test can be used to check whether there is significant correlations between environmental variables and distance matrix.
# 
# t1$cal_mantel(use_measure = "bray")
# # return t1$res_mantel
# head(t1$res_mantel)
# ## The result is stored in object$res_mantel ...
# Variables	Correlation coefficient	p.value	p.adjusted	Significance
# Temperature	0.452	0.001	0.002	**
#   Precipitation	0.2791	0.001	0.002	**
#   TOC	0.13	0.003	0.004	**
#   NH4	-0.05539	0.926	0.926	
# NO3	0.06758	0.05	0.05714	
# pH	0.4085	0.001	0.002	**
#   Conductivity	0.2643	0.001	0.002	**
#   TN	0.1321	0.002	0.0032	**
#   # mantel test for different groups
#   t1$cal_mantel(by_group = "Group", use_measure = "bray")
# # partial mantel test
# t1$cal_mantel(partial_mantel = TRUE)
# For the combination of mantel test and correlation heatmap, please see another example (https://chiliubio.github.io/microeco_tutorial/other-examples-1.html#mantel-test-correlation-heatmap).
#                                                                                         
#                                                                                         The correlations between environmental variables and taxa are important in analyzing and inferring the factors affecting community structure. Let’s first perform a correlation heatmap using relative abundance data at Genus level with the cal_cor function. The parameter p_adjust_type can control the p value adjustment type.
#                                                                                         
#                                                                                         t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
#                                                                                         ## Env data is stored in object$data_env ...
#                                                                                         # 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
#                                                                                         t1$cal_cor(use_data = "Genus", p_adjust_method = "fdr", p_adjust_type = "Env")
#                                                                                         ## The correlation result is stored in object$res_cor ...
#                                                                                         # return t1$res_cor
#                                                                                         Then, we can plot the correlation results using plot_cor function.
#                                                                                         
#                                                                                         # default ggplot2 method with clustering
#                                                                                         t1$plot_cor()
#                                                                                         There are too many genera. We can use the filter_feature parameter in plot_cor to filter some taxa that do not have any significance < 0.001.
#                                                                                         
#                                                                                         # filter genera that donot have at least one ***
#                                                                                         t1$plot_cor(filter_feature = c("", "*", "**"))
#                                                                                         Sometimes, if the user wants to do the correlation analysis between the environmental factors and some important taxa detected in the biomarker analysis, please use other_taxa parameter in cal_cor function.
#                                                                                         
#                                                                                         # first create trans_diff object as a demonstration
#                                                                                         t2 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", taxa_level = "Genus")
#                                                                                         # then create trans_env object
#                                                                                         t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
#                                                                                         # use other_taxa to select taxa you need
#                                                                                         t1$cal_cor(use_data = "other", p_adjust_method = "fdr", other_taxa = t2$res_diff$Taxa[1:40])
#                                                                                         t1$plot_cor()
#                                                                                         
#                                                                                         
#                                                                                         The pheatmap method is also available. Note that, besides the color_vector parameter, color_palette can also be used to control color palette with customized colors.
#                                                                                         
#                                                                                         # clustering heatmap; require pheatmap package
#                                                                                         # Let's take another color pallete
#                                                                                         t1$plot_cor(pheatmap = TRUE, color_palette = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))
#                                                                                         
#                                                                                         
#                                                                                         Sometimes, if it is needed to study the correlations between environmental variables and taxa for different groups, by_group parameter can be used for this goal.
#                                                                                         
#                                                                                         # calculate correlations for different groups using parameter by_group
#                                                                                         t1$cal_cor(by_group = "Group", use_data = "other", p_adjust_method = "fdr", other_taxa = t2$res_diff$Taxa[1:40])
#                                                                                         # return t1$res_cor
#                                                                                         t1$plot_cor()
#                                                                                         
#                                                                                         
#                                                                                         If the user is concerned with the relationship between environmental factors and alpha diversity, please use add_abund_table parameter in the cal_cor function.
#                                                                                         
#                                                                                         t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
#                                                                                         # use add_abund_table parameter to add the extra data table
#                                                                                         t1$cal_cor(add_abund_table = dataset$alpha_diversity)
#                                                                                         # try to use ggplot2 with clustering plot
#                                                                                         # require ggtree and aplot packages to be installed (https://chiliubio.github.io/microeco_tutorial/intro.html#dependence)
#                                                                                         t1$plot_cor(cluster_ggplot = "both")
#                                                                                         
#                                                                                         
#                                                                                         The function plot_scatterfit in trans_env class is designed for the scatter plot, adding the fitted line and statistics of correlation or regression.
#                                                                                         
#                                                                                         # use pH and bray-curtis distance
#                                                                                         # add correlation statistics
#                                                                                         t1$plot_scatterfit(
#                                                                                           x = "pH", 
#                                                                                           y = dataset$beta_diversity$bray[rownames(t1$data_env), rownames(t1$data_env)], 
#                                                                                           type = "cor",
#                                                                                           point_size = 3, point_alpha = 0.1, 
#                                                                                           label.x.npc = "center", label.y.npc = "bottom", 
#                                                                                           x_axis_title = "Euclidean distance of pH", 
#                                                                                           y_axis_title = "Bray-Curtis distance"
#                                                                                         )
#                                                                                         
#                                                                                         
#                                                                                         # regression with type = "lm", use group parameter for different groups
#                                                                                         t1$plot_scatterfit(
#                                                                                           x = dataset$beta_diversity$bray[rownames(t1$data_env), rownames(t1$data_env)],
#                                                                                           y = "pH",
#                                                                                           type = "lm", 
#                                                                                           group = "Group", 
#                                                                                           group_order = c("CW", "TW", "IW"),
#                                                                                           point_size = 3, point_alpha = 0.3, line_se = FALSE, line_size = 1.5, shape_values = c(16, 17, 7),
#                                                                                           y_axis_title = "Euclidean distance of pH", x_axis_title = "Bray-Curtis distance"
#                                                                                         ) + theme(axis.title = element_text(size = 17))
#                                                                                         
#                                                                                         
#                                                                                         Other examples.
#                                                                                         
#                                                                                         t1 <- trans_env$new(dataset = new_test, env_cols = 8:15)
#                                                                                         # with forward selection in RDA
#                                                                                         t1$cal_ordination(method = "dbRDA", feature_sel = TRUE)
#                                                                                         # CCA, canonical correspondence analysis
#                                                                                         t1$cal_ordination(method = "CCA", taxa_level = "Genus")
#                                                                                         t1$trans_ordination(adjust_arrow_length = TRUE)
#                                                                                         t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
#                                                                                         # correlation analysis without p adjustment
#                                                                                         t1$cal_cor(use_data = "Genus", p_adjust_method = "none", use_taxa_num = 30)
#                                                                                         # correlation heatmap with clustering based on the ggplot2 and aplot packages
#                                                                                         g1 <- t1$plot_cor(cluster_ggplot = "both")
#                                                                                         g1
#                                                                                         # clustering heatmap with ggplot2 depends on aplot package
#                                                                                         # to change the detail in the plot, please manipulate each element of g1
#                                                                                         g1[[1]]
#                                                                                         # standardize x axis text format
#                                                                                         g1[[1]] <- g1[[1]] + scale_x_discrete(labels = c(NH4 = expression(NH[4]^'+'-N), NO3 = expression(NO[3]^'-'-N)))
#                                                                                         g1[[1]]
#                                                                                         g1
#                                                                                         ggplot2::ggsave("test.pdf", g1, width = 8, height = 6)
#                                                                                         # For regression, lm_equation = FALSE can be applied to not display the equation.
#                                                                                         t1$plot_scatterfit(x = 1, y = 2, type = "lm", lm_equation = TRUE)
#                                                                                         t1$plot_scatterfit(x = 1, y = 2, type = "lm", lm_equation = FALSE)
#                                                                                         t1$plot_scatterfit(x = 1, y = 2, type = "lm", point_alpha = .3, line_se = FALSE)
#                                                                                         t1$plot_scatterfit(x = 1, y = 2, type = "lm", line_se_color = "grey90", label_sep = ",", label.x.npc = "center", label.y.npc = "bottom")
#                                                                                         t1$plot_scatterfit(x = 1, y = 2, line_se = FALSE, pvalue_trim = 3, cor_coef_trim = 3)
#                                                                                         t1$plot_scatterfit(x = "pH", y = "TOC", type = "lm", group = "Group", line_se = FALSE, label.x.npc = "center",
#                                                                                                            shape_values = 1:3, x_axis_title = "pH", y_axis_title = "TOC")
#                                                                                         # correlation between relative abundance of Genus-Arthrobacter and pH
#                                                                                         tmp <- unlist(dataset$taxa_abund$Genus["k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Micrococcales|f__Micrococcaceae|g__Arthrobacter", ])
#                                                                                         t1$plot_scatterfit(x = "pH", y = tmp, point_size = 3, point_alpha = 0.3, 
#                                                                                                            y_axis_title = "Arthrobacter", x_axis_title = "pH")
#                                                                                         7.1.2 Key points
#                                                                                         complete_na parameter in trans_env$new: used to fill the NA (missing value) of the environmental data based on the mice package.
#                                                                                         env_cols parameter in trans_env$new: select the variables from sample_table of your microtable object.
#                                                                                         add_abund_table parameter in cal_cor: other customized data can be also provided for the correlation analysis.
#                                                                                         use_cor parameter in plot_scatterfit: both the correlation and regression are available in this function.
#                                                                                         cal_mantel(): partial_mantel = TRUE can be used for partial mantel test.
#                                                                                         plot_ordination(): use plot_type parameter to select point types and env_nudge_x and taxa_nudge_x (also _y) to adjust the text positions.
#                                                                                         7.2 trans_func class
#                                                                                         Ecological researchers are usually interested in the the funtional profiles of microbial communities, because functional or metabolic data is powerful to explain the structure and dynamics of microbial communities. As metagenomic sequencing is complicated and expensive, using amplicon sequencing data to predict functional profiles is an alternative choice. Several software are often used for this goal, such as PICRUSt (Langille et al. 2013), Tax4Fun (Aßhauer et al. 2015) and FAPROTAX (Stilianos Louca et al. 2016; S. Louca, Parfrey, and Doebeli 2016). These tools are great to be used for the prediction of functional profiles based on the prokaryotic communities from sequencing results. In addition, it is also important to obtain the traits or functions for each taxa, not just the whole profile of communities. FAPROTAX database is a collection of the traits and functions of prokaryotes based on the known research results published in books and literatures. We match the taxonomic information of prokaryotes against this database to predict the traits of prokaryotes on biogeochemical roles. The NJC19 database (Lim et al. 2020) is also available for animal-associated prokaryotic data, such as human gut microbiota. We also implement the FUNGuild (Nguyen et al. 2016) and FungalTraits (Põlme et al. 2020) databases to predict the fungal traits. The idea identifying prokaryotic traits and functional redundancy was initially inspired by our another study (Liu et al. 2022).
#                                                                                         
#                                                                                         7.2.1 Example
#                                                                                         We first identify/predict traits of taxa with the prokaryotic example data.
#                                                                                         
#                                                                                         # create object of trans_func
#                                                                                         t2 <- trans_func$new(dataset)
#                                                                                         # mapping the taxonomy to the database
#                                                                                         # this can recognize prokaryotes or fungi automatically if the names of taxonomic levels are standard.
#                                                                                         # for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data
#                                                                                         # default database for prokaryotes is FAPROTAX database
#                                                                                         t2$cal_spe_func(prok_database = "FAPROTAX")
#                                                                                         ## FAPROTAX v1.2.6. Please also cite the original FAPROTAX paper: Louca et al. (2016).
#                                                                                         ## Decoupling function and taxonomy in the global ocean microbiome. Science, 353(6305), 1272.
#                                                                                         ## The functional binary table is stored in object$res_spe_func ...
#                                                                                         # return t2$res_spe_func, 1 represent trait exists, 0 represent no or cannot confirmed.
#                                                                                         t2$res_spe_func[1:5, 1:2]
#                                                                                         methanotrophy	acetoclastic_methanogenesis
#                                                                                         OTU_4272	0	0
#                                                                                         OTU_236	0	0
#                                                                                         OTU_399	0	0
#                                                                                         OTU_1556	0	0
#                                                                                         OTU_32	0	0
#                                                                                         The percentages of the OTUs having the same trait can reflect the functional redundancy of this function in the community.
#                                                                                         
#                                                                                         # calculate the percentages for communities
#                                                                                         # here do not consider the abundance
#                                                                                         t2$cal_spe_func_perc(abundance_weighted = FALSE)
#                                                                                         ## The result table is stored in object$res_spe_func_perc ...
#                                                                                         # t2$res_spe_func_perc[1:5, 1:2]
#                                                                                         methanotrophy	acetoclastic_methanogenesis
#                                                                                         S1	0.39	0.04
#                                                                                         S2	0.27	0
#                                                                                         S3	0.48	0
#                                                                                         S4	0.48	0
#                                                                                         S5	0.56	0
#                                                                                         Then we also take an example to show the percentages of the OTUs for each trait in network modules.
#                                                                                         
#                                                                                         # construct a network for the example
#                                                                                         network <- trans_network$new(dataset = dataset, cal_cor = "base", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")
#                                                                                         network$cal_network(p_thres = 0.01, COR_cut = 0.7)
#                                                                                         network$cal_module()
#                                                                                         # convert module info to microtable object
#                                                                                         meco_module <- network$trans_comm(use_col = "module")
#                                                                                         meco_module_func <- trans_func$new(meco_module)
#                                                                                         meco_module_func$cal_spe_func(prok_database = "FAPROTAX")
#                                                                                         meco_module_func$cal_spe_func_perc(abundance_weighted = FALSE)
#                                                                                         meco_module_func$plot_spe_func_perc(order_x = paste0("M", 1:10))
#                                                                                         
#                                                                                         
#                                                                                         # If you want to change the group list, reset the list t2$func_group_list
#                                                                                         t2$func_group_list
#                                                                                         # use show_prok_func to see the detailed information of prokaryotic traits
#                                                                                         t2$show_prok_func("methanotrophy")
#                                                                                         # then we try to correlate the res_spe_func_perc of communities to environmental variables
#                                                                                         t3 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
#                                                                                         t3$cal_cor(add_abund_table = t2$res_spe_func_perc, cor_method = "spearman")
#                                                                                         t3$plot_cor(pheatmap = TRUE)
#                                                                                         
#                                                                                         
#                                                                                         Tax4Fun (Aßhauer et al. 2015) requires a strict input file format associated with the taxonomic information. To analyze the trimmed or changed OTU data in R with Tax4Fun, we provide a link to the Tax4Fun functional prediction. Please check out the dependence part https://chiliubio.github.io/microeco_tutorial/intro.html#tax4fun for installing Tax4Fun package and download SILVA123 ref data.
#                                                                                         
#                                                                                         t1 <- trans_func$new(dataset)
#                                                                                         # https://chiliubio.github.io/microeco_tutorial/intro.html#tax4fun for the installation description
#                                                                                         # and provide the file path of SILVA123
#                                                                                         t1$cal_tax4fun(folderReferenceData = "./SILVA123")
#                                                                                         # return two files: t1$tax4fun_KO: KO file; t1$tax4fun_path: pathway file.
#                                                                                         # t1$tax4fun_KO$Tax4FunProfile[1:5, 1:2]
#                                                                                         K00001; alcohol dehydrogenase [EC:1.1.1.1]	K00002; alcohol dehydrogenase (NADP+) [EC:1.1.1.2]
#                                                                                         S1	0.0004823	5.942e-06
#                                                                                         S2	0.0005266	4.017e-06
#                                                                                         S3	0.0005054	6.168e-06
#                                                                                         S4	0.0005109	5.888e-06
#                                                                                         S5	0.0005083	5.547e-06
#                                                                                         We further analyze the abundance of predicted metabolic pathways.
#                                                                                         
#                                                                                         # must transpose to taxa row, sample column
#                                                                                         pathway_file <- t1$tax4fun_path$Tax4FunProfile %>% t %>% as.data.frame
#                                                                                         # filter rownames, only keep ko+number
#                                                                                         rownames(pathway_file) %<>% gsub("(^.*);\\s.*", "\\1", .)
#                                                                                         # load the pathway hierarchical metadata
#                                                                                         data(Tax4Fun2_KEGG)
#                                                                                         # further create a microtable object, familiar?
#                                                                                         func1 <- microtable$new(otu_table = pathway_file, tax_table = Tax4Fun2_KEGG$ptw_desc, sample_table = t1$sample_table)
#                                                                                         print(func1)
#                                                                                         ## microtable-class object:
#                                                                                         ## sample_table have 90 rows and 4 columns
#                                                                                         ## otu_table have 284 rows and 90 columns
#                                                                                         ## tax_table have 444 rows and 3 columns
#                                                                                         Now, we need to trim data and calculate abundance.
#                                                                                         
#                                                                                         func1$tidy_dataset()
#                                                                                         # calculate abundance automatically at three levels: Level.1, Level.2, Level.3
#                                                                                         func1$cal_abund()
#                                                                                         ## The result is stored in object$taxa_abund ...
#                                                                                         print(func1)
#                                                                                         ## microtable-class object:
#                                                                                         ## sample_table have 90 rows and 4 columns
#                                                                                         ## otu_table have 282 rows and 90 columns
#                                                                                         ## tax_table have 282 rows and 3 columns
#                                                                                         ## Taxa abundance: calculated for Level.1,Level.2,Level.3
#                                                                                         Then, we can plot the abundance.
#                                                                                         
#                                                                                         # bar plot at Level.1
#                                                                                         func2 <- trans_abund$new(func1, taxrank = "Level.1", groupmean = "Group")
#                                                                                         func2$plot_bar(legend_text_italic = FALSE)
#                                                                                         
#                                                                                         
#                                                                                         We can also do something else. For example, we can use lefse to test the differences of the abundances and find the important enriched pathways across groups.
#                                                                                         
#                                                                                         func2 <- trans_diff$new(dataset = func1, method = "lefse", group = "Group", alpha = 0.05, lefse_subgroup = NULL)
#                                                                                         func2$plot_diff_bar(threshold = 3, width = 0.8)
#                                                                                         
#                                                                                         
#                                                                                         Tax4Fun2 (Wemheuer et al. 2020) is another R package for the prediction of functional profiles of prokaryotic communities from 16S rRNA gene sequences. It also provides two indexes for the evaluation of functional gene redundancies. If the user want to use Tax4Fun2 method, the representative fasta file is necessary to be added in the microtable object. Please check out https://chiliubio.github.io/microeco_tutorial/intro.html#tax4fun2 to see how to read fasta file with read.fasta of seqinr package or readDNAStringSet of Biostrings package. Please also see https://chiliubio.github.io/microeco_tutorial/intro.html#tax4fun2 for downloading ncbi-blast and Ref99NR/Ref100NR. For windows system, ncbi-blast-2.5.0+ is recommended since other versions can not operate well.
#                                                                                         
#                                                                                         # first delete the dataset created before
#                                                                                         rm(dataset)
#                                                                                         # load the example dataset from microeco package as there is the rep_fasta object in it
#                                                                                         data(dataset)
#                                                                                         dataset
#                                                                                         
#                                                                                         t1 <- trans_func$new(dataset)
#                                                                                         # create a directory for result and log files
#                                                                                         dir.create("test_prediction")
#                                                                                         # https://chiliubio.github.io/microeco_tutorial/intro.html#tax4fun2 for installation
#                                                                                         # ignore blast_tool_path parameter if blast tools have been in path
#                                                                                         # the function can search whether blast tool directory is in the path, if not, automatically use provided blast_tool_path parameter
#                                                                                         t1$cal_tax4fun2(blast_tool_path = "ncbi-blast-2.5.0+/bin", path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
#                                                                                                         database_mode = "Ref99NR", path_to_temp_folder = "test_prediction")
#                                                                                         
#                                                                                         # prepare feature table and metadata
#                                                                                         data(Tax4Fun2_KEGG)
#                                                                                         # create a microtable object for pathways
#                                                                                         func2 <- microtable$new(otu_table = t1$res_tax4fun2_pathway, tax_table = Tax4Fun2_KEGG$ptw_desc, sample_table = dataset$sample_table)
#                                                                                         func2$tidy_dataset()
#                                                                                         func2$cal_abund()
#                                                                                         
#                                                                                         # calculate functional redundancies
#                                                                                         t1$cal_tax4fun2_FRI()
#                                                                                         7.2.2 Key points
#                                                                                         blast_tool_path parameter in cal_tax4fun2: if the blast tool has been in ‘environment variable’ of computer, it is ok to use blast_tool_path = NULL
#                                                                                         blast version: tax4fun2 require NCBI blast tool. However, some errors often come from the latest versions (https://www.biostars.org/p/413294/). An easy solution is to use previous version (such as v2.5.0).
#                                                                                         References
#                                                                                         Aßhauer, Kathrin P, Bernd Wemheuer, Rolf Daniel, and Peter Meinicke. 2015. “Tax4Fun: Predicting Functional Profiles from Metagenomic 16S rRNA Data.” Journal Article. Bioinformatics 31 (17): 2882–84.
#                                                                                         Langille, M. G., J. Zaneveld, J. G. Caporaso, D. McDonald, D. Knights, J. A. Reyes, J. C. Clemente, et al. 2013. “Predictive Functional Profiling of Microbial Communities Using 16S rRNA Marker Gene Sequences.” Journal Article. Nature Biotechnology 31 (9): 814–21. https://doi.org/10.1038/nbt.2676.
#                                                                                         Lim, Roktaek, Josephine Jill T. Cabatbat, Thomas L. P. Martin, Haneul Kim, Seunghyeon Kim, Jaeyun Sung, Cheol-Min Ghim, and Pan-Jun Kim. 2020. “Large-Scale Metabolic Interaction Network of the Mouse and Human Gut Microbiota.” Journal Article. Scientific Data 7 (1). https://doi.org/10.1038/s41597-020-0516-5.
#                                                                                         Liu, Chi, Xiangzhen Li, Felipe R. P. Mansoldo, Jiaxing An, Yongping Kou, Xiao Zhang, Junming Wang, Jianxiong Zeng, Alane B. Vermelho, and Minjie Yao. 2022. “Microbial Habitat Specificity Largely Affects Microbial Co-Occurrence Patterns and Functional Profiles in Wetland Soils.” Journal Article. Geoderma 418: 115866. https://doi.org/10.1016/j.geoderma.2022.115866.
#                                                                                         Louca, S, L. W. Parfrey, and M Doebeli. 2016. “Decoupling Function and Taxonomy in the Global Ocean Microbiome.” Journal Article. Science 353 (6305): 1272.
#                                                                                         Louca, Stilianos, Saulo MS Jacques, Aliny PF Pires, Juliana S Leal, Diane S Srivastava, Laura Wegener Parfrey, Vinicius F Farjalla, and Michael Doebeli. 2016. “High Taxonomic Variability Despite Stable Functional Structure Across Microbial Communities.” Journal Article. Nature Ecology & Evolution 1: 0015.
#                                                                                         Nguyen, Nhu H., Zewei Song, Scott T. Bates, Sara Branco, Leho Tedersoo, Jon Menke, Jonathan S. Schilling, and Peter G. Kennedy. 2016. “FUNGuild: An Open Annotation Tool for Parsing Fungal Community Datasets by Ecological Guild.” Journal Article. Fungal Ecology 20 (1): 241–48.
#                                                                                         Põlme, Sergei, Kessy Abarenkov, R. Henrik Nilsson, Björn D. Lindahl, Karina Engelbrecht Clemmensen, Havard Kauserud, Nhu Nguyen, et al. 2020. “FungalTraits: A User-Friendly Traits Database of Fungi and Fungus-Like Stramenopiles.” Journal Article. Fungal Diversity 105 (1): 1–16. https://doi.org/10.1007/s13225-020-00466-2.
#                                                                                         Wemheuer, Franziska, Jessica A. Taylor, Rolf Daniel, Emma Johnston, Peter Meinicke, Torsten Thomas, and Bernd Wemheuer. 2020. “Tax4Fun2: Prediction of Habitat-Specific Functional Profiles and Functional Redundancy Based on 16S rRNA Gene Sequences.” Journal Article. Environmental Microbiome 15 (1). https://doi.org/10.1186/s40793-020-00358-7.
# 
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         
# #####
#                                                                                         library(phyloseq)
#                                                                                         library(zCompositions)
#                                                                                         
#                                                                                         # Charger les données phyloseq
#                                                                                         
#                                                                                         physeq <- physeq_16S_nit
#                                                                                         
#                                                                                         # Prétraitement des données
#                                                                                         physeq <- filter_taxa(physeq, function(x) sum(x > 0) > 0, TRUE)
#                                                                                         
#                                                                                         # we first replace the zeros using
#                                                                                         # the Count Zero Multiplicative approach
#                                                                                         tmp <- zCompositions::cmultRepl(physeq@otu_table,
#                                                                                                                         method = "CZM",
#                                                                                                                         label = 0)
#                                                                                         
#                                                                                         # generate the centered log-ratio transformed. ASVs are in rows!!!!!
#                                                                                         physeq_clr_asv <- apply(tmp, 1, function(x) log(x) - mean(log(x)))
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         meta_RDA <- data.frame(sample_data(physeq))[, c(7,12:16,21:22)]
#                                                                                         
#                                                                                         
#                                                                                         #create a new phyloseq object with CLR tranformed counts
#                                                                                         physeq_clr <- transform(physeq,'clr')
#                                                                                         otu_table(physeq_clr) <- otu_table(t(physeq_clr_asv),
#                                                                                                                            taxa_are_rows = FALSE)
#                                                                                         data.frame(physeq_clr@otu_table@.Data[1:5,1:10])
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         # 8.1 Redundant analysis (RDA)
#                                                                                         # 8.1.1 Running the RDA
#                                                                                         # RDA is a method combining regression and principal component analysis (PCA). RDA computes axes that are linear combinations of the explanatory variables. In RDA, one can truly say that the axes explain or model (in the statistical sense) the variation of the dependent matrix.
#                                                                                         
#                                                                                         # RDA of the Aitchinson distance
#                                                                                         # constrained by all the environmental variables
#                                                                                         # contained in metadata
#                                                                                         #
#                                                                                         # Observe the shortcut formula
#                                                                                         spe_rda <- vegan::rda(t(physeq_clr_asv) ~ .,
#                                                                                                               metadata)
#                                                                                         head(summary(spe_rda))  # Scaling 2 (default)
#                                                                                         
#                                                                                         
#                                                                                         Call:
#                                                                                           rda(formula = t(physeq_clr_asv) ~ SiOH4 + NO2 + NO3 + NH4 + PO4 +      NT + PT + Chla + T + S + Sigma_t, data = metadata[, 11:21]) 
#                                                                                         
#                                                                                         Partitioning of variance:
#                                                                                           Inertia Proportion
#                                                                                         Total          328.56     1.0000
#                                                                                         Constrained    231.46     0.7044
#                                                                                         Unconstrained   97.11     0.2956
#                                                                                         
#                                                                                         Eigenvalues, and their contribution to the variance 
#                                                                                         
#                                                                                         Importance of components:
#                                                                                           RDA1     RDA2     RDA3     RDA4     RDA5     RDA6
#                                                                                         Eigenvalue            85.2928 30.29173 20.29415 18.85659 15.83909 12.98651
#                                                                                         Proportion Explained   0.2596  0.09219  0.06177  0.05739  0.04821  0.03952
#                                                                                         Cumulative Proportion  0.2596  0.35179  0.41355  0.47094  0.51915  0.55868
#                                                                                         RDA7     RDA8     RDA9   RDA10   RDA11      PC1
#                                                                                         Eigenvalue            11.78027 10.97738 10.18119 7.94385 7.01222 28.88564
#                                                                                         Proportion Explained   0.03585  0.03341  0.03099 0.02418 0.02134  0.08791
#                                                                                         Cumulative Proportion  0.59453  0.62794  0.65893 0.68310 0.70445  0.79236
#                                                                                         PC2     PC3      PC4      PC5     PC6
#                                                                                         Eigenvalue            16.45693 16.3958 15.58129 11.19715 8.59184
#                                                                                         Proportion Explained   0.05009  0.0499  0.04742  0.03408 0.02615
#                                                                                         Cumulative Proportion  0.84245  0.8923  0.93977  0.97385 1.00000
#                                                                                         
#                                                                                         Accumulated constrained eigenvalues
#                                                                                         Importance of components:
#                                                                                           RDA1    RDA2     RDA3     RDA4     RDA5     RDA6
#                                                                                         Eigenvalue            85.2928 30.2917 20.29415 18.85659 15.83909 12.98651
#                                                                                         Proportion Explained   0.3685  0.1309  0.08768  0.08147  0.06843  0.05611
#                                                                                         Cumulative Proportion  0.3685  0.4994  0.58706  0.66853  0.73696  0.79307
#                                                                                         RDA7     RDA8     RDA9   RDA10  RDA11
#                                                                                         Eigenvalue            11.7803 10.97738 10.18119 7.94385 7.0122
#                                                                                         Proportion Explained   0.0509  0.04743  0.04399 0.03432 0.0303
#                                                                                         Cumulative Proportion  0.8440  0.89139  0.93538 0.96970 1.0000
#                                                                                         
#                                                                                         Scaling 2 for species and site scores
#                                                                                         * Species are scaled proportional to eigenvalues
#                                                                                         * Sites are unscaled: weighted dispersion equal on all dimensions
#                                                                                         * General scaling constant of scores:  8.645047 
#                                                                                         
#                                                                                         
#                                                                                         Species scores
#                                                                                         
#                                                                                         RDA1      RDA2     RDA3
#                                                                                         ASV1_Cyanobiaceae_Synechococcus CC9902        -0.1033  0.108773  0.04666
#                                                                                         ASV2_Pseudoalteromonadaceae_Pseudoalteromonas -0.7807 -0.229145 -0.22860
#                                                                                         ASV3_Clade I_Clade Ia                         -0.2568  0.002182 -0.22536
#                                                                                         ASV4_NA_NA                                    -0.6996  0.193071  0.23547
#                                                                                         ASV5_Clade I_Clade Ia                          0.5264 -0.195773  0.23032
#                                                                                         ASV6_Clade II_NA                              -0.2542 -0.344583 -0.32380
#                                                                                         ....                                                                    
#                                                                                         RDA4     RDA5     RDA6
#                                                                                         ASV1_Cyanobiaceae_Synechococcus CC9902        -0.12535 -0.01552 -0.06487
#                                                                                         ASV2_Pseudoalteromonadaceae_Pseudoalteromonas  0.33352  0.13369 -0.08880
#                                                                                         ASV3_Clade I_Clade Ia                         -0.04191 -0.04528  0.11436
#                                                                                         ASV4_NA_NA                                     0.20648 -0.23531 -0.06807
#                                                                                         ASV5_Clade I_Clade Ia                         -0.05792  0.40196  0.22286
#                                                                                         ASV6_Clade II_NA                              -0.31352 -0.10920  0.06137
#                                                                                         ....                                                                    
#                                                                                         
#                                                                                         
#                                                                                         Site scores (weighted sums of species scores)
#                                                                                         
#                                                                                         RDA1     RDA2    RDA3    RDA4     RDA5    RDA6
#                                                                                         S11B -1.703 -1.23820  2.9437 -0.2362  1.13728  0.4405
#                                                                                         S1B   2.565 -0.13340 -0.7868 -5.7453  3.30268  3.3657
#                                                                                         S2B   3.022 -2.96571  0.4021 -0.9802 -3.09213 -0.9282
#                                                                                         S2S  -1.731 -1.82618  2.0707 -0.3281 -0.66853  1.6638
#                                                                                         S3B   3.624 -1.55655 -1.2829 -2.0701 -2.02586 -1.7347
#                                                                                         S3S   3.165 -0.08923  2.8998  2.0441 -0.08464 -2.0314
#                                                                                         ....                                                 
#                                                                                         
#                                                                                         
#                                                                                         Site constraints (linear combinations of constraining variables)
#                                                                                         
#                                                                                         RDA1    RDA2    RDA3    RDA4    RDA5    RDA6
#                                                                                         S11B -1.2105 -0.7764  3.0649 -0.2199  1.2569 -0.7586
#                                                                                         S1B   1.7387  0.3983 -0.3817 -5.4943  3.2411  2.7484
#                                                                                         S2B   2.0536 -3.3237  0.6260 -1.4897 -2.8936 -0.1774
#                                                                                         S2S   0.5936 -2.0609  1.1588 -0.1736 -0.8183  1.8069
#                                                                                         S3B   4.1498 -1.1569 -1.6837 -1.1942 -2.4216 -2.5295
#                                                                                         S3S   2.0704 -0.1285  3.6947  1.1733  0.3885 -1.8438
#                                                                                         ....                                                
#                                                                                         
#                                                                                         
#                                                                                         Biplot scores for constraining variables
#                                                                                         
#                                                                                         RDA1     RDA2     RDA3     RDA4     RDA5      RDA6
#                                                                                         SiOH4   -0.57424 -0.21106 -0.25450  0.25678 -0.02349  0.213981
#                                                                                         NO2     -0.51463 -0.10086 -0.08171 -0.34294  0.35340 -0.013696
#                                                                                         NO3      0.59878  0.05632 -0.04267  0.02065 -0.30772 -0.095439
#                                                                                         NH4     -0.63097 -0.49073 -0.01146  0.07457  0.25646 -0.259440
#                                                                                         PO4     -0.49369 -0.05367 -0.31521 -0.04459  0.19877 -0.304690
#                                                                                         NT       0.02778 -0.05873 -0.28198 -0.59590  0.14825  0.392684
#                                                                                         PT      -0.61634 -0.27995 -0.01129 -0.12013  0.07328  0.533916
#                                                                                         Chla    -0.47936 -0.07832 -0.06090  0.01293 -0.11376 -0.179421
#                                                                                         T       -0.57485  0.21879  0.26190 -0.53662 -0.42902 -0.007286
#                                                                                         S       -0.93622  0.00815 -0.06712 -0.05543  0.04078 -0.183950
#                                                                                         Sigma_t -0.52380 -0.20293 -0.31121  0.40702  0.43162 -0.205711
#                                                                                         The included environmental variables explain 70.44% of the variation in bacterial community composition across sites. 29.56 % of the variance is unexplained. However, we’ll see that the propotion of variance explained is much lower. The R2 from the summary measures the strength of the canonical relationship between the response variables (Y matrix, ASVs) and the explanatory variables (X matrix) by calculating the proportion of the variation of Y explained by the variables in X. However, this R2 is biased. We calculate an Adjusted R2, which also measures the strength of the relationship between Y and X, but applies a correction of the R2 to take into account the number of explanatory variables. This is the statistic that should be reported.
#                                                                                         
#                                                                                         # Unadjusted R^2 retrieved from the rda object
#                                                                                         R2 <- vegan::RsquareAdj(spe_rda)$r.squared
#                                                                                         R2
#                                                                                         
#                                                                                         [1] 0.7044457
#                                                                                         # Adjusted R^2 retrieved from the rda object
#                                                                                         R2adj <- vegan::RsquareAdj(spe_rda)$adj.r.squared
#                                                                                         R2adj
#                                                                                         
#                                                                                         [1] 0.1625961
#                                                                                         In reality, the proportion of variance explained dropped to 16.25 %. The numerical output shows that the first two canonical axes explain together 35.1% of the total variance of the data, the first axis alone explaining 25.9%. These are unadjusted values, however. Since R2 adj = 16.2 %, the percentages of accumulated constrained adj eigenvalues show that the first axis alone explains 0.162 * 0.368 = 0.059 or 5.9% variance. Because ecological data are generally quite noisy, one should never expect to obtain a very high value of R2 . Furthermore, the first unconstrained eigenvalue (PC1), the first unconstrained axe for the residuals, is comparatively high, which means that it does display an important residual structure of the response data that is not explain by the environmental parameters measure here.
#                                                                                         
#                                                                                         8.1.2 Significance testing
#                                                                                         The interpretation of the constrained ordination must be preceded by a test of statistical significance (see below). As in multiple regression, a non-significant result must not be interpreted and must be discarded.
#                                                                                         
#                                                                                         # Global test of the RDA result
#                                                                                         anova(spe_rda, step = 1000)
#                                                                                         
#                                                                                         Permutation test for rda under reduced model
#                                                                                         Permutation: free
#                                                                                         Number of permutations: 999
#                                                                                         
#                                                                                         Model: rda(formula = t(physeq_clr_asv) ~ SiOH4 + NO2 + NO3 + NH4 + PO4 + NT + PT + Chla + T + S + Sigma_t, data = metadata[, 11:21])
#                                                                                         Df Variance      F Pr(>F)  
#                                                                                         Model    11  231.456 1.3001  0.097 .
#                                                                                         Residual  6   97.109                
#                                                                                         ---
#                                                                                           Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#                                                                                         # Tests of all canonical axes
#                                                                                         anova(spe_rda, by = "axis", step = 1000)
#                                                                                         
#                                                                                         Permutation test for rda under reduced model
#                                                                                         Forward tests for axes
#                                                                                         Permutation: free
#                                                                                         Number of permutations: 999
#                                                                                         
#                                                                                         Model: rda(formula = t(physeq_clr_asv) ~ SiOH4 + NO2 + NO3 + NH4 + PO4 + NT + PT + Chla + T + S + Sigma_t, data = metadata[, 11:21])
#                                                                                         Df Variance      F Pr(>F)
#                                                                                         RDA1      1   85.293 5.2699  0.107
#                                                                                         RDA2      1   30.292 1.8716  0.615
#                                                                                         RDA3      1   20.294 1.2539  0.999
#                                                                                         RDA4      1   18.857 1.1651  1.000
#                                                                                         RDA5      1   15.839 0.9786  1.000
#                                                                                         RDA6      1   12.987 0.8024  1.000
#                                                                                         RDA7      1   11.780 0.7279  1.000
#                                                                                         RDA8      1   10.977 0.6783  1.000
#                                                                                         RDA9      1   10.181 0.6291  0.999
#                                                                                         RDA10     1    7.944 0.4908  0.993
#                                                                                         RDA11     1    7.012 0.4333  0.911
#                                                                                         Residual  6   97.109              
#                                                                                         Here we can see that ur full model is statistically non significant (p = 0.08), and every canonical axis resulting from the RDA are not either statistically significant (p > 0.05). This RDA model is not interpretable.
#                                                                                         
#                                                                                         Can you tell why?
#                                                                                           
#                                                                                           8.1.3 Selecting variables
#                                                                                         It happens sometimes that one wishes to reduce the number of explanatory variables. The reasons vary: search for parsimony, rich data set but poor a priori hypotheses and possible strong linear dependencies (correlations) among the explanatory variables in the RDA model, which could render the regression coefficients of the explanatory variables in the model unstable.
#                                                                                         
#                                                                                         A simple approach to identify collinearity among explanatory variables is the use of variance inflation factors (VIF). VIF calculations are straightforward and easily comprehensible; the higher the value, the higher the collinearity. VIF measure the proportion by which the variance of a regression coefficient is inflated in the presence of other explanatory variables. VIFs above 20 indicate strong collinearity. Ideally, VIFs above 10 should be at least examined, and avoided if possible.
#                                                                                         
#                                                                                         # Variance inflation factors (VIF)
#                                                                                         vegan::vif.cca(spe_rda)
#                                                                                         
#                                                                                         SiOH4         NO2         NO3         NH4         PO4          NT 
#                                                                                         4.066588    3.489186    3.634643   16.867288    8.819736    4.908553 
#                                                                                         PT        Chla           T           S     Sigma_t 
#                                                                                         6.835572    2.264012 5417.455601 8388.550079 6878.896122 
#                                                                                         Salinity, Temperature and Sigma.t have very hight VIFs wich confirm the collinearities observed earlier between explanatory variables (see the PERMANOVA section). A reduction of the number of explanatory variables is justified. In order to simplify this model, we can perform a forward selection (or backwards or stepwise). These types of selections help us select variables that are statistically important. However, it is important to note that selecting variables ecologically is much more important than performing selection in this way. If a variable of ecological interest is not selected, this does not mean it has to be removed from the RDA. Here, we will be performing forward selection on our 11 environmental variables. To do this, we can use the ordiR2step() function:
#                                                                                           
#                                                                                           # Forward selection of explanatory variables using vegan's ordiR2step()
#                                                                                           step_forward <- vegan::ordiR2step(vegan::rda(t(physeq_clr_asv) ~ 1,
#                                                                                                                                        data = metadata),
#                                                                                                                             scope = formula(spe_rda),
#                                                                                                                             direction = "forward",
#                                                                                                                             pstep = 1000)
#                                                                                         
#                                                                                         Step: R2.adj= 0 
#                                                                                         Call: t(physeq_clr_asv) ~ 1 
#                                                                                         
#                                                                                         R2.adjusted
#                                                                                         + S              0.18366030
#                                                                                         <All variables>  0.16259613
#                                                                                         + NH4            0.08392874
#                                                                                         + PT             0.07013415
#                                                                                         + T              0.06719602
#                                                                                         + NO3            0.05904665
#                                                                                         + SiOH4          0.05787026
#                                                                                         + Sigma_t        0.05002017
#                                                                                         + NO2            0.03846019
#                                                                                         + PO4            0.03190148
#                                                                                         + Chla           0.02451726
#                                                                                         <none>           0.00000000
#                                                                                         + NT            -0.01448677
#                                                                                         Here, we are essentially adding one variable at a time, and retaining it if it significantly increases the model’s adjusted R2. The forward selection show us that a model with only salinity has higher R2 adjust than with all variable and explain 18.4 % of the variance. Lets calculate this most parsimonious RDA and check its significance.
#                                                                                         
#                                                                                         # Parsimonious RDA
#                                                                                         spe_rda_pars <- vegan::rda(t(physeq_clr_asv) ~ S, data = metadata[, 11:21])
#                                                                                         anova(spe_rda_pars, step = 1000)
#                                                                                         
#                                                                                         Permutation test for rda under reduced model
#                                                                                         Permutation: free
#                                                                                         Number of permutations: 999
#                                                                                         
#                                                                                         Model: rda(formula = t(physeq_clr_asv) ~ S, data = metadata[, 11:21])
#                                                                                         Df Variance      F Pr(>F)    
#                                                                                         Model     1   76.122 4.8247  0.001 ***
#                                                                                           Residual 16  252.443                  
#                                                                                         ---
#                                                                                           Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#                                                                                         anova(spe_rda_pars, step = 1000, by = "axis")
#                                                                                         
#                                                                                         Permutation test for rda under reduced model
#                                                                                         Forward tests for axes
#                                                                                         Permutation: free
#                                                                                         Number of permutations: 999
#                                                                                         
#                                                                                         Model: rda(formula = t(physeq_clr_asv) ~ S, data = metadata[, 11:21])
#                                                                                         Df Variance      F Pr(>F)    
#                                                                                         RDA1      1   76.122 4.8247  0.001 ***
#                                                                                           Residual 16  252.443                  
#                                                                                         ---
#                                                                                           Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#                                                                                         R2adj_pars <- vegan::RsquareAdj(spe_rda_pars)$adj.r.squared
#                                                                                         
#                                                                                         # Compare variance inflation factors
#                                                                                         vegan::vif.cca(spe_rda)
#                                                                                         
#                                                                                         SiOH4         NO2         NO3         NH4         PO4          NT 
#                                                                                         4.066588    3.489186    3.634643   16.867288    8.819736    4.908553 
#                                                                                         PT        Chla           T           S     Sigma_t 
#                                                                                         6.835572    2.264012 5417.455601 8388.550079 6878.896122 
#                                                                                         vegan::vif.cca(spe_rda_pars)
#                                                                                         
#                                                                                         S 
#                                                                                         1 
#                                                                                         Now, both the model and the first canonical axis resulting from the RDA are statistically significant (p < 0.05). The VIF of salinity is only 1. This RDA model is interpretable. Lets plot it.
#                                                                                         
#                                                                                         8.1.4 RDA plot
#                                                                                         # Preparation of the data for the plot
#                                                                                         #
#                                                                                         # View analysis results
#                                                                                         ii <- summary(spe_rda_pars)
#                                                                                         
#                                                                                         # Depending on the drawing result
#                                                                                         # the drawing data can be enlarged or
#                                                                                         # reduced to a certain extent, as follows
#                                                                                         sp <- as.data.frame(ii$species[, 1:2]) * 2
#                                                                                         sp_top <- sp[order(abs(sp$RDA1), decreasing = TRUE), ][1:6, ]
#                                                                                         
#                                                                                         st <- as.data.frame(ii$sites[, 1:2])
#                                                                                         st <- merge(st,
#                                                                                                     metadata["Geo"],
#                                                                                                     by = "row.names")
#                                                                                         
#                                                                                         yz <- t(as.data.frame(ii$biplot[, 1:2]))
#                                                                                         row.names(yz) <- "Salinity"
#                                                                                         yz <- as.data.frame(yz)
#                                                                                         
#                                                                                         eigen_values <- format(100 *ii$cont[[1]][2,], digits=4)
#                                                                                         
#                                                                                         #plot
#                                                                                         ggplot() +
#                                                                                           geom_point(data = st, size = 4,
#                                                                                                      aes(x = RDA1, y = PC1,
#                                                                                                          shape = Geo, fill = Geo)) +
#                                                                                           scale_shape_manual(values = c(21:25)) +
#                                                                                           geom_segment(data = sp_top,
#                                                                                                        arrow = arrow(angle = 22.5,
#                                                                                                                      length = unit(0.35, "cm"),
#                                                                                                                      type = "closed"),
#                                                                                                        linetype = 1, size = 0.6, colour = "red",
#                                                                                                        aes(x = 0, y = 0, xend = RDA1, yend = PC1)) +
#                                                                                           ggrepel::geom_text_repel(data = sp_top,
#                                                                                                                    aes(x = RDA1, y = PC1, label = row.names(sp_top))) +
#                                                                                           geom_segment(data = yz,
#                                                                                                        arrow = arrow(angle = 22.5,
#                                                                                                                      length = unit(0.35,"cm"),
#                                                                                                                      type = "closed"),
#                                                                                                        linetype = 1, size = 0.6, colour = "blue",
#                                                                                                        aes(x = 0, y = 0, xend = RDA1, yend = PC1)) +
#                                                                                           ggrepel::geom_text_repel(data = yz, aes(RDA1, PC1, label=row.names(yz)))+
#                                                                                           labs(x = paste("RDA 1 (", eigen_values[1], "%)", sep = ""),
#                                                                                                y = paste("PC 1 (", eigen_values[2], "%)", sep = ""))+
#                                                                                           geom_hline(yintercept = 0,linetype = 3,size = 1) + 
#                                                                                           geom_vline(xintercept = 0,linetype = 3,size = 1)+
#                                                                                           guides(shape = guide_legend(title = NULL,
#                                                                                                                       color = "black"),
#                                                                                                  fill = guide_legend(title = NULL))+
#                                                                                           theme_bw() +
#                                                                                           theme(panel.grid = element_blank())
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         One of the most powerful aspects of RDA is the simultaneous visualization of your response and explanatory variables (i.e. species and environmental variables). From this ordination, we can really say now that salinity is the main environmental driver measured shaping bacterial communities. Among all the ASVs, some are more related to this gradient of salinity. This is the case of ASV 12 and 11 for which abundance increase when salinity decreases and ASV 7 which presents the opposite pattern. These differential abundance patterns can be explored with many kind of analyses (see next chapter) but what is really powerful with RDA is that you highlight gradient relationships not a difference of abundance between two conditions. However, a large part of the variance in the bacterial community remains unexplained. Variance in species communities can be explained by deterministic processes such as species sorting (influence of the environment as we’ve seen here) but also by stochastic processes such as dispersal which depend, among other things, of the distance between communities. Since we have this information, lets take a look at a very common pattern in community ecology: the distance-decay pattern.
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         #####################################################################
#                                                                                         
#                                                                                         remotes::install_github("gavinsimpson/ggvegan")
#                                                                                         library(ggvegan)
#                                                                                         
#                                                                                         
#                                                                                         install.packages(
#                                                                                           "microViz",
#                                                                                           repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
#                                                                                         )
#                                                                                         library(microViz)
#                                                                                         #ggvegan::
#                                                                                         library("ggord")
#                                                                                         
#                                                                                         autoplot(spe_rda, layers = c("sites","biplot"),
#                                                                                                  arrows = TRUE,
#                                                                                                  legend.position = "right")+theme_minimal()
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         #plotRDA(spe_rda)
#                                                                                         arrow_coords <- scores(spe_rda, display = "sites")
#                                                                                         
#                                                                                         sample_data <- cbind(sample_data(physeq_16S_nit), arrow_coords)
#                                                                                         
#                                                                                         
#                                                                                         summary(spe_rda)
#                                                                                         
#                                                                                         
#                                                                                         # Preparation of the data for the plot
#                                                                                         #
#                                                                                         # View analysis results
#                                                                                         ii <- summary(spe_rda)
#                                                                                         
#                                                                                         # Depending on the drawing result
#                                                                                         # the drawing data can be enlarged or
#                                                                                         # reduced to a certain extent, as follows
#                                                                                         sp <- as.data.frame(ii$species[, 1:2]) * 2
#                                                                                         sp_top <- sp[order(abs(sp$RDA1), decreasing = TRUE), ][1:6, ]
#                                                                                         
#                                                                                         st <- as.data.frame(ii$sites[, 1:2])
#                                                                                         st <- merge(st,
#                                                                                                     metadata["Geo"],
#                                                                                                     by = "row.names")
#                                                                                         
#                                                                                         yz <- t(as.data.frame(ii$biplot[, 1:2]))
#                                                                                         row.names(yz) <- "Salinity"
#                                                                                         yz <- as.data.frame(yz)
#                                                                                         
#                                                                                         eigen_values <- format(100 *ii$cont[[1]][2,], digits=4)
#                                                                                         
#                                                                                         #plot
#                                                                                         ggplot() +
#                                                                                           geom_point(data = st, size = 4,
#                                                                                                      aes(x = RDA1, y = PC1,
#                                                                                                          shape = Geo, fill = Geo)) +
#                                                                                           scale_shape_manual(values = c(21:25)) +
#                                                                                           geom_segment(data = sp_top,
#                                                                                                        arrow = arrow(angle = 22.5,
#                                                                                                                      length = unit(0.35, "cm"),
#                                                                                                                      type = "closed"),
#                                                                                                        linetype = 1, size = 0.6, colour = "red",
#                                                                                                        aes(x = 0, y = 0, xend = RDA1, yend = PC1)) +
#                                                                                           ggrepel::geom_text_repel(data = sp_top,
#                                                                                                                    aes(x = RDA1, y = PC1, label = row.names(sp_top))) +
#                                                                                           geom_segment(data = yz,
#                                                                                                        arrow = arrow(angle = 22.5,
#                                                                                                                      length = unit(0.35,"cm"),
#                                                                                                                      type = "closed"),
#                                                                                                        linetype = 1, size = 0.6, colour = "blue",
#                                                                                                        aes(x = 0, y = 0, xend = RDA1, yend = PC1)) +
#                                                                                           ggrepel::geom_text_repel(data = yz, aes(RDA1, PC1, label=row.names(yz)))+
#                                                                                           labs(x = paste("RDA 1 (", eigen_values[1], "%)", sep = ""),
#                                                                                                y = paste("PC 1 (", eigen_values[2], "%)", sep = ""))+
#                                                                                           geom_hline(yintercept = 0,linetype = 3,size = 1) + 
#                                                                                           geom_vline(xintercept = 0,linetype = 3,size = 1)+
#                                                                                           guides(shape = guide_legend(title = NULL,
#                                                                                                                       color = "black"),
#                                                                                                  fill = guide_legend(title = NULL))+
#                                                                                           theme_bw() +
#                                                                                           theme(panel.grid = element_blank())
#                                                                                         
#                                                                                         
#                                                                                         gglot()
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         ##
#                                                                                         ord_points(spe_rda, color_var = "group")
#                                                                                         ord_plot(spe_rda, axes = 1:2)
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         physeq_16S_nit %>%
#                                                                                           tax_transform("clr") %>%
#                                                                                           ord_calc(
#                                                                                             unconstraints = c("silicates","nitrites","phosphates","nitrates","ammonium","T.bassins","pH"),
#                                                                                             method = "RDA", # Note: you can specify RDA explicitly, and it is good practice to do so, but microViz can guess automatically that you want an RDA here (helpful if you don't remember the name?)
#                                                                                             scale_cc = FALSE # doesn't make a difference
#                                                                                           ) %>%
#                                                                                           ord_plot(
#                                                                                             
#                                                                                           )
# 
#                                                                                         
#                                                                                         
# ############################################################################################################
#                                                                                         
#                                                                                         physeq_clr
#                                                                                         
#                                                                                         #RDA redundan analysis----   
#                                                                                         
#                                                                                         
#                                                                                         library(phyloseq)
#                                                                                         library(zCompositions)
#                                                                                         
#                                                                                         # Charger les données phyloseq
#                                                                                         
#                                                                                         physeq <- physeq_16S_nit
#                                                                                         
#                                                                                         # Prétraitement des données
#                                                                                         physeq <- filter_taxa(physeq, function(x) sum(x > 0) > 0, TRUE)
#                                                                                         
#                                                                                         # we first replace the zeros using
#                                                                                         # the Count Zero Multiplicative approach
#                                                                                         tmp <- zCompositions::cmultRepl(physeq@otu_table,
#                                                                                                                         method = "CZM",
#                                                                                                                         label = 0)
#                                                                                         
#                                                                                         # generate the centered log-ratio transformed. ASVs are in rows!!!!!
#                                                                                         physeq_clr_asv <- apply(tmp, 1, function(x) log(x) - mean(log(x)))
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         meta_RDA <- data.frame(sample_data(physeq))[, c(7,12:16,21:22)]
#                                                                                         
#                                                                                         
#                                                                                         #create a new phyloseq object with CLR tranformed counts
#                                                                                         physeq_clr <- transform(physeq,'clr')
#                                                                                         otu_table(physeq_clr) <- otu_table(t(physeq_clr_asv),
#                                                                                                                            taxa_are_rows = FALSE)
#                                                                                         
#                                                                                         
#                                                                                         
#                                                                                         # 8.1 Redundant analysis (RDA)
#                                                                                         # 8.1.1 Running the RDA
#                                                                                         # RDA is a method combining regression and principal component analysis (PCA). RDA computes axes that are linear combinations of the explanatory variables. In RDA, one can truly say that the axes explain or model (in the statistical sense) the variation of the dependent matrix.
#                                                                                         
#                                                                                         # RDA of the Aitchinson distance
#                                                                                         # constrained by all the environmental variables
#                                                                                         # contained in metadata
#                                                                                         #
#                                                                                         # Observe the shortcut formula
#                                                                                         spe_rda <- vegan::rda(t(physeq_clr_asv) ~ .,
#                                                                                                               metadata)
#                                                                                         # RDA avec tous les variables environnementales
#                                                                                         
#                                                                                         
#                                                                                         # Sélection progressive des variables environnementales significatives
#                                                                                         fwd.sel <- ordiR2step(rda(mite.spe.hel ~ 1, data = mite.env),
#                                                                                                               scope = formula(mite.spe.rda),
#                                                                                                               direction = "forward",
#                                                                                                               R2scope = TRUE, pstep = 1000, trace = FALSE)
#                                                                                         fwd.sel$call
#                                                                                         # rda(formula = mite.spe.hel ~ WatrCont + Shrub + Substrate + Topo, 
#                                                                                         #     data = mite.env)
                                                                                        