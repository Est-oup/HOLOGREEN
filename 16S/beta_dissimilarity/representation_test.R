#beta_analysis


#loadings
library(ggplot2)
library(phyloseq)
library(zCompositions)

#normalisation of the data 

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









#PCA Treatment

# Installer les librairies si elles ne sont pas déjà installées
install.packages("phyloseq")
install.packages("ade4")
install.packages("ggplot2")

# Charger les librairies
library(phyloseq)
library(ade4)
library(ggplot2)

# Supposons que vous ayez déjà chargé votre objet phyloseq sous le nom physeq_clr

# Séparer les échantillons en fonction de la variable "type_enrich"
sample_data(physeq_clr)$type_enrich <- as.factor(sample_data(physeq_clr)$type_enrich)

# Calculer la distance de Aitchison entre les échantillons
dist_aitchison <- dist(otu_table(physeq_clr), method = "euclidean")

pca_result <- cmdscale(dist_aitchison, k = 2)

# Créer un plot de la PCA en utilisant ggplot2
pca_data <- as.data.frame(pca_result)
pca_data$type_enrich <- sample_data(physeq_clr)$type_enrich

pca_plot <- ggplot(pca_data, aes(x = pca_data[,1], y = pca_data[,2], color = type_enrich)) +
  geom_point() +
  labs(title = "PCA plot of Aitchison Transformed Data", x = "PC1", y = "PC2")

# Afficher le plot
print(pca_plot)


# Calculer les valeurs propres pour la variance expliquée
eigenvalues <- (pca_result$eig / sum(pca_result$eig)) * 100

# Créer un graphique en barre pour afficher la variance expliquée
variance_bar <- ggplot(NULL, aes(x = factor(1:2), y = eigenvalues)) +
  geom_bar(stat = "identity", fill = "dodgerblue") +
  labs(title = "Variance Explained by Principal Components",
       x = "Principal Component",
       y = "Percentage of Variance Explained") +
  scale_x_discrete(labels = c("PC1", "PC2"))

# Afficher le graphique de la variance expliquée
print(variance_bar)


# Calculer les eigenvalues pour l'éboulis des eigenvalues
eigenvalues <- pca_result$sdev^2

# Créer un graphique pour l'éboulis des eigenvalues
eigenvalues_plot <- ggplot(NULL, aes(x = 1:length(eigenvalues), y = eigenvalues)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(title = "Scree Plot of Eigenvalues",
       x = "Principal Component",
       y = "Eigenvalue") +
  scale_x_continuous(breaks = 1:length(eigenvalues))

# Afficher le graphique de l'éboulis des eigenvalues
print(eigenvalues_plot)



# Supposons que vous ayez déjà chargé votre objet phyloseq sous le nom physeq_clr

# Effectuer une PCA avec prcomp
pca_result <- prcomp(otu_table(physeq_clr), scale. = TRUE)

# Calculer la variance expliquée par chaque composante principale
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

# Créer un dataframe pour le plot
pca_data <- as.data.frame(pca_result$x)
pca_data$type_enrich <- sample_data(physeq_clr)$type_enrich

# Créer un plot de la PCA en utilisant ggplot2
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = type_enrich)) +
  geom_point() +
  labs(title = "PCA plot with Explained Variance", x = "PC1", y = "PC2") +
  theme_minimal()

# Afficher le plot
print(pca_plot)

# Afficher la variance expliquée
print(variance_explained)






########


# Charger les librairies
library(phyloseq)
library(composition)
library(ggplot2)

# Supposons que vous ayez déjà chargé votre objet phyloseq sous le nom physeq_clr

# Séparer les échantillons en fonction de la variable "type_enrich"
sample_data( )$type_enrich <- as.factor(sample_data(physeq_clr)$type_enrich)

# Effectuer une PCA avec princomp.acomp
pca_result <- princomp(physeq_clr_asv)

# Calculer la variance expliquée par chaque composante principale
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

# Créer un dataframe pour le plot
pca_data <- data.frame(pca_result$scores)
pca_data$type_enrich <- sample_data(physeq_clr)$type_enrich

# Créer un plot de la PCA en utilisant ggplot2
pca_plot <- ggplot(pca_data, aes(x = Comp.1, y = Comp.2, color = type_enrich)) +
  geom_point() +
  labs(title = "PCA plot of Aitchison Transformed Data", x = "PC1", y = "PC2")

# Afficher le plot
print(pca_plot)

# Afficher la variance expliquée
print(variance_explained)

# Créer un graphique pour l'éboulis des eigenvalues
eigenvalues_plot <- ggplot(NULL, aes(x = 1:length(pca_result$sdev), y = pca_result$sdev^2)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(title = "Scree Plot of Eigenvalues",
       x = "Principal Component",
       y = "Eigenvalue") +
  scale_x_continuous(breaks = 1:length(pca_result$sdev))

# Afficher le graphique de l'éboulis des eigenvalues
print(eigenvalues_plot)




#######





ordinate(physeq_clr, method="NMDS", distance="euclidien")






library(vegan)  # Cette librairie contient la fonction metaMDS pour la NMDS

# Supposons que vous ayez déjà chargé votre objet phyloseq sous le nom physeq_clr

# Calculer la matrice de distances euclidiennes entre les échantillons
dist_euclidean <- dist(otu_table(physeq_clr), method = "euclidean")

# Effectuer la NMDS en utilisant la distance euclidienne
nmds_result <- metaMDS(dist_euclidean)

# Visualiser le résultat de la NMDS en utilisant un plot
plot(nmds_result, type = "n")












# Effectuer une PCA sur les données de distance euclidienne
p <- PCAtools::pca(physeq_clr_asv, metadata =meta_RDA)
PCAtools::biplot(
  p,
  # loadings parameters
  showLoadings = TRUE,
  lengthLoadingsArrowsFactor = 1.5,
  sizeLoadingsNames = 3,
  colLoadingsNames = 'red4',
  ntopLoadings = 3,
  # other parameters
  #lab = p$metadata$Sample_ID,
  #colby = "type_enrich",
  hline = 0, vline = 0,
  legendPosition = "right"
)













gp <- ordinate(physeq_clr, "NMDS","euclidean")
gp <- ordinate(physeq_c
               lr, "PCoA","bray")

plot_ordination(physeq_clr, gp, type="samples", color="type_enrich", shape="type_enrich") +
  geom_point(size=2.5)+
  #scale_color_gradient(low = "blue", high = "red") +
  scale_shape_manual(values = c(15,24,18,8))+
  # new_scale_color()+
  # stat_ellipse(aes(color= type_enrich), show.legend = TRUE) +
  #scale_color_manual(values = c("biofilm_EDM"="#33CC33", "biofilm_ENR"="#99FF66", 
  #                              "filtre_EDM"="#3366FF", "filtre_ENR"="#0099CC"))+
  geom_text(aes(x=0.0,y=0.0, label = paste0("stress : ", round(gp$stress,3))), size= 4, family= "Arial", face= "plain" )+
  theme_bw()






####NMDS----
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

# North <- hull[hull$sample.Geo  == "North", ][chull(hull[hull$sample.Geo == 
#                                                                 "North", c("Axis.1", "Axis.2")]), ]  # hull values for North
# South <- hull[hull$sample.Geo == "South", ][chull(hull[hull$sample.Geo == 
#                                                                "South", c("Axis.1", "Axis.2")]), ]  # hull values for Jellyfishes  

# hull_data <- rbind(North, South)

# color <- rep("#a65628", length(hull_data$sample.Geo))
# color[hull_data$sample.Geo == "North"] <- "#1919ff"
# hull_data <- cbind(hull_data, color)


biofilm_SW_NMDS <- hull[hull$sample.type_enrich  == "biofilm_SW", ][chull(hull[hull$sample.type_enrich =="biofilm_SW", c("Axis.1", "Axis.2")]), ]  # hull values for biofilm_SW
biofilm_ENR_NMDS <- hull[hull$sample.type_enrich  == "biofilm_ENR", ][chull(hull[hull$sample.type_enrich =="biofilm_ENR", c("Axis.1", "Axis.2")]), ]  # hull values for biofilm_SW
filter_SW_NMDS <- hull[hull$sample.type_enrich  == "filter_SW", ][chull(hull[hull$sample.type_enrich =="filter_SW", c("Axis.1", "Axis.2")]), ]  # hull values for filter_SW
filter_ENR_NMDS <- hull[hull$sample.type_enrich  == "filter_ENR", ][chull(hull[hull$sample.type_enrich =="filter_ENR", c("Axis.1", "Axis.2")]), ]  # hull values for filter_SW


#Vector of color for hulls
hull_col <- c("#33CC33", "#99FF66", "#0099CC", "#33CCFF")
names(hull_col) <- c("biofilm_SW" , "biofilm_ENR", "filter_SW", "filter_ENR")



hull_data <- hull %>%
  dplyr::group_by(sample.type_enrich) %>%
  dplyr::slice(chull(Axis.1,Axis.2)) %>%
  dplyr::mutate(color = hull_col[sample.type_enrich])

#pdf(file="NMDS_Aitchison.pdf", wi = 7, he = 7)
ggplot(hull,aes(x = Axis.1, y = Axis.2)) +
  geom_hline(yintercept = 0, colour = "lightgrey", linetype = 2) + 
  geom_vline(xintercept = 0, colour = "lightgrey", linetype = 2) +
  geom_polygon(data = hull_data,
               aes(group = sample.type_enrich,
                   fill = sample.type_enrich),
               alpha = 0.3) + # add the convex hulls)
  scale_fill_manual(values =  c("#33CC33", "#99FF66", "#0099CC", "#33CCFF")) +
  geom_point(data = hull,
             aes(color = sample.type_enrich,
                 size = sample.date),
             alpha = 0.7) +
  scale_color_manual(values =  c("#33CC33", "#99FF66", "#0099CC", "#33CCFF")) +
  geom_text(data = hull_data,
            x = -50, y = 50,
            label = paste("Stress =", round(physeq_clr_nmds$stress, 2)),
            colour = "Black",
            size = 4)  +
  xlab(paste("MDS1")) +
  ylab(paste("MDS2")) +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=14), # remove x-axis labels
        axis.title.y = element_text(size=14), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

# Correlation with environmental data
data.frame(names(hull))


env <- hull

# The function envfit will add the environmental variables as vectors to the ordination plot
ef <- vegan::envfit(physeq_clr_nmds, env, permu = 1000,na.rm=TRUE)
ef


# The two last columns are of interest: the squared correlation coefficient and the associated p-value
# Plot the vectors of the significant correlations and interpret the plot
plot(physeq_clr_nmds, type = "t", display = "sites")
plot(ef, p.max = 0.05)




######PCA----
#6.1.1 Number of PCs to retain
#First, we will use a scree plot to examine the proportion of total variation explained by each PC.

#prepare the ASV table to add taxonomy
tax_CLR <-  as.data.frame(tax_table(physeq_clr)) #get taxnomic table
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


#Plotting the PCA
PCAtools::biplot(
  p,
  lab = p$metadata$Sample_ID,
  colby = "type_enrich",
  pointSize = 5,
  hline = 0, vline = 0,
  legendPosition = "right"
)

#6.1.3 Determine the variables that drive variation among each PC
PCAtools::biplot(
  p, 
  # # loadings parameters
  # showLoadings = TRUE,
  # lengthLoadingsArrowsFactor = 1.5,
  # sizeLoadingsNames = 3,
  # colLoadingsNames = 'red4',
  # ntopLoadings = 3,
  # # other parameters
  # lab = p$metadata$Sample_ID,
  # colby = "type_enrich",
  # hline = 0, vline = 0,
  # legendPosition = "right"
)

ggplot(data = p, x= p$rotating[,1], y=p$rotating[,2] )

#6.1.4 Correlate the principal components back to environmental data
PCAtools::eigencorplot(
  p,
  components = PCAtools::getComponents(p, 1:horn$n),
  metavars = c("silicates","nitrites","phosphates","nitrates","ammonium","Phaeo","t_ext","T.bassins","pH"),
  col = c('white', 'cornsilk1', 'gold',
                 'forestgreen', 'darkgreen'),
                 cexCorval = 1.2,
  fontCorval = 2,
  posLab = "all",
  rotLabX = 45,
  scale = TRUE,
  main = bquote(PC ~ Spearman ~ r^2 ~ environmental ~ correlates),
  plotRsquared = TRUE,
  corFUN = "spearman",
  corUSE = "pairwise.complete.obs",
  corMultipleTestCorrection = 'BH',
  signifSymbols = c("****", "***", "**", "*", ""),
  signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
)






######PCOA 

physeq_rar_jaccard <- phyloseq::distance(physeq_16S_nit_p,
                                         method = "jaccard",
                                         binary = TRUE)

# trick to avoid negative egein values in PCoA
# it recreates what ade4::dist.binary() does
physeq_rar_jaccard <- sqrt(physeq_rar_jaccard)


#Unweigted unifrac 
ape::is.rooted(physeq_16S_nit_p@phy_tree)

#phy_tree(physeq_rar) <- phangorn::midpoint(physeq_rar@phy_tree)
unifracs <- GUniFrac::GUniFrac(physeq_16S_nit_p@otu_table@.Data, physeq_16S_nit_p@phy_tree, alpha=c(0, 0.5, 1))$unifracs
physeq_rar_du <- unifracs[, , "d_UW"]   # Unweighted UniFrac

#Weigted unifrac 
physeq_rar_dw <- unifracs[, , "d_1"]   # Weighted UniFrac


#bray curtis
# physeq_rar_bray <- vegan::vegdist(physeq_rar@otu_table@.Data, method = "bray")
tmp <- transform_sample_counts(physeq_16S_nit_p,function(x) {x/sum(x)} )
physeq_rar_bray <- phyloseq::distance(tmp, method = "bray")


#Visualisation
dist_methods <- unlist(distanceMethodList)
data.frame(position = seq_along(dist_methods),
           dist_methods)

#Select the distances of interest
dist_methods <- dist_methods[c(1, 2, 10, 8)]
dist_methods

#Loop through each distance method, save each plot to a list, called plist.
plist <- vector("list")

for(i in dist_methods){
  # Calculate distance matrix
  iDist <- phyloseq::distance(physeq_16S_nit_p, method = i)
  # Calculate PCoA ordination
  iMDS <- ordinate(physeq_16S_nit_p, "MDS", distance = iDist)
  ## Make plot. Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(physeq_16S_nit_p, iMDS, color= "Geo")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to list
  plist[[i]] = p 
}

df <- plyr::ldply(plist, function(x) x$data)
head(df)

names(df)[1] <- "distance"

ggplot(df, aes(Axis.1, Axis.2, color = type_enrich)) +
  geom_point(size=3, alpha=0.5) +
  theme_bw() +
  facet_wrap(~distance, scales="free") +
  ggtitle("PCoA (MDS) on various distance metrics")
