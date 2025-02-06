#distance matrix calculation
physeq_clr_dist <- phyloseq::distance(physeq_clr, method = "euclidean")


#Different clustering methods----
####See diffent clustering obtained with the four aggregation criterion 

#Simple aggregation criterion
spe_single <- hclust(physeq_clr_dist, method = "single")

#Complete aggregation criterion
spe_complete <- hclust(physeq_clr_dist, method = "complete")

#Unweighted pair group method with arithmetic mean
spe_upgma <- hclust(physeq_clr_dist, method = "average")

#Ward criterion
spe_ward <- hclust(physeq_clr_dist, method = "ward.D")
spe_ward$labels <- sample_data(physeq_16S_nit)$Sample_ID_E


par(mfrow = c(2, 2))
plot(spe_single, main = "single")
plot(spe_complete, main = "complete")
plot(spe_upgma, main = "UPGMA")
plot(spe_ward, main = "ward")



#Cophenic correlation----
#to choose a method with the laximisation on distance 

#Cophenetic correlation
spe_single_coph <- cophenetic(spe_single)
cor(physeq_clr_dist, spe_single_coph)
spe_complete_coph <- cophenetic(spe_complete)
cor(physeq_clr_dist, spe_complete_coph)
spe_upgma_coph <- cophenetic(spe_upgma)
cor(physeq_clr_dist, spe_upgma_coph)
spe_ward_coph <- cophenetic(spe_ward)
cor(physeq_clr_dist, spe_ward_coph)


#let see graphically whats is it 
plot_coph_cor <- function(cophenetic_distance, hclust_type){
  
  # first calculate the correlation between
  # the cophenetic distance and the observed distance
  cor_res <- round(cor(physeq_clr_dist, cophenetic_distance),3)
  
  # generate a scatter plot to visualise
  # the relationship
  plot(x = physeq_clr_dist,
       y = cophenetic_distance,
       xlab = "Aitchison distance",
       ylab = "Cophenetic distance",
       xlim = c(10, 35), ylim = c(10, 35),
       main = c(hclust_type, paste("Cophenetic correlation ", cor_res)))
  abline(0, 1)
}

par(mfrow=c(2,2))

plot_coph_cor(cophenetic_distance = spe_complete_coph,
              hclust_type = "Single linkage")

plot_coph_cor(cophenetic_distance = spe_complete_coph,
              hclust_type = "Complete linkage")

plot_coph_cor(cophenetic_distance = spe_upgma_coph,
              hclust_type = "Average linkage")

plot_coph_cor(cophenetic_distance = spe_ward_coph,
              hclust_type = "Ward linkage")


#It seems clear that the UPGMA method give the most faithful representation of original distances.


#Cut clusters ----
#Looking for interpretable clusters
#Fusion level plot
par(mfrow = c(1, 1))

plot(x = spe_complete$height,
     y = phyloseq::nsamples(physeq_clr):2,
     type = "S",
     main = "Fusion levels - Aitchison - Average",
     ylab = "k (number of cluster)",
     xlab = "h (node height)")

text(x = spe_upgma$height,
     y = phyloseq::nsamples(physeq_clr):2,
     labels = phyloseq::nsamples(physeq_clr):2,
     col = "red",
     cex = 0.8)



install.packages("NbClust", lib = ".")
library("NbClust", lib.loc = ".")
nclust <- nb_clust_all(data = t(physeq_clr_asv), seed = 1000)



#suivant le nombre de cluster déterminéson coupe 


#Cut the dendrogram in order to obtain K groups and compare their compositionC
k <- 5 # Number of groups given by the fusion level plot

#Cut the dendrogram
spe_upgma_clust <- cutree(tree = spe_ward, k = k)
table(spe_upgma_clust)
spe_upgma_clust2 <- data.frame(UPGMA_clusters = spe_upgma_clust)
#x11()
# Plot dendrogram with group labels
plot(spe_ward,
     hang = -1,
     ylab = "Height",
     main="Aitchison distance - Ward")

rect.hclust(spe_ward,
            k = k,
            border = 2:6,
            cluster = spe_upgma_clust)

legend("topright",
       paste("Cluster", 1:k),
       pch = 22,
       col = 2:(k + 1),
       bty = "n")


#duun index 
cs <- fpc::cluster.stats(d = physeq_clr_dist,
                         clustering = spe_upgma_clust)

cs$dunn



#heigt of cut----
# Utilisez cutree pour obtenir les groupes
spe_upgma_clust <- cutree(tree = spe_ward, k = 5)

# Nombre total d'échantillons (ou de feuilles dans le dendrogramme)
n <- length(spe_upgma_clust)

# Trouvez la hauteur de la fusion qui donne k = 5 groupes
# Note : la longueur de spe_ward$height est (n - 1)
# La hauteur à laquelle le dendrogramme est coupé pour obtenir k groupes est stockée à l'indice (n - k) dans spe_ward$height
h_value <- spe_ward$height[n - 5]

# Affichez la valeur de la hauteur
print(paste("La hauteur à laquelle le dendrogramme est coupé pour obtenir 5 groupes est :", h_value))






#chosiir les cluster 
#Cophenetic correlation
spe_single_coph <- cophenetic(spe_single)
cor(physeq_clr_dist, spe_single_coph)
spe_complete_coph <- cophenetic(spe_complete)
cor(physeq_clr_dist, spe_complete_coph)
spe_upgma_coph <- cophenetic(spe_upgma)
cor(physeq_clr_dist, spe_upgma_coph)
spe_ward_coph <- cophenetic(spe_ward)
cor(physeq_clr_dist, spe_ward_coph)


plot_coph_cor <- function(cophenetic_distance, hclust_type){
  
  # first calculate the correlation between
  # the cophenetic distance and the observed distance
  cor_res <- round(cor(physeq_clr_dist, cophenetic_distance),3)
  
  # generate a scatter plot to visualise
  # the relationship
  plot(x = physeq_clr_dist,
       y = cophenetic_distance,
       xlab = "Aitchison distance",
       ylab = "Cophenetic distance",
       xlim = c(10, 35), ylim = c(10, 35),
       main = c(hclust_type, paste("Cophenetic correlation ", cor_res)))
  abline(0, 1)
}

par(mfrow=c(2,2))

plot_coph_cor(cophenetic_distance = spe_complete_coph,
              hclust_type = "Single linkage")

plot_coph_cor(cophenetic_distance = spe_complete_coph,
              hclust_type = "Complete linkage")

plot_coph_cor(cophenetic_distance = spe_upgma_coph,
              hclust_type = "Average linkage")

plot_coph_cor(cophenetic_distance = spe_ward_coph,
              hclust_type = "Ward linkage")




#Fusion level plot
par(mfrow = c(1, 1))

plot(x = spe_upgma$height,
     y = phyloseq::nsamples(physeq_clr):2,
     type = "S",
     main = "Fusion levels - Aitchison - Average",
     ylab = "k (number of cluster)",
     xlab = "h (node height)")

text(x = spe_upgma$height,
     y = phyloseq::nsamples(physeq_clr):2,
     labels = phyloseq::nsamples(physeq_clr):2,
     col = "red",
     cex = 0.8)


install.packages("NbClust", lib = ".")
library("NbClust", lib.loc = ".")
nclust <- NbClust(data = t(physeq_clr_asv), seed = 1000)
nclust <- NbClust(data = t(physeq_clr_asv),method = "ward.D", distance = 'euclidean', index = "all")
