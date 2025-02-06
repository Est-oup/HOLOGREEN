#Picrust analysis 

#loadingss
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("biomformat")

library(phyloseq)
library(Biostrings)
library(biomformat)
library(KEGGREST)



phylobiom <- physeq_16S_nit
phylobiom@otu_table <- phylobiom@otu_table[,colSums(phylobiom@otu_table) != 0]

#Write a fasta file
writeXStringSet(refseq(phylobiom), filepath = "F:/16S_HOLOGREEN/core_microbiome/physeq_16S_nit.fasta", format = "fasta")



#Write a biom file 
phylobiom <- as.matrix(otu_table(phylobiom))
write_biom(make_biom(phylobiom), "F:/16S_HOLOGREEN/core_microbiome/physeq_16S_nit.biom")
list.files(pattern = "*.biom")




#Output analysis----

View(read.table(file="F:/16S_HOLOGREEN/Picrust/Objets/picrust2_out_pipeline/EC_predicted.tsv"))
View(read.table(file="F:/16S_HOLOGREEN/Picrust/Objets/picrust2_out_pipeline/KO_predicted.tsv"))

#EC code pathway analysis
EC_tab <- read.table(file="F:/16S_HOLOGREEN/Picrust/Objets/picrust2_out_pipeline/EC_predicted.tsv")
colnames(EC_tab) <- EC_tab[1,]
rownames(EC_tab) <- EC_tab[,1]
EC_tab <- data.frame(EC_tab[-1,-1])


#récuperer les description des codes EC 
EC_description <- list()
for (i in 1:length(substring(colnames(EC_tab), 4))){
  EC_description[[i]] <- keggGet((substring(colnames(EC_tab),4)[-1])[i])
}


#ETduier fonction par fonction 

#NIFH 
View(prune_taxa( rownames(EC_tab[EC_tab$'EC.1.18.6.1' == 1, ]), physeq_16S_nit_p )@otu_table)

# Fixation de l'azote (Nitrogénase)
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.18.6.1' == 1, ]), ]

# Assimilation de l'ammoniac
# Glutamine synthétase
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.6.3.1.2' == 1, ]), ]
# Glutamate synthase (NADH)
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.4.1.13' == 1, ]), ]
# Glutamate synthase (NADPH)
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.4.1.14' == 1, ]), ]

# Nitrification
# Ammoniac monooxygénase
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.14.99.39' == 1, ]), ]
# Hydroxylamine oxydoréductase
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.7.2.6' == 1, ]), ]
# Nitrite oxydoréductase
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.7.2.1' == 1, ]), ]

# Dénitrification
# Nitrate réductase (anaérobie)
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.7.5.1' == 1, ]), ]
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.7.99.4' == 1, ]), ]
# Nitrite réductase (cytochrome)
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.7.2.2' == 1, ]), ]
# Nitrite réductase (ammoniaque)
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.7.1.4' == 1, ]), ]
# Oxyde nitreux réductase
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.7.2.4' == 1, ]), ]

# Uréase
diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.3.5.1.5' == 1, ]), ]



#KO code function analysis----
KO_tab <- read.table(file="F:/16S_HOLOGREEN/Picrust/Objets/picrust2_out_pipeline/KO_predicted.tsv")
colnames(KO_tab) <- KO_tab[1,]
rownames(KO_tab) <- KO_tab[,1]
KO_tab <- data.frame(KO_tab[-1,-1])



#Nitrogenase
View(prune_taxa( rownames(KO_tab[KO_tab$"K02588" == 1, ]), physeq_16S_nit_p )@otu_table)
rowSums(prune_taxa( rownames(KO_tab[KO_tab$"K02588" == 1, ]), physeq_16S_nit_p )@otu_table)
View(prune_taxa( rownames(KO_tab[KO_tab$"K02591" == 1, ]), physeq_16S_nit_p )@otu_table)





#presence in biofilm

KO_ASV_biofilm <- data.frame()
#Fixation de l'azote
  #nitrogenase 
KO_ASV_biofilm$fix <- row.names(diff_biofilm[row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K02588" == 1, ])==1,])
KO_ASV_biofilm$fix <- unique(c(KO_ASV_biofilm,row.names(diff_biofilm[row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K02591" == 1, ])==1,])))

#Dénitrification
  #nitrate reductase 
KO_ASV_biofilm$denit <- rownames(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K00370" == 1, ])==1,])
KO_ASV_biofilm$denit <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K00371" == 1, ])==1,])))
KO_ASV_biofilm$denit <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K00374" == 1, ])==1,])))
  #nitrite reductase
KO_ASV_biofilm$denit <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K03385" == 1, ])==1,])))
KO_ASV_biofilm$denit <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K00368" == 1, ])==1,])))
  #oxyde nitreux reductase 
KO_ASV_biofilm$oxidenitreuxred <- rownames(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K00376" == 1, ])==1,])

#Assimilation de l'ammoniac
  #Glutamine synthétase
KO_ASV_biofilm$assimilammonia <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K01915" == 1, ])==1,])))
  #Glutamate syntéthase
KO_ASV_biofilm$assimilammonia <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K00265" == 1, ])==1,])))
KO_ASV_biofilm$assimilammonia <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K00284" == 1, ])==1,])))

#Nitrification
  #Ammoniac monooxygenasae
KO_ASV_biofilm$nitrit <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K10944" == 1, ])==1,])))
KO_ASV_biofilm$nitrit <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K10945" == 1, ])==1,])))
KO_ASV_biofilm$nitrit <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K10946" == 1, ])==1,])))

  #Hydroxylamine oxydoréductase 
KO_ASV_biofilm$nitrit <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K10535" == 1, ])==1,])))

  #Nitrite oxydoréductase
KO_ASV_biofilm$nitrit <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K00362" == 1, ])==1,])))

#Uréase 
  #uréase
KO_ASV_biofilm$uree <- unique(c(KO_ASV_biofilm, row.names(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K01428" == 1, ])==1,])))




sum(prune_taxa(rownames(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K02588" == 1, ])==1,]), physeq_16S_nit_p@otu_table))
sum(prune_taxa(rownames(diff_biofilm [row.names(diff_biofilm ) %in%  rownames(KO_tab[KO_tab$"K02591" == 1, ])==1,]), physeq_16S_nit_p@otu_table))
sum(prune_taxa(rownames(diff_biofilm [row.names(diff_biofilm ) %in% rownames(EC_tab[EC_tab$'EC.1.18.6.1' == 1, ]), ]), physeq_16S_nit_p@otu_table))




#Liste récupération ASV----

##KO list-----
# Initialiser une liste principale
KO_ASV_biofilm <- list()

# Fixation de l'azote (Nitrogénase)
KO_ASV_biofilm$fixation_azote <- unique(c(
  row.names(t(t(diff_biofilm@otu_table)@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K02588" == 1, ]) == 1, ]),
  row.names(diff_biofilm[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K02591" == 1, ]) == 1, ])
))

# Dénitrification
KO_ASV_biofilm$denitrification <- unique(c(
  # Nitrate réductase
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K00370" == 1, ]) == 1, ]),
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K00371" == 1, ]) == 1, ]),
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K00374" == 1, ]) == 1, ]),
  # Nitrite réductase
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K03385" == 1, ]) == 1, ]),
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K00368" == 1, ]) == 1, ]),
  # Oxyde nitreux réductase
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K00376" == 1, ]) == 1, ])
))

# Assimilation de l'ammoniac
KO_ASV_biofilm$assimilation_ammoniac <- unique(c(
  # Glutamine synthétase
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(KO_tab[KO_tab$"K01915" == 1, ]) == 1, ]),
  # Glutamate synthétase
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(KO_tab[KO_tab$"K00265" == 1, ]) == 1, ]),
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(KO_tab[KO_tab$"K00284" == 1, ]) == 1, ])
))

# Nitrification
KO_ASV_biofilm$nitrification <- unique(c(
  # Ammoniac monooxygenase
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(KO_tab[KO_tab$"K10944" == 1, ]) == 1, ]),
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(KO_tab[KO_tab$"K10945" == 1, ]) == 1, ]),
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(KO_tab[KO_tab$"K10946" == 1, ]) == 1, ]),
  # Hydroxylamine oxydoréductase
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(KO_tab[KO_tab$"K10535" == 1, ]) == 1, ]),
  # Nitrite oxydoréductase
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(KO_tab[KO_tab$"K00362" == 1, ]) == 1, ])
))

# Uréase
KO_ASV_biofilm$urease <- unique(
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(KO_tab[KO_tab$"K01428" == 1, ]) == 1, ])
)

# Vous pouvez continuer de cette manière pour ajouter d'autres processus métaboliques à votre liste



##EC list----
# Initialiser une liste principale
EC_ASV_biofilm <- list()

# Fixation de l'azote (Nitrogénase)
EC_ASV_biofilm$fixation_azote <- unique(
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.18.6.1' == 1, ]) == 1, ])
)

# Assimilation de l'ammoniac
# Glutamine synthétase
EC_ASV_biofilm$assimilation_ammoniac <- unique(
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.6.3.1.2' == 1, ]) == 1, ])
)
# Glutamate synthase (NADH)
EC_ASV_biofilm$assimilation_ammoniac <- unique(c(
  EC_ASV_biofilm$assimilation_ammoniac,
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.4.1.13' == 1, ]) == 1, ])
))
# Glutamate synthase (NADPH)
EC_ASV_biofilm$assimilation_ammoniac <- unique(c(
  EC_ASV_biofilm$assimilation_ammoniac,
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.4.1.14' == 1, ]) == 1, ])
))

# Nitrification
# Ammoniac monooxygénase
EC_ASV_biofilm$nitrification <- unique(
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.14.99.39' == 1, ]) == 1, ])
)
# Hydroxylamine oxydoréductase
EC_ASV_biofilm$nitrification <- unique(c(
  EC_ASV_biofilm$nitrification,
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.7.2.6' == 1, ]) == 1, ])
))
# Nitrite oxydoréductase
EC_ASV_biofilm$nitrification <- unique(c(
  EC_ASV_biofilm$nitrification,
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.7.2.1' == 1, ]) == 1, ])
))

# Dénitrification
# Nitrate réductase (anaérobie)
EC_ASV_biofilm$denitrification <- unique(c(
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.7.5.1' == 1, ]) == 1, ]),
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.7.99.4' == 1, ]) == 1, ])
))
# Nitrite réductase (cytochrome et ammoniaque)
EC_ASV_biofilm$denitrification <- unique(c(
  EC_ASV_biofilm$denitrification,
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.7.2.2' == 1, ]) == 1, ]),
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.7.1.4' == 1, ]) == 1, ])
))
# Oxyde nitreux réductase
EC_ASV_biofilm$denitrification <- unique(c(
  EC_ASV_biofilm$denitrification,
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.1.7.2.4' == 1, ]) == 1, ])
))

# Uréase
EC_ASV_biofilm$urease <- unique(
  row.names(diff_biofilm[row.names(diff_biofilm) %in% rownames(EC_tab[EC_tab$'EC.3.5.1.5' == 1, ]) == 1, ])
)

# Vous pouvez continuer à ajouter d'autres processus métaboliques à votre liste



#Expression of metabolism 


colMeans(prune_taxa(KO_ASV_biofilm$fixation_azote ,otu_table(prune_taxa(row.names(test5), subset_samples(physeq_16S_nit_p, type=="biofilm")))))
colMeans(prune_taxa(EC_ASV_biofilm$denitrification ,otu_table(prune_taxa(row.names(test5), subset_samples(physeq_16S_nit_p, type=="biofilm")))))


phy_denitri_biof <- prune_taxa(KO_ASV_biofilm$denitrification ,prune_taxa(row.names(test5), subset_samples(physeq_16S_nit_p, type=="biofilm")))
phy_fixa_biof <- prune_taxa(KO_ASV_biofilm$fixation_azote ,prune_taxa(row.names(test5), subset_samples(physeq_16S_nit_p, type=="biofilm")))
phy_ammo_biof <- prune_taxa(KO_ASV_biofilm$assimilation_ammoniac ,prune_taxa(row.names(test5), subset_samples(physeq_16S_nit_p, type=="biofilm")))
phy_nitri_biof <- prune_taxa(KO_ASV_biofilm$nitrification ,prune_taxa(row.names(test5), subset_samples(physeq_16S_nit_p, type=="biofilm")))
phy_uree_biof <- prune_taxa(KO_ASV_biofilm$urease ,prune_taxa(row.names(test5), subset_samples(physeq_16S_nit_p, type=="biofilm")))


ggplot(psmelt(phy_uree_biof), aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Abondance des ASV dans le Biofilm pour la fixation N2",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))








##############################################
# Préparation des données
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ALDEx2)
library(pheatmap)
library(ANCOMBC)


KO_tab2 <- data.frame(lapply(KO_tab, as.numeric), row.names = rownames(KO_tab))
KO_tab2 <- KO_tab2[, colSums(KO_tab2) != 0] # Supprimer les colonnes vides

# Créer une liste de OTU tables pour chaque KO term
ko_otu_list <- list()
for (ko_term in colnames(KO_tab2)) {
  ko_otu_list[[ko_term]] <- prune_taxa(rownames(KO_tab2[KO_tab2[, ko_term, drop = FALSE] >= 1,]), physeq_16S_nit_p) %>% otu_table() %>% rowSums()
}

# Créer une liste de OTU tables pour chaque KO term
ko_otu_list_raw <- list()
for (ko_term in colnames(KO_tab2)) {
  ko_otu_list_raw[[ko_term]] <- prune_taxa(rownames(KO_tab2[KO_tab2[, ko_term, drop = FALSE] >= 1,]), physeq_16S_nit) %>% otu_table() %>% rowSums()
}

# Créer un objet phyloseq combiné
combined_otu <- t(do.call(rbind, ko_otu_list))
combined_otu_raw <- t(do.call(rbind, ko_otu_list_raw))
combined_otu2 <- t(do.call(rbind, ko_otu_list))


# donne les noms de fonction 
for (i in 1:length(colnames(KO_tab))){
  KO_description[[i]] <- tryCatch({
    keggGet(colnames(KO_tab)[i])
  }, error = function(e) NA)
}




# Parcourir chaque colonne de combined_otu
for (i in 1:ncol(combined_otu2)) {
  # Parcourir chaque élément de KO_description
  for (ko_item in KO_description) {
    if (!is.null(ko_item) && length(ko_item[[1]]) >= 3) {
      # Vérifier si le nom de la colonne correspond au nom dans KO_description
      if (colnames(combined_otu2)[i] == ko_item[[1]][[1]]) {
        # Mettre à jour le nom de la colonne
        colnames(combined_otu2)[i] <- paste0(colnames(combined_otu2)[i], " ", ko_item[[1]][[5]])
        break # Sortir de la boucle interne dès qu'une correspondance est trouvée
      }
    }
  }
}



# Parcourir chaque colonne de combined_otu
for (i in 1:ncol(combined_otu_raw)) {
  # Parcourir chaque élément de KO_description
  for (ko_item in KO_description) {
    if (!is.null(ko_item) && length(ko_item[[1]]) >= 3) {
      # Vérifier si le nom de la colonne correspond au nom dans KO_description
      if (colnames(combined_otu_raw)[i] == ko_item[[1]][[1]]) {
        # Mettre à jour le nom de la colonne
        colnames(combined_otu_raw)[i] <- paste0(colnames(combined_otu_raw)[i], " ", ko_item[[1]][[3]])
        break # Sortir de la boucle interne dès qu'une correspondance est trouvée
      }
    }
  }
}


phyloseq_metabo <- phyloseq(otu_table(combined_otu2, taxa_are_rows = FALSE), sample_data(physeq_16S_nit))
tax_table(phyloseq_metabo) <- 
  tax_table(matrix("NA", nrow = ncol(otu_table(phyloseq_metabo)), ncol = 7, dimnames = list(colnames(otu_table(phyloseq_metabo)), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))))


phyloseq_metabo_raw <- phyloseq(otu_table(combined_otu_raw, taxa_are_rows = FALSE), sample_data(physeq_16S_nit))
tax_table(phyloseq_metabo_raw) <- 
  tax_table(matrix("NA", nrow = ncol(otu_table(phyloseq_metabo_raw)), ncol = 7, dimnames = list(colnames(otu_table(phyloseq_metabo_raw)), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))))

# Analyse différentielle avec DESeq2
dds <- phyloseq_to_deseq2(phyloseq_metabo, ~ enrich)
dds <- DESeq(dds)
res_deseq2 <- results(dds)

# Afficher les résultats DESeq2
head(res_deseq2[order(res_deseq2$pvalue),])

# Analyse différentielle avec ALDEx2
aldex_res_pi <- aldex.clr(t(as(otu_table(phyloseq_metabo_raw), "matrix")), 
                          sample_data(phyloseq_metabo_raw)$enrich,
                          mc.samples=1000, denom="all", verbose=F)
aldex_kw_pi <- aldex.kw(aldex_res_pi)
write.table(aldex_kw_pi,"aldex_kw_pi",sep = "\t")    #je sauve le output au cas où il y ait un bug, pour ne pas avoir à refaire tourner
aldex_kw_pi_pval <- aldex_kw_pi[aldex_kw_pi$glm.eBH < 0.05,]  #garder seulement les ASV avec pvalue corrigée BH < 0.05

# Afficher les résultats ALDEx2
head(aldex_kw)

# Visualisation PCoA
pcoa_results <- ordinate(phyloseq_metabo, "PCoA", "bray")
pcoa_plot <- plot_ordination(phyloseq_metabo, pcoa_results, color = "enrich") + geom_point()
print(pcoa_plot)

# Heatmap des résultats DESeq2
top_taxa <- head(order(res_deseq2$pvalue), 20) # Sélectionner les 20 meilleurs taxons
data_heatmap <- otu_table(phyloseq_metabo)[,top_taxa ]
pheatmap(log2(data_heatmap + 1), cluster_rows = TRUE, cluster_cols = TRUE)

# Barplot des abondances
barplot_data <- as.data.frame(t(otu_table(phyloseq_metabo)))
ggplot(psmelt(barplot_data), aes(x = variable, y = value, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Abondance des ASV par KO term", x = "KO term", y = "Abondance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# Extraire et visualiser les résultats


#LEFSE
metabo_lefse <- microbiomeMarker::run_lefse(phyloseq_metabo, norm = "CPM",
                                            wilcoxon_cutoff = 0.01,
                                            group = "enrich",
                                            taxa_rank = "none",
                                            kw_cutoff = 0.01,
                                            multigrp_strat = TRUE,
                                            lda_cutoff = 4)
metabo_lefse_tab <- data.frame(metabo_lefse@marker_table)
metabo_lefse_tab

#ancomBC
metabo_ancombc <- run_ancombc_patched(
  phyloseq_metabo,
  group = "enrich",
  taxa_rank = "none",
  pvalue_cutoff = 0.001,
  p_adjust = "fdr"
)

metabo_ancombc_tab <- data.frame(metabo_ancombc@marker_table)
(metabo_ancombc_tab %>% 
  arrange(desc(ef_W)))[1:20,"feature"]
(metabo_ancombc_tab %>% 
  arrange(-desc(ef_W)))[1:20,"feature"]

KO_ASV_biofilm




###Dissimilation nitrate reduction nitrate => ammonia

#Kegg numbers
#K00370 K00371 k00374 K02567 K02568 K00362 K00363 K03385 K15876

dissimilation_asv__pathway <- unique(c(rownames(KO_tab[KO_tab$"K00370" == 1, ]),
rownames(KO_tab[KO_tab$"K00371" == 1, ]),
rownames(KO_tab[KO_tab$"k00374" == 1, ]),
rownames(KO_tab[KO_tab$"K02567" == 1, ]),
rownames(KO_tab[KO_tab$"K02568" == 1, ])))#,
# rownames(KO_tab[KO_tab$"K00362" == 1, ]),
# rownames(KO_tab[KO_tab$"K00363" == 1, ]),
# rownames(KO_tab[KO_tab$"K03385" == 1, ]),
# rownames(KO_tab[KO_tab$"K15876" == 1, ])))


dissimilation_asv__pathway <- unique(c(rownames(KO_tab[KO_tab$"K00370" == 1, ]),
                                       rownames(KO_tab[KO_tab$"K00371" == 1, ]),
                                       rownames(KO_tab[KO_tab$"k00374" == 1, ])))

dissimilation_asv__pathway <- unique(c(rownames(KO_tab[KO_tab$"K02567" == 1, ]),
                                       rownames(KO_tab[KO_tab$"K02568" == 1, ])))

prune_taxa(dissimilation_asv__pathway, diff_biofilm)
ggplot(psmelt(prune_taxa(dissimilation_asv__pathway, diff_biofilm)), aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Abondance des ASV dans le Biofilm pour la fixation N2",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(psmelt(prune_taxa(dissimilation_asv__pathway, diff_water)), aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Abondance des ASV dans le Biofilm pour la fixation N2",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



###Assimilatory nitrate reduction nitrate => ammonia

#Kegg numbers
#K00367 K10534 k00372 K00360 K00366 K17877 K26139 K26138 K00361

assimilation_asv__pathway <- unique(c(rownames(KO_tab[KO_tab$"K00367" == 1, ]),
                                       rownames(KO_tab[KO_tab$"K10534" == 1, ]),
                                       rownames(KO_tab[KO_tab$"k00372" == 1, ]),
                                       rownames(KO_tab[KO_tab$"K00360" == 1, ]),
                                       rownames(KO_tab[KO_tab$"K00366" == 1, ]),
                                       rownames(KO_tab[KO_tab$"K17877" == 1, ]),
                                       rownames(KO_tab[KO_tab$"K26139" == 1, ]),
                                       rownames(KO_tab[KO_tab$"K26138" == 1, ]),
                                       rownames(KO_tab[KO_tab$"K00361" == 1, ])))

prune_taxa(assimilation_asv__pathway, diff_biofilm)
ggplot(psmelt(prune_taxa(assimilation_asv__pathway, diff_biofilm)), aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Abondance des ASV dans le Biofilm pour la fixation N2",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))







##NarGHI_NapAB
NarGHI_NapAB <- c("K00370","K00371","K00374","K02567","K02568")

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NarGHI_NapAB, drop=F]) >= 1, ]), diff_biofilm)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans le Biofilm pour NarGHI_NapAB",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



##NirBD_NrfAH
NirBD_NrfAH <- c("K00362","K00363","K03385","K15876")

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NirBD_NrfAH, drop=F]) >= 1, ]), diff_biofilm)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans le Biofilm pour NirBD_NrfAH",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##NirK_NirS
NirK_NirS <- c("K00368","K15864")

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NirK_NirS, drop=F]) >= 1, ]), diff_biofilm)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans le Biofilm pour NirK_NirS",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##NorBC
NorBC <- c("K04561","K02305") 

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NorBC, drop=F]) >= 1, ]), diff_biofilm)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans le Biofilm pour NorBC",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##NosZ
NosZ <- c("K00376")

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NosZ, drop=F]) >= 1, ]), diff_biofilm)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans le Biofilm pour NosZ",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##NifDKH_AnfG_VnfDKGH
NifDKH_AnfG_VnfDKGH <- c("K02586", "K02591", "K02588", "K00531")#, "K22896", "K22897","K22898","K22899")

  ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NifDKH_AnfG_VnfDKGH, drop=F]) >= 1, ]), diff_biofilm)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans le Biofilm pour NifDKH_AnfG_VnfDKGH",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##Denitrification 
c(NarGHI_NapAB, NirK_NirS, NorBC, NosZ)


ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, c(NarGHI_NapAB, NirK_NirS, NorBC, NosZ), drop=F]) >= 1, ]), diff_biofilm)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance des ASV dans le Biofilm pour la fixation N2",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))









test <- unique(c(
  # Nitrate réductase
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K00370" == 1, ]) == 1, ]),
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K00371" == 1, ]) == 1, ]),
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K00374" == 1, ]) == 1, ]),
  # Nitrite réductase
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K03385" == 1, ]) == 1, ]),
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K00368" == 1, ]) == 1, ]),
  # Oxyde nitreux réductase
  row.names(t(diff_biofilm@otu_table)[row.names(t(diff_biofilm@otu_table)) %in% rownames(KO_tab[KO_tab$"K00376" == 1, ]) == 1, ])
))




ggplot(psmelt(prune_taxa(test, physeq_16S_nit_p)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance des ASV dans le Biofilm pour la fixation N2",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
















##NarGHI_NapAB
NarGHI_NapAB <- c("K00370","K00371","K00374","K02567","K02568")

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NarGHI_NapAB, drop=F]) >= 1, ]), diff_water)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans l'eau pour NarGHI_NapAB",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



##NirBD_NrfAH
NirBD_NrfAH <- c("K00362","K00363","K03385","K15876")

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NirBD_NrfAH, drop=F]) >= 1, ]), diff_water)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans le l'eau pour NirBD_NrfAH",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##NirK_NirS
NirK_NirS <- c("K00368","K15864")

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NirK_NirS, drop=F]) >= 1, ]), diff_water)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans l'eau pour NirK_NirS",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##NorBC
NorBC <- c("K04561","K02305") 

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NorBC, drop=F]) >= 1, ]), diff_water)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans l'eau pour NorBC",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##NosZ
NosZ <- c("K00376")

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NosZ, drop=F]) >= 1, ]), diff_water)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans l'eau pour NosZ",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##NifDKH_AnfG_VnfDKGH
NifDKH_AnfG_VnfDKGH <- c("K02586", "K02591", "K02588", "K00531")#, "K22896", "K22897","K22898","K22899")

ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, NifDKH_AnfG_VnfDKGH, drop=F]) >= 1, ]), diff_water)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance (diff) des ASV dans l'eau pour NifDKH_AnfG_VnfDKGH",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##Denitrification 
c(NarGHI_NapAB, NirK_NirS, NorBC, NosZ)


ggplot(psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, c(NarGHI_NapAB, NirK_NirS, NorBC, NosZ), drop=F]) >= 1, ]), diff_water)),
       aes(x = Sample_ID, y=Abundance, fill = enrich)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~type)+
  theme_minimal() +
  labs(title = "Abondance des ASV dans le Biofilm pour la fixation N2",
       x = "Échantillon",
       y = "Abondance",
       fill = "ASV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##plot combiné
#Assurez-vous que les bibliothèques nécessaires sont chargées
library(ggplot2)
library(patchwork)
library(phyloseq)
library(cowplot)
library(dplyr)


# Définition des groupes de gènes
gene_groups <- list(
  NarGHI_NapAB = c("K00370","K00371","K00374","K02567","K02568"),
  NirBD_NrfAH = c("K00362","K00363","K03385","K15876"),
  NirK_NirS = c("K00368","K15864"),
  NorBC = c("K04561","K02305"),
  NosZ = c("K00376"),
  NarB_NR_NasAB = c("K00367","K10534","K00372","K00360"),
  NITsix_NirA_NasBDE = c("K17877","K00366","K26139", "K26138","K00361")
  # NifDKH_AnfG_VnfDKGH = c("K02586", "K02591", "K02588", "K00531")
)

# Vérifier que chaque terme KO est dans le tableau KO_tab2
# et créer un nouveau gene_groups sans les termes manquants
valid_gene_groups <- lapply(gene_groups, function(ko_terms) {
  ko_terms[ko_terms %in% colnames(KO_tab2)]
})

# Filtrer les groupes vides
valid_gene_groups <- Filter(function(g) length(g) > 0, valid_gene_groups)

# Afficher les groupes de gènes valides
print(valid_gene_groups)
# Supposons que 'diff_biofilm' et 'diff_water' sont des subsets définis de vos données
# Supposons aussi que 'Sample_ID', 'Abundance', et 'enrich' sont des colonnes dans vos données fondues (melted)

# Fonction modifiée pour inclure le nom du groupe de gènes dans le titre
create_plot <- function(gene_group, gene_group_name, data_subset, title_suffix) {
  data_filtered <- prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, gene_group, drop=F]) >= 1, ]), data_subset)
  psmelted_data <- psmelt(data_filtered)
  
  ggplot(psmelted_data, aes(x = Sample_ID, y = Abundance, fill = enrich)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~type) +
    theme_minimal() +
    labs(title = paste("Abondance (diff) des ASV pour", gene_group_name, "dans", title_suffix),
         x = "Échantillon",
         y = "Abondance",
         fill = "ASV") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 5))
}

# Puis appelez cette fonction en passant le nom du groupe de gènes comme argument

# Générez tous les graphiques pour le biofilm
plots_biofilm <- lapply(list("NarGHI_NapAB", "NirBD_NrfAH", "NirK_NirS", "NorBC", "NosZ", "NifDKH_AnfG_VnfDKGH"), function(gene_group_name) {
  create_plot(get(gene_group_name), gene_group_name, data_subset = subset_samples(physeq_16S_nit, type == "biofilm"), title_suffix = "dans le Biofilm")
})

# Générez tous les graphiques pour l'eau
plots_water <- lapply(list("NarGHI_NapAB", "NirBD_NrfAH", "NirK_NirS", "NorBC", "NosZ"), function(gene_group_name) {
  create_plot(get(gene_group_name), gene_group_name, data_subset = subset_samples(physeq_16S_nit, type == "water"), title_suffix = "dans l'eau")
})

#dans l'eau on enlève fixation N2
# Ensuite, créez des paires de graphiques en vérifiant le type d'objet
paired_plots <- lapply(seq_along(plots_biofilm), function(i) {
  biofilm_plot <- plots_biofilm[[i]]
  water_plot <- if (i <= length(plots_water) && inherits(plots_water[[i]], "ggplot")) {
    plots_water[[i]]
  } else {
    plot_spacer() # Utilisez un espace vide si le graphique eau correspondant n'existe pas
  }
  
  if (inherits(biofilm_plot, "ggplot")) {
    plot_grid(biofilm_plot, water_plot, ncol = 2)
  } else {
    plot_spacer() # Utilisez un espace vide si le graphique biofilm n'est pas un ggplot
  }
})

# Organisez ces paires en pages
for (i in seq(1, length(paired_plots), by = 2)) {
  page_plots <- list(paired_plots[[i]])
  if (i + 1 <= length(paired_plots)) {
    page_plots <- c(page_plots, list(paired_plots[[i + 1]]))
  }
  
  combined_page_plot <- do.call(plot_grid, c(page_plots, ncol = 1, nrow = length(page_plots)))
  ggsave(paste("combined_plots_page_raw", ceiling(i / 2), ".png"), combined_page_plot, bg = "white", width = 20, height = 15, units = "cm")
}






diff_abund_metabo_bp <- function(physeq_obj, sample_type) {
  #subset of the pyloseq object
  # Créer un masque logique pour filtrer les échantillons
  mask <- rownames(physeq_obj@sam_data[physeq_obj@sam_data$type %in% sample_type,])
  
  # Sous-ensemble de l'objet phyloseq en utilisant le masque
  physeq <- prune_samples(mask, physeq_obj)
  
  
  # Liste des groupes de gènes
  gene_groups <- list(
    NarGHI_NapAB = c("K00370","K00371","K00374","K02567","K02568"),
    NirBD_NrfAH = c("K00362","K00363","K03385","K15876"),
    NirK_NirS = c("K00368","K15864"),
    NorBC = c("K04561","K02305"),
    NosZ = c("K00376"),
    NarB = c("K00367"),
    NR = c("K10534"),
    NasAB = c("K00372","K00360"),
    NITsix = c("K00366"),
    NasBDE = c("K26139", "K26138")
    
    #NifDKH_AnfG_VnfDKGH = c("K02586", "K02591", "K02588", "K00531")
  )
  
  calc_diff_abundance <- function(gene_group, gene_group_name, physeq, sample_type) {
    # Filtrer les échantillons en fonction du type et du traitement
    data_enr <- psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, gene_group, drop=F]) >= 1, ]), 
                                  subset_samples(physeq,enrich == "ENR"))) %>%
      group_by(OTU, date) %>%
      summarise(Abundance = mean(Abundance))
    
    data_sw <- psmelt(prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, gene_group, drop=F]) >= 1, ]), 
                                 subset_samples(physeq,enrich == "SW"))) %>%
      group_by(OTU, date) %>%
      summarise(Abundance = mean(Abundance))
    
    joined_data <- left_join(data_enr, data_sw, by = c("OTU", "date"), suffix = c("_ENR", "_SW")) %>%
      mutate(Gene_Group = gene_group_name)
    
    return(joined_data)
  }
  
  # Puis, modifiez l'appel à cette fonction dans `all_diffs` :
  all_diffs <- lapply(names(gene_groups), function(gene_group_name) {
    calc_diff_abundance(gene_groups[[gene_group_name]], gene_group_name, physeq, sample_type)
  })
  
  # Combinez tous les résultats
  all_diffs_combined <- do.call(rbind, all_diffs) %>%
    gather(key = "enrich", value = "Abundance", Abundance_ENR, Abundance_SW)
  
  
  all_diffs_combined2 <-all_diffs_combined%>% 
    group_by(date,enrich, Gene_Group)%>%
    summarise(Abundance = sum(Abundance))
  
  return(all_diffs_combined2)
  
}



# Créez le graphique
ggplot(diff_abund_metabo_bp(physeq_16S_nit_p, "water"), aes(x = Gene_Group, y = Abundance, fill = enrich)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = paste("Différence d'abondance entre ENR et SW pour chaque voie métabolique dans l'eau"),
       x = "",
       y = "") +
  theme_minimal()


# Créez le graphique
ggplot(diff_abund_metabo_bp(physeq_16S_nit_p, "biofilm"),aes(x=Gene_Group, y = Abundance, fill=enrich)) +
  geom_boxplot() +
  labs(title = "Différence d'abondance entre ENR et SW pour chaque voie métabolique dans le biofilm",
       x = "Date",
       y = "Différence d'Abondance")+
  theme_minimal()










#Plot Article ----
# Fonction pour calculer la différence d'abondance
diff_abund_metabo_diff <- function(physeq_obj, sample_type) {
  # Créer un masque logique pour filtrer les échantillons
  mask <- rownames(physeq_obj@sam_data[physeq_obj@sam_data$type %in% sample_type,])
  
  # Sous-ensemble de l'objet phyloseq en utilisant le masque
  physeq_subset <- prune_samples(mask, physeq_obj)
  
  # Initialiser un dataframe pour stocker les résultats
  results <- data.frame(Gene_Group = character(), Abundance = numeric(), stringsAsFactors = FALSE)
  
  # Calcul de la différence d'abondance pour chaque groupe de gènes
  for (gene_group_name in names(valid_gene_groups)) {
    gene_group <- valid_gene_groups[[gene_group_name]]
    asv_subset <- prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, gene_group, drop = FALSE]) >= 1, ]), physeq_subset)
    psmelt <- psmelt(asv_subset)
    
    if(sample_type == "biofilm"){
      psmelt <- psmelt %>% group_by(OTU, enrich, date) %>% summarise(Abundance =mean(Abundance))
    }
    
    asv_data <- psmelt %>% group_by(date, enrich) %>% summarise(Abundance =sum(Abundance))
    
    # Calculer la différence d'abondance
    enr_abundance <- data.frame(subset(asv_data, enrich == "ENR"))
    sw_abundance <- data.frame(subset(asv_data, enrich == "SW"))
    diff_abundance <- merge(sw_abundance, enr_abundance, by = "date")
    diff_abundance$Diff_Abundance <- diff_abundance$Abundance.y - diff_abundance$Abundance.x
    diff_abundance <- diff_abundance[, c("date", "Abundance.x", "Abundance.y", "Diff_Abundance")]
    colnames(diff_abundance) <- c("Date", "Abundance_SW", "Abundance_ENR", "Diff_Abundance")
    # Ajouter les résultats au dataframe
    results <- rbind(results, data.frame(Gene_Group = gene_group_name, date_ID = diff_abundance$Date, 
                                         Abundance_SW = diff_abundance$Abundance_SW, Abundance_ENR = diff_abundance$Abundance_ENR, 
                                         Diff_Abundance = diff_abundance$Diff_Abundance))
  }
  
  return(results)
}



# Exemple d'utilisation de la fonction
biofilm_diff_abund <- diff_abund_metabo_diff(physeq_16S_nit_p, "biofilm")
water_diff_abund <- diff_abund_metabo_diff(physeq_16S_nit_p, "water")
v <- c("NarGHI_NapAB", "NirBD_NrfAH", "NirK_NirS","NorBC","NosZ","NarB_NR_NasAB","NITsix_NirA_NasBDE")
# Créez le graphique pour le biofilm
library(patchwork)

# Créez le graphique pour le biofilm
biofilm_plot <- ggplot(biofilm_diff_abund, aes(x = factor(Gene_Group, levels = v), y = Diff_Abundance)) +
  geom_hline(yintercept = 0, color = "red", linetype = 'dashed') +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, fill = "#66e54c", alpha= 0.7) + 
  geom_jitter(position = position_dodge(0.8), size = 1) +
  labs(title = "", x = "", y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "bold", size = 9))

# Créez le graphique pour l'eau
water_plot <- ggplot(water_diff_abund, aes(x = factor(Gene_Group, levels = v), y = Diff_Abundance)) +
  geom_hline(yintercept = 0, color = "red", linetype = 'dashed') +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, fill = "#19b2e5", alpha = 0.7) +  
  geom_jitter(position = position_dodge(0.8), size = 1) +
  labs(title = "", x = "", y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "bold", size = 9))

# Empiler les deux graphiques avec la même largeur
combined_plot <- biofilm_plot / water_plot

# Afficher le graphique combiné
combined_plot



#Table article----
biofilm_diff_abund$compartment <- "biofilm"
water_diff_abund$compartment <- "water"

diff_abund_table  <- rbind(biofilm_diff_abund, water_diff_abund)

write.table(diff_abund_table, file = "F:/16S_HOLOGREEN/Picrust/Objets/diff_abund_table.txt", sep=";")

#Analyses article
mean(subset(water_diff_abund, Gene_Group == "NirBD_NrfAH")$Diff_Abundance)

mean(subset(biofilm_diff_abund, Gene_Group == "NITsix_NirA_NasBDE")$Abundance_ENR)
mean(subset(biofilm_diff_abund, Gene_Group == "NITsix_NirA_NasBDE")$Abundance_SW)


#Check si ASV dans KO----
# Liste des ASV d'intérêt (exemple)
test <- data.frame(physeq_16S_nit_p@otu_table[,colSums(subset_samples(physeq_16S_nit_p, type =="biofilm")@otu_table)>0])


colSums(KO_tab2[names(test),c("K00370", "K00371", "K00374", "K02567", "K02568","K00362", "K00363", "K03385", "K15876","K00368", "K15864","K04561", "K02305", "K00376","K00367", "K00372", "K00366")])




#Camembert plot----


# Fonction pour calculer la différence d'abondance
picrust_pie <- function(physeq_obj, sample_type, enrich_type) {
  # Créer un masque logique pour filtrer les échantillons
  mask <- rownames(physeq_obj@sam_data[physeq_obj@sam_data$type %in% sample_type,])
  
  # Sous-ensemble de l'objet phyloseq en utilisant le masque
  physeq_subset <- prune_samples(mask, physeq_obj)
  
  # Initialiser un dataframe pour stocker les résultats
  results <- data.frame(OTU  = character(),enrich = character(),Class = character(),
                        Genus = character(),Abundance = numeric(),stringsAsFactors = FALSE)
  
  # Calcul de la différence d'abondance pour chaque groupe de gènes
  for (gene_group_name in names(valid_gene_groups)) {
    gene_group <- valid_gene_groups[[gene_group_name]]
    asv_subset <- prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, gene_group, drop = FALSE]) >= 1, ]), physeq_subset)
    psmelt <- psmelt(asv_subset)
    results <- rbind(results,  psmelt)
  }
  # Remplacer les valeurs NA par "unassigned" et le niveau taxonomique supérieur
  for (i in 27:32) {
    results[, i] <- ifelse(is.na(results[, i]), 
                                     paste("unassigned", results[, (i-1)], sep = " "), 
                           results[, i])
  }
  
  
  results <- results %>% group_by(OTU, enrich, Class, Genus) %>% summarise(Abundance = mean(Abundance))
  results <- results %>% group_by(enrich, Class, Genus) %>% summarise(Abundance = sum(Abundance))
  
  for (c in unique(results$Class)) {
    tmp <- subset(results, Class ==c)
    other_genus <- tmp[, "Genus"]
    if (max(tmp$Abundance) < 4) {
      results[results[,"Class"]==c,][results[results[,"Class"]==c,]$Abundance <1,]$Genus <- "Others <1%"
    }
    if (max(tmp$Abundance) >= 4){
      results[results[,"Class"]==c,][results[results[,"Class"]==c,]$Abundance <1,]$Genus <- paste0("Others <1% ",c)
    }
  }
  
  results <- results %>% group_by(enrich, Class, Genus) %>% summarise(Abundance= sum(Abundance))
  
  return(results)
}



##partie water----
picrust_water_pie <- picrust_pie(physeq_16S_nit_p, "water")

picrust_water_pie_plot <- picrust_water_pie%>%
  group_by(enrich) %>%
  # Calculer le total de l'abondance pour chaque 'enrich'
  mutate(Total_Abundance = sum(Abundance)) %>%
  # Convertir l'abondance en pourcentage
  mutate(Percentage = (Abundance / Total_Abundance) * 100) %>%
  # Sélectionner les colonnes nécessaires
  select(enrich, Genus, Percentage)
picrust_water_pie_plot <- picrust_water_pie_plot %>% group_by(enrich, Genus) %>% summarise(Percentage= sum(Percentage))

# Définir l'ordre des genres
picrust_water_pie_order <- rev(as.factor(c("Marivita","Ponticoccus","Jannaschia","unassigned Rhodobacteraceae","Roseobacter","Yoonia-Loktanella", "Others <1% Alphaproteobacteria","Glaciecola", "Algoriphagus"  ,"Others <1%")))
picrust_water_pie_plot$Genus <- factor(picrust_water_pie_plot$Genus, levels = picrust_water_pie_order)

picrust_water_pie_col <- rev(c("#D9A46A","#C16A17", "#572B06","#FFDFDF","#B22222", "#8C450B" ,"#251800","#E7D4E8", "#0E509E", "#E0E0E0"))

ggplot() +
  geom_bar(data= picrust_water_pie_plot, aes(x = "", y = Percentage, fill = Genus),width = 1, stat = "identity") +
  facet_wrap(~enrich)+
  coord_polar("y", start = 0) +
  scale_fill_manual(values = picrust_water_pie_col) +
  labs(title = "Abondance relative des différents genres bactériens") +
  theme_minimal() +
  theme(legend.position = "bottom")


##partie biof----
picrust_biof_pie <- picrust_pie(physeq_16S_nit_p, "biofilm")

picrust_biof_pie_plot <- picrust_biof_pie%>%
  group_by(enrich) %>%
  # Calculer le total de l'abondance pour chaque 'enrich'
  mutate(Total_Abundance = sum(Abundance)) %>%
  # Convertir l'abondance en pourcentage
  mutate(Percentage = (Abundance / Total_Abundance) * 100) %>%
  # Sélectionner les colonnes nécessaires
  select(enrich, Genus, Percentage)
picrust_biof_pie_plot <- picrust_biof_pie_plot %>% group_by(enrich, Genus) %>% summarise(Percentage= sum(Percentage))


# Définir l'ordre des genres
picrust_biof_pie_order <- rev(as.factor(c("Agaribacterium", "Glaciecola","Others <1% Gammaproteobacteria","Jannaschia","unassigned Rhodobacteraceae","Roseobacter","unassigned Rhizobiaceae","Octadecabacter", "Others <1% Alphaproteobacteria","Croceitalea", "Maribacter", "Winogradskyella", "Others <1% Bacteroidia"  ,"Others <1%")))
picrust_biof_pie_plot$Genus <- factor(picrust_biof_pie_plot$Genus, levels = picrust_biof_pie_order)

picrust_biof_pie_col <- rev(c("#9970B5","#C2A5CF","#762A83","#F1E5D5","#D9A46A","#FFDFDF","#FF7F7F","#B22222","#251800","#EFF3FF","#C3DAF6","#96C1ED","#0C2A66","#E0E0E0") )



ggplot() +
  geom_bar(data= picrust_biof_pie_plot, aes(x = "", y = Percentage, fill = Genus),width = 1, stat = "identity") +
  facet_wrap(~enrich)+
  coord_polar("y", start = 0) +
  scale_fill_manual(values = picrust_biof_pie_col) +
  labs(title = "Abondance relative des différents genres bactériens") +
  theme_minimal() +
  theme(legend.position = "bottom")



##Supp table article pie----
picrust.pie.supp.article <- function(physeq_obj, type_enrich1) {
  # Créer un masque logique pour filtrer les échantillons
  mask <- rownames(physeq_obj@sam_data[physeq_obj@sam_data$type_enrich %in% type_enrich1,])
  
  # Sous-ensemble de l'objet phyloseq en utilisant le masque
  physeq_subset <- prune_samples(mask, physeq_obj)

    # Initialiser un dataframe pour stocker les résultats
  results <- data.frame(OTU  = character(),enrich = character(),Class = character(),
                        Genus = character(),Abundance = numeric(),stringsAsFactors = FALSE)
  
  # Calcul de la différence d'abondance pour chaque groupe de gènes
  for (gene_group_name in names(valid_gene_groups)) {
    gene_group <- valid_gene_groups[[gene_group_name]]
    asv_subset <- prune_taxa(rownames(KO_tab2[rowSums(KO_tab2[, gene_group, drop = FALSE]) >= 1, ]), physeq_subset)
    psmelt <- psmelt(asv_subset)
    results <- rbind(results,  psmelt)
  }
  # Remplacer les valeurs NA par "unassigned" et le niveau taxonomique supérieur
  for (i in 27:32) {
    results[, i] <- ifelse(is.na(results[, i]), 
                           paste("unassigned", results[, (i-1)], sep = " "), 
                           results[, i])
  }
  
  results$full_taxo <- paste(paste0("r___",results$Kingdom), paste0("p___",results$Phylum), paste0("c___",results$Class), paste0("o___",results$Order), paste0("f___",results$Family),paste0("g___",results$Genus) , sep = ";")
  
  results <- results %>% group_by(OTU, enrich, Class, Genus, full_taxo) %>% summarise(Abundance = mean(Abundance))
  results <- results %>% group_by(enrich, Class, Genus, full_taxo) %>% summarise(Abundance = sum(Abundance))
  
  results$Genus_graph <- results$Genus
  
  for (c in unique(results$Class)) {
    tmp <- subset(results, Class ==c)
    other_genus <- tmp[, "Genus"]
    if (max(tmp$Abundance) < 4) {
      results[results[,"Class"]==c,][results[results[,"Class"]==c,]$Abundance <1,]$Genus_graph <- "Others"
    }
    if (max(tmp$Abundance) >= 4){
      results[results[,"Class"]==c,][results[results[,"Class"]==c,]$Abundance <1,]$Genus_graph <- paste0("Others",c)
    }
  }
  
  results <- results %>% group_by(enrich, Class, Genus, Genus_graph, full_taxo) %>% summarise(Abundance= sum(Abundance))
  
  return(results)
}


picrust_biof_ENR_pie_supp_article <- picrust.pie.supp.article(physeq_obj = physeq_16S_nit_p, type_enrich1 = "biofilm_ENR")
picrust_biof_SW_pie_supp_article <- picrust.pie.supp.article(physeq_obj = physeq_16S_nit_p, type_enrich1 = "biofilm_SW")
picrust_water_ENR_pie_supp_article <- picrust.pie.supp.article(physeq_obj = physeq_16S_nit_p, type_enrich1 = "water_ENR")
picrust_water_SW_pie_supp_article <- picrust.pie.supp.article(physeq_obj = physeq_16S_nit_p,type_enrich1 = "water_SW")

picrust_biof_ENR_pie_supp_article$type <- "biofilm"
picrust_biof_SW_pie_supp_article$type <- "biofilm"
picrust_water_ENR_pie_supp_article$type <- "water"
picrust_water_SW_pie_supp_article$type <- "water"



write.table(rbind(picrust_biof_pie_supp_article,picrust_water_pie_supp_article), file ="F:/16S_HOLOGREEN/Picrust/Objets/picrust_pie_supp_article.txt", sep ="+")
write.csv2(rbind(picrust_biof_ENR_pie_supp_article,
                 picrust_biof_SW_pie_supp_article,
                 picrust_water_ENR_pie_supp_article,
                 picrust_water_SW_pie_supp_article), file ="F:/16S_HOLOGREEN/Picrust/Objets/picrust_pie_supp_article.csv")


#Combien les ASV représentent de % de séquences 
sum(picrust_biof_ENR_pie_supp_article$Abundance)
sum(picrust_biof_SW_pie_supp_article$Abundance)
sum(picrust_water_ENR_pie_supp_article$Abundance)
sum(picrust_water_SW_pie_supp_article$Abundance)

#nombre de genre associé a une condi 
sum(picrust_biof_SW_pie_supp_article$Abundance > 0)
sum(picrust_biof_ENR_pie_supp_article$Abundance > 0)
sum(picrust_water_SW_pie_supp_article$Abundance > 0)
sum(picrust_water_ENR_pie_supp_article$Abundance > 0)


#percentage in piechart
picrust_biof_pie_plot %>% 
  mutate(Percentage = round(Percentage, 1))

picrust_water_pie_plot %>% 
  mutate(Percentage = round(Percentage, 1))%>%
  View()

#Ecriture texte article----
round(sum(subset(picrust_biof_pie, Class=="Gammaproteobacteria"&enrich == "ENR")$Abundance),1)
round(sum(subset(picrust_biof_pie, Class=="Alphaproteobacteria"&enrich == "SW")$Abundance),1)
round(sum(subset(picrust_biof_pie, Class=="Bacteroidia"&enrich == "SW")$Abundance),1)

round(sum(subset(picrust_water_pie, Class=="Gammaproteobacteria"&enrich == "ENR")$Abundance),1)
round(sum(subset(picrust_water_pie, Class=="Alphaproteobacteria"&enrich == "SW")$Abundance),1)
round(sum(subset(picrust_water_pie, Class=="Bacteroidia"&enrich == "SW")$Abundance),1)

round(sum(subset(picrust_biof_pie, Genus=="Agaribacterium"&enrich == "SW")$Abundance),1)
round(sum(subset(picrust_biof_pie, Genus=="Glaciecola"&enrich == "ENR")$Abundance),1)
round(sum(subset(picrust_biof_pie, Genus=="Jannaschia"&enrich == "SW")$Abundance),1)
round(sum(subset(picrust_biof_pie, Genus=="unassigned Rhodobacteraceae"&enrich == "ENR")$Abundance),1)
round(sum(subset(picrust_biof_pie, Genus=="Croceitalea"&enrich == "ENR")$Abundance),1)
round(sum(subset(picrust_biof_pie, Genus=="Maribacter"&enrich == "ENR")$Abundance),1)

round(sum(subset(picrust_water_pie, Genus=="Ponticoccus"&enrich == "SW")$Abundance),1)
round(sum(subset(picrust_water_pie, Genus=="Marivita"&enrich == "ENR")$Abundance),1)
round(sum(subset(picrust_water_pie, Genus=="Glaciecola"&enrich == "ENR")$Abundance),1)
round(sum(subset(picrust_water_pie, Genus=="Algoriphagus"&enrich == "ENR")$Abundance),1)
