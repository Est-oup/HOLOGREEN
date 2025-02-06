#Aldex2 ----


#Loadings
#BiocManager::install("ALDEx2")

library(ALDEx2)
library(dplyr)
library(ggplot2)


physeq <- physeq_16S_nit

#Biofilm analysis 
#so
data_aldex_a <- subset_samples(physeq, type == "biofilm")
data_aldex_a <- physeq_biof_core
# #extraction des ASV d'intéret 
# data_aldex_a <- prune_taxa(taxa_names(data_aldex_a) %in% c(ASV_core_a_SW, ASV_core_a_ENR) ,data_aldex_a)


mm_aldex_a <- microbiomeMarker::run_aldex(data_aldex_a, group = "enrich",
                                          norm = "CPM",
                                          taxa_rank = "none",
                                          p_adjust = "fdr")

mm_aldex_table_a <- data.frame(mm_aldex_a@marker_table)
# mm_aldex_table_a

# mm_aldex_table_a$feature <- factor(unique(arrange(mm_aldex_table_a, ef_aldex)$feature))

# Ordonner le dataframe selon la colonne ef_aldex        non + faire la moyenne au niveau du genre
mm_aldex_table_a <- arrange(mm_aldex_table_a, ef_aldex)
mm_aldex_table_a <- merge(mm_aldex_table_a, cbind("feature" = noquote(row.names(data_aldex_a@tax_table)),as.data.frame(data_aldex_a@tax_table)), by= "feature")
#mm_aldex_table_a <- mm_aldex_table_a %>%dplyr::group_by(Genus, enrich_group)%>%summarize(mean_ef_aldex=mean(ef_aldex, na.rm=TRUE))

#etiquette des données 
names_diff_a <- paste(mm_aldex_table_a$Genus, mm_aldex_table_a$feature, sep = " ")





#Filtres analysis
data_aldex_f <- subset_samples(physeq_16S_nit, type == "filter")

#extraction des ASV d'intéret 
data_aldex_f <- prune_taxa(taxa_names(data_aldex_f) %in% c(ASV_core_f_SW, ASV_core_f_ENR) ,data_aldex_f)

mm_aldex_f <- microbiomeMarker::run_aldex(data_aldex_f, group = "enrich",
                                          norm = "CPM",
                                          taxa_rank = "none",
                                          p_adjust = "fdr")

mm_aldex_table_f <- data.frame(mm_aldex_f@marker_table)
mm_aldex_table_f

mm_aldex_table_f$feature <- factor(unique(arrange(mm_aldex_table_f, ef_aldex)$feature))

# Ordonner le dataframe selon la colonne ef_aldex + faire la moyenne au niveau du genre
mm_aldex_table_f <- arrange(mm_aldex_table_f, ef_aldex)
mm_aldex_table_f <- merge(mm_aldex_table_f, cbind("feature" = noquote(row.names(data_aldex_f@tax_table)),as.data.frame(data_aldex_f@tax_table)), by= "feature")
#mm_aldex_table_f <- mm_aldex_table_f %>%dplyr::group_by(Genus, enrich_group)%>%summarize(mean_ef_aldex=mean(ef_aldex, na.rm=TRUE))

#etiquette des données 
names_diff_f <- paste(mm_aldex_table_f$Genus, mm_aldex_table_f$feature, sep = " ")





#Plot général 



# Créer le graphique algue
aldex2_a_plot <- ggplot(mm_aldex_table_a , aes(x = ef_aldex, y = feature, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  xlab("ef_aldex") + ylab("Feature") +
  ggtitle("Differential abundance of genus in biofilm") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_discrete(limits=arrange(mm_aldex_table_a, ef_aldex)$feature, labels= names_diff_a)+
  scale_fill_manual(values = c("SW" = "#33CC33", "ENR" = "#99FF66"))


# Créer le graphique filtres
aldex2_f_plot <- ggplot(mm_aldex_table_f , aes(x = ef_aldex, y = feature, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  xlab("ef_aldex") + ylab("Feature") +
  ggtitle("Differential abundance of genus in filter") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_discrete(limits=arrange(mm_aldex_table_f, ef_aldex)$feature, labels = names_diff_f)+
  scale_fill_manual(values = c("SW" = "#0099CC", "ENR" = "#33CCFF"))

#graph groupé 


# Organiser les heatmaps dans une figure 2x2
aldex2_g_plot <- plot_grid(
  aldex2_f_plot,
  aldex2_a_plot,
  nrow = 1, ncol = 2,
  labels = c("", "", "", ""),
  align = "hv",
  rel_heights = c(1, 1),
  rel_widths = c(1, 1)
)

aldex2_g_plot









#Global analysis biofilm-filtre
data_aldex_1 <- Final1_16S

mm_aldex_1 <- microbiomeMarker::run_aldex(data_aldex_1, group = "type",
                                          norm = "CPM",
                                          taxa_rank = "none",
                                          p_adjust = "fdr")

mm_aldex_table_1 <- data.frame(mm_aldex_1@marker_table)
mm_aldex_table_1

mm_aldex_table_1$feature <- factor(unique(arrange(mm_aldex_table_1, ef_aldex)$feature))

# Ordonner le dataframe selon la colonne ef_aldex + faire la moyenne au niveau du genre
mm_aldex_table_1 <- arrange(mm_aldex_table_1, ef_aldex)
mm_aldex_table_1 <- merge(mm_aldex_table_1, cbind("feature" = noquote(row.names(Final2_16S_Ts@tax_table)),as.data.frame(Final2_16S_Ts@tax_table)), by= "feature")
mm_aldex_table_1 <- mm_aldex_table_1 %>%dplyr::group_by(Genus, enrich_group)%>%summarize(mean_ef_aldex=mean(ef_aldex, na.rm=TRUE))

mm_aldex_table_1$Genus <- as.factor(arrange(mm_aldex_table_1, mean_ef_aldex)$Genus)

# Créer le graphique
ggplot(mm_aldex_table_1 , aes(x = mean_ef_aldex, y = Genus, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  xlab("ef_aldex") + ylab("Feature") +
  ggtitle("Barplot avec ggplot") +
  theme_minimal() +
  scale_y_discrete(limits=arrange(mm_aldex_table_1, mean_ef_aldex)$Genus)




#Global analysis EDM-ENR
data_aldex_2 <- Final1_16S

mm_aldex_2 <- microbiomeMarker::run_aldex(data_aldex_2, group = "enrich",
                                          norm = "CPM",
                                          taxa_rank = "none",
                                          p_adjust = "fdr")

mm_aldex_table_2 <- data.frame(mm_aldex_2@marker_table)
mm_aldex_table_2

mm_aldex_table_2$feature <- factor(unique(arrange(mm_aldex_table_2, ef_aldex)$feature))

# Ordonner le dataframe selon la colonne ef_aldex + faire la moyenne au niveau du genre
mm_aldex_table_2 <- arrange(mm_aldex_table_2, ef_aldex)
mm_aldex_table_2 <- merge(mm_aldex_table_2, cbind("feature" = noquote(row.names(Final2_16S_Ts@tax_table)),as.data.frame(Final2_16S_Ts@tax_table)), by= "feature")
mm_aldex_table_2 <- mm_aldex_table_2 %>%dplyr::group_by(Genus, enrich_group)%>%summarize(mean_ef_aldex=mean(ef_aldex, na.rm=TRUE))

mm_aldex_table_2 <- data.frame(arrange(mm_aldex_table_2, mean_ef_aldex))

# Créer le graphique
ggplot(mm_aldex_table_2 , aes(y = mean_ef_aldex, x = Genus, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  xlab("ef_aldex") + ylab("Feature") +
  ggtitle("Barplot avec ggplot") +
  theme_minimal() +
  scale_x_discrete(limits=arrange(mm_aldex_table_2, mean_ef_aldex)$Genus)+
  coord_flip()





ggplot(mm_aldex_table_2 , aes(y = ef_aldex, x = feature, fill = Genus)) +
  geom_bar(stat = "identity") +
  xlab("ef_aldex") + ylab("feature") +
  ggtitle("Barplot avec ggplot") +
  theme_minimal() +
  scale_x_discrete(limits=arrange(mm_aldex_table_2, ef_aldex)$feature)+
  coord_flip()


microbiomeMarker::plot_ef_bar(mm_aldex_2)

cbind(subset(mm_aldex_table_2, enrich_group =="EDM"),subset(mm_aldex_table_2, enrich_group =="ENR"))





#ANF META------


#LEFSE
mm_lefse <- microbiomeMarker::run_lefse(physeq, norm = "CPM",
                                        wilcoxon_cutoff = 0.01,
                                        group = "enrich",
                                        taxa_rank = "none",
                                        kw_cutoff = 0.01,
                                        multigrp_strat = TRUE,
                                        lda_cutoff = 4)

mm_lefse_table <- data.frame(mm_lefse@marker_table)
mm_lefse_table

physeq@tax_table[mm_lefse_table$feature,]
physeq@otu_table[,mm_lefse_table$feature]


p_LDAsc <- microbiomeMarker::plot_ef_bar(mm_lefse)
p_abd <- microbiomeMarker::plot_abundance(mm_lefse, group = "enrich")
gridExtra::grid.arrange(p_LDAsc, p_abd, nrow = 1)


#ancomBC
mm_ancombc <- microbiomeMarker::run_ancombc(physeq, group = "enrich",
                                            taxa_rank = "all",
                                            pvalue_cutoff = 0.01,
                                            p_adjust = "fdr")

mm_ancombc_table <- data.frame(mm_ancombc@marker_table)
mm_ancombc_table







##########

library(mia)
library(patchwork)
library(tidySummarizedExperiment)
library(knitr)
library(tidyverse)
library(ALDEx2)
library(Maaslin2)
library(MicrobiomeStat)
library(ANCOMBC)
library(GUniFrac)



# set random seed because some tools can randomly vary and then produce 
# different results:
set.seed(13253)

# we use a demo dataset and restrict it to two geo locations
# for easy illustration
data(peerj13075)
tse0 <- peerj13075
tse0 <- tse0[ ,tse0$Geographical_location %in% c("Pune", "Nashik")]
# Let us make this a factor
tse0$Geographical_location <- factor(tse0$Geographical_location)

# how many observations do we have per group?
as.data.frame(colData(tse0)) %>% 
  count(Geographical_location) %>%
  kable()

####Prevalence filtering 


tse <- agglomerateByRank(tse0, rank = "genus") %>%
  transformCounts(assay.type = "counts",
                  method = "relabundance",
                  MARGIN = "samples") %>%
  # subset based on the relative abundance assay              
  subsetByPrevalentTaxa(detection = 0,
                        prevalence = 10/100,
                        assay.type = "relabundance")

# Add also clr abundances
tse <- transformCounts(tse, method="clr", pseudocount=1) # not bale to run

