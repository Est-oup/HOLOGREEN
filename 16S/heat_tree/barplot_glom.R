###Abundance >X%----

#Version améliorée avec les 4 variables 
glom <- tax_glom(merged_physeq2_16S_nit_p, taxrank = 'Genus', NArm = FALSE, bad_empty = c(""))

data_glom <- psmelt(glom)

for (i in (ncol(data_glom)-5):ncol(data_glom)) {  # Boucle de la colonne Genus à la colonne Kingdom
  data_glom[, i] <- ifelse(is.na(data_glom[, i]), 
                           paste("unknown", data_glom[, (i-1)], sep = "_"), 
                           data_glom[, i])
}

data_glom$Genus <- as.character(data_glom$Genus)
data_glom$Class <- as.character(data_glom$Class)
data_glom <-data.frame(data_glom %>% dplyr::group_by(OTU,Sample,Abundance,type_enrich,Kingdom,Phylum,Class,Order,Family,Genus)%>% summarize(Abundance = sum(Abundance)))


for (c in unique(data_glom$Class)) {
  tmp <- data_glom[data_glom$Class == c, ]
  other_genus <- tmp[, "Genus"]
  if (max(tmp$Abundance) < 2) {
    data_glom[rownames(data_glom) %in% rownames(unique(tmp)),][data_glom[rownames(data_glom) %in% rownames(unique(tmp)),]$Abundance < 2,]$Genus <- "Others <2%"
  }
  if (max(tmp$Abundance) >= 2){
    data_glom[rownames(data_glom) %in% rownames(unique(tmp)),][data_glom[rownames(data_glom) %in% rownames(unique(tmp)),]$Abundance < 2,]$Genus <- paste0("Others <2% ",c)
  }
}

# Créer une nouvelle colonne "title" en utilisant les conditions demandées
data_glom$title <- ifelse(data_glom$Genus != "Others <2%", data_glom$Class, "Others")


# Réordonner le tableau en fonction de l'ordre obtenu
data_glom_ordered <- data_glom

data_glom_ordered <- data_glom_ordered %>% 
  group_by(type_enrich ,Class, Genus) %>% 
  summarize(Abundance = sum(Abundance)) %>% 
  ungroup()

data_glom_ordered <- data_glom_ordered[order(data_glom_ordered$Class), ]


genus_levels <- factor(c(unique(data_glom_ordered$Genus)[-3], "Others <2%"))
genus_levels <- factor(c("Others <2%","Others <2% Acidimicrobiia", "unknown_Microtrichaceae", "Others <2% Deinococci","Truepera", "Others <2% Alphaproteobacteria" ,"Jannaschia","unknown_Rhodobacteraceae","Litorimonas","Roseobacter","unknown_Hyphomonadaceae","Nereida","Ponticoccus","Marivita","Sulfitobacter","Yoonia-Loktanella",
                         "Others <2% Bacteroidia","Croceitalea","Maribacter","Winogradskyella","Aurantivirga","NS3a marine group","Polaribacter","unknown_Cryomorphaceae","Others <2% Gammaproteobacteria","Agaribacterium","Glaciecola","Granulosicoccus"))




pal1 <- colorRampPalette(c("#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "#084594", "#08306B", "#002040"))
# pal2 <- colorRampPalette(c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704", "#651C00", "#540A00"))
# pal2 <- colorRampPalette(c("#F9D71C", "#F3B71C", "#E4A502", "#C78502", "#AD6E00", "#855C00", "#5E4D00", "#3F3A00", "#292200"))
# pal2 <- colorRampPalette(c("#FFFBE6", "#FFF5CC", "#FFF0B3", "#FFE699", "#FFD966", "#FFC61A", "#FFB300", "#E6A600", "#CCA300", "#B38F00", "#997300"))
# pal2 <- colorRampPalette(c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704", "#651C00", "#4F1300"))
# pal2 <- colorRampPalette(c("#C3A562", "#D5C29F", "#DDBB85", "#E1CFA5", "#E3D7B3", "#EAE0C5", "#EFD8A7", "#F1E1B6", "#F6E9C7", "#F9ECCF", "#FBF0D8"))
pal2 <- colorRampPalette(c("#F1E5D5", "#E5C29F", "#D9A46A", "#CD8740", "#C16A17", "#A95611", "#8C450B", "#703508", "#572B06", "#3E2103", "#251800"))
# pal2 <- colorRampPalette(c("#F4C430", "#E7B416", "#D9A200", "#C49400", "#A67700", "#8A5A00", "#6E3D00", "#532000", "#3D1500", "#280B00", "#120500"))





pal3 <- colorRampPalette(c("#E7D4E8", "#C2A5CF", "#9970B5", "#762A83"))
pal4 <- colorRampPalette(c("#FEE5D9", "#FCAE91"))
pal5 <- colorRampPalette(c("#CCEBC5", "#78C679"))
pal6 <- c("#E0E0E0")


pal_fan <- c(pal6,pal5(2),pal4(2), pal2(11), pal1(8) ,pal3(4)) 



data_glom$Genus <- factor(data_glom$Genus, levels = genus_levels)

ggplot(data=data_glom, aes(x=type_enrich, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual( values = pal_fan, breaks = genus_levels) +
  guides(fill=guide_legend(ncol=1))+theme_minimal()





