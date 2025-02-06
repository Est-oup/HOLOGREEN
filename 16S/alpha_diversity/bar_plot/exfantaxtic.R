###Fantaxtic package----   

# Get the top taxa
top_level <- "Class"
nested_level <- "Genus"
sample_order <- NULL



top_asv<- top_taxa(Final2_16S_Ts_P , n_taxa = 25)
top_asv <- top_taxa(Final2_16S_Ts_P, n_taxa = 25)
top_asv <- top_taxa(physeq2_16S_nit_p, n_taxa = 25)
top_asv <- top_taxa(subset_samples(physeq2_16S_nit_p , type == "biofilm"), n_taxa = 45)


# Generate a palette based  on the phyloseq object
#  pal <- taxon_colours(top_asv$ps_obj,tax_level = top_level)
pal= c(Alphaproteobacteria= "green" , Gammaproteobacteria= "hotpink", Bacteroidia = "lightcyan", Other="gray90", Deinococci="#FF9900", Phycisphaerae="#993300")
# Create names for NA taxa
ps_tmp <- top_asv$ps_obj %>%
  name_na_taxa( na_label = "Unidentified <tax> (<rank>)") 

# Add labels to taxa with the same names
ps_tmp <- ps_tmp %>%
  label_duplicate_taxa(tax_level = nested_level)

# Convert physeq to df
psdf <- psmelt(ps_tmp)

# Move the merged labels to the appropriate positions in the plot:
# Top merged labels need to be at the top of the plot,
# nested merged labels at the bottom of each group
psdf <- move_label(psdf = psdf,
                   col_name = top_level,
                   label = "Other",
                   pos = 0)

psdf <- move_nested_labels(psdf,
                           top_level = top_level,
                           nested_level = nested_level,
                           top_merged_label = "Other",
                           nested_label = "Other",
                           pos = Inf)

# Move the other label to the start, and Bacteroidetes to the end
psdf <- move_label(psdf, col_name = "Class", label =  "Other", pos = 0)
psdf <- move_label(psdf, col_name = "Class", label =  "Alphaproteobacteria", pos = 1)
psdf <- move_label(psdf, col_name = "Class", label =  "Gammaproteobacteria", pos = 2)
psdf <- move_label(psdf, col_name = "Class", label =  "Bacteroidia", pos =3)

levels(psdf$Phylum)
etiquettes_date <- sample_data(subset_samples( Final2_16S_Ts_P , date <= "210601"))$date
psdf$date <- factor(psdf$date)
# Generate a base plot
p <- ggnested(psdf,
              aes_string(main_group = top_level,
                         sub_group = nested_level,
                         x = "date2",
                         y = "Abundance"), )+ #  main_palette = pal) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_nested(theme_light) +
  theme(axis.text.x = element_text(angle=90, hjust=0)) + 
  xlab("Echantillon")+ ylab("Abondance relative %") 
# scale_x_discrete(labels = c(rep(subset_samples(subset_samples( Final2_16S_Ts_P , date <= "210601") , type_enrich == "biofilm_EDM")@sam_data$date,2),rep(subset_samples(subset_samples( Final2_16S_Ts_P , date <= "210601") , type_enrich == "filtre_ENR")@sam_data$date,2)))
#  theme(legend.position="bottom" ) 
#   guides(fill=guide_legend(ncol = 3, nrow= 10 ,byrow = F))

# Add relative abundances
p <- p + geom_col(position = position_fill())

#Faire les facets 
type_enrichlabs<- c("Biofilm eau de mer", "Biofilm engrais","Filtre eau de mer", "Filtre engrais")
names(type_enrichlabs) <- c("biofilm_EDM", "biofilm_ENR", "filtre_EDM", "filtre_ENR")
p <- p + facet_wrap(~type_enrich, scales = "free_x", labeller = labeller(type_enrich=type_enrichlabs))+
  theme(strip.text.x= element_text(size = 13, color = "black", face="bold"))

p



#test pour le ettiquettes de noms 
order(subset_samples( Final2_16S_Ts_P , date <= "210601")@sam_data$Sample_ID , x=c( algae_EDM[1:21] , algae_ENR[1:21], filter_EDM[1:12], filter_ENR[1:12]))
c( algae_EDM[1:21] , algae_ENR[1:21], filter_EDM[1:16], filter_ENR[1:16])
c(rep(subset_samples(subset_samples( Final2_16S_Ts_P , date <= "210601") , type_enrich == "biofilm_EDM")@sam_data$date,2),rep(subset_samples(subset_samples( Final2_16S_Ts_P , date <= "210601") , type_enrich == "filtre_EDM")@sam_data$date,2)) 









# To identify and plot the top 3 most abundant Phyla, and the top 3 most abundant species within those Phyla, run:

data_fantaxtic <- pyseq_merged_type_enrich
data_fantaxtic <- Final2_16S_Ts_P
data_fantaxtic <- merged_physeq2_16S_nit_p
data_fantaxtic <- merged_physeq2_16S_nitur_p
data_fantaxtic <- physeq2_16S_nitur_p
data_fantaxtic <- physeq2_16S_nit_p
data_fantaxtic <- subset_samples(physeq2_16S_nit_p , type == "biofilm")


top_nested <- nested_top_taxa(data_fantaxtic,
                              top_tax_level = "Class",
                              nested_tax_level = "Genus",
                              n_top_taxa = 5, 
                              n_nested_taxa = 4)

plot_nested_bar(top_nested$ps_obj,
                top_level = "Class",
                nested_level = "Genus",
                nested_merged_label = "NA and other",
                palette = c(Actinobacteriota = "cornflowerblue",Bacteroidota = "darkgoldenrod", Planctomycetota = "green",  Proteobacteria = "hotpink", Deinococcota = "#FF9900"),
)+
  xlab("Sample's category")+ ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size=9, face="bold", color = "black"))+
  facet_wrap(~type )


top_nested <- nested_top_taxa(data_fantaxtic,
                              top_tax_level = "Class",
                              nested_tax_level = "Genus",
                              n_top_taxa = 10, 
                              n_nested_taxa = )

plot_nested_bar(top_nested$ps_obj,
                top_level = "Class",
                nested_level = "Genus",
                nested_merged_label = "NA and other",
                palette = c(Actinobacteriota = "cornflowerblue",Bacteroidota = "darkgoldenrod", Planctomycetota = "green",  Proteobacteria = "hotpink", Deinococcota = "#FF9900"),)+
  xlab("Sample's category")+ ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size=9, face="bold", color = "black"))+
  facet_wrap(~type_enrich )

