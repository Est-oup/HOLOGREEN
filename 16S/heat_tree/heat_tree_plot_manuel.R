  library("phyloseq")
  library("metacoder")
  library("exactRankTests")
  library(phyloseq)
  library(zCompositions)
  
 
  
  #input 
  #objets
  
  
  #normalisation 

  
  # # Charger les données phyloseq
  # 
  # physeq <- physeq_16S_nit
  # 
  # #Enelever les colonnes avec des 0 
  # physeq@otu_table <- physeq@otu_table[,colSums(physeq@otu_table) != 0]
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
  # #create a new phyloseq object with CLR tranformed counts
  # physeq_clr <- physeq
  # phyloseq::otu_table(physeq_clr) <- phyloseq::otu_table(t(physeq_clr_asv),
  #                                    taxa_are_rows = FALSE)
  # 
  
  
  
  physeq <- physeq_clr
  #physeq <- microbiome::transform(physeq_16S_nit, 'clr')

  #
  samp <- "type"
  varsamp <- "biofilm"
  sepenv <- "enrich"
  nbsamp <- 100
  
  
  ## Sous dimensionnement de lobjet phyloseq ####
  ### Choix des groupes ####
#subset selon la var env 
  physeq <- phyloseq::subset_samples(physeq, get(samp) == varsamp)
  

  
  
  ### Choix du nombre de taxons ####
    physeq <- prune_taxa(names(sort(taxa_sums(physeq), TRUE)[1:nbsamp]), physeq)
  
  
  ## Creation de hmp_otu et hmp_samples ####
  
  ### Creation de la taxonomie speciale ####

  taxTab_physeq <- tax_table(physeq)
  
  taxTab_format_special <- rep(NA,length(rownames(tax_table(physeq))))

    for (i in 1:length(rownames(tax_table(physeq)))) {
      
      taxTab_format_special[i] <- paste0(
        "r__", taxTab_physeq[i, 1],
        ";p__", taxTab_physeq[i, 2],
        ";c__", taxTab_physeq[i, 3],
        ";o__", taxTab_physeq[i, 4],
        ";f__", taxTab_physeq[i, 5],
        ";g__", taxTab_physeq[i, 6],
        ";s__", taxTab_physeq[i, 7]
      )
    }

  classif <- sub(".{2}__NA.*", "", taxTab_format_special)
  
  
  ### Creation de lobjet hmp_otus ####
  hmp_otus <- t(otu_table(physeq))
  hmp_otus <- data.frame(hmp_otus)
  
  hmp_otus <- cbind(classif, hmp_otus)
  hmp_otus <- data.frame(hmp_otus)
  colnames(hmp_otus) <- c("classif",rownames(sample_data(physeq)))
  
  hmp_otus <- hmp_otus[!grepl("r__NA", hmp_otus$classif), ]
  hmp_otus <- hmp_otus[!grepl("r__unclassified sequences", hmp_otus$classif), ]
  hmp_otus <- hmp_otus[!grepl("^r__Bacteria$", hmp_otus$classif), ]
  hmp_otus <- hmp_otus[!grepl("^r__Archaea", hmp_otus$classif), ]
  
  ### Creation de lobjet hmp_sample ####    
  
  hmp_samples <- data.frame(sample_data(physeq))
  rownames(hmp_samples) <- rownames(sample_data(physeq))
  
  
  
  ## Creation de lobjet metacoder ####
    obj_metacoder <- parse_tax_data(hmp_otus,
                                  class_cols = "classif", # the column that contains taxonomic information
                                  class_sep = ";", # The character used to separate taxa in the classification
                                  class_regex = "^(.+)__(.+)$",
                                  class_key = c(tax_rank = "info", # A key describing each regex capture group
                                                tax_name = "taxon_name"))

  
  
  
  
  ### Detection de la colonne correspondant a la variable dinteret

  env_var_number <- which(colnames(hmp_samples) == sepenv)

  
  obj_metacoder$data$tax_data <- calc_obs_props(obj_metacoder, "tax_data")
  
  
  obj_metacoder$data$tax_abund <- calc_taxon_abund(obj_metacoder, "tax_data",
                                                   cols = rownames(hmp_samples))
  
  
  obj_metacoder$data$tax_occ <- calc_n_samples(obj_metacoder, "tax_abund", groups = hmp_samples[,env_var_number], cols = rownames(hmp_samples))
  
  
  # Plotting read depth:
  #  To plot read depth, you first need to add up the number of reads per taobj_metacoderon.
  #  The function `calc_taobj_metacoderon_abund` is good for this. 
  obj_metacoder$data$taxon_counts <- calc_taxon_abund(obj_metacoder, data = "tax_data")
  obj_metacoder$data$taxon_counts$total <- rowMeans(obj_metacoder$data$taxon_counts[, -1]) # -1 = taxon_id column
  
 

  
  ### Choix de la figure a realiser ####
 
    
 
        
        obj_metacoder$data$diff_table <- compare_groups(obj_metacoder,
                                                        data = "tax_abund",
                                                        cols = rownames(hmp_samples), # What columns of sample data to use
                                                        groups = hmp_samples[,env_var_number])# What category each sample is assigned to)
        
        
        obj_metacoder$data$diff_table$wilcox_p_value <- p.adjust(obj_metacoder$data$diff_table$wilcox_p_value,
                                                                 method = "fdr")
        paste0("Le range de p-value est : ",range(obj_metacoder$data$diff_table$wilcox_p_value, finite = TRUE))
        
        obj_metacoder$data$diff_table$log2_median_ratio[obj_metacoder$data$diff_table$wilcox_p_value > 0.05] <- 0
        

        if (varsamp == "biofilm"){
          #pour mettre les bonnes couleurs sur sw et enr il faut inverser les treatment 
          #changer le nom 
          colnames(obj_metacoder$data$diff_table) <- c("taxon_id","treatment_2","treatment_1","log2_median_ratio","median_diff","mean_diff","wilcox_p_value")
          
          #changer l'ordre + inverser la valeur du log 
          
          obj_metacoder$data$diff_table <- cbind(obj_metacoder$data$diff_table[,1],obj_metacoder$data$diff_table[,3],obj_metacoder$data$diff_table[,2],-obj_metacoder$data$diff_table[,4], obj_metacoder$data$diff_table[,5:7])
          
        }
       
        
       res_biofilm_heat <-  heat_tree_matrix(obj_metacoder,
                         data = "diff_table",
                         node_size = total, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                         node_label = taxon_names,
                         node_color = log2_median_ratio, # A column from `obj_metacoder$data$diff_table`
                         node_color_range = diverging_palette(), # The built-in palette for diverging data
                         node_color_trans = "linear", # The default is scaled by circle area
                         node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                         edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                         node_size_axis_label = "Abundance",
                         node_color_axis_label = "Log2 ratio median proportions",
                         layout = "davidson-harel", # The primary layout algorithm
                         initial_layout = "reingold-tilford")
                         
        
       
       
       
       
       
       #######POUR LES FILTRES
       
       
       
       physeq <- physeq_clr
       #physeq <- microbiome::transform(physeq_16S_nit, 'clr')
       
       #
       samp <- "type"
       varsamp <- "filter"
       sepenv <- "enrich"
       nbsamp <- 100
       
       
       ## Sous dimensionnement de lobjet phyloseq ####
       ### Choix des groupes ####
       #subset selon la var env 
       physeq <- phyloseq::subset_samples(physeq, get(samp) == varsamp)
       
       
       
       
       ### Choix du nombre de taxons ####
       physeq <- prune_taxa(names(sort(taxa_sums(physeq), TRUE)[1:nbsamp]), physeq)
       
       
       ## Creation de hmp_otu et hmp_samples ####
       
       ### Creation de la taxonomie speciale ####
       
       taxTab_physeq <- tax_table(physeq)
       
       taxTab_format_special <- rep(NA,length(rownames(tax_table(physeq))))
       
       for (i in 1:length(rownames(tax_table(physeq)))) {
         
         taxTab_format_special[i] <- paste0(
           "r__", taxTab_physeq[i, 1],
           ";p__", taxTab_physeq[i, 2],
           ";c__", taxTab_physeq[i, 3],
           ";o__", taxTab_physeq[i, 4],
           ";f__", taxTab_physeq[i, 5],
           ";g__", taxTab_physeq[i, 6],
           ";s__", taxTab_physeq[i, 7]
         )
       }
       
       classif <- sub(".{2}__NA.*", "", taxTab_format_special)
       
       
       ### Creation de lobjet hmp_otus ####
       hmp_otus <- t(otu_table(physeq))
       hmp_otus <- data.frame(hmp_otus)
       
       hmp_otus <- cbind(classif, hmp_otus)
       hmp_otus <- data.frame(hmp_otus)
       colnames(hmp_otus) <- c("classif",rownames(sample_data(physeq)))
       
       hmp_otus <- hmp_otus[!grepl("r__NA", hmp_otus$classif), ]
       hmp_otus <- hmp_otus[!grepl("r__unclassified sequences", hmp_otus$classif), ]
       hmp_otus <- hmp_otus[!grepl("^r__Bacteria$", hmp_otus$classif), ]
       hmp_otus <- hmp_otus[!grepl("^r__Archaea", hmp_otus$classif), ]
       
       ### Creation de lobjet hmp_sample ####    
       
       hmp_samples <- data.frame(sample_data(physeq))
       rownames(hmp_samples) <- rownames(sample_data(physeq))
       
       
       
       ## Creation de lobjet metacoder ####
       obj_metacoder <- parse_tax_data(hmp_otus,
                                       class_cols = "classif", # the column that contains taxonomic information
                                       class_sep = ";", # The character used to separate taxa in the classification
                                       class_regex = "^(.+)__(.+)$",
                                       class_key = c(tax_rank = "info", # A key describing each regex capture group
                                                     tax_name = "taxon_name"))
       
       
       
       
       
       ### Detection de la colonne correspondant a la variable dinteret
       
       env_var_number <- which(colnames(hmp_samples) == sepenv)
       
       
       obj_metacoder$data$tax_data <- calc_obs_props(obj_metacoder, "tax_data")
       
       
       obj_metacoder$data$tax_abund <- calc_taxon_abund(obj_metacoder, "tax_data",
                                                        cols = rownames(hmp_samples))
       
       
       obj_metacoder$data$tax_occ <- calc_n_samples(obj_metacoder, "tax_abund", groups = hmp_samples[,env_var_number], cols = rownames(hmp_samples))
       
       
       # Plotting read depth:
       #  To plot read depth, you first need to add up the number of reads per taobj_metacoderon.
       #  The function `calc_taobj_metacoderon_abund` is good for this. 
       obj_metacoder$data$taxon_counts <- calc_taxon_abund(obj_metacoder, data = "tax_data")
       obj_metacoder$data$taxon_counts$total <- rowMeans(obj_metacoder$data$taxon_counts[, -1]) # -1 = taxon_id column
       
       
       
       
       ### Choix de la figure a realiser ####
       
       
       
       
       obj_metacoder$data$diff_table <- compare_groups(obj_metacoder,
                                                       data = "tax_abund",
                                                       cols = rownames(hmp_samples), # What columns of sample data to use
                                                       groups = hmp_samples[,env_var_number])# What category each sample is assigned to)
       
       
       obj_metacoder$data$diff_table$wilcox_p_value <- p.adjust(obj_metacoder$data$diff_table$wilcox_p_value,
                                                                method = "fdr")
       paste0("Le range de p-value est : ",range(obj_metacoder$data$diff_table$wilcox_p_value, finite = TRUE))
       
       obj_metacoder$data$diff_table$log2_median_ratio[obj_metacoder$data$diff_table$wilcox_p_value > 0.05] <- 0
       
       
       res_filter_heat<-  heat_tree_matrix(obj_metacoder,
                                        data = "diff_table",
                                        node_size = total, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                                        node_label = taxon_names,
                                        node_color = log2_median_ratio, # A column from `obj_metacoder$data$diff_table`
                                        node_color_range = diverging_palette(), # The built-in palette for diverging data
                                        node_color_trans = "linear", # The default is scaled by circle area
                                        node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                                        edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                                        node_size_axis_label = "Abundance",
                                        node_color_axis_label = "Log2 ratio median proportions",
                                        layout = "davidson-harel", # The primary layout algorithm
                                        initial_layout = "reingold-tilford")
        
       
       
       #Plot général 
       library(gridExtra)
       library(ggplot2)
       
       
       #plot 
       grid.arrange(res_filter_heat,
                    res_biofilm_heat,
                         nrow = 1)
           
       