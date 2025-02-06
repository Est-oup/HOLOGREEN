

heat_tree_fonction <- function (physeq){
  
  library("phyloseq")
  library("metacoder")
  library("exactRankTests")
  
  # Récupérer la date et l'heure actuelles
  date_time <- Sys.time()
  
  # Formater la date et l'heure dans le format souhaité
  formatted_date_time <- format(date_time, "%Y-%m-%d_%H-%M-%S")
  
  
  ## Sous dimensionnement de lobjet phyloseq ####
  ### Choix des groupes ####
  print("1.0")
  sam_answer <- readline("Do you want to subset a sepcial sample group ? yes or no :")
  if (sam_answer == "yes"){
    env_var <- readline("Which variable ?")
    env_val <- readline("Which value ?")
    eval_expr_1 <- paste0("subset_samples(", "physeq, ", "as.factor(sample_data(physeq)$", env_var , ") %in% '", env_val, "')")
    physeq <- eval(parse(text = eval_expr_1))
    print(physeq)
  }
  
  ### Choix des taxons ####
  print("1.1")
  taxon_answer <- readline("Do you want to subset a sepcial group from the taxonomy ? yes or no: ")
  if (taxon_answer == "yes"){
    taxon_level <- readline("Enter the taxon level : ")
    taxon_name <- readline("Enter the taxon name: ")
    eval_expr_2 <- paste0("subset_taxa(physeq, ", taxon_level, " == '", taxon_name, "')")
    physeq <- eval(parse(text = eval_expr_2)) 
    print(physeq)
  }
  
  
  
  ### Choix du nombre de taxons ####
  print("1.2")
  taxon_number_answer <- readline("Do you want to subset a special number of taxon ? yes or no: ")
  if (taxon_number_answer == "yes"){
    taxon_number <- readline("Enter the number of taxa to be analysed : ")
    taxon_number <- as.numeric(taxon_number)
    physeq <- prune_taxa(names(sort(taxa_sums(physeq), TRUE)[1:taxon_number]), physeq)
    print(physeq)
  }
  
  
  ### Donner le chemin du dossier où vous voulez enregistrer le fichier ####
  
  nom_repertoire_figure <- readline("Give the directory path ne pas oublier le type de slash windows ou linux  à la fin : ")
  
  print("1.3")
  
  
  
  
  
  assign("physeq", physeq, envir = .GlobalEnv)
  
  
  
  
  
  ## Creation de hmp_otu et hmp_samples ####
  
  ### Creation de la taxonomie speciale ####
  
  print("1.4")
  taxo_number_level <- length(colnames(tax_table(physeq)))
  taxTab_physeq <- tax_table(physeq)
  
  taxTab_format_special <- rep(NA,length(rownames(tax_table(physeq))))
  
  print("1.5")
  
  if (taxo_number_level == 7) {
    
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
  }
  
  
  
  if (taxo_number_level == 9) {
    
    for (i in 1:length(rownames(tax_table(physeq)))) {
      
      
      taxTab_format_special[i] <- paste0(
        "r__", taxTab_physeq[i, 1],
        ";p__", taxTab_physeq[i, 2],
        ";c__", taxTab_physeq[i, 3],
        ";o__", taxTab_physeq[i, 4],
        ";f__", taxTab_physeq[i, 5],
        ";g__", taxTab_physeq[i, 6],
        ";s__", taxTab_physeq[i, 7],
        ";t__", taxTab_physeq[i, 8],
        ";u__", taxTab_physeq[i, 9]
      )
    }
  }
  
  print("1.6")
  classif <- sub(".{2}__NA.*", "", taxTab_format_special)
  
  
  ### Creation de lobjet hmp_otus ####
  
  print("1.7")
  
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
  
  print("1.8")
  hmp_samples <- data.frame(sample_data(physeq))
  rownames(hmp_samples) <- rownames(sample_data(physeq))
  
  
  
  ## Creation de lobjet metacoder ####
  
  print("1.9")      
  obj_metacoder <- parse_tax_data(hmp_otus,
                                  class_cols = "classif", # the column that contains taxonomic information
                                  class_sep = ";", # The character used to separate taxa in the classification
                                  class_regex = "^(.+)__(.+)$",
                                  class_key = c(tax_rank = "info", # A key describing each regex capture group
                                                tax_name = "taxon_name"))
  print(obj_metacoder)
  
  
  
  # Suppression des ASV avec un faible nombre de read (optionel)####
  print("1.10")
  low_abun_ASV <- readline("Do you want to remove the low abundance ASV or not ? yes or no: ")
  
  
  
  if (low_abun_ASV == "yes"){
    
    count_mini_ASV <- readline("What threshold limit do you want to apply (number of reads) : ")
    count_mini_ASV <- as.numeric(count_mini_ASV)
    
    obj_metacoder$data$tax_data <- zero_low_counts(obj_metacoder, data = "tax_data", min_count = count_mini_ASV)
    no_reads <- rowSums(obj_metacoder$data$tax_data[, rownames(hmp_samples)]) == 0
    print(no_reads)
    
    print(sum(no_reads))
    
    obj_metacoder <- filter_obs(obj_metacoder, data = "tax_data", ! no_reads, drop_taxa = TRUE)
    print(obj_metacoder)
    
  }
  
  ### Detection de la colonne correspondant a la variable dinteret
  env_var <- readline("Veuillez rentrer la variable dinteret pour construire la matrice?")
  env_var_number <- which(colnames(hmp_samples) == env_var)
  print(env_var_number)
  
  print("1.11")
  obj_metacoder$data$tax_data <- calc_obs_props(obj_metacoder, "tax_data")
  
  
  print(colnames(obj_metacoder$data$tax_data))
  print("1.12")    
  obj_metacoder$data$tax_abund <- calc_taxon_abund(obj_metacoder, "tax_data",
                                                   cols = rownames(hmp_samples))
  
  
  print("1.13")
  obj_metacoder$data$tax_occ <- calc_n_samples(obj_metacoder, "tax_abund", groups = hmp_samples[,env_var_number], cols = rownames(hmp_samples))
  
  
  # Plotting read depth:
  #  To plot read depth, you first need to add up the number of reads per taobj_metacoderon.
  #  The function `calc_taobj_metacoderon_abund` is good for this. 
  obj_metacoder$data$taxon_counts <- calc_taxon_abund(obj_metacoder, data = "tax_data")
  obj_metacoder$data$taxon_counts$total <- rowSums(obj_metacoder$data$taxon_counts[, -1]) # -1 = taxon_id column
  
  
  assign("objmetacoder",obj_metacoder, envir = .GlobalEnv)
  
  
  ### Choix de la figure a realiser ####
  
  fig <- TRUE
  
  while (fig == TRUE) {
    
    fig <- readline(" Quelle figure voulez vous faire (1,2,3,4) : ")
    fig <- as.numeric(fig)
    
    
    if (fig == 1){
      
      
      if(sam_answer =="no"){
        print("test")
        vect_requis <- rep("vect_requis",length(rownames(hmp_samples)))
      }
      
      set.seed(1) # This makes the plot appear the same each time it is run 
      eval_expr_2 <- paste0("heat_tree(obj_metacoder, ",
                            "node_label = taxon_names, ",
                            "node_size = n_obs, ",
                            "node_color = ", env_val, ", ",
                            "node_size_axis_label = 'OTU count', ",
                            "node_color_axis_label = 'Samples with reads', ",
                            "layout = 'davidson-harel', ",
                            "initial_layout = 'reingold-tilford', ",
                            "output_file = paste0(nom_repertoire_figure, 'heat_tree_', formatted_date_time, '.pdf'))")
      
      eval(parse(text = eval_expr_2))
      
      
      
    }
    
    
    if (fig == 2){
      
      obj_metacoder$data$diff_table <- compare_groups(obj_metacoder,
                                                      data = "tax_abund",
                                                      cols = rownames(hmp_samples), # What columns of sample data to use
                                                      groups = hmp_samples[,env_var_number],# What category each sample is assigned to
                                                      func = function(abund_1, abund_2) {
                                                        log_ratio <- log2(median(abund_1) / median(abund_2))
                                                        if (is.nan(log_ratio)) {
                                                          log_ratio <- 0
                                                        }
                                                        list(log2_median_ratio = log_ratio,
                                                             median_diff = median(abund_1) - median(abund_2),
                                                             mean_diff = mean(abund_1) - mean(abund_2),
                                                             wilcox_p_value = wilcox.exact(abund_1, abund_2)$p.value)})
      
      
      obj_metacoder$data$diff_table$wilcox_p_value <- p.adjust(obj_metacoder$data$diff_table$wilcox_p_value,
                                                               method = "fdr")
      paste0("Le range de p-value est : ",range(obj_metacoder$data$diff_table$wilcox_p_value, finite = TRUE))
      
      obj_metacoder$data$diff_table$log2_median_ratio[obj_metacoder$data$diff_table$wilcox_p_value > 0.05] <- 0
      
      heat_tree(obj_metacoder, 
                node_label = taxon_names,
                node_size = total, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                node_color = log2_median_ratio, # A column from `obj_metacoder$data$diff_table`
                node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
                node_color_range = c("#75b8d1", "gray92", "#FFD700"), # The color palette used
                node_size_axis_label = "Abundance",
                node_color_axis_label = "Log 2 ratio of median proportions",
                layout = "davidson-harel", # The primary layout algorithm
                initial_layout = "reingold-tilford",
                overlap_avoidance=1.5,
                output_file = paste0(nom_repertoire_figure, 'heat_tree_', formatted_date_time, '.pdf'))
      
      
    }
    
    
    if (fig == 3){
      
      if (nlevels(as.factor(hmp_samples[,2])) < 2 ){
        
        print("Le jeu de donnees que vous avez selectionner ne contient pas assez de groupes differents pour efffectuer la figure 3, minimum 2 groupes)")
        
        
      }else {
        
        
        obj_metacoder$data$diff_table <- compare_groups(obj_metacoder,
                                                        data = "tax_abund",
                                                        cols = rownames(hmp_samples), # What columns of sample data to use
                                                        groups = hmp_samples[,env_var_number],# What category each sample is assigned to
                                                        func = function(abund_1, abund_2) {
                                                          log_ratio <- log2(median(abund_1) / median(abund_2))
                                                          if (is.nan(log_ratio)) {
                                                            log_ratio <- 0
                                                          }
                                                          list(log2_median_ratio = log_ratio,
                                                               median_diff = median(abund_1) - median(abund_2),
                                                               mean_diff = mean(abund_1) - mean(abund_2),
                                                               wilcox_p_value = wilcox.exact(abund_1, abund_2)$p.value)})
        
        
        obj_metacoder$data$diff_table$wilcox_p_value <- p.adjust(obj_metacoder$data$diff_table$wilcox_p_value,
                                                                 method = "fdr")
        paste0("Le range de p-value est : ",range(obj_metacoder$data$diff_table$wilcox_p_value, finite = TRUE))
        
        obj_metacoder$data$diff_table$log2_median_ratio[obj_metacoder$data$diff_table$wilcox_p_value > 0.05] <- 0
        
        
        
        
        heat_tree_matrix(obj_metacoder,
                         data = "diff_table",
                         node_size = total, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                         node_label = taxon_names,
                         node_color = log2_median_ratio, # A column from `obj_metacoder$data$diff_table`
                         node_color_range = diverging_palette(), # The built-in palette for diverging data
                         node_color_trans = "linear", # The default is scaled by circle area
                         node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                         edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                         node_size_axis_label = "Abundance of OTUs",
                         node_color_axis_label = "Log2 ratio median proportions",
                         layout = "davidson-harel", # The primary layout algorithm
                         initial_layout = "reingold-tilford",
                         output_file = paste0(nom_repertoire_figure, 'heat_tree_', formatted_date_time, '.pdf'))
      }
      
      
      
      
    }  
    
    
    
    
    
    
    
    ### Continuation ou fin de la boucle ####
    
    fig <- readline(" Voulez vous faire une autre figure, si yes laquelle (1,2,3,4) sinon remplissez par no : ")   
    
    if (fig == "no") {
      fig <- FALSE
    } else {
      fig <- TRUE
    }
    
    
  }
  
  
}




# # Figure 1 seul groupe 
# 
# print("1.11")
# print(group_var)
# 
# heat_tree(obj_metacoder, 
#           node_label = taxon_names,
#           node_size = n_obs,
#           node_color = get(env_val),
#           node_size_axis_label = "OTU count",
#           node_color_axis_label = "Samples with reads",
#           layout = "davidson-harel",
#           initial_layout = "reingold-tilford",
#           output_file = paste0(nom_repertoire_figure,"heat_tree_",formatted_date_time,".pdf"))
# 
# 
# print("1.12")
# }
# 
# 
# 
# 
# 
# 
# 
#     ### Statistiques sur l'objet dans le cas de comparaisons entre 2 groupes ###
#     
#     
#     # Comparaison des groupes
#     
#     print(hmp_samples[,num_col_groupes])
#     
#     library("plyr")
#     
#     if (nlevels(factor(hmp_samples[,num_col_groupes])) >= 2) {
#       
#       library(exactRankTests)
#       obj_metacoder$data$diff_table <- compare_groups(obj_metacoder,
#                                             data = "tax_abund",
#                                             cols = hmp_samples$Stations, # What columns of sample data to use
#                                             groups = hmp_samples[,num_col_groupes],# What category each sample is assigned to
#                                             func = function(abund_1, abund_2) {
#                                               log_ratio <- log2(median(abund_1) / median(abund_2))
#                                               if (is.nan(log_ratio)) {
#                                                 log_ratio <- 0
#                                               }
#                                               list(log2_median_ratio = log_ratio,
#                                                    median_diff = median(abund_1) - median(abund_2),
#                                                    mean_diff = mean(abund_1) - mean(abund_2),
#                                                    wilcox_p_value = wilcox.exact(abund_1, abund_2)$p.value)})
#       
#     }
#     
#     print("1.7")
#     
#     # Stats sur les groupes
#     
#     obj_metacoder$data$diff_table$wilcox_p_value <- p.adjust(obj_metacoder$data$diff_table$wilcox_p_value,
#                                                    method = "fdr")
#     
#     print("1.8")
#     
#     range(obj_metacoder$data$diff_table$wilcox_p_value, finite = TRUE)
#     
#     print("1.9")
#     
#     obj_metacoder$data$diff_table$log2_median_ratio[obj_metacoder$data$diff_table$wilcox_p_value > 0.05] <- 0
#     
#     print("1.10")
#     print(head(obj_metacoder$data$diff))
#     
#     
#     if (type_fig == 1){
#       # Figure 1 seul groupe 
#       
#       print("1.11")
#       print(group_var)
#       
#       X11()
#       
#       set.seed(1) # This makes the plot appear the same each time it is run 
#       heat_tree(obj_metacoder, 
#                 node_label = taxon_names,
#                 node_size = n_obs,
#                 node_color = group_var,
#                 node_size_axis_label = "OTU count",
#                 node_color_axis_label = "Samples with reads",
#                 layout = "davidson-harel", # The primary layout algorithm
#                 initial_layout = "reingold-tilford",
#                 output_file = paste0(nom_repertoire_figure,"heat_tree_",group_var,"_",formatted_date_time,".pdf")) # The layout algorithm that initializes node locations
#       
#       print("1.12")
#     }
#     
#     
#     
#     if (type_fig == 2){
#       
#       # Comparaison 2 groupes en fonction des stats calculées avant
#       
#       print("1.2")
#       
#       x11()
#       
#       print("1.21")
#       
#       set.seed(999)
#       heat_tree(obj_metacoder, 
#                 node_label = taxon_names,
#                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
#                 node_color = log2_median_ratio, # A column from `obj_metacoder$data$diff_table`
#                 node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
#                 node_color_range = c("#75b8d1", "gray92", "#FFD700"), # The color palette used
#                 node_size_axis_label = "OTU count",
#                 node_color_axis_label = "Log 2 ratio of median proportions",
#                 layout = "davidson-harel", # The primary layout algorithm
#                 initial_layout = "reingold-tilford",
#                 overlap_avoidance=1.5,
#                 output_file = paste0(nom_repertoire_figure,"differential_heat_tree",group_var,"_",formatted_date_time,".pdf")) # The layout algorithm that initializes node locations
#       
#       print("1.22")
#     }
#     
#     
#     
#     
#     
#     if (type_fig == 3){
#       
#       # n Comparaison entre n groupes, 1 figure de descriptions + toutes les autres combinaisons
#       
#       print("1.31")
#       
#       
#       x11()
#       set.seed(1)
#       heat_tree_matrix(obj_metacoder,
#                        data = "diff_table",
#                        node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
#                        node_label = taxon_names,
#                        node_color = log2_median_ratio, # A column from `obj_metacoder$data$diff_table`
#                        node_color_range = diverging_palette(), # The built-in palette for diverging data
#                        node_color_trans = "linear", # The default is scaled by circle area
#                        node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
#                        edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
#                        node_size_axis_label = "Number of OTUs",
#                        node_color_axis_label = "Log2 ratio median proportions",
#                        layout = "davidson-harel", # The primary layout algorithm
#                        initial_layout = "reingold-tilford",
#                        output_file = paste0(nom_repertoire_figure,"matrix_comparison_heat_tree",group_var,"_",formatted_date_time,".pdf"))
#       print("1.32")
#     }
#     
#   }
#   
#   print("2.0")
#   
#   if (data_particuliere == FALSE){
#     
#     
#     print("2.1")  
#     
#     hmp_otus <- t(otu_table(phylo_obj))
#     hmp_otus <- data.frame(hmp_otus)
#     hmp_samples <- data.frame(sample_data(phylo_obj))
#     
#     classif <- sub(".{2}__NA.*", "", taxa_speciale)
#     
#     hmp_otus <- cbind(classif, hmp_otus)
#     hmp_otus <- data.frame(hmp_otus)
#     colnames(hmp_otus) <- c("classif",stations)
#     
#     hmp_otus <- hmp_otus[!grepl("__NA", hmp_otus$classif), ]
#     hmp_otus <- hmp_otus[!grepl("r__unclassified sequences", hmp_otus$classif), ]
#     hmp_otus <- hmp_otus[!grepl("^r__Bacteria$", hmp_otus$classif), ]
#     hmp_otus <- hmp_otus[!grepl("^r__Archaea", hmp_otus$classif), ]
#     
#     
#     print("2.2") 
#     
#     
#     #### Creation de l'objet metacoder ####
#     
#     obj <- parse_tax_data(hmp_otus,
#                           class_cols = "classif", # the column that contains taxonomic information
#                           class_sep = ";", # The character used to separate taxa in the classification
#                           class_regex = "^(.+)__(.+)$",
#                           class_key = c(tax_rank = "info", # A key describing each regex capture group
#                                         tax_name = "taxon_name"))
#     
#     print("2.3")
#     
#     #obj$data$tax_data <- zero_low_counts(obj, data = "tax_data", min_count = 5)
#     
#     #no_reads <- rowSums(obj$data$tax_data[, hmp_samples$Stations]) == 0
#     #sum(no_reads)
#     
#     # <- filter_obs(obj, data = "tax_data", ! no_reads, drop_taxa = TRUE)
#     
#     obj$data$tax_data <- calc_obs_props(obj, "tax_data")
#     
#     print("2.4")
#     
#     obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
#                                            cols = hmp_samples$Stations)
#     
#     print("2.5")
#     
#     obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = hmp_samples$groupes, cols = hmp_samples$Stations)
#     
#     
#     print("2.6")
#     
#     ### Statistiques sur l'objet dans le cas de comparaisons entre 2 groupes ###
#     
#     
#     # Comparaison des groupes
#     
#     library(exactRankTests)
#     obj$data$diff_table <- compare_groups(obj,
#                                           data = "tax_abund",
#                                           cols = hmp_samples$Stations, # What columns of sample data to use
#                                           groups = hmp_samples[,num_col_groupes],# What category each sample is assigned to
#                                           func = function(abund_1, abund_2) {
#                                             log_ratio <- log2(median(abund_1) / median(abund_2))
#                                             if (is.nan(log_ratio)) {
#                                               log_ratio <- 0
#                                             }
#                                             list(log2_median_ratio = log_ratio,
#                                                  median_diff = median(abund_1) - median(abund_2),
#                                                  mean_diff = mean(abund_1) - mean(abund_2),
#                                                  wilcox_p_value = wilcox.exact(abund_1, abund_2)$p.value)})
#     
#     
#     print("2.7")
#     
#     # Stats sur les groupes
#     
#     obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value,
#                                                    method = "fdr")
#     
#     print("2.8")
#     
#     range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
#     
#     print("2.9")
#     
#     obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0
#     
#     print("2.10")
#     print(head(obj$data$diff))
#     
#     
#     if (type_fig == 1){
#       # Figure 1 seul groupe 
#       
#       print("2.11")
#       print(as.character(groupe_seul))
#       X11()
#       
#       set.seed(1) # This makes the plot appear the same each time it is run 
#       heat_tree(obj, 
#                 node_label = taxon_names,
#                 node_size = n_obs,
#                 node_color = as.character(groupe_seul),
#                 node_size_axis_label = "OTU count",
#                 node_color_axis_label = "Samples with reads",
#                 layout = "davidson-harel", # The primary layout algorithm
#                 initial_layout = "reingold-tilford",
#                 output_file = paste0(nom_repertoire_figure,"heat_tree_",group_var,"_",formatted_date_time,".pdf")) # The layout algorithm that initializes node locations
#       
#       print("2.13")
#     }
#     
#     
#     
#     if (type_fig == 2){
#       
#       # Comparaison 2 groupes en fonction des stats calculées avant
#       
#       print("2.12")
#       
#       x11()
#       
#       print("2.121")
#       
#       set.seed(999)
#       heat_tree(obj, 
#                 node_label = taxon_names,
#                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
#                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
#                 node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
#                 node_color_range = c("#75b8d1", "gray92", "#FFD700"), # The color palette used
#                 node_size_axis_label = "OTU count",
#                 node_color_axis_label = "Log 2 ratio of median proportions",
#                 layout = "davidson-harel", # The primary layout algorithm
#                 initial_layout = "reingold-tilford",
#                 overlap_avoidance=1.5,
#                 output_file = paste0(nom_repertoire_figure,"differential_heat_tree",group_var,"_",formatted_date_time,".pdf")) # The layout algorithm that initializes node locations
#       
#       print("2.122")
#     }
#     
#     
#     
#     
#     
#     if (type_fig == 3){
#       
#       # n Comparaison entre n groupes, 1 figure de descriptions + toutes les autres combinaisons
#       
#       print("2.131")
#       
#       
#       x11()
#       set.seed(1)
#       heat_tree_matrix(obj,
#                        data = "diff_table",
#                        node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
#                        node_label = taxon_names,
#                        node_color = log2_median_ratio, # A column from `obj$data$diff_table`
#                        node_color_range = diverging_palette(), # The built-in palette for diverging data
#                        node_color_trans = "linear", # The default is scaled by circle area
#                        node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
#                        edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
#                        node_size_axis_label = "Number of OTUs",
#                        node_color_axis_label = "Log2 ratio median proportions",
#                        layout = "davidson-harel", # The primary layout algorithm
#                        initial_layout = "reingold-tilford",
#                        output_file = paste0(nom_repertoire_figure,"matrix_comparison_heat_tree",group_var,"_",formatted_date_time,".pdf"))
#       print("2.132")
#     }
#     
#   }
#   
# }





# #aldex2 results on heat_tree
# 
# library(ALDEx2)
# library(microbiomeMarker)
# 
# #Biofilm analysis 
# data_aldex_a <- subset_samples(physeq_16S_nit, type == "biofilm")
# 
# test<- t(data.frame(physeq_16S_nit@otu_table)[,1:400])
# 
# testna1 <- colnames(test)[1:42]
# testna2 <- colnames(test)[43:64]
# 
# x.all <- aldex(test, c(testna1, testna2), mc.samples=16, test="t", effect =TRUE, include.sample.summary = F, denom="all", verbose=T, paired.test=F)
# x.all <- aldex.clr(test, c(testna1, testna2), mc.samples=100, denom="all", verbose=T)
# 
# 
# par(mfrow=c(1,2))
# aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance",
#            ylab="Difference")
# aldex.plot(x.all, type="MW", test="welch", xlab="Dispersion",
#            ylab="Difference")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# mm_aldex_a <- run_aldex(data_aldex_a, group = "enrich",
#                                           norm = "CPM",
#                                           taxa_rank = "none",
#                                           p_adjust = "fdr")
# 
# mm_aldex_table_a <- data.frame(mm_aldex_a@marker_table)
# # mm_aldex_table_a
# 
# 
# # Ordonner le dataframe selon la colonne ef_aldex + faire la moyenne au niveau du genre
# mm_aldex_table_a <- arrange(mm_aldex_table_a, ef_aldex)
# mm_aldex_table_a <- merge(mm_aldex_table_a, cbind("feature" = noquote(row.names(data_aldex_a@tax_table)),as.data.frame(data_aldex_a@tax_table)), by= "feature")
# mm_aldex_table_a <- mm_aldex_table_a %>%dplyr::group_by(Genus, enrich_group)%>%summarize(mean_ef_aldex=mean(ef_aldex, na.rm=TRUE))
# 
# 
# 
# 
# #metacoder 
# 
# 
# 
# 



## Creation de hmp_otu et hmp_samples ####

### Creation de la taxonomie speciale ####

taxo_number_level <- length(colnames(tax_table(data_aldex_a)))
taxTab_physeq <- tax_table(data_aldex_a)

taxTab_format_special <- rep(NA,length(rownames(tax_table(data_aldex_a))))

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

print("1.7")

hmp_otus <- t(otu_table(data_aldex_a))
hmp_otus <- data.frame(hmp_otus)

hmp_otus <- cbind(classif, hmp_otus)
hmp_otus <- data.frame(hmp_otus)
colnames(hmp_otus) <- c("classif",rownames(sample_data(data_aldex_a)))

hmp_otus <- hmp_otus[!grepl("r__NA", hmp_otus$classif), ]
hmp_otus <- hmp_otus[!grepl("r__unclassified sequences", hmp_otus$classif), ]
hmp_otus <- hmp_otus[!grepl("^r__Bacteria$", hmp_otus$classif), ]
hmp_otus <- hmp_otus[!grepl("^r__Archaea", hmp_otus$classif), ]

### Creation de lobjet hmp_sample ####    

print("1.8")
hmp_samples <- data.frame(sample_data(data_aldex_a))
rownames(hmp_samples) <- rownames(sample_data(data_aldex_a))



## Creation de lobjet metacoder ####

print("1.9")      
obj_metacoder <- parse_tax_data(hmp_otus,
                                class_cols = "classif", # the column that contains taxonomic information
                                class_sep = ";", # The character used to separate taxa in the classification
                                class_regex = "^(.+)__(.+)$",
                                class_key = c(tax_rank = "info", # A key describing each regex capture group
                                              tax_name = "taxon_name"))
print(obj_metacoder)


#ajout de l'abondance 
obj_metacoder$data$tax_data <- zero_low_counts(obj_metacoder, data = "tax_data", min_count = 500)
no_reads <- rowSums(obj_metacoder$data$tax_data[, rownames(hmp_samples)]) == 0

obj_metacoder <- filter_obs(obj_metacoder, data = "tax_data", ! no_reads, drop_taxa = TRUE)
print(obj_metacoder)


obj_metacoder$data$tax_abund <- calc_taxon_abund(obj_metacoder, "tax_data",
                                                 cols = rownames(hmp_samples))
obj_metacoder$data$tax_occ <- calc_n_samples(obj_metacoder, "tax_abund", groups = hmp_samples[,"type"], cols = rownames(hmp_samples))






obj_metacoder$data$diff_table <- compare_groups(obj_metacoder,
                                                data = "tax_abund",
                                                cols = rownames(hmp_samples), # What columns of sample data to use
                                                groups = hmp_samples[,"enrich"],# What category each sample is assigned to
                                                func = function(abund_1, abund_2) {
                                                  log_ratio <- log2(median(abund_1) / median(abund_2))
                                                  if (is.nan(log_ratio)) {
                                                    log_ratio <- 0
                                                  }
                                                  list(log2_median_ratio = log_ratio,
                                                       median_diff = median(abund_1) - median(abund_2),
                                                       mean_diff = mean(abund_1) - mean(abund_2),
                                                       wilcox_p_value = wilcox.exact(abund_1, abund_2)$p.value)})


obj_metacoder$data$diff_table$wilcox_p_value <- p.adjust(obj_metacoder$data$diff_table$wilcox_p_value,
                                                         method = "fdr")
paste0("Le range de p-value est : ",range(obj_metacoder$data$diff_table$wilcox_p_value, finite = TRUE))

obj_metacoder$data$diff_table$log2_median_ratio[obj_metacoder$data$diff_table$wilcox_p_value > 0.05] <- 0


obj_metacoder$data$taxon_counts <- calc_taxon_abund(obj_metacoder, data = "tax_data")
obj_metacoder$data$taxon_counts$total <- rowSums(obj_metacoder$data$taxon_counts[, -1]) # -1 = taxon_id column





heat_tree(obj_metacoder, 
          node_label = taxon_names,
          node_size = total, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
          node_color = log2_median_ratio, # A column from `obj_metacoder$data$diff_table`
          node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
          node_color_range = c("#75b8d1", "gray92", "#FFD700"), # The color palette used
          node_size_axis_label = "Abundance",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford",
          overlap_avoidance=1.5,
          output_file = paste0(nom_repertoire_figure, 'heat_tree_', formatted_date_time, '.pdf'))


}

