



## Creation de lobjet metacoder ####

print("1.9")      
obj_metacoder <- parse_tax_data(hmp_otus,
                                class_cols = "classif", # the column that contains taxonomic information
                                class_sep = ";", # The character used to separate taxa in the classification
                                class_regex = "^(.+)__(.+)$",
                                class_key = c(tax_rank = "info", # A key describing each regex capture group
                                              tax_name = "taxon_name"))
print(obj_metacoder)


obj_metacoder$data$tax_data

obj_metacoder$data$class_data




mm_aldex_a <- run_aldex(data_aldex_a, group = "enrich",
                        norm = "CPM",
                        taxa_rank = "none",
                        p_adjust = "fdr")

mm_aldex_table_a <- data.frame(mm_aldex_a@marker_table)
# mm_aldex_table_a


# Ordonner le dataframe selon la colonne ef_aldex + faire la moyenne au niveau du genre
mm_aldex_table_a <- arrange(mm_aldex_table_a, ef_aldex)
mm_aldex_table_a <- merge(mm_aldex_table_a, cbind("feature" = noquote(row.names(data_aldex_a@tax_table)),as.data.frame(data_aldex_a@tax_table)), by= "feature")
mm_aldex_table_a <- mm_aldex_table_a %>%dplyr::group_by(Genus, enrich_group)%>%summarize(mean_ef_aldex=mean(ef_aldex, na.rm=TRUE))





obj_metacoder$data$diff_table <- as_tibble(mm_aldex_table_a)






####adaptation aldex

phyloseq_obj <- physeq_16S_nit


tax_levels <- c("Kingdom" ,"Phylum" , "Class"  , "Order" ,  "Family" , "Genus"  , "Species")  # Niveaux taxonomiques souhaités

aldex_results <- list()

for (level in tax_levels) {
  aldex_result <- run_aldex(phyloseq_obj, group = "enrich", norm = "CPM",
                            taxa_rank = level, p_adjust = "fdr")
  aldex_results[[level]] <- aldex_result
}




