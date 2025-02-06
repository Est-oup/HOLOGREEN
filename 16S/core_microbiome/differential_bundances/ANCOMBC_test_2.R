

#http://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html


library(ANCOMBC)
library(tidyverse)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))


#Run the ancombc function

out = ancombc2(data = NULL, assay_name = NULL,
                            tax_level = "Genus", phyloseq = physeq_16S_nit,
                            formula = "type_enrich",
                            p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
                            group = "biofilm_SW", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                            max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                            n_cl = 1, verbose = TRUE)

res = out$res
res_global = out$res_global

# ancombc also supports importing data in phyloseq format
# tse_alt = agglomerateByRank(tse, "Family")
# pseq = makePhyloseqFromTreeSummarizedExperiment(tse_alt)
# out = ancombc(data = NULL, assay_name = NULL,
#               tax_level = "Family", phyloseq = pseq,
#               formula = "age + region + bmi",
#               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#               group = "region", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
#               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
#               n_cl = 1, verbose = TRUE)



View(data(atlas1006, package = "microbiome"))
data(atlas1006, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)

ancombc2(physeq_16S_nit,tax_level = "Genus", group = ~ type_enrich)



set.seed(123)
output = ancombc2(data = physeq_16S_nit, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "type_enrich", rand_formula = NULL,
                  p_adj_method = "holm", 
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "cat_cov", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)

res_prim = output$res
res_pair = output$res_pair