#' Plot rarefaction curves
#'
#' The function plot rarefaction curves...
#'
#' @return Print ggplot rarefaction curves
#' @param physeq A phyloseq class object, from which abundance
#' data are extracted
#' @param step Step size for sample size in rarefaction curves
#' @param label Default `NULL`. Character string.
#' The name of the variable to map to text labels on the plot.
#' Similar to color option but for plotting text.
#' @param color (Optional). Default ‘NULL’. Character string.
#' The name of the variable to map to colors in the plot.
#' This can be a sample variable (among the set returned by
#' ‘sample_variables(physeq)’ ) or taxonomic rank (among the set returned by
#' ‘rank_names(physeq)’).
#' The color scheme is chosen automatically by
#' ‘link{ggplot}’, but it can be modified afterward with an
#' additional layer using ‘scale_color_manual’.
#' @param plot Logical, should the graphic be plotted.
#' @param parallel should rarefaction be parallelized (using parallel framework)
#' @param se Default TRUE. Logical. Should standard errors be computed.

ggrare <- function(physeq, step = 10, label = NULL, color = NULL,
                   plot = TRUE, parallel = FALSE, se = TRUE) {
  
  require("ggplot2")
  
  x <- as(phyloseq::otu_table(physeq), "matrix")
  if (phyloseq::taxa_are_rows(physeq)) x <- t(x)
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(phyloseq::sample_data(physeq, FALSE))) {
    sdf <- as(phyloseq::sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if (length(color) > 1) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if (length(label) > 1) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data,
              aes_string(x = "Size", y = ".S",
                         group = "Sample", color = color)) +
    labs(x = "Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + geom_text(data = labels,
                       aes_string(x = "x", y = "y",
                                  label = label, color = color),
                       size = 4, hjust = 0)
  }
  
  p <- p + geom_line()
  
  if (se) { ## add standard error if available
    p <- p +
      geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se",
                             color = NULL, fill = color),
                  alpha = 0.2)
  }
  
  if (plot) {
    plot(p)
  }
  
  invisible(p)
  
}

#' Check for normal distribution of alpha diversity measurements
#'
#'
#' @return Print ggplot rarefaction curves
#' @param rich A dataframe with alpha diversity indices as columns
#' and samples as rows as the output of pyloseq::estimate_richness()
#' @param file Output file name

indices_normality <- function(rich, nrow, ncol) {
  
  ### p-value < 0.05 means data failed normality test
  
  par(mfrow = c(nrow, ncol))
  
  for (i in names(rich)) {
    shap <- shapiro.test(rich[, i])
    qqnorm(rich[, i], main = i, sub = shap$p.value)
    qqline(rich[, i])
  }
  
  par(mfrow = c(1, 1))
}





####Fonction de détermination de différences significatives entres deux vateurs pour de la diversité alpha

# Récap de la fonction
#Cette fonction est constituée de la manière suivante :

#  Si le nombre de valeur du vecteur 1 et/ou du vecteur 2 est strictement inférieur à 15, je veux que la fonction réalise un test  non paramétrique de Wilcoxon MannWhitney "wilcox.test()" et informe l'utilisateur 
# de l'existence d'une différence ou pas entre les 2 vecteurs pour un interval de confiance de 0.05
# Si le nombre de valeur du vecteur 1 et du vecteur 2 est supérieur ou égale à 15, ou que le vecteur 1 et/ou le vecteur 2 soit inférieur ou égale à 30, je veux que la fonction teste la normalité des distributions 
# de chaque vecteur via le test "shapiro.test( )". Si la normalité des distributions n'est pas obtenue alors la fonction effectuera un test non paramétrique de Wilcoxon Mann-Whitney via le test "wilcoxon.test()". 
# Si la normalité est obtenue, la fonction réalisera un test de fisher via le test "var.test()". Dans les deux cas, un test non paramétrique de Wilcoxon MannWhitney "wilcox.test( )" sera ensuite réalisé. 
# Si les variances sont égales, la fonction effectuera un test paramétrique Student t test "t.test( …, var.equal=T)" et indiquera à l'utilisateur s'il y a une différence significative entre les deux vecteurs.
# Si les variances ne sont égales, la fonction effectuera un test paramétrique Student-Welch "t.test()" et indiquera à l'utilisateur s'il y a une différence significative entre les deux vecteurs.
# Si le nombre de valeur du vecteur 1 et du vecteur 2 est strictement supérieur à 30, la fonction effectuera un  test paramétrique Student-Welch "t.test()" et indiquera à l'utilisateur 
# s'il y a une différence significative entre les deux vecteurs.



significant.difference <- function(vec1, vec2,inter_conf) {
  
  
  n1 <- length(na.omit(vec1))
  n2 <- length(na.omit(vec2))
  
  print("La longueur du premier vecteur est")
  print(n1)
  print("La longueur du second vecteur est")
  print(n2)
  print("Vous avez choisi un interval de confiance de")
  print(inter_conf)
  
  
  
  if ((n1 < 15 | n2 < 15)  | (n1 < 15 & n2 < 15))
  {
    print("Vérifiez que la longueur des vecteurs 1 et 2 correspond bien à l'option : (n1 < 15 | n2 < 15)  | (n1 < 15 & n2 < 15)")
    wilcox_test <- wilcox.test(vec1, vec2, conf.int = inter_conf, correct = FALSE)
    if (wilcox_test$p.value < 0.05) {
      print("Il y a une différence statistiquement significative entre les deux vecteurs (Test non paramétrique de Wilcoxon-Mann-Whitney)")
      print("La p_value est de")
      print(wilcox_test$p.value)
      type_test <- "wilcox.test"
      assign("type_test",type_test, envir = .GlobalEnv)
      assign("p_value_fonction_significante_difference",wilcox_test$p.value, envir = .GlobalEnv)
    } 
    else {
      print("Il n'y a pas de différence statistiquement significative entre les deux vecteurs (Test non paramétrique de Wilcoxon-Mann-Whitney)")
      print("La p_value est de")
      print(wilcox_test$p.value)
      type_test <- "wilcox.test"
      assign("type_test",type_test, envir = .GlobalEnv)
      assign("p_value_fonction_significante_difference",wilcox_test$p.value, envir = .GlobalEnv)
      
    }
  }
  
  else if (((n1 >= 15 & n2 >= 15) & ((n1 <= 30 & n2 <= 30) | (n1 <= 30 | n2 <= 30))))
  {
    print("Vérifiez que la longueur des vecteurs 1 et 2 correspond bien à l'option : ((n1 >= 15 & n2 >= 15) & ((n1 <= 30 & n2 <= 30) | (n1 <= 30 | n2 <= 30)))")
    shapiro_vec1 <- shapiro.test(vec1)
    shapiro_vec2 <- shapiro.test(vec2)
    if (shapiro_vec1$p.value < 0.05 | shapiro_vec2$p.value < 0.05) {
      wilcox_test <- wilcox.test(vec1, vec2, conf.int = inter_conf, correct = FALSE)
      if (wilcox_test$p.value < 0.05) {
        print("vec 1 ou vec2 n'est pas distribué normalement ")
        print("Le résultat du test de chapiro pour le vec1 est")
        print(shapiro_vec1$p.value)
        print("Le résultat du test de chapiro pour le vec2 est")
        print(shapiro_vec2$p.value)
        print("Le test de Wilcoxon-Mann-Whitney a été réalisé et il y a une différence statistiquement significative entre les deux vecteurs")
        print(wilcox_test$p.value)
        type_test <- "wilcox.test"
        assign("type_test",type_test, envir = .GlobalEnv)
        assign("p_value_fonction_significante_difference",wilcox_test$p.value, envir = .GlobalEnv)
      } else {
        print("vec 1 ou vec2 n'est pas distribuées normalement ")
        print("Le résultat du test de chapiro pour le vec1 est")
        print(shapiro_vec1$p.value)
        print("Le résultat du test de chapiro pour le vec2 est")
        print(shapiro_vec2$p.value)
        print("Le test de Wilcoxon-Mann-Whitney a été réalisé et il n'y a pas de différence statistiquement significative entre les deux vecteurs")
        print(wilcox_test$p.value)
        type_test <- "wilcox.test"
        assign("type_test",type_test, envir = .GlobalEnv)
        assign("p_value_fonction_significante_difference",wilcox_test$p.value, envir = .GlobalEnv)
      }
      
    } else {
      
      print("Les données sont distribué normalement")
      print("Le résultat du test de chapiro pour le vec1 est")
      print(shapiro_vec1$p.value)
      print("Les données sont distribuées normalement")
      print("Le résultat du test de chapiro pour le vec2 est")
      print(shapiro_vec2$p.value)
      
      var_test <- var.test(vec1, vec2,conf.level = inter_conf)
      if (var_test$p.value < 0.05) {
        
        print("Le test de Fisher a montré que les variances étaient significativement différentes")
        print(var_test$p.value)
        t_test <- t.test(vec1, vec2, var.equal = F, conf.level = inter_conf)
        type_test <- "t.test Welch"
        assign("type_test",type_test, envir = .GlobalEnv)
        assign("p_value_fonction_significante_difference",t_test$p.value, envir = .GlobalEnv)
      } 
      else {
        print("Le test de Fisher a montré que les variances n'étaient pas significativement différentes")
        print(var_test$p.value)
        t_test <- t.test(vec1, vec2, var.equal = T, conf.level = inter_conf)
        type_test <- "t.test_(var=T)"
        assign("type_test",type_test, envir = .GlobalEnv)
        assign("p_value_fonction_significante_difference",t_test$p.value, envir = .GlobalEnv)
      }
      
      if (t_test$p.value < 0.05) {
        print("Le test paramétrique de student a montré une différence significative entre les deux vecteurs")
        print(t_test$p.value)
        
        wilcox_test <- wilcox.test(vec1, vec2, conf.int = inter_conf, correct = FALSE)
        
        if (wilcox_test$p.value < 0.05) {
          print("Le test de wilcoxon a été réalisé en complément et à également montré des différences significatives entre les deux vecteurs")
          print(wilcox_test$p.value)
        } 
        else {
          print("Le test de wilcoxon a été réalisé en complément et n'a en revenche pas montré de différences significatives entre les deux vecteurs")
          print(wilcox_test$p.value)
        }}
      
      else {
        print("Il n'y a pas de différence statistiquement significative entre les deux vecteurs")
        print(t_test$p.value)
        type_test <- "t.test"
        assign("type_test",type_test, envir = .GlobalEnv)
        assign("p_value_fonction_significante_difference",t_test$p.value, envir = .GlobalEnv)
        
        wilcox_test <- wilcox.test(vec1, vec2, conf.int = inter_conf, correct = FALSE)
        if (wilcox_test$p.value < 0.05) {
          print("En revanche le test de Wilcoxon montre une différence statistiquement significative entre les deux vecteurs")
          print(wilcox_test$p.value)
          type_test <- "wilcox.test"
          assign("type_test",type_test, envir = .GlobalEnv)
          assign("p_value_fonction_significante_difference",wilcox_test$p.value, envir = .GlobalEnv)
        } 
        else {
          print("Le test de Wilcoxon ne montre pas non plus de différence statistiquement significative entre les deux vecteurs")
          print(wilcox_test$p.value)
          type_test <- "wilcox.test"
          assign("type_test",type_test, envir = .GlobalEnv)
          assign("p_value_fonction_significante_difference",wilcox_test$p.value, envir = .GlobalEnv)
        }
      }
    }
  }
  if (n1 > 30 & n2 > 30) {
    print("Vérifiez que la longueur des vecteurs 1 et 2 correspond bien à l'option : (n1 > 30 & n2 > 30)")
    t_test <- t.test(vec1, vec2, var.equal = FALSE,conf.level = inter_conf)
    if (t_test$p.value < 0.05) {
      print("Il y a une différence statistiquement significative entre les deux vecteurs (Test paramétrique de Student-Welch)")
      print(t_test$p.value)
      type_test <- "t.test Welch"
      assign("type_test",type_test, envir = .GlobalEnv)
      assign("p_value_fonction_significante_difference",t_test$p.value, envir = .GlobalEnv)
    } else {
      print("Il n'y a pas de différence statistiquement significative entre les deux vecteurs (Test paramétrique de Student-Welch)")
      print(t_test$p.value)
      type_test <- "t.test Welch"
      assign("type_test",type_test, envir = .GlobalEnv)
      assign("p_value_fonction_significante_difference",t_test$p.value, envir = .GlobalEnv)
    }
  }
}



####fonction automatique pour les test ou il faut lui donner un objet phyloseq et la colonne du samdata qui sépare les groupe et lui dire quels test il faut réaliser  

test.alpha <- function(phyloseq_obj, variable_group, indices){
  mat_alpha <- cbind(alpha(phyloseq_obj), phyloseq_obj@sam_data[,variable_group])
  groups_alpha <- unique(mat_alpha[, ncol(mat_alpha)])
  results_alpha <- list()
  
  for (n in indices){
    res_alpha <- data.frame(matrix(data = NA, nrow = length(groups_alpha), ncol= length(groups_alpha)))
    colnames(res_alpha) <- groups_alpha
    rownames(res_alpha) <- groups_alpha
    
    for (i in combn(groups_alpha, 2, simplify=F)){
      var1alpha <- subset(mat_alpha, mat_alpha[, ncol(mat_alpha)] == i[1])
      var2alpha <- subset(mat_alpha, mat_alpha[, ncol(mat_alpha)] == i[2])
      print(paste0("Index: ", n, " Group 1: ", i[1], " Group 2: ", i[2]))
      significant.difference(var1alpha[[n]], var2alpha[[n]], 0.05)
      res_alpha[i[1],i[2]] <-  paste0(type_test , "  p-value: ", round(p_value_fonction_significante_difference,4))
      print("#######################TEST SUIVANT############################")
    }
    results_alpha[[paste0("res_alpha", n)]] <- res_alpha
  }
  assign("results_alpha",results_alpha, envir = .GlobalEnv)
  print(results_alpha)
}

