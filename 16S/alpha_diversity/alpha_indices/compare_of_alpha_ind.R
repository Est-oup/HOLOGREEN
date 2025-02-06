
sum(c(0.98,rep(x = (1-0.98)/2000, 2000)))
sum(c(0.60,rep(x = (1-0.60)/2000, 2000)))

test_pop1 <- c(0.98,rep(x = (1-0.98)/2000, 2000))
test_pop2 <- c(0.60,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.03, rep(x = (1-0.98)/2000, 2000))


# Calcul de l'indice de Shannon
test_shannon1 <- vegan::diversity(test_pop1, index =  "shannon")
test_shannon2 <- vegan::diversity(test_pop2, index =  "shannon")

# Calcul de l'indice de Simpson
test_simpson1 <- vegan::diversity(test_pop1, index = "simpson")
test_simpson2 <- vegan::diversity(test_pop2, index = "simpson")

