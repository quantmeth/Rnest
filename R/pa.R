pa <- function(data, n = nrow(data), p = ncol(data), nreps = 1000, alpha = .05){
  eig <- eigen(cor(data), only.values = TRUE)$values
  E <- replicate(nreps,
                 eigen(cor(matrix(rnorm(n * p),
                                  ncol = p,
                                  nrow = n)), 
                       only.values = TRUE)$values)
  crit <- apply(E, 
                MARGIN = 1, 
                FUN = quantile, 
                probs = (1-alpha))
  
  nfactors = max(c(0, which(eig >= crit)))
  return(list(nfactors = nfactors))
}
