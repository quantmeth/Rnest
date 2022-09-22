#' Parallel analysis
#'
#' @param data data.frame
#' @param n number of subjects
#' @param p number of variables
#' @param nreps number of replications
#' @param alpha type I error rate
#' @param ... Other arguments
#'
#' @return nfactors (if data is supplied) and sampled eigenvalues
#' @export
#'
#' @examples
#' pa(ex_2factors)
#' E <- pa(n = 10, p = 2, nrep = 5)
pa <- function(data = NULL, n = nrow(data), p = ncol(data), nreps = 1000, alpha = .05, ...){
  
  E <- replicate(nreps,
                 eigen(cor(matrix(rnorm(n * p),
                                  ncol = p,
                                  nrow = n)), 
                       only.values = TRUE)$values)
  
  if(is.null(data)){
    nfactors = NULL
  } else {
    eig <- eigen(cor(data), only.values = TRUE)$values
    
    crit <- apply(E, 
                  MARGIN = 1, 
                  FUN = quantile, 
                  probs = (1-alpha))
    
    nfactors = max(c(0, which(eig >= crit)))
  }
  return(list(nfactors = nfactors,
              sample.eig = E))
}
