#' Parallel analysis
#'
#' @param data data.frame.
#' @param n number of subjects.
#' @param p number of variables.
#' @param nrep number of replications.
#' @param alpha type I error rate.
#' @param crit Critical values to compare the eigenvalues.
#' @param ... Other arguments
#'
#' @return nfactors (if data is supplied) and sampled eigenvalues
#' @export
#'
#' @examples
#' pa(ex_2factors, n = 42)
#' E <- pa(n = 10, p = 2, nrep = 5)
pa <- function(data = NULL, n = NULL, p = NULL, nrep = 1000, alpha = .05, crit = NULL, ...){
  
  eig <- NULL
  
  if(inherits(data, "nest")){
    crit <- data$Eig[[1]]
    alpha <- data$alpha
    eig <- data$values
    data <- NULL
  }
  
  nf <- .nf(alpha)
  
  if(!is.null(data)){
    R <- prepare.nest(data, n = n)
    eig <- R$values
    n <- R$n
    p <- length(eig)
    #nfactors = nf$nfactors
  }
  
  if(is.null(crit)){
    
    E <- replicate(nrep,
                   eigen(cor(matrix(rnorm(n * p),
                                    ncol = p,
                                    nrow = n)), 
                         only.values = TRUE)$values)
    
    crit <- apply(E, 
                  MARGIN = 1, 
                  FUN = quantile, 
                  probs = (1-nf$alpha))
  } else {
    
    E <- NULL
    
  }
  
  
  if(!is.null(eig)){
    
    crit <- matrix(crit, nrow = length(alpha)); rownames(crit) <- nf$CI
    nf$nfactors[] <- apply(crit, MARGIN = 1, function(crit) {max(c(0, which(eig >= crit)))})
    
  } else {
    
    nfactors <- NULL
    
  }
  
  return(structure(list(nfactors = nf$nfactors,
                        values = eig,
                        Eig = list(crit),
                        sample.eig = E,
                        alpha = nf$alpha,
                        stopping.rule = "Parallel Analysis"), class = "nest"))
  
}
