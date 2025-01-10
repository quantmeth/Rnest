#' Parallel analysis
#'
#' @param data a data.frame.
#' @param n the number of subjects.
#' @param nv the number of variables.
#' @param nreps the number of replications.
#' @param alpha type I error rate.
#' @param crit critical values to compare the eigenvalues.
#' @param ... other arguments
#'
#' @aliases PA
#' 
#' @references 
#' Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. \emph{Psychometrika}, \emph{30}(2), 179–185. \doi{10.1007/BF02289447}
#' 
#' @return nfactors (if data is supplied) and sampled eigenvalues
#' @export
#'
#' @examples
#' pa(ex_2factors, n = 42)
#' pa(n = 10, nv = 2, nreps = 100)
pa <- function(data = NULL, n = NULL, nv = NULL, nreps = 1000, alpha = .05, crit = NULL, ...){
  
  eig <- NULL
 # cl <- "stoppingrules"
  
  if(inherits(data, "nest")){
    crit <- data$Eig[[1]]
    #alpha <- data$alpha
    eig <- data$values
    data <- NULL
    #cl <- "nest"
  }
  
  #nf <- .nf(alpha)
  
  if(!is.null(data)){
    if(!(is.matrix(data) || is.data.frame(data) || is.array(data))){
      ls <- data
      if(!is.null(ls$n)) n <- ls$n
      if(!is.null(ls$covmat)) {data <- ls$covmat
      } else {
        data <- ls$.data
      }
    }
    
    R <- prepare.nest(data, n = n, ...)
    eig <- R$values
    n <- R$n
    nv <- length(eig)
    #nfactors = nf$nfactors
  }
  
  if(is.null(crit)){
    
    E <- replicate(nreps,
                   eigen(cor(matrix(rnorm(n * nv),
                                    ncol = nv,
                                    nrow = n)), 
                         only.values = TRUE)$values)
    
    crit <- apply(E, 
                  MARGIN = 1, 
                  FUN = quantile, 
                  probs = (1-alpha))
  } else {
    
    E <- NULL
    
  }
  
  
  if(!is.null(eig)){
    
    crit <- matrix(crit, nrow = length(alpha))#; rownames(crit) <- nf$CI
    nfactors <- apply(crit, MARGIN = 1, function(crit) {max(c(0, which(eig >= crit)))})
    
  } else {
    
    nfactors <- NULL
    
  }
  
  return(structure(list(nfactors = nfactors,
                        values = eig,
                        Eig = crit,
                        #sample.eig = E,
                        alpha = alpha,
                        stopping.rule = "Parallel Analysis"), class =  "stoppingrules"))
}

#' @export
PA <- pa
