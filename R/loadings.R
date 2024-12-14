#' Print Loadings in NEST
#'
#' @param x an object of class "nest".
#' @param nfactors the number of factors to retains.
#' @param method a method used to compute loadings and uniquenesses.
#' @param ... further arguments to methods in "nest" or the \code{stats::loadings} function.
#' 
#' @return A \eqn{p \times k} matrix containing loadings where \eqn{p} is the number of variables and \eqn{k} is the number of factors (\code{nfactors}).
#' 
#' @note See \code{stats::loadings} for the original documentation.
#' @export
#'
#' @examples 
#' results <- nest(ex_2factors, n = 100)
#' loadings(results)
loadings <- function(x, nfactors = x$nfactors, method = x$method, ...){
  if(inherits(x, "nest")){
    
    if(any(nfactors == 0)) stop("The number of factor is 0")
  
    L <- list()
    
    for(i in 1:length(x$alpha)){
    M <- do.call(method[[1]],
                 list(covmat = x$cor,
                      n = x$n,
                      factors = nfactors[i]))
    L[[i]] <- M$loadings
    
    }
    names(L) <- rownames(nfactors)
    
    return(loadings = L)
    
  } else {
    
    stats::loadings(x, ...)
    
  }
}
