#' Print Loadings in NEST
#'
#' @param x An object of class "nest".
#' @param nfactors The number of factors to retains.
#' @param method A method used to compute loadings and uniquenesses.
#' @param ... Further arguments to methods in "nest" or the \code{stats::loadings} function.
#'
#' @note See \code{stats::loadings} for the original documentation.
#' @export
#'
#' @examples 
#' \dontrun{
#' results <- nest(ex_2factors, n = 100)
#' loadings(results)
#' }
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
