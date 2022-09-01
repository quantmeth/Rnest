#' Title
#'
#' @param x an object of class "nest".
#' @param nfactors The number of factors to retains.
#' @param method A method used to compute loadings and uniquenesses.
#' @param ... further arguments to methods in "nest" or the \code{stats::loadings} function.
#'
#' @note See \code{stats::loadings} for the original documentation.
#' @export
#'
#' @examples 
#' results <- nest(ex_2factors, n = 100)
#' loadings(results)
loadings <- function(x, nfactors = x$nfactors, method = x$method, ...){
  if(class(x) == "nest"){
    
    if(any(nfactors == 0)) stop("The number of factor is 0")
    if(length(nfactors) > 1) stop("Choose a number of factors to extract loadings.")
    M <- do.call(method[[1]],
                 list(covmat = x$cor,
                      n = x$n,
                      factors = nfactors, 
                      ...))
    
    return(loadings = M$loadings)
    
  } else {
    
    stats::loadings(x, ...)
    
  }
}
