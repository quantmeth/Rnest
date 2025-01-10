#' Ledermann bound.
#'
#' @description Returns the maximum number of latent factors in a factor analysis model.
#'
#' @param p The number of variables.
#'
#' @author 
#' André Achim (Matlab)
#' 
#' P.-O. Caron (R)
#' 
#' @references 
#' Ledermann, W. (1937). On the rank of reduced correlation matrices in multiple factor analysis. \emph{Psychometrika}, {2}, 85–93.
#' 
#' @return The Ledermann bound.
#' @export
#'
#' @examples
#' Ledermann(ncol(ex_4factors_corr))
Ledermann <- function(p){
  max.fact <- floor((2*p+1-sqrt(8*+p+1))/2)
  return(max.fact = max.fact)
}
.max.fact <- function(nv){
  # CARON
  dof <- function(nv, i = 0:nv){
    (nv-i)*(nv-i-1)/2 - i
  }
  max(which(dof(nv,0:nv)>(nv/2)))-1
}

