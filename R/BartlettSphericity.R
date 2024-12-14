#' Bartlett Sphericity Test
#'
#' @description \code{BartlettSphericity} tests if variables are orthogonal.
#'
#' @param R the correlation matrix.
#' @param n the sample size.
#'
#' @author 
#' André Achim (Matlab)
#' P.-O. Caron (R)
#'
#' @references
#' Bartlett, M. S. (1937). Properties of sufficiency and statistical tests. \emph{Proceedings of the Royal Statistical Society, Series A}, \emph{160}, 268–282
#'
#' @return The \eqn{\chi^2} test of the correlation matrix \code{R} with sample size \code{n}.
#' 
#' @import stats
#' @export
#'
#' @examples
#' BartlettSphericity(ex_4factors_corr, 42)
BartlettSphericity <- function(R, n){
  p <- ncol(R)
  chisq <- -((n-1)-(2*p-5)/6)*log(det(R))
  df <- p*(p-1)/2
  p <- pchisq(chisq, df, lower.tail = FALSE)
  return(stat = list(chisq = chisq,
                     df = df,
                     p = p))
}
