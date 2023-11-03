#' Simplify the the generation from a Multivariate Normal Distributions
#'
#' @description Speeds up the use of MASS::mvrnorm
#'
#' @param n the number of samples required.
#' @param R a positive-definite symmetric matrix specifying the covariance matrix of the variables.
#' @param mean an optinal vector giving the means of the variables. Default is 0.
#' @param ... Arguments for \code{MASS::mvrnorm()}, such as \code{tol}, \code{empirical}, and \code{EISPACK}.
#'
#' @return A data frame of size n by ncol(R).
#' @import MASS
#' @export
#'
#' @examples
#' set.seed(19)
#' R <- caron2016$mat1
#' mydata <- genr8(n = nrow(R)+1, R = R, empirical = TRUE)
#' round(mydata, 2)
#' round(cov(mydata), 2)
genr8 <- function(n = 1, R = diag(10), mean = rep(0, ncol(R)), ...){
  as.data.frame(MASS::mvrnorm(n = n,
                              mu = mean,
                              Sigma = R,
                              ...))
}
