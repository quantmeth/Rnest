#' Nest Eigenvalue Sufficiency Test (NEST)
#'
#' @description \code{nest} is used to identify the number of factors to retain in exploratory factor analysis.
#'
#' @param data A data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param n The number of cases (subjects, participants, or units) if a covariance matrix is supplied in \code{data}.
#' @param nrep The number of replications to simulate. Default is 1000.
#' @param alpha A vector of type I error rates or \code{(1-alpha)*100\%} confidence intervals. Default is .05.
#' @param max.fact An optional maximum number of factor to extract. Default is \code{max.fact = ncol(data)}.
#' @param method A method used to compute loadings and uniquenesses. Four methods are implemented in \code{Rnest} : maximum likelihood \code{method = "ml"} (default), regularized common factor analysis \code{method = "rcfa"}, minimum rank factor analysis \code{method = "mrfa"}, and principal axis factoring \code{method = "paf"}. See details for custom methods.
#' @param ... Arguments for \code{method} that can be supplied. See details.
#'
#' @details 
#' The Next Eigenvalues Sufficiency Test (NEST) is an extension of parallel analysis by adding a sequential hypothesis testing procedure for every \eqn{k = 1, ..., p} factor until the hypothesis is not rejected. 
#' 
#' At \eqn{k = 1}, NEST and parallel analysis are identical. Both use an Identity matrix as the correlation matrix. Once the first hypothesis is rejected, NEST uses a correlation matrix based on the loadings and uniquenesses of the \eqn{k^{th}} factorial structure. NEST then resamples the eigenvalues of this new correlation matrix. NEST stops when the $k_1^2$ eigenvalues is within the confidence interval.  
#' 
#' There is two \code{method} already implemented in \code{nest} to extract loadings and uniquenesses: maximum likelihood (\code{"ml"}; default), principal axis factoring (\code{"paf"}), and minimum rank factor analysis (\code{"mrfa"}). The functions use as arguments: \code{covmat}, \code{n}, \code{factors}, and \code{...} (supplementary arguments passed by \code{nest}). They return \code{loadings} and \code{uniquenesses}. Any other user-defined functions can be used as long as it is programmed likewise.
#'
#' @return \code{nest} returns an object of class code{nest}. The functions \code{summary} and \code{plot} are used to obtain and show a summary of the results.
#' 
#' An object of class \code{nest} is a list containing the following components:
#' 
#' \itemize{
#'   \item \code{nfactors} - The number of factors to retains (one by \code{alpha}).
#'   \item \code{cor} - The supplied correlation matrix.
#'   \item \code{n} - The number of cases (subjects, participants, or units).
#'   \item \code{values} - The eigenvalues of the supplied correlation matrix.
#'   \item \code{alpha} - The type I error rate.
#'   \item \code{method} - The method used to compute loadings and uniquenesses.
#'   \item \code{Eig} - A list of simulated eigenvalues.
#' }
#'
#' @section Generic function:
#'
#' \code{plot.nest} Scree plot of the eigenvalues and the simulated confidence intervals for \code{alpha}.
#'
#' \code{loadings} Extract loadings. It does not overwrite \code{stat::loadings}.
#'
#' @author
#' P.-O. Caron
#'
#' @references
#' Achim, A. (2017). Testing the number of required dimensions in exploratory factor analysis. \emph{The Quantitative Methods for Psychology}, \emph{13}(1), 64-74. \url{https://doi.org/10.20982/tqmp.13.1.p064}
#'
#' @import stats
#' @import EFA.MRFA
#' @import fungible
#' @export  
#'
#' @examples
#'\dontrun{
#'nest(ex_2factors, n = 100)
#'nest(mtcars)
#'}
nest <- function(data, n = NULL, nrep = 1000, alpha = .05, max.fact = ncol(data), method = "ml", ...){
  
  R <- prepare.nest(data, n = n)

  R$alpha <- alpha
  R$method <- method
  R$Eig <- list()
  test.eig <- rep(TRUE, length(R$alpha))
  
  nf <- .nf(alpha)
  nfactors <- nf$nfactors
  CI <- nf$CI
  R$alpha <- nf$alpha
  
  for (i in 0:max.fact){
    if(all(!test.eig)) {
      
      break
      
    } else {
      
      if (i == 0) {
        
        M <- diag(length(R$values))
        
      } else {
        
        M <- do.call(method[[1]],
                     list(covmat = R$cor,
                          n = R$n,
                          factors = i, 
                          ...))
        M <- cbind(M$loadings, diag(sqrt(M$uniquenesses)))
        
      }
      
      Rep <- as.matrix(replicate(n = nrep,
                                 expr = .reig(n = R$n,
                                              M = M)))
      
      R$Eig[[i+1]] <- matrix(apply(X = Rep,
                                   MARGIN = 1,
                                   FUN = quantile,
                                   probs = 1-R$alpha),
                             nrow = length(R$alpha),
                             dimnames = list(CI))
      
      test.eig <- as.logical((R$Eig[[i+1]][,i+1] <= R$values[i+1]) * test.eig)
      nfactors <- nfactors + test.eig
    }
    
    
  }
  
  return(structure(c(list(nfactors = nfactors), R, list(stopping.rule = "Nest Eigenvalue Sufficiency Test (NEST)")), class = "nest"))
  
}

# prepare.nest ####
prepare.nest <- function(data, n = NULL){
  
  data <- as.matrix(data)
  out <- list()
  
  if(isSymmetric.matrix(data)){
    
    if(all(diag(data) == 1)){
      out$cor <- data
    }else{
      out$cor <-  cov2cor(data)
    }
    
    if(is.null(n)){
      stop("argument \"n\" is missing with covariance matrix")
    } else {
      out$n <- n
    }
  } else {
    out$cor <- cor(data)
    out$n <- nrow(data)
  }
  
  p <- ncol(data)
  out$values <- t(as.matrix(eigen(out$cor, symmetric = TRUE)$values))
  
  if((length(Re(out$values)) != p) && (sum(out$values) != p)){
    stop("Correlation matrix is not positive semi definite")
  }
  
  return(out)
}


# .reig ####
.reig <- function(M, n){
  d <- ncol(M)
  D <- M %*% matrix(rnorm(d * n), nrow = d)
  E <- eigen(cov(t(D)), symmetric = TRUE, only.values = TRUE)$values
  return(E)
}

# .paf ####
.paf <- function(covmat, factors, convergence = 1e-7, maxit = 500, ...){
  
  for (jj in 1:maxit){
    
    res <- eigen(covmat, symmetric = TRUE)
    ld <- res$vectors[,1:factors] %*% diag(sqrt(res$values[1:factors]), ncol = ncol(covmat))
    co <- rowSums(ld^2)
    
    # Check communalities
    
    diff <- diag(covmat) - co
    
    # Check if convergence is met
    if (all(abs(diff) < convergence)) break
    
    covmat <- covmat - diag(diff)
  }
  
  if(any(abs(diff) > convergence)) warning("Convergence not met.")
  return(list(loadings = ld, uniquenesses = 1-co))
  
}

# paf ####
paf <- function(covmat, factors, ...){
  fa <- .paf(covmat = covmat, factors = factors, ...)
  list(loadings = fa$loadings, uniquenesses = fa$uniquenesses)
}

# ml ####
ml <- function(covmat, n, factors, ...){
  fa <- factanal(covmat = covmat, n.obs = n, factors = c(factors), ...)
  list(loadings = fa$loadings[], uniquenesses = fa$uniquenesses)
}

mrfa <- function(covmat, n, factors, ...){
  fa <- EFA.MRFA::mrfa(SIGMA = covmat, dimensionality = factors, ...)
  list(loadings = fa$A, uniquenesses = 1-fa$gam)
}

rcfa <- function(covmat, n, factors, ...){
  fa <- fungible::fareg(R = covmat, numFactors =  factors)
  list(loadings = fa$loadings, uniquenesses = 1-fa$h2)
}

.nf <- function(alpha){
  alpha <- sort(alpha)
  CI <- paste0((1 - alpha) * 100,"%")
  nfactors <- t(setNames(data.frame(matrix(0,
                                           ncol = length(alpha),
                                           nrow = 1)),
                         nm = CI))
  colnames(nfactors) <- "nfactors"
  out <- list(nfactors = nfactors,  CI = CI, alpha = alpha)
  return(out)
}

