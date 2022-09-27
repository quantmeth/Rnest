#' Super Nest Eigenvalue Sufficiency Test (SNEST)
#'
#' @description \code{nest} is used to identify the number of factors to retain in exploratory factor analysis.
#'
#' @param data A data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param n.obs The number of cases (subjects, participants, or units) if a covariance matrix is supplied in \code{data}.
#' @param nreps The number of replications to simulate. Default is 1000.
#' @param alpha A vector of type I error rates or \code{(1-alpha)*100\%} confidence intervals. Default is .05.
#' @param max.fact An optional maximum number of factor to extract. Default is \code{max.factors = ncol(data)}.
#' @param method A method used to compute loadings and uniquenesses. Two methods are implemented in \code{Rnest} : maximum likelihood \code{method = "lm"} (default) and principal axis factoring \code{method = "paf"}. See details for custom methods.
#' @param ... Arguments for \code{method} that can be supplied. See details.
#'
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
#' @export  
#'
#' @examples
#' \dontrun{
#' snest(ex_2factors, n = 100)
#' snest(mtcars)
#' }
snest <- function(data, 
                  n.obs = NULL, 
                  nreps = 1000, 
                  alpha = .05, 
                  max.fact = max(which(dof(nv)>(nv/2)))-1, 
                  method = "ml", ...){
  
  out <- list()
  
  if(isSymmetric.matrix(data)){
    
    if(all(diag(data) == 1)){
      RR <- data
    }else{
      RR <- cov2cor(data)
    }
    
    if(is.null(n.obs)){
      stop("argument \"n.obs\" is missing with covariance matrix")
    }
  } else {
    RR <- cor(data)
    n.obs <- nrow(data)
  }
  
  nv <- ncol(data)
  out$values <- t(as.matrix(eigen(RR, symmetric = TRUE, only.values = TRUE)$values))
  
  if((length(Re(out$values)) != nv) && (sum(out$values) != nv)){
    stop("Correlation matrix is not positive semi definite")
  }
  
  out$alpha <- sort(alpha)
  
  # Prepare YASR
  ndl <- dof(nv,0)
  Rt <- RR[lower.tri(RR)]
  rr <- replicate(nreps, gen_R(RR, n.obs))
  za <- (atanh(rr)-atanh(Rt)) * sqrt(n.obs-3)
  zs <- colSums(za^2)/ndl
  crit <- quantile(zs, 1-out$alpha) #1.582; 1.57
  
  #R <- prepare.nest(data, n = n)
  CI <- paste0((1 - out$alpha) * 100,"%")
  out$Eig <- list()
  test.eig <- rep(TRUE, length(out$alpha))
  test.cor <- rep(TRUE, length(out$alpha))
  
  # nfactors <- t(setNames(data.frame(matrix(0,
  #                                          ncol = length(out$alpha),
  #                                          nrow = 1)),
  #                        nm = CI))
  # colnames(nfactors) <- "nfactors"
  nfactors <- data.frame(matrix(0, ncol = 2, nrow = length(out$alpha),
                                dimnames = list(CI,c("eig","cor"))))
  
  
  for (i in 0:max.fact){
    if(all(!test.eig) && all(!test.cor)) {
      
      break
      
    } else {
      
      if (i == 0) {
        
        M <- diag(nv)
        
      } else {
        
        M <- do.call(method[[1]],
                     list(covmat = RR,
                          n = n.obs,
                          factors = i)) # ...
        M <- cbind(M$loadings, diag(sqrt(M$uniquenesses)))
        
      }
      
      if(any(test.eig)){
        Rep <- as.matrix(replicate(n = nreps,
                                   expr = .reig(n = n.obs,
                                                M = M)))
        
        out$Eig[[i+1]] <- matrix(apply(X = Rep,
                                       MARGIN = 1,
                                       FUN = quantile,
                                       probs = 1-out$alpha),
                                 nrow = length(out$alpha),
                                 dimnames = list(CI))
        test.eig <- as.logical((out$Eig[[i+1]][,i+1] <= out$values[i+1]) * test.eig)
        nfactors$eig <- nfactors$eig + test.eig
      }
      
      if(any(test.cor)){
        R.test <- (M%*%t(M))[lower.tri(RR)]
        z <- (atanh(Rt)-atanh(R.test))*sqrt(n.obs-3)
        stat <- sum(z^2)/ndl
        test.cor <- as.logical((stat > crit) * test.cor)
        nfactors$cor <- nfactors$cor + test.cor
      }
      
      nfactors
      
    }
    
    
  }
  
  structure(c(list(nfactors = nfactors), out), class = "nest")
  
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
    co = rowSums(ld^2)
    
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

#multivariate case
gen_R <- function(R, n){
  S <- cor(gen_data(R, n))
  S[lower.tri(S)]
}

dof <- function(nv, i = 0:nv){
  (nv-i)*(nv-i-1)/2 - i
}

gen_data <- function(R, n){
  MASS::mvrnorm(n = n, Sigma = R, mu = rep(0,ncol(R)))
}

