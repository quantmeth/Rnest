#' Next Eigenvalue Sufficiency Test (NEST)
#'
#' @description \code{nest} is used to identify the number of factors to retain in exploratory factor analysis.
#'
#' @param .data a data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param n the number of cases (subjects, participants, or units) if a covariance matrix is supplied in \code{.data}.
#' @param nreps the number of replications to derive the empirical probability distribution of each eigenvalue. Default is 1000.
#' @param alpha a vector of type I error rates or \code{(1-alpha)*100\%} confidence intervals. Default is .05.
#' @param max.fact an optional maximum number of factor to extract. Default is \code{NULL}, so the maximum number possible.
#' @param method a method used to compute loadings and uniquenesses. Four methods are implemented in \code{Rnest} : maximum likelihood \code{method = "ml"} (default), regularized common factor analysis \code{method = "rcfa"}, minimum rank factor analysis \code{method = "mrfa"}, and principal axis factoring \code{method = "paf"}. See details for custom methods.
#' @param na.action how should missing data be removed. \code{"na.omit"} removes complete rows with at least one single missing data. \code{"fiml"} uses full information maximum likelihood to compute the correlation matrix. Other options are \code{"everything"}, \code{"all.obs"}, \code{"complete.obs"}, \code{"na.or.complete"}, or \code{"pairwise.complete.obs"}. Default is \code{"fiml"}.
#' @param ... arguments for \code{method} that can be supplied. See details.
#'
#' @details 
#' The Next Eigenvalues Sufficiency Test (NEST) is an extension of parallel analysis by adding a sequential hypothesis testing procedure for every \eqn{k = 0, ..., \code{max.fact}} factor until the hypothesis is not rejected. 
#' 
#' At \eqn{k = 0}, NEST and parallel analysis are identical. Both use an identity matrix as the correlation matrix. Once the first hypothesis is rejected, NEST uses a correlation matrix based on the loadings and uniquenesses of the \eqn{k^{th}} factorial structure. NEST then resamples \code{nreps} times the \eqn{k^{th}} eigenvalue of this new correlation matrix. NEST stops when the \eqn{k^{th}} eigenvalues is below the \eqn{1-\alpha}*100%} quantile threshold.  
#' 
#' There is four \code{method} already implemented in \code{nest} to estimate loadings and uniquenesses: maximum likelihood (\code{"ml"}; default), principal axis factoring (\code{"paf"}), regularized common factor analysis \code{method = "rcfa"}, and minimum rank factor analysis (\code{"mrfa"}). These functions use as arguments: \code{covmat}, \code{n}, \code{factors}, and \code{...} (supplementary arguments passed by \code{nest}). They return \code{loadings} and \code{uniquenesses}. Any other user-defined functions can be used as long as it is programmed likewise.
#'
#' The \code{method = "paf"} is the same as Achim's (2017) NESTip.
#'
#' @return \code{nest()} returns an object of class \code{nest}. The functions \code{summary} and \code{plot} are used to obtain and show a summary of the results.
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
#'   \item \code{nreps} - The number of replications used.
#'   \item \code{prob} - Probabilities of each factor.
#'   \item \code{Eig} - A list of simulated eigenvalues.
#' }
#'
#' @section Generic function:
#'
#' \code{plot.nest} Scree plot of the eigenvalues and the simulated confidence intervals for \code{alpha}.
#'
#' \code{loadings} Extract loadings. It does not overwrite \code{stat::loadings}.
#' 
#' \code{summary.nest} Summary statistics for the number of factors.
#' 
#' @author
#' P.-O. Caron
#'
#' @references
#' Achim, A. (2017). Testing the number of required dimensions in exploratory factor analysis. \emph{The Quantitative Methods for Psychology}, \emph{13}(1), 64-74. \doi{10.20982/tqmp.13.1.p064}
#'
#' @import stats
#' @import EFA.MRFA
#' @import fungible
#' @export  
#' 
#' @aliases NEST
#'
#' @examples
#' nest(ex_2factors, n = 100)
#' nest(mtcars)
nest <- function(.data, ..., n = NULL, nreps = 1000, alpha = .05, max.fact = NULL, method = "ml", na.action = "fiml"){
  
  if(!(is.matrix(.data) || is.data.frame(.data) || is.array(.data))){
    ls <- .data
    if(!is.null(ls$n)) n <- ls$n
    if(!is.null(ls$covmat)) {.data <- ls$covmat
    } else {
      .data <- ls$.data
    }
  }
  
  R <- prepare.nest(.data, n = n, na.action = na.action, ...)
  
  R$alpha <- alpha
  R$method <- method
  R$na.action <- na.action 
  R$nreps <- nreps
  R$Eig <- list()
  R$prob <- numeric()
  test.eig <- rep(TRUE, length(R$alpha))
  
  nf <- .nf(alpha)
  nfactors <- nf$nfactors
  CI <- nf$CI
  R$alpha <- nf$alpha
  R$convergence <- TRUE
  if(is.null(max.fact)) max.fact <- .max.fact(ncol(R$cor))
  
  for (i in 0:max.fact){
    if(all(!test.eig)) {
      
      break
      
    } else {
      
      if (i == 0) {
        
        M <- diag(length(R$values))
        
      } else {
        
        tryCatch({
          M <- do.call(method[[1]],
                       list(covmat = R$cor,
                            n = R$n,
                            factors = i,
                            ...))
        }, error = function(e){
          R$convergence <- FALSE
        })
        
        if(!R$convergence) break
        
        M <- cbind(M$loadings, diag(sqrt(M$uniquenesses)))
        
      }
      
      Rep <- as.matrix(replicate(n = nreps,
                                 expr = .reig(n = R$n,
                                              M = M)))
      R$prob[i+1] <- sum(R$values[i+1] < Rep[i+1,]) / nreps 
      
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
  
  return(structure(c(list(nfactors = nfactors), R, list(stopping.rule = "Next Eigenvalue Sufficiency Test (NEST)")), class = "nest"))
  
}

# prepare.nest ####
prepare.nest <- function(data, n = NULL, na.action = "fiml", ...){
  
  data <- as.matrix(data)
  out <- list()
  
  if(isSymmetric(data, check.attributes = FALSE) ){
    
    if(all(diag(data) == 1)){
      out$cor <- data
    }else{
      out$cor <-  cov2cor(data)
    }
    
    if(is.null(n)){
      stop("Argument \"n\" is missing with covariance matrix.")
    } else {
      out$n <- n
    }
  } else {
    if(anyNA(data)){
      # TODO
      # to opt
      if(na.action == "fiml") {
        out$cor <- cor_nest(.data = data, ...)$covmat
      } else if (na.action == "na.omit"){
        out$cor <- cor(na.omit(data))
      } else {
        out$cor <- cor(data, use = na.action)
      }
    } else {
      out$cor <- cor(data)
    }
    
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

