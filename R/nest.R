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
#' @param missing how should missing data be removed. \code{"fiml"} uses full information maximum likelihood to compute the correlation matrix. Other options are \code{"ml"}, \code{"pairwise"}, \code{"listwise"}. Default is \code{"fiml"}.
#' @param cluster a (single) variable name in the data frame defining the clusters in a two-level dataset.
#' @param ordered a character vector to treat the variables as ordered (ordinal) variables. If TRUE, all observed endogenous variables are treated as ordered (ordinal).
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
#' @importFrom fungible fareg
#' @export  
#' 
#' @aliases NEST
#'
#' @examples
#' nest(ex_2factors, n = 100)
#' nest(mtcars)
nest <- function(.data, ..., n = NULL, nreps = 1000, alpha = .05, max.fact = NULL, method = "ml", missing = "fiml", cluster = NULL, ordered = NULL){
  
  R <- list()
  
  if(!(is.matrix(.data) || is.data.frame(.data) || is.array(.data))){
    ls <- .data
    if(!is.null(ls$n)) n <- ls$n
    if(!is.null(ls$covmat)) {.data <- ls$covmat
    } else {
      .data <- ls$.data
    }
  }
  
  if(!(cluster %in% colnames(.data)) && !is.null(cluster)){
    cluster <- NULL
    warning("Cluster variable ", cluster," is not found in .data. Cluster is ignored.")
    cluster <- NULL
  }
  if((!any(ordered %in% colnames(.data))) && ((!is.null(ordered)) * (!is.logical(ordered)))){
    if(length(ordered) == 1){
      warning("The ordered variable ", ordered ," was not found in .data. Ordered is ignored.")
    } else {
      warning("The ordered variables ", ordered ," were not found in .data. Ordered is ignored.")
    }
    ordered <- NULL
  }
  
  if(nrow(.data) == ncol(.data)){
    
    if(!is.null(cluster)){
      .d2 <- as.matrix(.data[which(!(colnames(.data) %in% cluster)), which(!(colnames(.data) %in% cluster))])
    } else {
      .d2 <- as.matrix(.data)
    }
    
    if(isSymmetric(.d2, check.attributes = FALSE)){
      
      if(!is.null(cluster)) warning("cluster is ignored with covariance matrix.")
      if(!is.null(cluster)) warning("ordered is ignored with covariance matrix.")
      
      if(is.null(n)){
        stop("Argument \"n\" is missing with covariance matrix.")
      } else {
        R$n <- n
      }
      if(all(diag(.data) == 1)){
        R$cor <- .d2
      } else {
        R$cor <-  cov2cor(.d2)
      }
    }
  }
  
  if(anyNA(.data)){
    if(missing %in% c("ml", "fiml", "two.stage", "robust.two.stage","listwise","pairwise","direct",
                      "ml.x","fiml.x","direct.x","available.cases","doubly.robust")){
      warning("The value of missing does not match a possible argument. It is overwritten to \"listwise\".")
      missing <- "listwise"
    }
  }
  
  if(is.null(n) || is.null(R$n)){
    R$cor <- cor_nest(.data = .data, 
                      ordered = ordered, 
                      missing = missing, 
                      cluster = cluster, 
                      ...)$covmat
    
    R$n <- nrow(.data)
    if(R$n != n && !is.null(n)) warning("The value of n does not match the number of rows in .data. It is overwritten to ",R$n,".")
    
  } else {
    
    R$n <- n
    
  }
  
  nv <- ncol(R$cor)
  R$values <- t(as.matrix(eigen(R$cor, symmetric = TRUE)$values))
  
  if((length(Re(R$values)) != nv) && (sum(R$values) != nv) && all(R$values > 0)){
    stop("Correlation matrix is not positive semi definite.")
  }
  
  #R <- prepare.nest(.data, n = n, missing = missing, ...)
  
  R$alpha <- alpha
  R$method <- method
  R$na.action <- missing
  R$nreps <- nreps
  R$Eig <- list()
  R$prob <- numeric()
  test.eig <- rep(TRUE, length(R$alpha))
  
  nf <- .nf(alpha)
  nfactors <- nf$nfactors
  CI <- nf$CI
  R$alpha <- nf$alpha
  R$convergence <- TRUE
  if(is.null(max.fact)) {
    max.fact <- .max.fact(ncol(R$cor))
  } else {
    if(max.fact > .max.fact(ncol(R$cor))){
      max.fact <- .max.fact(ncol(R$cor))
      warning("max.fact is too high. It is overwritten to ", max.fact,".")
    }
  }
  
  for (i in 0:max.fact){
    #cat(i)
    if(all(!test.eig)) {
      
      break
      
    } else {
      
      if (i == 0) {
        
        M <- diag(length(R$values))
        
      } else {
        
        R$convergence <- tryCatch({
          M <- do.call(method[[1]],
                       list(covmat = R$cor,
                            n = R$n,
                            factors = i,
                            ...
                       ))
          TRUE
        }, error = function(e){
          FALSE
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
# prepare.nest <- function(data, n = NULL, na.action = "fiml", ...){
#   
#   data <- as.matrix(data)
#   out <- list()
#   
#   if(isSymmetric(data, check.attributes = FALSE) ){
#     
#     if(all(diag(data) == 1)){
#       out$cor <- data
#     }else{
#       out$cor <-  cov2cor(data)
#     }
#     
#     if(is.null(n)){
#       stop("Argument \"n\" is missing with covariance matrix.")
#     } else {
#       out$n <- n
#     }
#   } else {
#     if(anyNA(data)){
#       # TODO
#       # to opt
#       if(na.action == "fiml") {
#         out$cor <- cor_nest(.data = data, ...)$covmat
#       } else if (na.action == "na.omit"){
#         out$cor <- cor(na.omit(data))
#       } else {
#         out$cor <- cor(data, use = na.action)
#       }
#     } else {
#       out$cor <- cor(data)
#     }
#     
#     out$n <- nrow(data)
#   }
#   
#   p <- ncol(data)
#   out$values <- t(as.matrix(eigen(out$cor, symmetric = TRUE)$values))
#   
#   if((length(Re(out$values)) != p) && (sum(out$values) != p)){
#     stop("Correlation matrix is not positive semi definite")
#   }
#   
#   return(out)
# }


# .reig ####
.reig <- function(M, n){
  d <- ncol(M)
  D <- M %*% matrix(rnorm(d * n), nrow = d)
  E <- eigen(cov(t(D)), symmetric = TRUE, only.values = TRUE)$values
  return(E)
}

# .paf ####
.paf <- function(covmat, factors, convergence = 1e-7, maxit = 500, ...){
  
  S <- covmat
  
  for (jj in 1:maxit){
    
    res <- eigen(S, symmetric = TRUE)
    ld <- res$vectors[,1:factors] %*% diag(sqrt(res$values[1:factors]), ncol = factors)
    co <- rowSums(ld^2)
    
    # Check communalities
    
    diff <- diag(S) - co
    
    # Check if convergence is met
    if (all(abs(diff) < convergence) && all(1-co >= convergence)) break
    
    S <- S - diag(diff)
  }
  
  if(any(abs(diff) > convergence)) stop("Convergence not met in the ",factors," factor model.\n")
  if(any(1-co <= convergence)) stop("Communality overflow in the ",factors," factor model.\n")
  if(maxit == jj) stop("Convergence not met in the ",factors," factor model.\n")
  
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

