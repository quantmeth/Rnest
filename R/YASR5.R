#' Yet Another Stopping Rule (.5)
#'
#' @param data A data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param n.obs The number of cases (subjects, participants, or units) if a covariance matrix is supplied in \code{data}.
#' @param alpha Type I error rate.
#' @param max.fact  An optional maximum number of factor to extract.
#' @param ... Arguments for methods.
#'
#' @return \code{YASR} returns the number of factors to retain.
#' @export
#'
#' @examples
#' YASR5(ex_4factors_corr, n.obs = 42)
YASR5 <- function(data, 
                 n.obs = NULL,
                 alpha = .05, 
                 max.fact = max(which(dof(ncol(data))>0))-2,
                 ...){
  # CHECK max.fact's dof < 0
  
  if(isSymmetric.matrix(data)){
    
    if(all(diag(data) == 1)){
      R <- data
    }else{
      R <- cov2cor(data)
    }
    
    if(is.null(n.obs)){
      stop("argument \"n.obs\" is missing with covariance matrix")
    }
  } else {
    R <- cor(data)
    n.obs <- nrow(data)
  }
  
  nv <- ncol(data)
  eig <- t(as.matrix(eigen(R, symmetric = TRUE, only.values = TRUE)$values))
  
  if((length(Re(eig)) != nv) && (sum(eig) != nv)){
    stop("Correlation matrix is not positive semi definite")
  }
  
  # add method
  res <- as.data.frame(matrix(0, ncol = 5,
                              nrow = max.fact + 1,
                              dimnames = list(0:max.fact,
                                              c("-2ll", "dof", "diff", "diff_dof", "pv"))))
  
  Rt <- R[lower.tri(R)]
  
  for(i in 0:max.fact){
    
    res$dof[i+1] <- dof(nv, i)
    
    if(i == 0){
      R.test <- diag(diag(R))
    } else {
      fa <- factanal(covmat = R, n.obs = n.obs, factors = i, ...) #,...
      ld <- cbind(fa$loadings, diag(sqrt(fa$uniquenesses)))
      R.test <- ld %*% t(ld)
    }
    
    R.test <- R.test[lower.tri(R.test)]
    z <- (log((1+Rt)/(1-Rt))/2 - log((1+R.test)/(1-R.test))/2) / sqrt(.5/(n.obs-3))
    res$`-2ll`[i+1] <- -2*sum(log(dnorm(z)))
    if(i > 0){
      res$diff <- res$`-2ll`[i] - res$`-2ll`[i+1]
      res$diff_dof[i+1] <- res$dof[i] - res$dof[i+1]
      res$pv[i+1] <- round(1-pchisq(q = res$diff,
                                    df = res$diff[i+1]), 4)
      
      if(res$pv[i+1]>alpha) break
      
    }
  }
  nfactors = i - 1 
  results <- res[1:(i+1),]
  
  # res$diff[-1] <- res$dof[-(max.fact+1)] - res$dof[-1]
  # res$pv[-1] <- round(1-pchisq(q = res$`-2ll`[-(max.fact+1)] - res$`-2ll`[-1], 
  #                              df = res$diff[-1]), 4)
  #results <- res
  #nfactors <- min(which(res$pv >= .05)) - 2
  #results <- res#list(YASR = res)#,
  #cormat = R,
  #values = eig,
  #loadings = fa$loadings[], # How to add loadings for nfactors (case = 0 or else)
  #alpha = alpha,
  #method = "method")
  
  return(list(nfactors = nfactors, results = results))
}

dof <- function(nv, i = 0:nv){
  (nv-i)*(nv-i-1)/2 - i
}
