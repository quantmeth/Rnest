#' Yet Another Stopping Rule (new) chi2(atanh)
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
#' YASRnew(ex_4factors_corr, n.obs = 42)
YASRnew <- function(data, 
                 n.obs = NULL,
                 alpha = .05, 
                 max.fact = max(which(dof(nv)>(nv/2)))-1,
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
  
  RR <- cor(data)
  Rt <- RR[lower.tri(RR)]
  res <- as.data.frame(matrix(0, ncol = 1,
                              nrow = max.fact + 1,
                              dimnames = list(0:max.fact,
                                              c("stat"))))
  for(i in 1:(max.fact)){
    if(i == 1){
      RRR <- diag(nv)
    } else{
      fa <- factanal(covmat = RR, n.obs = n.obs, factors = i-1,...) #,...
      ld <- cbind(fa$loadings, diag(sqrt(fa$uniquenesses)))
      RRR <- ld %*% t(ld)
    }
    R.test <- RRR[lower.tri(RRR)]
    z <- (atanh(Rt)-atanh(R.test))*sqrt(n.obs-3)
    res$stat[i] <- sum(z^2)
    res$dof[i] <- dof(nv,i-1)
    res$crit[i] <- qchisq(1-alpha,res$dof[i])
    res$pv[i] <- round(pchisq(res$stat[i], res$dof[i], lower = FALSE), 5)
    res$sig[i] <- (res$stat[i] >= res$crit[i])
    if(!res$sig[i]) break
  }
  res <- res[1:i,]
  # res$dof <- dof(nv,0:max.fact)
  # res$crit <- qchisq(1-alpgha,res$dof)
  # res$pv <- round(pchisq(res$stat, res$dof, lower = FALSE), 5)
  # res$sig <- res$stat >= res$crit 
  nfactors <- min(which(!res$sig))-1
  return(list(nfactors = nfactors, results = res))
}

dof <- function(nv, i = 0:nv){
  (nv-i)*(nv-i-1)/2 - i
}

gen_data <- function(R, n = 20*ncol(R)+3){
  MASS::mvrnorm(n = n, Sigma = R, mu = rep(0, ncol(R)))
}

#multivariate case
gen_R <- function(R, n){
  S <- cor(gen_data(R, n))
  S[lower.tri(S)]
}
