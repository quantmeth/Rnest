#' YASR (new - 29-09-2022)
#'
#' @param data A data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param n.obs The number of cases (subjects, participants, or units) if a covariance matrix is supplied in \code{data}.
#' @param nreps The number of replications to simulate. Default is 1000.
#' @param alpha Type I error rates or \code{(1-alpha)*100\%} confidence intervals. Default is .10.
#' @param max.fact An optional maximum number of factor to extract.
#' @param ...  Arguments for \code{method} that can be supplied. See details.
#'
#' @return The number of factors
#' 
#' @import MASS
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' new(ex_4factors_corr, n.obs = 84)
#' }
new <- function(data, 
                 n.obs = NULL, 
                 nreps = 1000, 
                 alpha = .05, 
                 max.fact = max(which(dof(nv,0:nv)>(nv/2)))-1, 
                 ...){
  
  if(isSymmetric(data)){
    
    if(all(diag(data) == 1)){
      R <- data
    }else{
      R <-  cov2cor(data)
    }
    
    if(is.null(n.obs)){
      stop("argument \"n.obs\" is missing with covariance matrix")
    }
    
  } else {
    
    R <- cor(data)
    n.obs <- nrow(data)
    
  }
  
  nv <- ncol(data)
  eig <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  
  if((length(Re(eig)) != nv) && (sum(eig) != nv)){
    stop("Correlation matrix is not positive semi definite")
  }
  
  dl <- dof(nv, 0) # choose(ncol(R), 2)
  #cl <- makeCluster(detectCores()[1]-1)
  #registerDoParallel(cl)
  ztest0 <- list()
  
  res <- data.frame(p = rep(0, max.fact),
                    a = 0,
                    z = 0,
                    pz = 0)
  #res <- data.frame(mean = rep(0, max.fact),
  #                  aire = rep(0, max.fact))

  
  for(i in 1:(max.fact+1)){
    
    #T1 <- Sys.time();
    # rtest <- foreach(l = 1:nreps, .combine = cbind
    #                  ,.export = c("dof","fa","tril","diffR","gen_data")) %dopar% { #,.export = c("dof","fa","tril","diffR","gen_data")
    #                  diffR(RHO = R, fact = i-1, n = n.obs)
    #                }
    #;T2 <- Sys.time(); 
    rtest <- replicate(nreps, diffR(RHO = R, fact = i-1, n = n.obs, ...))
    #;T3 <- Sys.time()
    #T3-T2
    #T2-T1
    ztest0[[i]] <- colSums(rtest^2)/dl
    if(i > 1){
      res$z[i-1] <- (mean(ztest0[[i-1]]) - mean(ztest0[[i]])) / (sd(ztest0[[i-1]]) + sd(ztest0[[i-1]]))
      
      res$pz[i-1] <- round(pnorm(res$z[i-1], lower.tail = FALSE),4)
      
      #res$mean[i-1] <- as.logical(quantile(ztest0[[i-1]], alpha) < mean(ztest0[[i]]))
      
      res$p[i-1] <- round(sum(ztest0[[i]] > quantile(ztest0[[i-1]],alpha))/nreps,4)
        #sum(ztest0[[i-1]] < quantile(ztest0[[i]], 1-alpha)) / nreps
      res$a[i-1] <- alpha / 2 * (i-1)
      if(res$p[i-1] > res$a[i-1]) break
    }
  }
  
  res <- res[1:(i-1),]
  #list(nfactors = i-2, results = res)
  #stopCluster(cl)
  return(list(nfactors = i-2, results = res))
}

dof <- function(nv, i = 0:nv){
  (nv-i)*(nv-i-1)/2 - i
}

fa <- function(x, fact, ...){
  if(fact == 0){
    diag(ncol(x))
  } else {
    fanal <- factanal(covmat = x, n.obs = nrow(x), factors = fact, ...)
    FS <- cbind(fanal$loadings, sqrt(diag(fanal$uniquenesses)))
    FS %*% t(FS)
  }
}

tril <- function(x){x[lower.tri(x)]}

diffR <- function(RHO, fact, n = 20*ncol(RHO)+3, ...){ # n = 20*ncol(R)+3
   D <- gen_data(RHO, n)
  # add bootstrap eventually
  # D <- RHO[sample(nrow(RHO)),]
  if(fact == ncol(RHO)){
    W <- cor(D)
    Rs <- RHO
    denom = 1
  } else {
    W <- fa(cor(D), fact, ...)
    Rs <- cor(D) #cor(RHO)#cor(D)
    denom = 1
  }
  (atanh(tril(W))-atanh(tril(Rs))) * sqrt(n-3) / denom
}

gen_data <- function(R, n = 20*ncol(R)+3){
  MASS::mvrnorm(n = n, Sigma = R, mu = rep(0, ncol(R)))
}

