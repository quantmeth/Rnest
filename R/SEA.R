#' Significant Eigenvalues Assessment
#'
#' @param data A data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param n.obs The number of cases (subjects, participants, or units) if a covariance matrix is supplied in \code{data}.
#' @param alpha Type I error rate
#' @param nreps The number of replications to simulate. Default is 1000.
#'
#' @return \code{SEA} returns the number of factors to retain.
#' @export
#' 
#' @import utils
#' @import MASS
#'
#' @examples
#' SEA(ex_2factors, n.obs = 50)
SEA <- function(data, n.obs, alpha = .05, nreps = 5000) {
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
  
  Rr <- reduit(R = R, k = 0)$RR
  Rb <- blanchit(Rr)
  eig <- eigen(Rb, only.values = TRUE)$values
  
  res <- matrix(eig, nrow = nv, ncol = nreps)
  
  # parallel analysis ####
  tbl <- replicate(n = nreps,
                   expr = eigen(cor(MASS::mvrnorm(n = n.obs,
                                                  Sigma = diag(nv),
                                                  mu = rep(0,nv))), 
                                only.values = TRUE)$values)
  ####
  
  pr <- (1 + rowSums(eig <= tbl)) / (nreps + 1)
  nf <- sum(pr < alpha)
  # nf <- which(pr > alpha)[1]-1 (?)
  
  # return(list(nf = nf,
  #             pr = pr,
  #             eig = eig,
  #             R = Rb))
  return(nfactors = nf)
}

critreduitA <- function(a, Rr, k){
  eig <- 0
  RR <- Rr + diag(a)
  d <- diag(RR)
  if(max(d) >= .999){
    crit <- 9999
  } else {
    if(k >= 0){
      d <- sqrt(1-d)
      d <- diag(1/d)
      RR <- d %*% RR %*% d
    } else {
      k <-  -k
    }
    eig <- eigen(RR, only.values = TRUE)$values
    if(tail(eig, 1) < 1e-9){
      crit <- 999
    }  else {
      crit <- max(c(0,sum(eig[-c(1:(k-1))])))
    }
  }
  return(crit = crit)
}

blanchit <- function(Rr){
  d <- diag(Rr)
  d <- sqrt(1-d)
  d <- diag(1/d)
  d %*% Rr %*% d
}

reduit <- function(R, k){
  
  # Étape 0
  if(k == 0){
    # reduitM
    Rr <- R - diag(1/diag(solve(R))) 
  } else {
    # reduitR
    fa <- factanal(covmat = R, factors = k)
    Rr <- fa$loadings %*% t(fa$loadings)
  }
  
  # Étape 1
  res <- eigen(Rr)
  eig <- res$values; vec <- res$vectors
  j <- which(eig < 0)
  lim <- 2000
  # iter <- matrix(0, ncol = v, 2)
  for(i in j){
    
    crit <- eig[i]
    
    for(it in 1:lim){
      
      if(is.complex(eig[i]) || (crit >= -1e-10)) break
      
      d <- sqrt(abs(eig[i])) * vec[,i]
      Rr <- apply(Rr+diag(d), 
                  MARGIN = c(1, 2), 
                  function(x) min(c(x, 1)))
      
      res <- eigen(Rr)
      eig <- res$values; vec <- res$vectors
      
      if(eig[i]<crit){
        Rr <- Rr-diag(2*d)
        res <- eigen(Rr)
        eig <- res$values; vec <- res$vectors
        if(is.complex(eig)) break
      }
      
      crit <- eig[i]
      
      if(it == lim){
        # Cas extrême
        # reduitA1
        # Message
        warnings("Limit was achieved. Convergence is not guaranteed. Results may be inaccurate.")
      }
    }
  }
  
  # Étape 2
  if (k >= 0){
    k1 <- k + 1
    k1s <- k1
  } else {
    k1 <- -k + 1
    k1s <- -k1
  }
  
  delta <- rep(0, ncol(R))
  crit1 <- critreduitA(delta, Rr, k1s)
  
  res.opt <- optim(delta,  
                   fn = critreduitA,
                   Rr = Rr,
                   k = k1s)
  RR <- Rr + diag(delta)
  return(list(RR = RR,
              eig = eigen(RR, only.values = TRUE)$values,
              delta = res.opt$par,
              crit = c(crit1, crit)))
}