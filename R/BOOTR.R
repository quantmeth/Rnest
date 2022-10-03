#' BOOTR
#'
#' @param D data
#' @param alpha Type I error
#' @param nreps Number of replications
#' @param max.fact max number of factors
#'
#' @return The number of factors
#' @export
#'
#' @examples
#' \dontrun{
#' BOOT(ex_2factors)
#' }
BOOTR <- function(D, alpha = .05, nreps = 1000, max.fact = max.fact <- max(which(dof(nv,0:nv)>(nv/2)))-1){
  nv <- ncol(D)
  z <- list()
  res <- data.frame(mean = rep(0, max.fact + 1),
                    crit = 0,
                    crit2 = 0,
                    q = 0,
                    pv = 0,
                    sig = 0)
  
  for(h in 0:max.fact){
    if(h == 0){
      bts0 <- replicate(nreps, mean(df(tril(cor(gen_data(R = cor(D),
                                                                n = nrow(D)))),
                                              tril(cor(D)),
                                              nrow(D))^2))
      
      crit <- quantile(bts0, .95)
    }  
      
      W <- fa(D, h)
      bts <- replicate(nreps, bootR(D, W))
      z[[h+1]] <- data.frame(test = paste0("H",h),
                             z = bts)

    
    res[h+1, ] <- c(mean(bts),
                    crit, 
                    1+2/sqrt(nrow(D)-3)*1.64, 
                    quantile(bts, alpha),
                    sum(bts<1) / (nreps + 1),
                    quantile(bts, alpha) < 1)
    
    if(res$sig[h+1]) break
    
  }
  res <- res[1:(h+1),]
  return(list(nfactor = h, results = res, data = z))
}

df <- function(Rt, RR, n) (atanh(tril(Rt))-atanh(tril(RR))) * sqrt(n-3)

between <- function(x, lu){
  (x >= lu[1]) & (x <= lu[2])
}

fa <- function(x, fact,...){
  if(fact == 0){
    diag(ncol(x))
  } else {
    fa <- factanal(covmat = cor(x), n.obs = nrow(x), factors = fact,...)
    FS <- cbind(fa$loadings, sqrt(diag(fa$uniquenesses)))
    FS %*% t(FS)
  }
}

dof <- function(nv, i = 0:nv){
  (nv-i)*(nv-i-1)/2 - i
}

gen_data <- function(R, n = 20*ncol(R)+3){
  MASS::mvrnorm(n = n, Sigma = R, mu = rep(0, ncol(R)))
}

tril <- function(x){x[lower.tri(x)]}

bootR <- function(D,W){
  D.boot <- D[sample(nrow(D), replace  = TRUE),]
  n <- nrow(D)
  R <- cor(D.boot)
  mean(df(R, W, n)^2)
}  

