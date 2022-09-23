#' Yet Another Stopping Rule - Steiger
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
#' YASR(ex_2factors, n.obs = 50)
steiger <- function(data, 
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
  res <- as.data.frame(matrix(0, ncol = 3,
                              nrow = max.fact + 1,
                              dimnames = list(0:max.fact,
                                              c("stat", "dof", "pv"))))
  
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
    
    resu <- cortest.normal(R1 = R.test,
                           R2 = R, 
                           n1 = n.obs, 
                           n2 = n.obs)
    
    res[i+1] <- c(resu$stat, resu$dof, resu$pv)
    
    if(res$pv[i+1] > alpha) break
  }
  nfactors = i
  results <- res[1:(i+1),]
  
  return(list(nfactors = nfactors, results = results))
}

dof <- function(nv, i = 0:nv){
  (nv-i)*(nv-i-1)/2 - i
}
#psych::cortest.norma
cortest.normal <- function(R1,R2=NULL, n1=NULL,n2=NULL,fisher=TRUE) {
  cl <- match.call()
  if (dim(R1)[1] != dim(R1)[2]) {n1 <- dim(R1)[1] 
  message("R1 was not square, finding R from data")
  R1 <- cor(R1,use="pairwise")}
  
  if(!is.matrix(R1) ) R1 <- as.matrix(R1)  #converts data.frames to matrices if needed
  
  p <- dim(R1)[2]
  if(is.null(n1)) {n1 <- 100 
  warning("n not specified, 100 used") }
  if(is.null(R2)) { if(fisher) {R <- 0.5*log((1+R1)/(1-R1))
  R <- R*R} else {R <- R1*R1}
    diag(R) <- 0
    E <- (sum(R*lower.tri(R)))
    chisq <- E *(n1-3)
    df <- p*(p-1)/2
    p.val <- pchisq(chisq,df,lower.tail=FALSE)
  } else {         #end of 1 matrix test
    if (dim(R2)[1] != dim(R2)[2]) {n2 <- dim(R2)[1] 
    message("R2 was not square, finding R from data")
    R2 <- cor(R2,use="pairwise")}
    if(!is.matrix(R2) ) R2 <- as.matrix(R2)
    
    
    if(fisher) { 
      R1 <- 0.5*log((1+R1)/(1-R1)) 
      R2 <-  0.5*log((1+R2)/(1-R2))
      diag(R1) <- 0
      diag(R2) <- 0 }
    R <-  R1 -R2   #direct difference 
    R <- R*R
    if(is.null(n2)) n2 <- n1
    n <- (n1*n2)/(n1+n2)    #why do I do this? should it be 2 * (n1*n2)/(n1+n2)   or 
    #n <- harmonic.mean(c(n1,n2))  #no, actually this gives the right results
    E <- (sum(R*lower.tri(R)))
    chisq <- E *(n-3)
    df <- p*(p-1)/2
    p.val <- pchisq(chisq,df,lower.tail=FALSE)
  }
  result <- list(stat=chisq,pv=p.val,dof=df)
  return(result)
}
#version of 2008
#commented 2018


