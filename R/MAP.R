#' Minimum average partial correlation (MAP)
#'
#' @param .data a data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param ... further argument for \code{cor_nest()}.
#'
#' @return The number of factors to retain.
#' @aliases map
#' @export
#'
#' @references
#' Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations. \emph{Psychometrika}, \emph{41}(3), 321-327. \doi{10.1007/BF02293557}
#'
#' @examples
#' D <- genr8(n = 42, R = ex_4factors_corr)
#' MAP(D)
MAP <- function(.data, ...){
  
  if(!(is.matrix(.data) || is.data.frame(.data) || is.array(.data))){
    ls <- .data
    #if(!is.null(ls$n)) n <- ls$n
    if(!is.null(ls$covmat)) {.data <- ls$covmat
    } else {
      .data <- ls$.data
    }
  }
  if(isSymmetric(as.matrix(.data))){
    #if(is.null(n)) stop("Argument \"n\" is missing with covariance matrix.")
    if(!all(diag(as.matrix(.data) == 1))) {R <- cov2cor(.data)} else {R <- .data}
  } else {
    R <- cor_nest(.data, ...)$covmat
    #n <- nrow(.data)
  }
  
  nv <- ncol(R)
  
  E <- eigen(R, symmetric = TRUE)
  loadings <- E$vectors %*% sqrt(diag(E$values))
  
  fm <- cbind(seq(1, nv), seq(1, nv))
  fm[1, 2] <- (sum(R ^ 2) - nv) / (nv * (nv - 1))
  
  for (m in 1:(nv - 1)){
    A <- loadings[,1:m]
    partcov <- R - (A %*% t(A))
    d <- diag ((1 / sqrt(diag(partcov))))
    pr <- d %*% partcov %*% d
    fm[m+1, 2]  = (sum(pr^2) - nv) / (nv * (nv - 1))
  }
  
  minfm <-  fm[1, 2]
  nfact <-  0
  for (s in 1:nv){
    fm[s, 1] <- s - 1
    if(fm[s, 2] < minfm){
      minfm <- fm[s, 2]
      nfact <-  s-1}
  }   
  
  structure(list(nfactors = nfact, stopping.rule = "Minimum average partial correlation (MAP)"), class = "stoppingrules")
  #return(OUT=list(name = "MAP", nfact = nfact))
}
