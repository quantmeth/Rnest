#' Empirical Kaiser Criterion (EKC)
#'
#' @param .data a data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param n the number of cases (subjects, participants, or units) if a covariance matrix is supplied in \code{.data}.
#' @param nv the number of variables if the critical values are required.
#' @param lowest.eig minimal eigenvalues to retain. Default is Kaiser's suggestion of 1.
#' @param ... further argument for \code{cor_nest()}.
#'
#' @return The number of factors to retain or the crititical eigenvalues.
#' @export
#' 
#' @note As Rnest version >= 1.2, a correction to EKC was done which was reported by Marcel van Assen (personnal communication, june 2025), which was found in Rnest and other packages as well. There was a confusion in the sample and critical eigenvalues in equation 2 (Braeken & van Assen, 2017, p. 454).
#' 
#' @aliases ekc
#'
#' @references 
#' Braeken, J., & van Assen, M. A. L. M. (2017). An empirical Kaiser criterion. \emph{Psychological Methods}, \emph{22}(3), 450–466. \doi{10.1037/met0000074}
#' 
#' @examples
#' EKC(ex_4factors_corr, n = 42)
EKC <- function(.data = NULL, n = NULL, nv = NULL, lowest.eig = 1, ...){
  
  if(!is.null(.data)){
    if(!(is.matrix(.data) || is.data.frame(.data) || is.array(.data))){
      ls <- .data
      if(!is.null(ls$n)) n <- ls$n
      if(!is.null(ls$covmat)) {.data <- ls$covmat
      } else {
        .data <- ls$.data
      }
    }
    if(isSymmetric(as.matrix(.data))){
      if(is.null(n)) stop("Argument \"n\" is missing with covariance matrix.")
      if(!all(diag(as.matrix(.data) == 1))) {R <- cov2cor(.data)} else {R <- .data}
    } else {
      R <- cor_nest(.data, ...)$covmat
      n <- nrow(.data)
    }
    
    E <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    nv <- length(E)
  }
  
  if(is.null(n) & is.null(nv)) stop("Please supplied both \"n\" and \"nv\" or a \".data\" argument.")

  crit <- as.numeric()
  for(i in 1:nv){
    #crit[i] <- max(((1 + sqrt(nv / n)) ^ 2) * (nv - sum(crit)) / (nv - i + 1), lowest.eig)
    # which one?
    crit[i] <- max(((1 + sqrt(nv / n)) ^ 2) * (nv - sum(E[0:(i-1)])) / (nv - i + 1), lowest.eig)
    }
  if(!is.null(.data)){
    nfact <- min(which(E < crit))-1
  } else {
    nfact <- NULL
  }
  structure(list(nfactors = nfact, 
                 Eig = crit,
                 stopping.rule = "Empirical Kaiser Criterion (EKC)"), class = "stoppingrules")
  #return(OUT=list(name = "EKC", nfactors = nfact))
}
#structure(c(list(nfactors = nfactors), R, list(stopping.rule = "Next Eigenvalue Sufficiency Test (NEST)")), class = "nest")

#' @export
#' @importFrom crayon blue
print.stoppingrules <- function(x, ...){
  if(!is.null(x$nfactors)){
    cat(x$stopping.rule, " suggests ", crayon::blue(x$nfactors, .s(x$nfactors, "factor")), ". \n", sep = "")
  } else {
    cat("Critical values are : \n")
    cat(round(x$Eig, 3),"\n")
    # print(x$Eig, digit = 3, row.names = FALSE)
  }
}

#' @export
ekc <- EKC