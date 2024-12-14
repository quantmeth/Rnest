#' Remove unique variables
#'
#' @param .data a data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param ... further arguments for \code{unique_variable()} and \code{cor_nest()}.
#' @param alpha type I error rate.
#'
#' @return A list containing the unique variables and a data frame containing their probabilities and the \code{.data} with the unique variable removed.
#' @export
#'
#' @examples
#' remove_unique(ex_3factors_doub_unique, n = 420)
remove_unique <- function(.data, ..., alpha = .05){
  
  if(!(is.matrix(.data) || is.data.frame(.data) || is.array(.data))){
    ls <- c(.data, list(...))
    #if(!is.null(ls$n)) n <- ls$n
    if(!is.null(ls$covmat)) {D <- ls$covmat
    } else {
      D <- ls$.data
    }
  } else {
    D <- .data
    ls <- list()
    #ls$.data <- D <- .data
  }
  
  .col <- apply(D, 2, is.numeric)
  ls$.data <- D[,.col]
  
  out <- unique_variable(ls, ...)
  pval <- out$Results["p"]
  
  .rcol <- names(.col)[which(pval > alpha)]
  
  if(isSymmetric(as.matrix(D))){
    D <- D[!(rownames(D) %in% .rcol), 
           !(colnames(D) %in% .rcol)]
  } else {
    D[,.rcol] <- NULL
  }
  out$Variables <- .rcol
  out$.data <- D
  out
}

