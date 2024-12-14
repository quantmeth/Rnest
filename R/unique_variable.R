#' Probability of unique variables
#'
#' @param .data a data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param n the number of cases (subjects, participants, or units) if a covariance matrix is supplied in \code{.data}.
#' @param ... further arguments for \code{cor_nest()}.
#' 
#' @return A data frame containing the F-values and probabilities of the variable to be an unique variable. 
#' @export
#'
#' @author 
#' P.-O. Caron (R)
#' André Achim (Matlab)
#' 
#' @examples
#' exData <- genr8(n = 420, R = ex_3factors_doub_unique)
#' unique_variable(exData)
unique_variable <- function(.data, n = NULL, ...){

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
    if(!all(diag(as.matrix(.data) == 1))) {.R <- cov2cor(.data)} else {.R <- .data}
  } else {
    .R <- cor_nest(.data, ...)$covmat
    n <- nrow(.data)
  }
  
  p <- ncol(.R)
  db <- n - p - 1
  
  B <- 1 / diag(solve(.R))
  R2 <- 1 - B
  FF <- R2 * db / (B * p)
  pval <- pf(FF, p, db, lower.tail = FALSE)
  
  out <- data.frame(Fvalue = FF,
                    df1 = p,
                    df2 = db,
                    p = pval)
  
  return(structure(list(Results = out,
                        .data = .data,
                        n = n),
                   class = "puniquevar"))
}

#' @export
print.puniquevar <- function(x, alpha = .05, ...){
  out <- round(x$Results, 3)
  out$p[out$p == 0] <- paste("     < 0.001  ")
  out$p[out$p > alpha] <- paste0("", out$p[out$p > .05]," *")
  colnames(out) <- c("F-value","df1",'df2',"Pr(F>)")
  print(out)
  if(!is.null(x$Variables)){
  cat("\n")
  cat(.s(length(x$Variables), "Variable"), crayon::blue(x$Variables),.ve(length(x$Variables)),"been removed.\n")
  }
}
