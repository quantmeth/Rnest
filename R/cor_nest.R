#' Compute covariance or correlation matrix with treatments for clusters and missing values
#'
#' @param .data a data frame, a numeric matrix.
#' @param ... further arguments.
#' @param cluster a variable name defining the clusters in a two-level dataset in the data frame.
#' @param missing treatment to deal with missing values. Options are \code{"listwise"} or \code{"pairwise"}. Default if \code{"fiml"}.
#' @param pvalue an argument to indicate if \eqn{p}-values are required.
#'
#' @usage cor_nest(.data, ..., cluster = NULL, missing = "fiml", pvalue = FALSE)
#' @usage cov_nest(.data, ..., cluster = NULL, missing = "fiml", pvalue = FALSE)
#' 
#' @aliases cor_nest
#' 
#' @return A list of class "covnest"
#' @importFrom lavaan sem lavInspect 
#' @importFrom utils combn
#' @export
#'
#' @examples
#' cov_nest(airquality)
cov_nest <- function(.data, ..., cluster = NULL, missing = "fiml", pvalue = FALSE){
  
  .data <- as.data.frame(.data)
  nom <- colnames(.data)
  if(is.null(nom)) colnames(.data) <- nom <- paste0("v",1:ncol(.data))
  nom <- nom[!(nom %in% cluster)]
  #colnames(.data) <- nv <- paste0("v",1:ncol(.data))
  nv <- combn(nom,2)
  
  mod <- paste0(paste0(nv[1,],"~~",nv[2,]), collapse = "\n")
  fit <- lavaan::sem(mod, .data, missing = missing, cluster = cluster, warn = FALSE, ...)
  S <- lavaan::lavInspect(fit, c("cov.ov"))[]
  n <- lavaan::lavInspect(fit, c("nobs"))
  
  colnames(S) <- row.names(S) <- nom
  
  out <- list(covmat = S,
       n = n)
  
  if(pvalue){
    r <- cov2cor(S)
    df <- t(!is.na(.data)) %*% (!is.na(.data)) -2
    vt <-  r  * sqrt(df) / sqrt(1 -  r  ^ 2)
    vp <- (1 - pt(abs(vt), df = df)) * 2
    out <- c(out,
             list(df = df,
                  tvalue = vt,
                  pvalue = vp))
  }
  class(out) <- "covnest"
  out
}

#' @export
cor_nest <- function(.data, ..., cluster = NULL, missing = "fiml", pvalue = FALSE){
  out <- cov_nest(.data, ..., cluster = cluster, missing = missing, pvalue = pvalue)
  out$covmat <- cov2cor(out$covmat)
  out
}

#' @export
print.covnest <- function(x, digit = 3, alpha = .05, ...){
  S <- x$covmat
  S[] <- sprintf(paste0("%.", digit, "f"),S)
  #S <- round(S, digit = digit)
  if(all(diag(S) == "1.000")){
    diag(S) <- "  -"
  }
  if(!is.null(x$pvalue)){
    S[lower.tri(S)] <- paste0(S[lower.tri(S)], ifelse(x$pvalue[lower.tri(S)] < alpha,"*"," "))
  }
  S[upper.tri(S)] <- ""
  print(S[-1, -ncol(S)], quote = FALSE)
}
