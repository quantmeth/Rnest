#' Summary results of NEST
#'
#' @description summary method for class "nest".
#' @param object an object of class "nest".
#' @param ... further arguments for other methods, ignored for "nest".
#'
#' @return No returned value, called for side effects.
#'
#' @importFrom crayon blue
#' @export
#'
#' @examples 
#' results <- nest(ex_2factors, n = 100)
#' summary(results)
summary.nest <- function(object, ...){
  conv <- ifelse(object$convergence, "ended normally", paste0("failed at testing k = ", length(object$Eig)-1))
  cat("\n")
  cat("nest",  paste0(unlist(packageVersion("Rnest")), collapse = "."),conv," \n \n")
  cat("   Estimator                      ", toupper(object$method),"\n")
  cat("   Missing data treatment         ", toupper(object$na.action),"\n")
  cat("   Number of model parameters     ", (ncol(object$cor))*(ncol(object$cor)-1)/2, "\n")
  cat("   Resampling                     ", object$nrep,"\n")
  cat("   Sample size                    ", object$n, "\n")
  cat("   Stopped at                     ", length(object$prob),"\n \n \n")
  
  cat("Test that k factors are sufficient \n")
  cat("\n")
  df <- data.frame(`k factor` = paste0("k = ", 0:(length(object$Eig)-1)),
                   `NextEig` = round(object$values[1:length(object$Eig)],3),
                   CritEig = round(sapply(1:length(object$Eig), function(i) object$Eig[[i]][1,i]),3),
                   `Prob` = ifelse(object$prob < .001," < .001",paste("  ", numformat(object$prob))),
                   check.names = FALSE)

  
print(df, row.names = FALSE, justify = "centre")

  cat("\n \n")
  print(object)
  cat("\n")
  cat("Try", crayon::blue("plot(nest())"), "to see a graphical representation of the results. \n \n")
  
}

numformat <- function(val) {sub("^(-?)0.", "\\1.", sprintf("%.3f", val))}
