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
  cat("\n")
  cat("nest",  paste0(unlist(packageVersion("Rnest")), collapse = "."),"ended normally \n \n")
  cat("   Estimator                      ", toupper(object$method),"\n")
  cat("   Number of model parameters     ", (ncol(object$cor))*(ncol(object$cor)-1)/2, "\n")
  cat("   Resampling                     ", object$nrep,"\n")
  cat("   Sample size                    ", object$n, "\n")
  cat("   Stopped at                     ", length(object$prob),"\n \n \n")
  
  cat("Probabilities of factors \n")
  cat("  Factor     Eigenvalue     Prob \n")
  for(i in 1:length(object$Eig)){
    cat("   ",paste0("F", i),"       ", sprintf("%.3f", object$values[i]), "     ", 
        
        
        ifelse(object$prob[i] < .001,"< .001",paste(" ", numformat(object$prob[i]))),"\n")
  }
  cat("\n \n")
  print(object)
  cat("Try", crayon::blue("plot(nest())"), "to see a graphical representation of the results. \n \n")
 
}
numformat <- function(val) {sub("^(-?)0.", "\\1.", sprintf("%.3f", val))}
