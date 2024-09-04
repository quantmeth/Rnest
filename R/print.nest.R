#' Print results of NEST
#'
#' @description Print the number of factors to retain according to confidence levels.
#' @param x An object of class "nest".
#' @param ... Further arguments for other methods, ignored for "nest".
#'
#' @importFrom crayon blue
#' @importFrom utils packageVersion
#' @export
#' 
#' @return No return value, called for side effects.
#'
#' @examples 
#' results <- nest(ex_2factors, N = 100)
#' print(results)
print.nest <- function(x, ...){
  for(i in 1:length(x$alpha)){
    al <- paste0("At ",rownames(x$nfactors)[i]," confidence", sep = "")
    cat(al, ", ", x$stopping.rule, " suggests ", crayon::blue(x$nfactors[i,], .s(x$nfactors[i,], "factor")), ". \n", sep = "")
  }
}

# .s ####
.s <- function(x, w = NULL){
  paste0(w, c("s")[x>1])
}
