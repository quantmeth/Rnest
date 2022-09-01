#' Print results of NEST
#'
#' @description Print the number of factors to retain according to confidence levels.
#' @param x an object of class "nest".
#' @param ... further arguments for other methods, ignored for "nest".
#'
#' @importFrom crayon red blue green
#' @export
#'
#' @examples 
#' results <- nest(ex_2factors, n = 100)
#' print(results)
print.nest <- 
  function(x, ...){
  for(i in 1:length(x$alpha)){
    al <- paste0("At ",rownames(x$nfactors)[i]," confidence", sep = "")
    cat(al, ", NEST suggests ", crayon::blue(x$nfactors[i,], .s(x$nfactors[i,], "factor")), ". \n", sep = "")
  }
  }

# .s ####
.s <- function(x, w = NULL){
  paste0(w, c("s")[x>1])
}
