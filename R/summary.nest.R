summary.nest <- function(x, ...){
  ld <- loadings(x)
  ld <- round(ld, 3)
  ld[ld < .1] <- ""
  cat("\n")
  print(x)
  cat("\n")
  cat("Loadings: \n")  
  cat("\n")
  noquote(ld)
}
