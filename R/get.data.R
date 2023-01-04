#' Get RDS files from Rnest package
#'
#' @description Get RDS files from Rnest package (compatibility issues)
#'
#' @param x can be \itemize{
#' \item{"cormat.l"}
#' \item{"cormat"}
#' \item{"ex_2factors"}
#' \item{"ex_3factors_doub_unique"}
#' \item{"ex_4factors_corr"}
#' }
#'
#' @return A data file
#' @export
#'
#' @examples \dontrun{get.data(ex_2factors)}
get.data <- function(x = list("cormat.l", 
                              "cormat",
                              "ex_2factors",
                              "ex_3factors_doub_unique",
                              "ex_4factors_corr")){
  sapply(x, function(x){
    path <- system.file(file.path("extdata", paste0(x,".rds")), 
                        package = "Rnest")
    assign(x, readRDS(path))
  })
}
