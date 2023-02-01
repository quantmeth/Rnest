#' Get rds files from Rnest into the environment
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
#' @param pkg is the package
#'
#' @return A data file
#' @export
#'
#' @examples \dontrun{get.data(ex_2factors)}
get.data <- function(x = list("cormat.l", 
                              "cormat",
                              "ex_2factors",
                              "ex_3factors_doub_unique",
                              "ex_4factors_corr",
                              "eig_pa",
                              "vcritproto"), pkg = "Rnest"){
  for(i in 1:length(x)){
    path <- system.file(file.path("extdata", paste0(x[[i]],".rds")), 
                        package = pkg)
    #assign(x[[i]], readRDS(path), envir = as.environment(1L))
    (function(key, val, pos) assign(key, 
                                    val, 
                                    envir =
                                      as.environment(pos)))(x[[i]], readRDS(path), 1L) 
  }
}
# get.data()
# global env set hack (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(myKey, myVal, 1L) 