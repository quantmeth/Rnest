#' @importFrom cli style_hyperlink
#' @importFrom crayon blue
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is ", pkgname," ", packageVersion(pkgname),"\n", 
                        pkgname," is FREE software!\n",
                        "Please report any bugs to ", crayon::blue(cli::style_hyperlink(
                          text = "https://github.com/quantmeth/Rnest",
                          url = "https://github.com/quantmeth/Rnest"
                        )))
}

