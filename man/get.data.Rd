% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.data.R
\name{get.data}
\alias{get.data}
\title{Get files from packages into the environment}
\usage{
get.data(
  x = list("cormat.l", "cormat", "ex_2factors", "ex_3factors_doub_unique",
    "ex_4factors_corr"),
  pkg = "Rnest"
)
}
\arguments{
\item{x}{can be \itemize{
\item{"cormat.l"}
\item{"cormat"}
\item{"ex_2factors"}
\item{"ex_3factors_doub_unique"}
\item{"ex_4factors_corr"}
}}

\item{pkg}{is the package}
}
\value{
A data file
}
\description{
Get RDS files from Rnest package (compatibility issues)
}
\examples{
\dontrun{get.data(ex_2factors)}
}
