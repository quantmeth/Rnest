% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_nest.R
\name{cov_nest}
\alias{cov_nest}
\alias{cor_nest}
\title{Compute covariance or correlation matrix with treatments for clusters and missing values}
\usage{
cor_nest(.data, ..., cluster = NULL, missing = "fiml", pvalue = FALSE)

cov_nest(.data, ..., cluster = NULL, missing = "fiml", pvalue = FALSE)
}
\arguments{
\item{.data}{a data frame, a numeric matrix.}

\item{...}{further arguments.}

\item{cluster}{a variable name defining the clusters in a two-level dataset in the data frame.}

\item{missing}{treatment to deal with missing values. Options are \code{"listwise"} or \code{"pairwise"}. Default if \code{"fiml"}.}

\item{pvalue}{an argument to indicate if \eqn{p}-values are required.}
}
\value{
A list of class "covnest"
}
\description{
Compute covariance or correlation matrix with treatments for clusters and missing values
}
\examples{
cov_nest(airquality)
}
