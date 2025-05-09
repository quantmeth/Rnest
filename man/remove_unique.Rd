% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_unique.R
\name{remove_unique}
\alias{remove_unique}
\title{Remove unique variables}
\usage{
remove_unique(.data, ..., alpha = 0.05)
}
\arguments{
\item{.data}{a data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.}

\item{...}{further arguments for \code{unique_variable()} and \code{cor_nest()}.}

\item{alpha}{type I error rate.}
}
\value{
A list containing the unique variables and a data frame containing their probabilities and the \code{.data} with the unique variable removed.
}
\description{
Remove unique variables
}
\examples{
remove_unique(ex_3factors_doub_unique, n = 420)
}
