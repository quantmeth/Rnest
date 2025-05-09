% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shem.R
\name{shem}
\alias{shem}
\title{Split-Half Eigenvector Matching (SHEM)}
\usage{
shem(.data, nIts = 30)
}
\arguments{
\item{.data}{a data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.}

\item{nIts}{number of iterations.}
}
\value{
\code{shem} returns a list containing the number of components, \code{nfactors}, whether the additional step in case of zero true latent components was carried, \code{zeroComponents}, the \code{eigenvalues} and the \code{eigenvectors} of the solution.
}
\description{
\code{shem} estimates the number of principal components via Split-Half Eigenvector Matching (SHEM).
}
\examples{
jd <- genr8(n = 404, R = ex_4factors_corr)
shem(jd)
}
\references{
Galdwin, T. E. (2023) Estimating the number of principal components via Split-Half Eigenvector Matching (SHEM). \emph{MethodsX}, \emph{11}, 102286. \doi{10.1016/j.mex.2023.102286}
}
