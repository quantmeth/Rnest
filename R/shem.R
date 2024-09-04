#' Split-Half Eigenvector Matching (SHEM)
#'
#' @description \code{shem} estimates the number of principal components via Split-Half Eigenvector Matching (SHEM).
#'
#' @param data A data frame, a numeric matrix, covariance matrix or correlation matrix from which to determine the number of factors.
#' @param nIts Number of iterations.
#'
#' @return \code{shem} returns a list containing the number of components, \code{nfactors}, whether the additional step in case of zero true latent components was carried, \code{zeroComponents}, the \code{eigenvalues} and the \code{eigenvectors} of the solution. 
#' 
#' @import MASS
#' 
#' @export
#'
#' @references 
#' Galdwin, T. E. (2023) Estimating the number of principal components via Split-Half Eigenvector Matching (SHEM). \emph{MethodsX}, \emph{11], 102286. \doi{10.1016/j.mex.2023.102286}
#'
#' @examples
#' jd <- genr8(n = 404, R = ex_4factors_corr)
#' shem(jd)
shem <- function(data, nIts = 30){
  get_n_components(X = data, nIts)
}

# Similarity of eigenvectors over split-halves
run_PCA <- function(X) {
  X <- scale(X, center = TRUE, scale = FALSE)
  C <- cov(X)
  eigenvalues <- eigen(C)$values
  eigenvectors <- eigen(C)$vectors
  return(list(eigenvalues = eigenvalues, eigenvectors = eigenvectors))
}

get_eigenvector_sims <- function(X1, X2) {
  pca1 <- run_PCA(X1)
  pca2 <- run_PCA(X2)
  Sims <- abs(t(pca1$eigenvectors) %*% pca2$eigenvectors)
  sims <- numeric()
  for (iL in 1:ncol(Sims)) {
    most_sim_ind <- which.max(Sims[, iL])
    s <- Sims[most_sim_ind, iL]
    sims <- c(sims, s)
    Sims[most_sim_ind, iL] <- 0
  }
  return(sims)
}

sim_per_component <- function(X, nIts = 20) {
  nObs <- nrow(X)
  nVar <- ncol(X)
  Sims <- matrix(0, nIts, nVar)
  for (iIt in 1:nIts) {
    indices <- sample(1:nObs, nObs, replace = FALSE)
    half <- floor(length(indices) / 2)
    X1 <- X[indices[1:half], ]
    X2 <- X[indices[(half + 1):nObs], ]
    sims <- get_eigenvector_sims(X1, X2)
    Sims[iIt, ] <- sims
  }
  m <- colMeans(Sims)
  return(m)
}

sims_split <- function(sims) {
  if (length(sims) < 2) {
    return(0)
  }
  k_v <- numeric()
  scores <- numeric()
  scores_adjusted <- numeric()
  for (k in 1:(length(sims) - 1)) {
    lhs <- sims[1:k]
    rhs <- sims[(k + 1):length(sims)]
    ws_left <- var(lhs)
    ws_right <- var(rhs)
    m_lhs <- mean(lhs)
    m_rhs <- mean(rhs)
    betw_v <- c(rep(m_lhs, length(lhs)), rep(m_rhs, length(rhs)))
    bs <- var(betw_v)
    score <- bs / (ws_left + ws_right)
    scores <- c(scores, score)
    score_adjusted <- score * (1 - (k - 1) / length(sims))
    scores_adjusted <- c(scores_adjusted, score_adjusted)
    k_v <- c(k_v, k)
  }
  nComp <- k_v[which.max(scores)]
  zeroComp <- FALSE
  if (nComp > 0 && scores_adjusted[nComp] < 1) {
    zeroComp <- TRUE
  }
  return(list(nComponents = nComp, zeroComponents = zeroComp))
}

get_n_components <- function(X, nIts = 30) {
  pca <- run_PCA(X)
  sims <- sim_per_component(X, nIts = nIts)
  res <- sims_split(sims)
  return(list(nfactors = res$nComponents, zeroComponents = res$zeroComponents, eigenvalues = pca$eigenvalues, eigenvectors = pca$eigenvectors))
}