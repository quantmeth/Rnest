#' Full Information Maximum Likelihood (FIML) correlation or covariance matrix
#'
#' @param data a data frame of rdata matrix.
#' @param tol tolerance.
#' @param maxiter maximum number of iterations.
#' @param pvalue an argument to indicate if \eqn{p}-values are required.
#'
#' @usage covFIML(data, tol = 1e-6, maxiter = 1000, pvalue = FALSE)
#' @usage corFIML(data, tol = 1e-6, maxiter = 1000, pvalue = FALSE)
#' 
#' @aliases covFIML corFIML
#' 
#' @note A not so efficient function. See \code{?cor_nest} instead.
#'
#' @return A list containing the means, th correlation or covariance matrix, and optionnaly the degree of freedom and the p-values.
#' @export
#' @importFrom mvtnorm dmvnorm
#'
#' @examples
#' covFIML(airquality)
covFIML <- function(data, tol = 1e-6, maxiter = 1000, pvalue = FALSE){

  # mod <- paste0(paste0(colnames(dta),"~~", colnames(dta)), collapse ="\n")
  # fit <- sem(mod, dta, missing = "fiml")
  # S <- lavInspect(fit, "cov.ov")
  
  mu <- colMeans(data, na.rm = TRUE)
  sigma <- cov(data, use = "pairwise.complete.obs")
  dl <- vp <- NULL
  # EM algorithm
  n <- nrow(data)
  p <- ncol(data)
  
  data <- as.matrix(data)
  
  iter <- 0
  prev_ll <- -Inf
  repeat {
    iter <- iter + 1
    mu_new <- rep(0, p)
    sigma_new <- matrix(0, p, p)
    # E-step
    for (i in 1:n) {
      obs <- !is.na(data[i, ])
      mis <- is.na(data[i, ])
      y_obs <- data[i, obs]
      if (sum(mis) > 0) {
        mu_mis <- mu[mis] + sigma[mis, obs] %*% solve(sigma[obs, obs]) %*% (y_obs - mu[obs])
        sigma_mis <- sigma[mis, mis] - sigma[mis, obs] %*% solve(sigma[obs, obs]) %*% sigma[obs, mis]
      } else {
        mu_mis <- numeric(0)
        sigma_mis <- matrix(0, 0, 0)
      }
      y_complete <- data[i, ]
      y_complete[mis] <- mu_mis
      mu_new <- mu_new + y_complete
      sigma_new <- sigma_new + outer(y_complete - mu, y_complete - mu) #+ sigma_mis
      if (length(sigma_mis) > 0) {
        sigma_new[mis, mis] <- sigma_new[mis, mis] + sigma_mis
      }
    }
    mu_new <- mu_new / n
    sigma_new <- sigma_new / n
    # M-step
    ll <- log_likelihood(mu_new, sigma_new, data, n)
    if (abs(ll - prev_ll) < tol || iter >= maxiter) break
    mu <- mu_new
    sigma <- sigma_new
    prev_ll <- ll
  }
  
  if(pvalue){
    r <- cov2cor(sigma)
    dl <- t(!is.na(data)) %*% (!is.na(data)) -2
    vt <-  r  * sqrt(n - 2) / sqrt(1 -  r  ^ 2)
    vp <- (1 - pt(abs(vt), df = dl)) * 2
  }
  
  return(list(mu = mu, sigma = sigma, df = dl, pvalue = vp))
}

log_likelihood <- function(mu, sigma, data, n) {
  ll <- 0
  for (i in 1:n) {
    obs <- !is.na(data[i, ])
    y_obs <- data[i, obs]
    mu_obs <- mu[obs]
    sigma_obs <- sigma[obs, obs]
    ll <- ll + mvtnorm::dmvnorm(y_obs, mean = mu_obs, sigma = sigma_obs, log = TRUE)
  }
  return(ll)
}

#' @export
corFIML <- function(data, tol = 1e-6, maxiter = 1000, pvalue = FALSE){
  S <- covFIML(data, tol, maxiter, pvalue)
  S$sigma <- cov2cor(S$sigma)
  S
}