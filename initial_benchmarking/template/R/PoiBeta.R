#' PoiBetaMMFit
#'
#' Estimate a,b, and k of Poisson-Beta distribution using mean,
#' variance and skewness of x.
#'
#' @param x \code{integer} vector of sample values
#'
#' @return list of estimates for a, b, k
#'
#' @export
#'
#' @importFrom moments moment
#'
PoiBetaMMFit <- function(x) {
  #Calculate the central moments
  mu <- mean(x)
  mu2 <- var(x)
  mu3 <- moments::moment(x, order = 3, central = T)
  #Estimate the parameters
  k <- 2 * mu - (mu3 * (mu^2 + mu2)) / (- 2 * mu2^2 + mu * mu3)
  a <- (2 * mu * (mu^2 * mu2 + mu3 * mu - mu2^2)) /
       (-mu3 * mu^2 + 4 * mu * mu2^2 + mu3 * mu2)
  b <- -(2 * mu2 * (mu3 + 2 * mu * mu2) * (mu^2 * mu2 + mu3 * mu - mu2^2)) /
       ((- 2 * mu2^2 + mu * mu3) * (- mu3 * mu^2 + 4 * mu * mu2^2 + mu3 * mu2))
  return(list("k" = k, "a" = a, "b" = b))
}

#' PoiBetaLogLikelihood
#'
#' Estimate Log-Likelihood for a Poisson-Beta distribution parameter estimation
#' by the 1st 3 moments.
#'
#'
#' @export
#'
#' @importFrom hypergeo genhypergeo
#' @importFrom orthopolynom lpochhammer
#'
PoiBetaLogLikelihood <- function(x, k = NaN, a = NaN, b = NaN,
                                 nMC = 1e3, mcLim = 1e2) {
  if (is.nan(k) | is.nan(a) | is.nan(b)) {
    res <- PoiBetaMMFit(x)
    if (is.nan(k)) { k <- res$k }
    if (is.nan(a)) { a <- res$a }
    if (is.nan(b)) { b <- res$b }
  }
  n <- length(x)
  logLike <- rep(n, 0)
  for (i in 1:n) {
    if (x[i] < mcLim & x[i] > 0) {
      logLike[i] <- x[i]*log(k) -
                    k +
                    orthopolynom::lpochhammer(a, x[i]) -
                    lgamma(x[i]) -
                    orthopolynom::lpochhammer(a + b, x[i]) +
                    log(hypergeo::genhypergeo(a, a + b + x[i], k))
    } else {
      q <- rbeta(nMC, a, b)
      p <- dpois(x[i], k*q)
      logLike[i] <- log(mean(p))
    }
  }
  return(sum(logLike))
}
