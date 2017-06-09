#' GetModes
#' 
#' Estimate modes by first estimating density using \code{density}, then computing maxima (including border Maxima).
#' \code{density} uses a \emph{Gauss} kernel and \emph{Silverman} rule of thumb bandwith estimation.
#'
#' @param x  \code{num} array to estimate modes
#' @param bw either \code{chr} or \code{num} describing bandwith (see \code{density} function)
#' @param w  \code{num} vector of weights (same length as \code{x}) or \code{NULL} (equal weights)
#'
#' @return \code{num} matrix with x and y positions of modes in density estimation or \code{NULL}
#' @export
#'
GetModes <- function( x, bw = "nrd0", w = NULL ) 
{
  # can't calc dens with less than 2
  if( length(x) < 2 ) return(NULL)
  
  # calc dens
  dens <- density(x, bw, kernel = "gaussian", n = 512, weights = w)
  
  # get local max without borders
  ms <- which(diff(sign(diff(dens$y))) == -2) + 1
  
  # border maxima
  if( dens$y[1] > dens$y[2] ) ms <- c(1, ms)
  if( dens$y[512] > dens$y[511] ) ms <- c(ms, 512)
  
  # monotonous case
  if( identical(ms, numeric()) ) return(NULL) 
  
  # return coords
  out <- cbind(dens$x[ms], dens$y[ms])
  colnames(out) <- c("x", "y")
  return(out)
}
