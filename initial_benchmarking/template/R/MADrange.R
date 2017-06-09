#' MADrange
#' 
#' Computes upper and lower limit of an array for a given number of absolute deviations from the median.
#' With option \code{log} this calculation is done on the log10 transformed array.
#'
#' @param arr     \code{num} array for which to compute limits
#' @param nmads   \code{int} scalar for the number of MADs to use
#' @param log     if \code{TRUE} limits are computed based on log10 of the arr
#'
#' @return \code{num} vector of length 2 with lower and upper limit
#' @export
#'
MADrange <- function( arr, nmads = 5, log = FALSE ) 
{
  if(log) arr <- log10(arr)
  cur.med <- median(arr)
  cur.mad <- mad(arr, center = cur.med)
  upper.limit <- cur.med + nmads * cur.mad
  lower.limit <- cur.med - nmads * cur.mad
  out <- c(lower.limit, upper.limit)
  if(log) out <- 10^out
  return(out)
}