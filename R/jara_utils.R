#' rqpois()
#'
#' Random generator for quasi-poisson variables    
#' @param n number of random observations
#' @param lambda vector of (non-negative) expected means
#' @param phi over-dispersion (variance) parameter 
#' @return random quasi-poisson observations
#' @export
rqpois = function(n, lambda, phi) {
  mu = lambda
  k = mu/phi/(1-1/phi)
  r = stats::rnbinom(n, mu = mu, size = k)
  return(r)
}



