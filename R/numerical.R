# numerical stuff on the R-side


#' Numeric integration with Gauss-Hermite quadrature.
#'
#' Weight function is \code{exp(-x^2)} and it is an improper integral from -Inf to Inf.
#' E.g. for standard normal density integration, use f=1 constant and get sqrt(pi)-times the integral value.
#'
#' The given function can either be a constant or a function to be integrated  after factoring out the weight function w = exp(-x^2).
#' The function must be vectorized in its variable of integration.
#'
#' @export
#' @param f constant or function to be integrated
#' @param ord order of quadrature rule
#' @param ... additional arguments passed down to integrand function f
#' @return numeric, approximation of integral
int_gh <- compiler::cmpfun(function(f, ord=8L, ...){
  ind <- which(ghQuadRule$order == ord)
  stopifnot( ord >=2L, length(ind) == ord )

  xi <- ghQuadRule$abscissa[ind]
  wi <- ghQuadRule$weight[ind]

  if (is.numeric(f)) f[1L] * sum(wi) else
    sum(wi * f(xi, ...))
})



