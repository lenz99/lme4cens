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




#' Efficient calculation of `log(1-exp(x))` for negative `x`.
#'
#' It combines two robust ways for the calculation. It is vectorized.
#' Useful for some likelihood calculations, e.g. exponential distribution with left-censoring.
#' @param x numeric, non-positive values
#' @seealso Martin Maechler
#' @export
log1mexp <- function(x){
  stopifnot( is.numeric(x), all(x <= 0L) )

  log1pBranch <- x < -log(2)
  x[log1pBranch] <- log1p(-exp(x[log1pBranch]))
  x[!log1pBranch] <- log(-expm1(x[!log1pBranch]))

  x
}

#' Calculates `log(x+y)` where the (small positive) arguments are given on log-scale.
#'
#' This implementation avoids underflow when numbers are small. The trick is to do the summation on log-space.
#' The implementation is vectorized.
#' @param lx 1st argument on log-scale
#' @param ly 2nd argument on log-scale
#' @seealso http://stackoverflow.com/questions/5802592/dealing-with-very-small-numbers-in-r
logxpy <- function(lx, ly) { stopifnot( length(lx) == length(ly) ); lxy_max <- lx; lxy_max[lx < ly] <- ly[lx < ly]; lxy_max + log1p(exp(-abs(lx-ly))) }

#' Calculates `log(x-y)` where `x > y` holds
#' @param lx 1st argument on log-scale
#' @param ly 2nd argument on log-scale
logxmy <- function(lx, ly) { stopifnot( length(lx) == length(ly), all(lx > ly) ); lx + log1mexp(ly-lx) }



