# mkuhn, 2017-03-13
# package stuff that has no better place than here

## Fixed by using lme4cens::ghQuadRule (see http://r.789695.n4.nabble.com/no-visible-binding-for-global-variable-for-data-sets-in-a-package-td4696053.html)
# .onLoad <- function(libname = find.package("lme4cens"), pkgname = "lme4cens"){
#
#   # CRAN Note avoidance
#   if(getRversion() >= "2.15.1")
#     utils::globalVariables(names = c("ghQuadRule"))
#
#   invisible()
# }






#' Efficient calculation of \code{log(1-exp(x))}.
#'
#' useful for some likelihood calculations, e.g. exponential distribution with left-censoring.
#' @param x numeric, non-positive values
#' @seealso Martin Maechler
#' @export
log1mexp <- function(x){
  stopifnot( all(x <= 0L), is.numeric(x) )

  log1pBranch <- x < -log(2)
  x[log1pBranch] <- log1p(-exp(x[log1pBranch]))
  x[!log1pBranch] <- log(-expm1(x[!log1pBranch]))

  x
}

#' Calculates log(x+y)
#' Avoid underflow when numbers are small.
#' @param lx 1st argument on log-scale
#' @param ly 2nd argument on log-scale
#' @seealso http://stackoverflow.com/questions/5802592/dealing-with-very-small-numbers-in-r
logxpy <- function(lx,ly) max(lx,ly) + log1p(exp(-abs(lx-ly)))
