# mkuhn, 2017-03-13
# package stuff that has no better place than here

## warning of no visible binding
## Fixed by using lme4cens::ghQuadRule (see http://r.789695.n4.nabble.com/no-visible-binding-for-global-variable-for-data-sets-in-a-package-td4696053.html)
# .onLoad <- function(libname = find.package("lme4cens"), pkgname = "lme4cens"){
#
#   # CRAN Note avoidance
#   if(getRversion() >= "2.15.1")
#     utils::globalVariables(names = c("ghQuadRule"))
#
#   invisible()
# }


#' @importFrom survival Surv
#' @export
survival::Surv




# mkuhn, 2017-11-23:
# add some S3-convenience functions for q&d lmercens

#' @export
print.lmercens <- function(obj){
  cat("Linear Mixed Model with Censored Observations\n")

  cat("\nCall:\n", paste(deparse(obj$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("\nCoefficients:\nFixed coefs: ", fixef(obj))
  cat("\nRandom effect coefs: log(S_betw) = ", obj$par[1L], "\t log(S_within) = ", obj$par[2L])
}

#' Extract residual standard deviation 'sigma'
#' @export
sigma.lmercens <- function(obj, ...){
  exp(obj$par[2L])
}

#' Extract the between subject standard deviation 'sigma_betw'
#' @export
sigmaBetw <- function(obj){
  exp(obj$par[1L])
}

#' @export
fixef.lmercens <- function(obj){
  obj$fixef
}

#' @export
summary.lmercens <- function(obj){
  print.lmercens(obj)
}
