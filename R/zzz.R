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





# mkuhn, 2017-11-23:
# add some S3-convenience functions for q&d lmercens

#' @export
print.lmercens <- function(obj){
  cat("Linear Mixed Model with Censored Observations\n")

  cat("\nCoefficients:\nFixed coefs: ", fixef(obj))
  cat("\nRandom effect coefs: ", obj$par[1L:2L])
}

#' @export
fixef.lmercens <- function(obj){
  obj$fixef
}
