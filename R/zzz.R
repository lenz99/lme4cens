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
  cat("\nCoefficients:\nFixed coefs: ", fixef.lmercens(obj))
  cat("\nRandom effect coefs: log(S_betw) = ", log(sigma(obj, which = "between")),
      "\t log(S_within) = ", log(sigma(obj, which = "residual")), "\n")
}

#' Extract variance estimates on standard-deviation scale
#' @param which which type of variance estimate
#' @export
sigma.lmercens <- function(obj, which = c("residual", "between"), ...){
  which <- match.arg(which)

  switch (which,
          between= exp(obj$par[1L]),
          residual=,
          exp(obj$par[2L])
  )
}


#' @export
fixef.lmercens <- function(obj){
  setNames(obj$fixef, nm = colnames(obj$ingredients$X))
}

#' @export
summary.lmercens <- function(obj){
  if (!inherits(obj, "lmercens")) warning("calling summary.lmercens(<fake-lmercens-object>) ...")
  ans <- obj[c("call", "ingredients", "fval")]
  ans$sigmas <- c(between = as.numeric(sigma(obj, which = "between")),
                  residual = as.numeric(sigma(obj, which = "residual")))

  est <- fixef(obj)
  stdError <- sqrt(diag(vcov(obj)))
  tval <- est / stdError
  ans$coefs <- cbind(Estimate = est, `Std. Error` = stdError,
                 `z value` = tval,
                 `Pr(>|z|)` = 2 * pnorm(abs(tval), lower.tail = FALSE))

  class(ans) <- "summary.lmercens"

  ans
}

#' @export
print.summary.lmercens <- function(obj, ...) {
  cat("Linear Mixed Model with Censored Observations\n")

  cat("\nCall:\n", paste(deparse(obj$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("\nCoefficients:\nFixed coefs:\n ")
  coefs <- obj$coefs
  stats::printCoefmat(coefs, ...)
  cat("\nRandom effects:\n")
  cat("S_betw = ", obj$sigmas["between"], "\t S_within = ", obj$sigmas["residual"], "\n")
  cat("Log-likelihood: ", round(-obj$fval, 2), "\n")
}

#' Variance-covariance matrix for fixed effect coefficients.
#' It uses the observed Fisher information matrix at the ML parameter estimate.
#' Therefore, it is only asymptotically usefule.
#' @param obj a fitted `lmercens`-object
#' @return estimated variance-covariance matrix for fixed effect coefficients
#' @export
vcov.lmercens <- function(obj, ...){
  stopifnot( inherits(obj, what = "lmercens") )

  retMat <- NULL

  if (!is.null(obj$hess) && is.matrix(obj$hess) && isSymmetric(obj$hess)){
    retMat <- solve(obj$hess)
    dimnames(retMat) <- list(colnames(obj$ingredients$X), colnames(obj$ingredients$X))
  }

  retMat
}


#' @export
confint.lmercens <- function(obj, parm, level = 0.95, ...){
  stopifnot( inherits(obj, "lmercens"), is.numeric(level), level > 0L, level < 1L )

  cf <- fixef(obj)
  pnames <- names(cf)
  ses <- sqrt(diag(vcov(obj)))

  if (missing(parm)) parm <- pnames else if (is.numeric(parm)) parm <- pnames[parm]

  a <- (1L - level)/2L
  a <- c(a, 1L - a)
  qnt <- qnorm(a)

  ci <- array(NA_real_,
              dim = c(length(parm), 2L), dimnames = list(parm, stats:::format.perc(a, 3)))
  ci[] <- cf[parm] + tcrossprod(ses[parm],qnt)

  ci
}
