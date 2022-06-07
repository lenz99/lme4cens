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
print.lmercens <- function(x, ...){
  cat("Linear Mixed Model with Censored Observations\n")

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("\nCoefficients:\nFixed coefs: ", fixef.lmercens(x))
  cat("\nRandom effect coefs: log(S_betw) = ", log(sigma(x, which = "between")),
      "\t log(S_within) = ", log(sigma(x, which = "residual")), "\n")
}

#' Extract variance estimates on standard-deviation scale
#' @param which which type of variance estimate
#' @export
sigma.lmercens <- function(object, which = c("residual", "between"), ...){
  which <- match.arg(which)

  switch(which,
          between  = exp(object$par[[1L]]),
          residual = ,
          exp(object$par[[2L]])
  )
}


#' Fixed effect coefficient estimates.
#' @export
fixef.lmercens <- function(object, ...){
  stats::setNames(object$fixef, nm = colnames(object$ingredients$X))
}

#' Random effect predictions for the random intercept.
#' @export
ranef.lmercens <- function(object, ...){
  stopifnot( inherits(object, what = "lmercens") )

  # use only observed cases
  respMatrix <- as.matrix(object$ingredients$fr[[1L]])
  stopifnot( "status" %in% colnames(respMatrix) )
  obsInd <- respMatrix[,"status", drop=TRUE] == 1L  ## eventually could use also interval-censored (=> take mean of interval)
  fixef_resids <- respMatrix[obsInd, 1L, drop=FALSE] - object$ingredients$X[obsInd,, drop=FALSE] %*% fixef.lmercens(object)

  # cf. Demidenko, section 3.7
  Zt <- object$ingredients$reTrms$Zt
  D <- diag(exp(2 * (object$par[[1L]] - object$par[[2L]])), nrow = NROW(Zt))

  # with random intercept only, it is enough to give out a vector (not a named list)
  setNames(as.numeric(D %*% Matrix::solve(a = diag(1, nrow = NROW(Zt)) + Matrix::tcrossprod(Zt) %*% D,
                                          b = Zt[, obsInd] %*% (fixef_resids))),
           nm = levels(object$ingredients$reTrms$flist[[1]]))
}

#' @export
logLik.lmercens <- function(object, ...){
  myNobs <- NROW(object$ingredients$X)
  structure(-object$fval,
            nobs = myNobs,
            nall = myNobs,
            df = length(object$par),
            class = 'logLik')
}

#' Predict method for `lmercens` objects.
#'
#' @param object `lmercens` model object
#' @param newdata dataframe with covariate values for prediction. Default is `NULL` which is to fall back to the training data.
#' @param re.form right-side formula for random effects to condition on in case of prediction on training data. If `NULL`, include all random effects; if `~0` or `NA`, include no random effects.
#' @export
predict.lmercens <- function(object, newdata = NULL, re.form = NULL, ...){
  stopifnot( inherits(object, what = "lmercens") )

  ret <- NULL

  if ( is.null(newdata) ){ # use training data itself for predictions

    ret <- object$ingredients$X %*% fixef.lmercens(object)

    if ( is.null(re.form) || re.form == ~1 ){ ##!isTRUE(is.na(re.form)) && !re.form == ~0 ){
      # expect simple random effects model (only a single random intercept)
      stopifnot( length(object$ingredients$reTrms$cnms) == 1L && object$ingredients$reTrms$cnms[[1L]] == '(Intercept)' )

      ret <- ret + Matrix::crossprod(object$ingredients$reTrms$Zt, y = ranef(object))
    }

  } else { # new data used without random effect BLUPs
    X <- model.matrix(lme4:::getFixedFormula(object$call$formula[-2]), data = newdata[,-1])
    ret <- X %*% fixef.lmercens(object)
  }

  as.vector(ret)
}


#' @export
fitted.lmercens <- function(object, ...){
  predict.lmercens(object, ...)
}

#' @export
residuals.lmercens <- function(object, ...){
  stopifnot( inherits(object, what = "lmercens") )

  y_pred <- predict(object)
  y_obs  <- as.matrix(object$ingredients$fr[[1L]])
  y_obs_status <- y_obs[, "status"]

  # NA as default (for all censored outcomes)
  res <- rep(NA_real_, NROW(y_obs))
  # res as obs - pred for observed outcomes (i.e. no censorings)
  res[y_obs_status == 1L] <- y_obs[y_obs_status == 1L, 1L] - y_pred[y_obs_status == 1L]

  res
}

#' @export
summary.lmercens <- function(object, ...){
  if (!inherits(object, "lmercens")) warning("calling summary.lmercens with fake lmercens-object!", call. = FALSE)
  ans <- object[c("call", "ingredients")]
  ans$sigmas <- c(between = as.numeric(sigma(object, which = "between")),
                  residual = as.numeric(sigma(object, which = "residual")))
  ans$logLik <- logLik(object, ...)

  est <- fixef(object)
  stdError <- sqrt(diag(vcov(object)))
  tval <- est / stdError
  ans$coefs <- cbind(Estimate = est,
                     `Std. Error` = stdError,
                     `z value` = tval,
                     `Pr(>|z|)` = 2L * pnorm(abs(tval), lower.tail = FALSE))

  class(ans) <- "summary.lmercens"

  ans
}


#' @export
print.summary.lmercens <- function(x, ...) {
  cat("Linear Mixed Model with Censored Observations\n")

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("\nCoefficients:\nFixed coefs:\n ")
  coefs <- x$coefs
  stats::printCoefmat(coefs, ...)
  cat("\nRandom effects:\n")
  cat("S_betw = ", x$sigmas["between"], "\t S_within = ", x$sigmas["residual"], "\n")
  cat("Log-likelihood: ", round(as.numeric(x$logLik), 2L), "\n")
}

#' Variance-covariance matrix for fixed effect coefficients.
#' It uses the observed Fisher information matrix at the ML parameter estimate.
#' Therefore, it is only useful asymptotically.
#' @param object a fitted `lmercens`-object
#' @return estimated variance-covariance matrix for fixed effect coefficients
#' @export
vcov.lmercens <- function(object, ...){
  stopifnot( inherits(object, what = "lmercens") )

  retMat <- NULL

  if (!is.null(object$hess) && is.matrix(object$hess) && isSymmetric(object$hess)){
    retMat <- solve(object$hess)
    dimnames(retMat) <- list(colnames(object$ingredients$X), colnames(object$ingredients$X))
  }

  retMat
}


#' @export
confint.lmercens <- function(object, parm, level = 0.95, ...){
  stopifnot( inherits(object, "lmercens"), is.numeric(level), level > 0L, level < 1L )

  cf <- fixef(object)
  pnames <- names(cf)
  ses <- sqrt(diag(vcov(object)))

  if (missing(parm)) parm <- pnames else if (is.numeric(parm)) parm <- pnames[parm]

  a <- (1L - level)/2L
  a <- c(a, 1L - a)
  qnt <- qnorm(a)

  ci <- array(NA_real_,
              dim = c(length(parm), 2L), dimnames = list(parm, stats:::format.perc(a, 3)))
  ci[] <- cf[parm] + tcrossprod(ses[parm],qnt)

  ci
}
