# Censored observations in linear models


#' Fit a linear model with censored observations.
#'
#' Optimization via [stats::optim].
#' Residuals are not implemented, yet. Which type of residuals are best for censored observations?
#' @seealso [stats::lm]
#' @param start numeric vector of start parameters. If `NULL` use ordinary linear regression for start values
#' @param offset offset vector that is subtracted from the response variable (in
#'   case of interval-censoring both boundaries are adapted)
#' @param method optimization method used by [stats::optim]. Defaults to BFGS (as we have analytical gradient).
#' @param ... further arguments passed to [stats::optim].
#' @return list with regression fit stuff (e.g. coefficients, fitted.values
#'   effects, rank, ..)
#' @export
lmcens <- function(formula, data, subset, weights, contrasts = NULL, offset = NULL, start = NULL,
                   method = c("BFGS", "L-BFGS-B", "Nelder-Mead", "SANN", "CG"), ...){

  method <- match.arg(method)
  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data", "subset", "weights", "offset"),
             table = names(mf), nomatch = 0L)

  # model.frame and model.matrix -----
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, envir = parent.frame())

  mt <- attr(mf, "terms")

  if (is.empty.model(mt)){
    stop("Empty model provided!")
  }

  x <- stats::model.matrix(mt, mf, contrasts)
  p <- NCOL(x)



  # response variable ------
  y <- prepSurvResp(stats::model.response(mf, type = "any"))

  stopifnot( inherits(y, what = "Surv"), attr(y, which = "type") == "interval" )

  yMat <- as.matrix(y)
  stopifnot( identical(colnames(yMat), c("time1", "time2", "status")) )

  n <- NROW(yMat)
  yTime1 <- yMat[, "time1"]
  yTime2 <- yMat[, "time2"]
  yStat <- yMat[, "status"]


  # weights and offset ------
  w <- stats::model.weights(mf)
  if (! is.null(w) && ! is.numeric(w))
    stop("'weights' must translate into a numeric vector!")

  offset <- as.vector(stats::model.offset(mf))
  if (!is.null(offset) && length(offset) != NROW(y) )
    stop(gettextf("number of offsets is %d, should equal number of obs %d",
                  length(offset), NROW(y), domain = NA))





  # objective function -----
  negLogLikFun <- lmcens.objFun(x, yTime1, yTime2, yStat, w, offset)
  negLogLikGradFun <- attr(negLogLikFun, "grad")


  # start values ----
  start <- if (is.null(start) || ! length(start) %in% c(p, p+1L)){

    # quick en dirty start values:
    # setNames(c(mean(yTime1), rep(0L, p)),
    #          c(colnames(x), "lSigma"))

    lmy <- flattenResponse(yTime1 = yTime1, yTime2 = yTime2, yStat = yStat)

    lmfit <- if (is.null(w))
      lm.fit(x, lmy, offset = offset, singular.ok = TRUE, method = "qr")
    else lm.wfit(x, lmy, w, offset = offset, singular.ok = TRUE, method = "qr")

    setNames(c(lmfit[["coefficients"]], log(sd(lmy) + 0.001)),
             c(colnames(x), "lSigma"))

  } else {
    if (length(start) == p) start <- c(start, 0L)
    start
  }


  # optimization -----
  # optimize the negative log-likelihood function
  res_optim <- optim(par = start, fn = negLogLikFun, gr = negLogLikGradFun, hessian = TRUE, method = method, ...)


  # return value ------
  fit_coef <- res_optim$par
  fit_vals <- x %*% fit_coef[1:p]
  if (!is.null(offset))
    fit_vals <- fit_vals + offset



  z <- list(coefficients = fit_coef, logLik=-res_optim$value,
            start = start,
            logLik.contribs = negLogLikFun(fit_coef),
            # hessian relates to the negative log-likelihood function
            hess=res_optim$hessian,
            fitted.values = fit_vals, y=y,
            rank = p, qr = qr(crossprod(x)),
            df.residual = NROW(x) - (p + 1L),
            negLogLikFun = negLogLikFun)

  if (! is.null(w)) z <- append(z, values = list(weights = w))


  class(z) <- c("lmcens", "lm")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- match.call()
  z$terms <- mt

  z$model <- x

  z
}


#' Objective function for linear regression with censored observations
#'
#' Factory method to create the negative log-likelihood function. The
#' (analytical) gradient is provided for improved performance for optimization.
#'
#' @param x model matrix nxp
#' @param yTime1 first response time
#' @param yTime2 second response time
#' @param yStat status variable in interval style
#' @param w weights vector
#' @param offset offset vector, can be `NULL`
#' @return negative log-likelihood as objective function. Attribute \code{grad} contains the analytical gradient.
#'
#' @export
lmcens.objFun <- function(x, yTime1, yTime2, yStat, w, offset){

  n <- length(yTime1)
  p <- NCOL(x)


  # weights -----
  weighting <- ! missing(w) && ! is.null(w)
  if ( weighting ){
    if (any(w < 0 | is.na(w))) stop("missing or negative weights not allowed")

  } else {
    w <- rep(1L, n)
  }
  # normalize weights to sum to n
  w <- w/sum(w) * n
  stopifnot( length(w) == n )


  # offset ---
  if (!is.null(offset)){
    yTime1 <- yTime1 - offset
    yTime2 <- yTime2 - offset
  }

  # negative log-likelihood function to minimize.
  # weights are applied to the log-likelihood contributions.
  # returns negative log-likelihood for given parameter vector for data sample, with positive likelihood contribution of individual observations as attribute
  negLogLikFun <- function(paramVect){
    stopifnot( is.numeric(paramVect), length(paramVect) == p+1L ) # residual log(σ) as extra parameter
    linPred <- x %*% paramVect[1L:p]

    resSD <- exp(paramVect[p+1L])


    logLiks <- numeric(length = length(yStat))

    # point observations
    logLiks[yStat == 1L] <- w[yStat == 1L] * dnorm(x = yTime1[yStat == 1L], mean = linPred[yStat == 1L], sd = resSD, log = TRUE)
    # right cens
    logLiks[yStat == 0L] <- w[yStat == 0L] * pnorm(q = yTime1[yStat == 0L], mean = linPred[yStat == 0L], sd = resSD, lower.tail = FALSE, log.p = TRUE)
    # left cens
    logLiks[yStat == 2L] <- w[yStat == 2L] * pnorm(q = yTime1[yStat == 2L], mean = linPred[yStat == 2L], sd = resSD, lower.tail = TRUE, log.p = TRUE)
    # interval cens
    ln_t2_intCens <- pnorm(q = yTime2[yStat == 3L], mean = linPred[yStat == 3L], sd = resSD, log.p = TRUE)
    ln_t1_intCens <- pnorm(q = yTime1[yStat == 3L], mean = linPred[yStat == 3L], sd = resSD, log.p = TRUE)
    logLiks[yStat == 3L] <- w[yStat == 3L] * (ln_t2_intCens + log1mexp(ln_t1_intCens - ln_t2_intCens))


    retVal <-   - sum(logLiks)
    attr(retVal, "loglik.contribs") <- logLiks

    retVal
  }

  # gradient of negative log-likelihood function
  # weights are on the log-likelihood scale
  # returns gradient numeric vector of length equal to number of parameters
  negLogLikGradFun <- function(paramVect){
    stopifnot( is.numeric(paramVect), length(paramVect) == p+1L ) # residual log(σ) as extra parameter
    linPred <- x %*% paramVect[1:p] # nx1 vector

    resSD <- exp(paramVect[p+1L])


    # the partial derivatives sum up from the individual observations
    # we build up a matrix with n rows and p+1 columns (the order of observations [=rows] is not retained, though)
    # we return the column sums (with negative sign because it is negative log-lik function)

    # point obs
    # weights are not fused here with the resid_obs-vector as it is used as sum-of-square and would be squared
    resid_obs <- (yTime1[yStat == 1L] - linPred[yStat == 1L])
    contrib_obs <- c(crossprod(x[yStat == 1L,], w[yStat == 1L] * resid_obs)[,1L] / resSD^2, crossprod(w[yStat == 1L] * resid_obs, resid_obs) / resSD^2 - sum(w[yStat == 1L]))
    # right cens
    factor_right <- w[yStat == 0L] * dnorm(x = yTime1[yStat == 0L], mean = linPred[yStat == 0L], sd = resSD) / pnorm(q = linPred[yStat == 0L], mean = yTime1[yStat == 0L], sd = resSD, lower.tail = TRUE)
    contrib_right <- c( crossprod(x[yStat == 0L,], factor_right), crossprod(yTime1[yStat == 0L] - linPred[yStat == 0L], factor_right) )
    # left cens
    factor_left <- - w[yStat == 2L] * dnorm(x = yTime1[yStat == 2L], mean = linPred[yStat == 2L], sd = resSD) / pnorm(q = yTime1[yStat == 2L], mean = linPred[yStat == 2L], sd = resSD, lower.tail = TRUE)
    contrib_left <- c( crossprod(x[yStat == 2L,], factor_left), crossprod(yTime1[yStat == 2L] - linPred[yStat == 2L], factor_left) )
    # interval cens
    factor_int <- w[yStat == 3L] / (pnorm(q = yTime2[yStat == 3L], mean = linPred[yStat == 3L], sd = resSD, lower.tail = TRUE) - pnorm(q = yTime1[yStat == 3L], mean = linPred[yStat == 3L], sd = resSD, lower.tail = TRUE))
    contrib_int <- c( crossprod(-x[yStat == 3L,], (dnorm(x = yTime2[yStat == 3L], mean = linPred[yStat == 3L], sd = resSD) - dnorm(x = yTime1[yStat == 3L], mean = linPred[yStat == 3L], sd = resSD)) * factor_int),
                      crossprod( -dnorm(x = yTime2[yStat == 3L], mean = linPred[yStat == 3L], sd = resSD), (yTime2[yStat == 3L] - linPred[yStat == 3L])*factor_int) +
                        crossprod(dnorm(x = yTime1[yStat == 3L], mean = linPred[yStat == 3L], sd = resSD), (yTime1[yStat == 3L] - linPred[yStat == 3L])*factor_int) )

    - colSums(rbind(contrib_obs, contrib_right, contrib_left, contrib_int), na.rm = TRUE)
  }

  attr(negLogLikFun, "grad") <- negLogLikGradFun

  negLogLikFun
}


#' @export
summary.lmcens <- function(object){
  z <- object

  lSigmaCoef <- z$coefficients[length(z$coefficients)]

  est <- z$coefficients[-length(z$coefficients)]
  # hessian is for the negative log-likelihood function
  se <- sqrt(diag(vcov(object))[-length(z$coefficients)])
  tval <- est/se
  rdf <- z$df.residual

  ans <- z[c("call", "terms", "logLik", if (!is.null(z$weights)) "weights")] #weights last
  ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval),
                                                  rdf, lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(names(est),
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

  ans$sigma <- exp(lSigmaCoef)
  ans$df <- c(z$rank, z$df.residual, z$qr$rank)  ## think here: what df values are needed?!

  class(ans) <- "summary.lmcens"

  ans
}


#' @export
print.summary.lmcens <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = FALSE, ...){
  df <- x$df
  rdf <- df[2L]

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("\nCoefficients:\n")
  coefs <- x$coefficients
  stats::printCoefmat(coefs, digits = digits, na.print = "NA", signif.stars = signif.stars, ...)

  cat("\nResidual standard error:", format(round(x$sigma, digits)), "(",format(round(log(x$sigma), digits)), "on log-scale)",
      "on", rdf, "degrees of freedom\n")
  cat("Log-Likelihood: ", format(round(x$logLik), digits))
  cat("\n")

  cat("\n")
}

#' @export
logLik.lmcens <- function(obj) obj$logLik


#' Variance covariance function for [lmcens]-models.
#'
#' Includes the residual standard error as own parameter
#' @return asymptotic variance-covariance matrix for the parameters, including residual std error as last parameter
#' @export
vcov.lmcens <- function(obj) solve(obj$hess)
