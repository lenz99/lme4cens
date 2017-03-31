# Censored observations in linear models


#' Prepare the survival response.
#' It uses the interval-coding that supports left-, right- and interval-censoring.
#' @param y0 response vector of class [survival::Surv].
#' @return response vector, normalized to interval-coding.
prepSurvResp <- function(y0) {
  stopifnot( inherits(y0, what = "Surv") )

  survType <- attr(y0, which = "type")

  if (survType == 'interval') y0 else {
    yMat <- as.matrix(y0)
    yTime <- yMat[, "time"]
    yStat <- yMat[, "status"]
    censLevel <- if (max(yStat) == 2) 1 else 0
    isCens <- (yStat == censLevel)

    switch(survType,
           right={Surv(time  = yTime,
                       time2 = ifelse(isCens, yes = NA_real_, no = yTime),
                       type = "interval2")},
           left={ Surv(time  = ifelse(isCens, yes = NA_real_, no = yTime),
                       time2 = yTime,
                       type = "interval2")},
           stop("this type of censoring is not supported!")
    )
  }
}


#' Fit a linear model with censored observations.
#'
#' Residuals are not implemented, yet. What is it for censored observations?
#' @seealso [stats::lm]
#' @export
lmcens <- function(formula, data, subset, weights, contrasts = NULL, offset){

  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data", "subset", "weights", "offset"),
             table = names(mf), nomatch = 0L)

  # model.frame -----
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, envir = parent.frame())

  mt <- attr(mf, "terms")

  if (is.empty.model(mt)){
    stop("Empty model provided!")
  }

  # response variable ------
  y <- prepSurvResp(stats::model.response(mf, type = "any"))


  w <- stats::model.weights(mf)
  if (! is.null(w) && !is.numeric(w))
    stop("'weights' must translate into a numeric vector!")

  offset <- as.vector(stats::model.offset(mf))
  if (!is.null(offset) && length(offset) != NROW(y) )
    stop(gettextf("number of offsets is %d, should equal number of obs %d",
                  length(offset), NROW(y), domain = NA))


  # model matrix ------
  x <- stats::model.matrix(mt, mf, contrasts)


  z <- lmcens.fit(x, y, w, offset = offset)


  # return value ------
  class(z) <- c("lmcens", "lm")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- match.call()
  z$terms <- mt

  z$model <- x

  z
}


#' Internal workhorse function for linear regression with censored observations.
#'
#' It uses maximum likelihood famework to come up with estimates.
#' Optimization via [stats::optim()].
#'
#' @param x model matrix
#' @param y response that uses the interval-type censoring
#' @param w weights vector
#' @param offset offset vector that is subtracted from the time variable (in
#'   case of interval-censoring both boundaries are adapted)
#' @return list with regression fit stuff (e.g. coefficients, fitted.values
#'   effects, rank, ..)
lmcens.fit <- function(x, y, w, offset = NULL, tol = 1e-07){

  stopifnot( inherits(y, what = "Surv"), attr(y, which = "type") == "interval" )

  yMat <- as.matrix(y)
  stopifnot( identical(colnames(yMat), c("time1", "time2", "status")) )


  n <- NROW(yMat)
  yTime1 <- yMat[, "time1"]
  yTime2 <- yMat[, "time2"]
  yStat <- yMat[, "status"]

  # apply offset to y-variable
  if (!is.null(offset)){
    yTime1 <- yTime1 - offset
    yTime2 <- yTime2 - offset
  }


  # weights -----
  weighting <- ! missing(w) && ! is.null(w)
  if ( weighting ){
    if (any(w < 0 | is.na(w)))
      stop("missing or negative weights not allowed")

  } else {
    w <- rep(1L, n)
  }
  # normalize weights to sum to n
  w <- w/sum(w) * n
  stopifnot( length(w) == n )

  p <- NCOL(x)
  qr_x <- qr(crossprod(x))


  #' negative log-likelihood
  negLogLikFun <- function(beta){
    stopifnot( length(beta) == p+1L ) # residual log(σ) as extra parameter
    linPred <- x %*% beta[1:p]

    resSD <- exp(beta[p+1L])

    - sum(
      # point observations
      w[yStat == 1] * dnorm(x = yTime1[yStat == 1], mean = linPred[yStat == 1], sd = resSD, log = TRUE),
      # right cens
      w[yStat == 0] * pnorm(q = yTime1[yStat == 0], mean = linPred[yStat == 0], sd = resSD, lower.tail = FALSE, log.p = TRUE),
      # left cens
      w[yStat == 2] * pnorm(q = yTime1[yStat == 2], mean = linPred[yStat == 2], sd = resSD, lower.tail = TRUE, log.p = TRUE),
      # interval cens
      w[yStat == 3] * (log(pnorm(q = yTime2[yStat == 3], mean = linPred[yStat == 3], sd = resSD) -
            pnorm(q = yTime1[yStat == 3], mean = linPred[yStat == 3], sd = resSD)))
    )
  }


  # hessian relates to the negative log-likelihood function
  beta_init <- setNames(c(median(yTime1), rep(0, p)),
                        c(colnames(x), "lSigma"))
  optimRes <- optim(par = beta_init, fn = negLogLikFun, hessian = TRUE)

  fit_coef <- optimRes$par
  fit_vals <- x %*% fit_coef[1:p]
  if (!is.null(offset))
    fit_vals <- fit_vals + offset

  retVal <- list(coefficients = fit_coef, logLik=-optimRes$value, hess=optimRes$hessian,
       fitted.values = fit_vals, y=y,
       rank = p, qr = qr_x,
       df.residual = NROW(x) - p - 1L)

  if (weighting) retVal <- append(retVal, values = list(weights = w))

  retVal
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

  cat("\nResidual standard error:", format(signif(x$sigma, digits)),
      "on", rdf, "degrees of freedom\n")
  cat("Log-Likelihood: ", format(signif(x$logLik, digits)))
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
