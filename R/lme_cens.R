
#' Factory for objective function closure for simple LMER-models with censored response for minimization.
#'
#' This implements the objective function for simple scalar random effect models with censoring.
#' It is implemented in R. The variance parameter and the fixed effect parameter is not profiled out because of the censoring.
#' It uses ML, REML not implemented.
#'
#' @export
#' @param fr model frame with response variable in its first column
#' @param reTrms list with random effect terms
#' @param quadrature type of numeric integration method
#' @return objective function which maps the parameters to the corresponding negative log-likelihood
mkLmerCensDevfun_rInt_R <- function(fr, X, reTrms, REML = FALSE, verbose = 0, quadrature = c("gh", "stats"), ...){

  quadrature <- match.arg(quadrature)

  stopifnot( is.data.frame(fr), NROW(fr) > 0L )
  stopifnot( is.matrix(X), NROW(X) > 0L )

  if (isTRUE(REML)) stop("Restricted maximum likelihood not implemented, -- only ML.")

  y <- prepSurvResp(stats::model.response(fr, type = "any"))
  stopifnot( inherits(y, "Surv"), attr(y, which = "type") == "interval" )

  yMat <- as.matrix(y)
  stopifnot( identical(colnames(yMat), c("time1", "time2", "status")) )


  n <- NROW(yMat)
  yTime1 <- yMat[, "time1"]
  yTime2 <- yMat[, "time2"]
  yStat <- yMat[, "status"]

  offset <- as.vector(stats::model.offset(fr))

  # apply offset to y-variable
  #ZZZ delegate this into the deviance function! (like for `lmcens`)
  if (!is.null(offset)){
    yTime1 <- yTime1 - offset
    yTime2 <- yTime2 - offset
  }


  Zt <- reTrms[["Zt"]]
  Lambdat <- reTrms[["Lambdat"]]

  theta <- reTrms[["theta"]]
  mapping <- function(theta) theta[reTrms[["Lind"]]]

  p <- NCOL(X)
  q <- NROW(Zt)
  m <- length(theta)

  if ( m > 1L ) stop("only simple-random intercept models are implemented!")
  stopifnot( NROW(X) == n, NCOL(Zt) == n)



  # weights -----

  w <- stats::model.weights(fr)
  # sqrtW <- if (! is.null(w) && is.numeric(w))
  #   Matrix::Diagonal(n=n, x=sqrt(w)) else Matrix::Diagonal(n=n)

  # W <- W / sum(W)

  # # normalize weights (necessary?)
  w <- if (is.null(w) || ! is.numeric(w)) rep(1L, n) else w / sum(w)

  #return(list(y=y, X=X, Zt=Zt, Lambdat=Lambdat, reTrms=reTrms, sqrtW=sqrtW))


  # negative log-likelihood ------
  negLogLikFun <- function(param){
    stopifnot( length(param) == p+2L ) # betw SD and residual SD (on log-scale) as extra parameter

    beta <- param[1L:p]
    # std. deviation parameters are on log-scale
    betwSD <- exp(param[p+1L])
    resSD <- exp(param[p+2L])

    linPred <- X %*% beta


    # integrand function for likelihood contribution for a subject
    intFun <- Vectorize(function(mu, Ztrow, betwSD, resSD) {
      subjInd <- as.vector(Zt[Ztrow,]) > 0L

      ## QQQ /!\ weighting happens here on likelihood scale (for linear models it was on log-likelihood scale). Think over!
      prod(
        # point obs  #> factor 1/resSD is part of dnorm!
        w[yStat == 1 & subjInd] * dnorm(x = yTime1[yStat == 1 & subjInd], mean = linPred[yStat == 1 & subjInd] + mu, sd = resSD),
        # right cens
        w[yStat == 0 & subjInd] * pnorm(q = yTime1[yStat == 0 & subjInd], mean = linPred[yStat == 0 & subjInd] + mu, sd = resSD, lower.tail = FALSE),
        # left cens
        w[yStat == 2 & subjInd] * pnorm(q = yTime1[yStat == 2 & subjInd], mean = linPred[yStat == 2 & subjInd] + mu, sd = resSD, lower.tail = TRUE),
        # interval cens
        w[yStat == 3 & subjInd] * (pnorm(q = yTime2[yStat == 3 & subjInd], mean = linPred[yStat == 3 & subjInd] + mu, sd = resSD) -
                                     pnorm(q = yTime1[yStat == 3 & subjInd], mean = linPred[yStat == 3 & subjInd] + mu, sd = resSD)),
        if (quadrature == 'stats') dnorm(x=mu, sd = betwSD) else 1L
      )

    }, vectorize.args = "mu")

    ##-sum(sapply(1:q, function(i) log(integrate(intFun, lower=-Inf, upper = Inf, Ztrow = i, betwSD = betwSD, resSD = resSD)$value)))
    Li <- vector("numeric", length = q)

    switch(quadrature,
           gh = {
              for (i in 1:q){
                Li[i] <- 1/sqrt(pi) * int_gh(f = function(mu, Ztrow, betwSD, resSD) intFun(mu = sqrt(2) * betwSD * mu, Ztrow, resSD = resSD),
                                             Ztrow = i, betwSD = betwSD, resSD = resSD)
              }
           },
           stats = {
             for (i in 1:q){
               Li[i] <- integrate(intFun, lower=-Inf, upper = Inf, Ztrow = i, betwSD = betwSD, resSD = resSD)$value
             }
           },
           stop("This quadrature is not implemented yet!")

    )


    retVal <- -sum(log(Li+.Machine$double.xmin))
    attr(retVal, "lik.contribs") <- Li

    retVal
  }

  negLogLikFun
}




#' Simple random intercept mixed models with censored response.
#'
#' This function is modelled after the official function `lmer`.
#' @export
lmercens <- function (formula, data = NULL, REML = TRUE, control = lmerControl(),
                      start = NULL, verbose = 0L, subset, weights, na.action, offset,
                      contrasts = NULL, devFunOnly = FALSE, quadrature = c("gh", "stats"), ...)   {

  quadrature <- match.arg(quadrature)

  mcout <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  missCtrl <- missing(control)
  if (!missCtrl && !inherits(control, "lmerControl")) {
    if (!is.list(control))
      stop("'control' is not a list; use lmerControl()")
    warning("passing control as list is deprecated: please use lmerControl() instead",
            immediate. = TRUE)
    control <- do.call(lmerControl, control)
  }
  if (!is.null(list(...)[["family"]])) {
    stop("calling lmercens with 'family' is not supported!")
    ## original code from `lmer`:
    # mc[[1]] <- quote(lme4::glmer)
    # if (missCtrl)
    #   mc$control <- glmerControl()
    # return(eval(mc, parent.frame(1L)))
  }
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset", "REML"), names(mc), nomatch = 0L)
  mc <- mc[c(1L, m)]
  mc$control <- control
  mc[[1]] <- quote(lme4::lFormula)
  lmod <- eval(mc, parent.frame(1L))
  mcout$formula <- lmod$formula
  lmod$formula <- NULL

  p <- NCOL(lmod$X)

  # objective function ------
  stopifnot( ! is.null(lmod$reTrms), ! is.null(lmod$reTrms$cnms) )

  if (length(lmod$reTrms$cnms) > 1L || lmod$reTrms$cnms[[1L]] != "(Intercept)" ){
    stop("Only simple random intercept models are supported!")
  }

  devfun <- do.call(mkLmerCensDevfun_rInt_R, c(lmod, list(quadrature = quadrature,
                                                          verbose = verbose, control = control)))
  if (devFunOnly)
    return(devfun)



  # start value ------
  start <- if (is.null(start) || ! length(start) %in% c(p, p+2L)){

    y <- prepSurvResp(stats::model.response(lmod$fr, type = "any"))
    stopifnot( inherits(y, "Surv"), attr(y, which = "type") == "interval" )

    yMat <- as.matrix(y)
    stopifnot( identical(colnames(yMat), c("time1", "time2", "status")) )



    # operate with call to update formula is complicated..
    # mc[[1]] <- quote(lme4::lmer)
    data_start <- data.frame(eval(mc$data))
    # approximate y-vector for regular lmer
    data_start$y_start <- flattenResponse(yTime1 = yMat[, "time1"], yTime2 = yMat[, "time2"], yStat = yMat[, "status"])
    # mc$data <- quote(data_start)
    # mc$formula <- quote(update(formula(mc$formula), y_start ~ .))
    ##lmerObj <- eval(mc, parent.frame(1L))

    # new call to lmer for start values
    lmerStart <- do.call(what = lme4::lmer, args = list(data = data_start, formula = update(formula, y_start ~ .)))
    c(fixef(lmerStart), log(as.data.frame(lme4::VarCorr(lmerStart))[, "sdcor"]+.001))

  } else {
    if (length(start) == p) start <- c(start, 0, 0)
    start
  }

  stopifnot( is.numeric(start) )


  # optimization -----
  if (identical(control$optimizer, "none"))
    stop("deprecated use of optimizer=='none'; use NULL instead")

  opt <- if (length(control$optimizer) == 0L) {
    stop("start values are required if no optimization")
    #s <- getStart(start, environment(devfun)$lower, environment(devfun)$pp)
    list(par = s, fixef = s[1L:p], fval = devfun(s), conv = 1000, message = "no optimization")
  }
  else {
    res_optim <- optim(par = start, fn = devfun)
    list(par = res_optim$par, fixef = res_optim$par[1L:p], fval = res_optim$value, conv = res_optim$convergence, message = "call to optim", start = start)
    # optimizeLmer(devfun, optimizer = control$optimizer, restart_edge = control$restart_edge,
    #              boundary.tol = control$boundary.tol, control = control$optCtrl,
    #              verbose = verbose, start = start, calc.derivs = control$calc.derivs,
    #              use.last.params = control$use.last.params)
  }
  cc <- NULL
  # cc <- checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,
  #                 lbound = environment(devfun)$lower)


  # return value -----
  ## ZZZ fix the rho$pp thing of class "merPredD"
  # mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
  #          mc = mcout, lme4conv = cc)
  opt
}



