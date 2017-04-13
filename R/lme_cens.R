
#' Objective function closure for LMER-models with censored response.
#'
#' This implements the objective function for simple scalar models with censoring.
#' It is implemented in R. The variance parameter and the fixed effect parameter is not profiled out because of the censoring.
#' It uses ML, REML not implemented.
#'
#' @export
#' @param fr dataframe with response variable in its first column
#' @param reTrms list with random effect terms
#' @return objective function which maps the parameters to the corresponding negative log-likelihood
mkLmerCensDevfun_rInt_R <- function(fr, X, reTrms, REML = FALSE, start = NULL, verbose = 0, quadrature = c("gh", "stats"), ...){

  quadrature <- match.arg(quadrature)

  stopifnot( is.data.frame(fr), NROW(fr) > 0L )
  stopifnot( is.matrix(X), NROW(X) > 0L )

  if (isTRUE(REML)) stop("Restricted maximum likelihood not implemented, -- only ML.")

  y <- prepSurvResp(stats::model.response(fr, type = "any"))
  stopifnot( inherits(y, "Surv") )

  yMat <- as.matrix(y)
  stopifnot( identical(colnames(yMat), c("time1", "time2", "status")) )


  n <- NROW(yMat)
  yTime1 <- yMat[, "time1"]
  yTime2 <- yMat[, "time2"]
  yStat <- yMat[, "status"]

  offset <- as.vector(stats::model.offset(fr))

  # apply offset to y-variable
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

  if ( m > 1 ) stop("only simple-random intercept models are implemented!")
  stopifnot( NROW(X) == n, NCOL(Zt) == n)

  # raw weights matrix
  w <- stats::model.weights(fr)
  # sqrtW <- if (! is.null(w) && is.numeric(w))
  #   Matrix::Diagonal(n=n, x=sqrt(w)) else Matrix::Diagonal(n=n)

  # W <- W / sum(W)

  # # normalize weights (necessary?)
  w <- if (is.null(w) || ! is.numeric(w)) rep(1L, n) else w / sum(w)

  #return(list(y=y, X=X, Zt=Zt, Lambdat=Lambdat, reTrms=reTrms, sqrtW=sqrtW))


  # negative log-likelihood
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

      ## /!\ weighting happens here on likelihood scale (for linear models it was on log-likelihood scale). Think over!
      prod(
        # point obs  ### 1/resSD *  <-- but this is part of dnorm!
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
                Li[i] <- 1/sqrt(pi) * int_gh(f = function(mu, Ztrow, betwSD, resSD) intFun(mu = sqrt(2) * betwSD * mu, Ztrow, betwSD=betwSD, resSD = resSD),
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


    -sum(log(Li+.Machine$double.xmin))
  }

  negLogLikFun
}

