
#' Factory for objective function closures for simple LMER-models with censored response.
#'
#' This implements the objective function -- i.e. the negative log-likelihood -- for simple scalar random effect models with censoring.
#' It is implemented in R. The variance parameters and the fixed effect parameter are not profiled out of the log-likelihood
#' because of the censoring.
#' The variance parameters are on log-scale (as to avoid a constrained optimization that variances are non-negative).
#' It uses ML, REML currently not implemented.
#'
#' @param fr model frame with response variable in its first column
#' @param X design matrix of fixed effects
#' @param reTrms list with random effect terms
#' @param control an object of class `lmerControl`
#' @param formula model formula passed from `lmerCens` to get start values
#' @param start list of start values
#' @return objective function which maps the parameters to the corresponding negative log-likelihood
#' @export
mkLmerCensDevfun_rInt_R <- function(fr, X, reTrms, REML = FALSE, verbose = 0, control,
                                    formula = stop("provide model formula"), start, ...){

  stopifnot( is.data.frame(fr), NROW(fr) > 0L )
  stopifnot( is.matrix(X), NROW(X) > 0L )
  stopifnot( inherits(control, "lmerControl"), all(c("quadrature", "quadrature_ord") %in% names(control)) )

  quadrature <- control[["quadrature"]]
  gh_ord <- control[["quadrature_ord"]]

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
  if (!is.null(offset)) {
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

  if ( m > 1L ) stop("only simple random intercept models are implemented!")
  stopifnot( NROW(X) == n, NCOL(Zt) == n)



  # weights -----

  w <- stats::model.weights(fr)
  # sqrtW <- if (! is.null(w) && is.numeric(w))
  #   Matrix::Diagonal(n=n, x=sqrt(w)) else Matrix::Diagonal(n=n)

  # normalize weights (necessary?)
  w <- if (is.null(w) || !is.numeric(w)) rep(1L, n) else w / sum(w)

  #return(list(y=y, X=X, Zt=Zt, Lambdat=Lambdat, reTrms=reTrms, sqrtW=sqrtW))


  # pp quick & start values -----
  lower <- rep(-Inf, 2L)  ## was: p+2L, but it is applied only to the variance parameters later in `optimizeLmerCens`
  pp <- list()


  if (missing(start) || is.null(start)) {
    ##data_start <- data.frame(eval(mc$data))
    # approximate y-vector for regular lmer
    #ZZZ 2017-09-20: add bluntly a column named y_start!! I should make that safe. (e.g. check that name is not used yet)
    #ZZZ problem: lmer has deviance function as a function of theta only, while for my lmercens the deviance is fucntion of beta and theta
    fr$y_start <- flattenResponse(yTime1 = yMat[, "time1"], yTime2 = yMat[, "time2"], yStat = yMat[, "status"])
    # mc$data <- quote(data_start)
    # mc$formula <- quote(update(formula(mc$formula), y_start ~ .))
    ##lmerObj <- eval(mc, parent.frame(1L))

    # new call to lmer for start values
    lmerStart <- do.call(what = lme4::lmer, args = list(data = fr, formula = update(formula, y_start ~ .)))

    pp$theta <- log(as.data.frame(lme4::VarCorr(lmerStart))[, "sdcor"] + .001)
    pp$delb <- lme4::fixef(lmerStart)
  } else {
    stopifnot( is.list(start), all(c("theta", "fixef") %in% names(start)) )
    pp$theta <- start[["theta"]]
    pp$delb <- start[["fixef"]]
  }



  # negative log-likelihood ------

  #' @param param parameter vector, first the variance parameters
  negLogLikFun <- function(param){
    stopifnot( is.numeric(param), length(param) == p + 2L ) # betw SD and residual SD (on log-scale) as extra first two parameter

    # std. deviation parameters are on log-scale
    betwSD <- exp(param[1L])
    resSD <- exp(param[2L])
    beta <- param[-c(1L,2L)]
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
                                            Ztrow = i, betwSD = betwSD, resSD = resSD, ord = gh_ord)
             }
           },
           stats = {
             for (i in 1:q){
               Li[i] <- integrate(intFun, lower=-Inf, upper = Inf, Ztrow = i, betwSD = betwSD, resSD = resSD)$value
             }
           },
           stop("This quadrature is not implemented yet!")

    )

    # fix numeric instabilities
    Li[Li <= 0L] <- .Machine$double.xmin


    retVal <- -sum(log(Li))
    #if (isTRUE(logLik.contribs))
    attr(retVal, "lik.contribs") <- Li

    retVal
  }




  # gradient ----------------------------------------------------------------

  # gradient for neg. log-likelihood function
  #+when Gauß-Hermite quadrature is used
  if (quadrature == 'gh'){

    # given a vector as point in parameter space
    # returns a vector of partial derivatives
    negLogLikGradFun <- function(paramVect){
      stopifnot( length(paramVect) == p + 2L)

      # std. deviation parameters are on log-scale
      betwSD <- exp(paramVect[1L])
      resSD <- exp(paramVect[2L])
      beta <- paramVect[-c(1L, 2L)]
      linPred <- X %*% beta

      # weights and abscissa for GH-quadrature of order gh_ord
      ind <- which(ghQuadRule$order == gh_ord)
      stopifnot( gh_ord >= 2L, length(ind) == gh_ord )

      xi <- ghQuadRule$abscissa[ind]
      wi <- ghQuadRule$weight[ind]

      Li <- attr(negLogLikFun(paramVect), "lik.contribs")
      stopifnot( length(Li) > 0L ) # number of subjects

      # the partial derivatives are a sum of the contributions from the different subjects (e.g. levels of the random effects factor).
      # But per subject, it's not a sum of contributions from individual observations any more because of the (numerical) integration.
      # GH transforms the integration into a sum with `gh_ord` summand terms.
      # We build up an array (dim: gh_ord  x  params  x  subjects), i.e. per subject (q) a matrix with gh_ord rows and p+2 columns
      # we return the sums over the indices of gh_ord and over the individuals (with negative sign because it's negative log-likelihood)

      res <- array(data = 0L,
                   dim = c(length(ind), length(paramVect), length(Li)),
                   dimnames = list(gh_ind = seq_along(ind), params = names(paramVect), subj = seq_along(Li)))


      # calculate gradient
      # per subject i
      for (i in seq_along(Li)) {

        subjInd <- as.vector(Zt[i,]) > 0L
        # Xi <- X[subjInd,]

        for (h in seq_along(ind)) {

          sbpsi <- sqrt(2L) * betwSD * xi[h]

          log_diff_pnorm_int <- logxmy(pnorm(q = yTime2[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD, log.p = TRUE),
                                       pnorm(q = yTime1[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD, log.p = TRUE))

          log_dnorm_int2 <- dnorm(x = yTime2[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD, log = TRUE)
          log_dnorm_int1 <- dnorm(x = yTime1[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD, log = TRUE)

          ## ZZZ weights?
          factor1_ih <- prod(
            #obs
            dnorm(x = yTime1[yStat == 1L & subjInd], mean = linPred[yStat == 1L & subjInd] + sbpsi, sd = resSD),
            #right
            pnorm(q = yTime1[yStat == 0L & subjInd], mean = linPred[yStat == 0L & subjInd] + sbpsi, sd = resSD, lower.tail = FALSE),
            # left
            pnorm(q = yTime1[yStat == 2L & subjInd], mean = linPred[yStat == 2L & subjInd] + sbpsi, sd = resSD, lower.tail = TRUE),
            # interval
            # direct naive implementation
            # pnorm(q = yTime2[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD) -
            #   pnorm(q = yTime1[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD)
            #
            # numerically more robust is:
            exp(sum(log_diff_pnorm_int))
          )

          vector2_ih <- c(
            # par: log(betwSD)
            sbpsi * sum(
              # obs
              ( yTime1[yStat == 1L & subjInd] - linPred[yStat == 1L & subjInd] - sbpsi ) / resSD^2,
              # right
              + dnorm(x = yTime1[yStat == 0L & subjInd], mean = linPred[yStat == 0L & subjInd] + sbpsi, sd = resSD) /
                pnorm(q = yTime1[yStat == 0L & subjInd], mean = linPred[yStat == 0L & subjInd] + sbpsi, sd = resSD, lower.tail = FALSE),
              # left
              - dnorm(x = yTime1[yStat == 2L & subjInd], mean = linPred[yStat == 2L & subjInd] + sbpsi, sd = resSD) /
                pnorm(q = yTime1[yStat == 2L & subjInd], mean = linPred[yStat == 2L & subjInd] + sbpsi, sd = resSD, lower.tail = TRUE),
              # interval: see for interval observations above at par betas
              - exp( log_dnorm_int2 - log_diff_pnorm_int) + exp( log_dnorm_int1 - log_diff_pnorm_int)
            ),

            # par: log(resSD)
            sum(
              # obs
              (yTime1[yStat == 1L & subjInd] - linPred[yStat == 1L & subjInd] - sbpsi)^2 / resSD^2 -1L, # -1L gets recycled and yields - |D|
              # right
              + dnorm(x = yTime1[yStat == 0L & subjInd], mean = linPred[yStat == 0L & subjInd] + sbpsi, sd = resSD) /
                pnorm(q = yTime1[yStat == 0L & subjInd], mean = linPred[yStat == 0L & subjInd] + sbpsi, sd = resSD, lower.tail = FALSE) *
                (yTime1[yStat == 0L & subjInd] - linPred[yStat == 0L & subjInd] - sbpsi),
              # left
              - dnorm(x = yTime1[yStat == 2L & subjInd], mean = linPred[yStat == 2L & subjInd] + sbpsi, sd = resSD) /
                pnorm(q = yTime1[yStat == 2L & subjInd], mean = linPred[yStat == 2L & subjInd] + sbpsi, sd = resSD, lower.tail = TRUE) *
                (yTime1[yStat == 2L & subjInd] - linPred[yStat == 2L & subjInd] - sbpsi),
              # interval
              # naive implementation
              # - ( dnorm(x = yTime2[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD) *
              #       (yTime2[yStat == 3L & subjInd] - linPred[yStat == 3L & subjInd] - sbpsi) -
              #       dnorm(x = yTime1[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD) *
              #       (yTime1[yStat == 3L & subjInd] - linPred[yStat == 3L & subjInd] - sbpsi) ) /
              #   ( pnorm(q = yTime2[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD) -
              #       pnorm(q = yTime1[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD) )

              - crossprod(exp( log_dnorm_int2 - log_diff_pnorm_int), yTime2[yStat == 3L & subjInd] - linPred[yStat == 3L & subjInd] - sbpsi),
              + crossprod(exp( log_dnorm_int1 - log_diff_pnorm_int), yTime1[yStat == 3L & subjInd] - linPred[yStat == 3L & subjInd] - sbpsi)
            ),

            # par: betas
            crossprod(
              c(
                # obs
                (yTime1[yStat == 1L & subjInd] - linPred[yStat == 1L & subjInd] - sbpsi) / resSD^2,
                # right
                + dnorm(x = yTime1[yStat == 0L & subjInd], mean = linPred[yStat == 0L & subjInd] + sbpsi, sd = resSD) /
                  pnorm(q = yTime1[yStat == 0L & subjInd], mean = linPred[yStat == 0L & subjInd] + sbpsi, sd = resSD, lower.tail = FALSE),
                # left
                - dnorm(x = yTime1[yStat == 2L & subjInd], mean = linPred[yStat == 2L & subjInd] + sbpsi, sd = resSD) /
                  pnorm(q = yTime1[yStat == 2L & subjInd], mean = linPred[yStat == 2L & subjInd] + sbpsi, sd = resSD, lower.tail = TRUE),
                # interval
                # direct naive implementation
                # - ( dnorm(x = yTime2[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD) -
                #       dnorm(x = yTime1[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD) ) /
                #   ( pnorm(q = yTime2[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD) -
                #       pnorm(q = yTime1[yStat == 3L & subjInd], mean = linPred[yStat == 3L & subjInd] + sbpsi, sd = resSD) )
                - exp( log_dnorm_int2 - log_diff_pnorm_int) + exp( log_dnorm_int1 - log_diff_pnorm_int)
              ),
              X[c(which(yStat == 1L & subjInd), which(yStat == 0L & subjInd), which(yStat == 2L & subjInd), which(yStat == 3L & subjInd)),]
            )
          )
          stopifnot( length(vector2_ih) == length(paramVect) )

          res[h,,i] <- wi[h] * factor1_ih * vector2_ih
        }# rof h (Gauss-Hermite order)

      }# rof subject i

      # numerical integration as sum:
      #+sum over the Gauß-Hermite indices within each subject
      #+results in matrix: subjects x params
      #+ 1/Li is there because we derive the *log* of likelihood
      #+
      subjMat <- diag(1/(sqrt(pi) * Li)) %*% t(colSums(res))

      # sum up over the subjects
      # negative sign because of *neg* log-likelihood
      - colSums(subjMat)

    } # fun negLogLikGradFun

    # set gradient as attribute
    attr(negLogLikFun, "grad") <- negLogLikGradFun

  }# fi quadrature == 'gh'

  negLogLikFun
}





#' Simple random intercept mixed models with censored response.
#'
#' This function is modelled like the main function \code{\link[lme4]{lmer}}.
#' @param control list-like control object of class \code{\link[lme4]{lmerControl}}. Defaults to use \code{optimx}'s BFGS optimizer.
#' @return lmercens object
#' @export
lmercens <- function(formula, data = NULL, REML, control = lme4::lmerControl(optimizer = "optimx", optCtrl = list(method = "L-BFGS-B")),
                      start = NULL, verbose = 0L, subset, weights, na.action, offset,
                      contrasts = NULL, devFunOnly = FALSE, ...)   {


  stopifnot( !missing(REML) )

  #ZZZ default arguments are not taken into account by match.call
  mcout <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  missCtrl <- missing(control)
  if (!missCtrl && !inherits(control, "lmerControl")) {
    if (!is.list(control))
      stop("'control' is not a list; use lmerControl()")
    warning("passing control as list is deprecated: please use lmerControl() instead",
            immediate. = TRUE)
    control <- do.call(lme4::lmerControl, control)
  }

  # mkuhn, 2017-09-18:
  # handle own control parameter default values
  # Idea: use own sub-class from merControl?
  if (!"quadrature" %in% names(control)) {
    control[["quadrature"]] <- "gh"
  }

  if (!"quadrature_ord" %in% names(control)){
    control[["quadrature_ord"]] <- 8L
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
  # save formula in mcout
  mcout$formula <- lmod$formula
  ##lmod$formula <- NULL  #mkuhn, 20180215: keep it in lmod because I pass the complete lmod to lmercens.fit

  opt <- lmercens.fit(lmod, start = start, verbose = verbose, control = control, devFunOnly = devFunOnly)
  ##names(opt$par) <- c("ln_SD_betw", "ln_SD_res", colnames(lmod$X)) ## name like `(Intercept)` leads to problem with numerical Hessian

  # return value -----
  if (isTRUE(devFunOnly)) return(opt)

  ## ZZZ fix the rho$pp thing of class "merPredD"
  # mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
  #          mc = mcout, lme4conv = cc)


  thetaParInd <- c(1L, 2L)

  # Hessian matrix for negative log-likelihood
  myHess <- estimHessian(opt, thetaParInd)

  # ZZZ for the time being, return own list with optimization information
  structure(  list(call = mcout, ingredients = lmod, par = opt$par, fixef = opt$par[-thetaParInd],
                   fval = opt$fval, hess = myHess,
                   conv = opt$convergence, message = sprintf("call to %s w/ method %s. Resulting code: %d and msg: %s",
                                                             control$optimizer[1L], control$optCtrl[["method"]], opt$ierr, opt$msg),
                   ##"with", paste(c("fn", "gr"), res_optim$counts, sep = ": ", collapse = " - "), "evaluations"),
                   start = start, control = control, negLogLikFun = opt$devfun),
              class = "lmercens"
  )
}

#' Estimate ('observed') Hessian matrix at MLE
#' @param opt optimization result
#' @return Hessian matrix at MLE for fixed-effect coefficients
estimHessian <- function(opt, thetaParInd){

  coef_theta <- opt$par[thetaParInd]
  coef_fixef <- opt$par[-thetaParInd]

  myHess <- attr(opt, "derivs")[["hessian"]]
  # use numeric 2nd derivative at ('observed') MLE if not provided by optimization routine
  if ( is.null(myHess)) {

    negScoreFun <- attr(opt$devfun, "grad", exact = TRUE)
    if (!is.null(negScoreFun) && is.function(negScoreFun)) {
      # for numeric derivation I need individual parameters
      #+I take only the fixef-coefficients and use the observed MLE for the theta parameters
      negScoreFun2 <- function(){ negScoreFun(unlist(c(coef_theta, mget(x = names(coef_fixef))))) }
      formals(negScoreFun2) <- as.list(coef_fixef) #set parameters of the function
      myenv <- new.env() # evaluation environment
      for (fpn in names(coef_fixef)) assign(fpn, as.numeric(coef_fixef[fpn]), envir = myenv)
      myHess <- attr(stats::numericDeriv(str2lang(paste0("negScoreFun2(", paste(names(coef_fixef), collapse = ", "), ")")),
                                         theta = names(coef_fixef), rho = myenv), which = "gradient", exact = TRUE)[-seq_along(coef_theta),]

      # ensure symmetric matrix
      myHess <- (myHess + t(myHess)) / 2L
    }#fi negScoreFun
  }#fi myHess

  myHess
}




#' Refit a lmercens model with a different response.
#'
#' Basically, open up the existing lmercens-fit and inject the new response and call \code{\link{lmercens}} again.
#'
#' @param newresp a \code{\link[survival]{Surv}}-object carrying a censored response
#' @return a lmercens-object from the new censored respone variable
#' @export
refit.lmercens <- function(object, newresp, control = NULL, devFunOnly = NULL, start = NULL, verbose = NULL, ...){
  octrl <- if (! is.null(control)) control else object$control
  ostart <- if (! is.null(start)) start else object$start
  ocall <- object$call
  odevFunOnly <- if (! is.null(devFunOnly)) devFunOnly else  if ("devFunOnly" %in% names(ocall)) eval(ocall$devFunOnly) else FALSE
  overbose <- if (! is.null(verbose)) verbose else if ("verbose" %in% names(ocall)) eval(ocall$verbose) else 0L
  lmod <- object$ingredients

  # update response in model
  lmod$formula <- update(lmod$formula, ysim ~ .)
  attr(lmod$fr, "formula") <- lmod$formula

  ocall$formula <- lmod$formula

  modFrame <- lmod$fr
  modFrame[[1L]] <- newresp
  names(modFrame)[1L] <- "ysim"
  lmod$fr <- modFrame

  opt <- lmercens.fit(lmod, start = ostart, control = octrl, # complex arguments (i.e. with pre-processing in function lmercens) are stored in the object
                      devFunOnly = odevFunOnly, verbose = overbose)

  thetaParInd <- c(1L, 2L)
  myHess <- estimHessian(opt, thetaParInd)
  # return value -----
  structure(  list(call = ocall, ingredients = lmod, par = opt$par, fixef = opt$par[-thetaParInd],
                   fval = opt$fval, hess = myHess,
                   conv = opt$convergence, message = sprintf("call to %s w/ method %s. Resulting code: %d and msg: %s",
                                                             octrl$optimizer[1L], octrl$optCtrl[["method"]], opt$ierr, opt$message),
                   ##"with", paste(c("fn", "gr"), res_optim$counts, sep = ": ", collapse = " - "), "evaluations"),
                   start = ostart, control = octrl, negLogLikFun = opt$devfun),
              class = "lmercens"
  )
}




#' Internal lmercens-fitting routine.
#' @param lmod the description of the model as prepared by \code{\link[lme4]{lFormula}}
#' @param start list of start values with entries `theta` and `fixef`. Or NULL
#' @param control a list that carries control-settings for the optimization, see \code{\link[lme4]{lmerControl}}
#' @param devFunOnly flag if we want to have only the deviance function
#' @return result of optimization call (or deviance function)
lmercens.fit <- function(lmod, start = NULL, verbose = 0L, control = lme4::lmerControl(), devFunOnly){
  # objective function ------
  p <- NCOL(lmod$X)
  stopifnot( !is.null(lmod$reTrms), !is.null(lmod$reTrms$cnms) )

  if (length(lmod$reTrms$cnms) > 1L || lmod$reTrms$cnms[[1L]] != "(Intercept)" ) {
    stop("Only simple random intercept models are supported here (for now?)!")
  }

  ##mkuhn, 20170920: provide the formula for mkLmerCensDevfun_rInt_R to have the start-value code work (quick-pp)
  devfun <- do.call(mkLmerCensDevfun_rInt_R, c(lmod, list(verbose = verbose, control = control, start = start))) #mcout$

  if (isTRUE(devFunOnly)) return(devfun)

  #negLogLikGradFun <- attr(devfun, "grad")
  ##if ( ! is.null(negLogLikGradFun)) cat("\nWe have a gradient!\n")



  #mkuhn, 20170920: start value code is for now in the deviance function


  # optimization -----
  if (identical(control$optimizer, "none"))
    stop("deprecated use of optimizer='none'; use NULL instead")

  opt <- if (length(control$optimizer) == 0L) {
    #stop("start values are required if no optimization")
    #s <- getStart(start, environment(devfun)$lower, environment(devfun)$pp)
    list(par = start, fixef = start[["fixef"]], fval = devfun(start), conv = 1000, message = "no optimization", negLogLikFun = devfun)
  } else {
    # delegate optimization
    optimizeLmerCens(devfun, optimizer = control$optimizer, restart_edge = control$restart_edge,
                     boundary.tol = control$boundary.tol, control = control$optCtrl,
                     verbose = verbose, start = start, calc.derivs = control$calc.derivs)

    # ZZZ adapt this old commented code  e.g. control$optimizer == 'optimx' and use instead: method = control$optCtr$method
    # stopifnot( length(control$optimizer) >= 1L, is.character(control$optimizer),
    #            pmatch(x = control$optimizer[1L], table = eval(formals(stats::optim)$method), nomatch = -1L) > 0L )
    #
    # res_optim <- optim(par = start, fn = devfun, gr = negLogLikGradFun, hessian = TRUE, method = control$optimizer[1L], control = list(trace = 1L))
    # list(par = res_optim$par, fixef = res_optim$par[1L:p], fval = res_optim$value, hess=res_optim$hessian,
    #      conv = res_optim$convergence, message = paste("call to optim w/ method", control$optimizer[1L], res_optim$message,
    #                                                    "with", paste(c("fn", "gr"), res_optim$counts, sep = ": ", collapse = " - "), "evaluations"),
    #      start = start, negLogLikFun = devfun)
  }
  cc <- NULL
  # cc <- checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,
  #                 lbound = environment(devfun)$lower)
  opt[["devfun"]] <- devfun

  opt
}

