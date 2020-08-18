# mkuhn, 2017-10-16
# optimization stuff


#' Optimization function for linear mixed models with censoring.
#'
#' @seealso [lme4::optimizeLmer]
#' @param devfun a deviance function
#' @return Results of an optimization.
#' @export
optimizeLmerCens <- function(devfun,
                             optimizer=     formals(lmerControl)$optimizer,
                             restart_edge=  formals(lmerControl)$restart_edge,
                             boundary.tol = formals(lmerControl)$boundary.tol,
                             start = NULL,
                             verbose = 0L,
                             control = list(),
                             ...) {
  verbose <- as.integer(verbose)
  rho <- environment(devfun)


  opt <- optwrapCens(optimizer,
                     devfun,
                     lme4:::getStart(start, pred=rho$pp, returnVal = "all"),
                     lower=rho$lower,
                     control=control,
                     verbose=verbose,
                     ...)


  if (restart_edge) {
    ## FIXME: should we be looking at rho$pp$theta or opt$par
    ##  at this point???  in koller example (for getData(13)) we have
    ##   rho$pp$theta=0, opt$par=0.08
    if (length(bvals <- which(rho$pp$theta==rho$lower)) > 0) {
      ## *don't* use numDeriv -- cruder but fewer dependencies, no worries
      ##  about keeping to the interior of the allowed space
      theta0 <- new("numeric",rho$pp$theta) ## 'deep' copy ...
      d0 <- devfun(theta0)
      btol <- 1e-5  ## FIXME: make user-settable?
      bgrad <- sapply(bvals,
                      function(i) {
                        bndval <- rho$lower[i]
                        theta <- theta0
                        theta[i] <- bndval+btol
                        (devfun(theta)-d0)/btol
                      })
      ## what do I need to do to reset rho$pp$theta to original value???
      devfun(theta0) ## reset rho$pp$theta after tests
      ## FIXME: allow user to specify ALWAYS restart if on boundary?
      if (any(bgrad < 0)) {
        if (verbose) message("some theta parameters on the boundary, restarting")
        opt <- optwrapCens(optimizer,
                           devfun,
                           opt$par,
                           lower=rho$lower, control=control,
                           adj=FALSE, verbose=verbose,
                           ...)
      }
    }
  }
  if (boundary.tol > 0)
    lme4:::check.boundary(rho, opt, devfun, boundary.tol)
  else
    opt
}




#' @title Get the optimizer function and perform minimal checks
#' Internal utility, only used in [optwrapCens]
#' @param optimizer character string ( = function name) *or* function
#' @return optimizer function
#' @seealso [lme4::optwrap]
getOptfun <- function(optimizer) {
  optfun <-
    if (((is.character(optimizer) && optimizer == "optimx") ||
         deparse(substitute(optimizer)) == "optimx")) {
      if (!requireNamespace("optimx")) {
        stop(shQuote("optimx")," package must be installed order to ",
             "use ",shQuote('optimizer="optimx"'))
      }
      optimx::optimx
    } else if (is.character(optimizer)) {
      tryCatch(get(optimizer), error = function(e) NULL)
    } else optimizer

  if (is.null(optfun)) stop("couldn't find optimizer function ",optimizer)
  if (!is.function(optfun)) stop("non-function specified as optimizer")
  needArgs <- c("fn","par","lower","control")
  if (anyNA(match(needArgs, names(formals(optfun)))))
    stop("optimizer function must use (at least) formal parameters ",
         paste(sQuote(needArgs), collapse = ", "))
  optfun
}



#' optwrap for censored responses.
#'
#' Started as a copy of [lme4:::optwrap]
#' but adapted to the censoring situation (e.g. optimization is for all parameters) and removed some options.
#' @param optimizer specifies optimization method, either directly as function object or as function name (as character)
#' @seealso [lme4::optwrap]
optwrapCens <- function(optimizer, fn, par, lower = -Inf, upper = Inf,
                        control = list(), calc.derivs = TRUE,
                        verbose = 0L)
{


  ## control must be specified if adj==TRUE;
  ##  otherwise this is a fairly simple wrapper
  optfun <- getOptfun(optimizer)
  optName <- if (is.character(optimizer)) optimizer
  else ## "good try":
    deparse(substitute(optimizer))[[1L]]

  lower <- rep(lower, length.out = length(par))
  upper <- rep(upper, length.out = length(par))


  switch(optName,
         "bobyqa" = {
           if (all(par == 0)) par[] <- 0.001  ## minor kludge
           if (!is.numeric(control$iprint)) control$iprint <- min(verbose, 3L)
         },
         "Nelder_Mead" = control$verbose <- verbose,
         "nloptwrap" = control$print_level <- min(as.numeric(verbose),3L),
         # mkuhn, pass verbose to optim
         "optim" = control$trace <- verbose,
         ## otherwise:
         if(verbose) warning(gettextf(
           "'verbose' not yet passed to optimizer '%s'; consider fixing optwrapCens()",
           optName), domain = NA)
  )

  arglist <- list(
    # we drop the attributes of the return value of the function [needed for optimx]
    fn = function(param) as.vector(fn(param)),
    par = par, lower = lower, upper = upper, control = control
  )

  ## optimx: must pass method in control (?) because 'method' was previously
  ## used in lme4 to specify REML vs ML
  if (optName == "optimx") {
    if (is.null(method <- control$method))
      stop("must specify 'method' explicitly for optimx")
    arglist$control$method <- NULL

    # mkuhn, add gradient for optimx (if available)
    arglist <- c(arglist, list(method = method, gr = attr(fn, "grad")))
  }


  curWarnings <- list()
  opt <- withCallingHandlers(do.call(optfun, arglist),
                             warning = function(w) {
                               curWarnings <<- append(curWarnings,list(w$message))
                             })

  ## post-fit tweaking
  if (optName == "bobyqa") {
    opt$convergence <- opt$ierr
  }
  else if (optName == "optimx") {
    opt <- list(par = coef(opt)[1,],
                fvalues = opt$value[1],
                method = method,
                conv = opt$convcode[1],
                feval = opt$fevals + opt$gevals,
                message = attr(opt,"details")[,"message"][[1]])
  }
  if ((optconv <- getConv(opt)) != 0) {
    wmsg <- paste("convergence code",optconv,"from",optName)
    if (!is.null(opt$msg)) wmsg <- paste0(wmsg,": ",opt$msg)
    warning(wmsg)
    curWarnings <<- append(curWarnings,list(wmsg))
  }
  ## pp_before <- environment(fn)$pp
  ## save(pp_before,file="pp_before.RData")

  if (calc.derivs) {

    if (verbose > 10) cat("computing derivatives\n")
    derivs <- lme4:::deriv12(fn,opt$par,fx = opt$value)
  } else derivs <- NULL

  ## run one more evaluation of the function at the optimized
  ##  value, to reset the internal/environment variables in devfun ...
  fn(opt$par)

  structure(opt, ## store all auxiliary information
            optimizer = optimizer,
            control   = control,
            warnings  = curWarnings,
            derivs    = derivs)
}



#' get convergence status
#'
#' internal helper function
#' @seealso copy from [lme4:::getConv]
getConv <- function(x) {
  if (!is.null(x[["conv"]])) {
    x[["conv"]]
  } else x[["convergence"]]
}

